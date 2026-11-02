unit RTheory.RatSparse;

{ Sparse exact linear algebra over Q (TRational from
  RTheory.ExactField), plus dense small-matrix helpers over Q.

  Sparse rows are sorted lists of (Col, Val) with no zero entries.
  The Kronecker-system builder appends the 256 rows of
  kron(P, I) - kron(I, Q) for dense 16x16 P, Q, in Sage's
  tensor_product convention: row index (i,j) -> i*16+j, column index
  (k,l) -> k*16+l, (A x B)[(i,j),(k,l)] = A[i,k]*B[j,l]; the unknown
  vector is vec(X) row-major, component (r,c) at index r*16+c.

  SparseRREF does full Gauss-Jordan (forward elimination +
  back-substitution, pivot column eliminated from ALL rows), ending
  in true reduced row echelon form. SparseKernelBasis returns one
  dense vector per free column in the standard RREF form, matching
  Sage right_kernel().basis() semantics (the basis matrix whose rows
  are in RREF). Everything is exact rational arithmetic -- no floating
  point anywhere. }

{$mode objfpc}{$H+}

interface

uses
  SysUtils, RTheory.ExactField;

type
  TRatMatrix = array of array of TRational;
  TRatMatrices = array of TRatMatrix;
  TRatVector = array of TRational;
  TRatVectors = array of TRatVector;

  TSparseEntry = record
    Col: Integer;
    Val: TRational;
  end;
  TSparseRow = array of TSparseEntry; { sorted by Col, no zero Vals }
  TSparseRows = array of TSparseRow;

{ --- sparse rows --- }

function SparseFind(const Row: TSparseRow; Col: Integer): Integer;
procedure SparseRowAddEntry(var Row: TSparseRow; Col: Integer;
  const Val: TRational);
{ Append all N*N rows of kron(P, I_N) - kron(I_N, Q) to Rows, where
  N = Length(P) (16 in every use here). Built generically over the
  nonzeros of P and Q; for the monomial P, Q of the Clifford family
  every row keeps <= 2 entries. }
procedure AppendKroneckerSystem(var Rows: TSparseRows;
  const P, Q: TRatMatrix);
{ Full reduced row echelon form of a sparse row system. }
procedure SparseRREF(const ARows: TSparseRows; NumCols: Integer;
  out RREFRows: TSparseRows; out PivotCols: TIntArray);
{ One kernel basis vector (dense, length NumCols) per free column,
  in standard RREF form. }
function SparseKernelBasis(const ARows: TSparseRows;
  NumCols: Integer): TRatVectors;
function SparseKernelBasisRank(const ARows: TSparseRows;
  NumCols: Integer; out Rank: Integer): TRatVectors;

{ --- exact rational-square test (integer sqrt on TBigInt limbs) --- }

{ True iff R is the square of a rational; then Root is that rational
  (the nonnegative root, normalized). }
function RatIsSquare(const R: TRational; out Root: TRational): Boolean;

{ --- dense small-matrix helpers over Q --- }

function RMatZero(Rows, Cols: Integer): TRatMatrix;
function RMatIdent(N: Integer): TRatMatrix;
function RMatAdd(const A, B: TRatMatrix): TRatMatrix;
function RMatSub(const A, B: TRatMatrix): TRatMatrix;
function RMatNeg(const A: TRatMatrix): TRatMatrix;
function RMatScale(const A: TRatMatrix; const S: TRational): TRatMatrix;
function RMatMul(const A, B: TRatMatrix): TRatMatrix;
function RMatTranspose(const A: TRatMatrix): TRatMatrix;
function RMatTrace(const A: TRatMatrix): TRational;
function RMatEq(const A, B: TRatMatrix): Boolean;
function RMatIsZero(const A: TRatMatrix): Boolean;
function RMatRank(const A: TRatMatrix): Integer;

implementation

{ ================= sparse rows ================= }

function SparseFind(const Row: TSparseRow; Col: Integer): Integer;
var
  Lo, Hi, Mid: Integer;
begin
  Lo := 0;
  Hi := High(Row);
  while Lo <= Hi do
  begin
    Mid := (Lo + Hi) div 2;
    if Row[Mid].Col = Col then
      Exit(Mid);
    if Row[Mid].Col < Col then
      Lo := Mid + 1
    else
      Hi := Mid - 1;
  end;
  Result := -1;
end;

procedure SparseRowAddEntry(var Row: TSparseRow; Col: Integer;
  const Val: TRational);
var
  Lo, Hi, Mid, Idx, I: Integer;
  S: TRational;
begin
  if RatIsZero(Val) then
    Exit;
  Lo := 0;
  Hi := High(Row);
  Idx := Length(Row);
  while Lo <= Hi do
  begin
    Mid := (Lo + Hi) div 2;
    if Row[Mid].Col = Col then
    begin
      S := RatAdd(Row[Mid].Val, Val);
      if RatIsZero(S) then
      begin
        for I := Mid to High(Row) - 1 do
          Row[I] := Row[I + 1];
        SetLength(Row, Length(Row) - 1);
      end
      else
        Row[Mid].Val := S;
      Exit;
    end;
    if Row[Mid].Col < Col then
      Lo := Mid + 1
    else
      Hi := Mid - 1;
  end;
  Idx := Lo;
  SetLength(Row, Length(Row) + 1);
  for I := High(Row) - 1 downto Idx do
    Row[I + 1] := Row[I];
  Row[Idx].Col := Col;
  Row[Idx].Val := Val;
end;

function SparseRowScale(const Row: TSparseRow;
  const S: TRational): TSparseRow;
var
  I: Integer;
begin
  if RatIsZero(S) then
  begin
    SetLength(Result, 0);
    Exit;
  end;
  SetLength(Result, Length(Row));
  for I := 0 to High(Row) do
  begin
    Result[I].Col := Row[I].Col;
    Result[I].Val := RatMul(Row[I].Val, S);
  end;
end;

{ Result = A + Fac * B (Fac <> 0), sparse merge, zeros dropped. }
function SparseRowAddScaled(const A, B: TSparseRow;
  const Fac: TRational): TSparseRow;
var
  I, J, K: Integer;
  V: TRational;
begin
  SetLength(Result, Length(A) + Length(B));
  I := 0;
  J := 0;
  K := 0;
  while (I <= High(A)) or (J <= High(B)) do
  begin
    if (J > High(B)) or ((I <= High(A)) and (A[I].Col < B[J].Col)) then
    begin
      Result[K] := A[I];
      Inc(K);
      Inc(I);
    end
    else if (I > High(A)) or (B[J].Col < A[I].Col) then
    begin
      V := RatMul(Fac, B[J].Val);
      if not RatIsZero(V) then
      begin
        Result[K].Col := B[J].Col;
        Result[K].Val := V;
        Inc(K);
      end;
      Inc(J);
    end
    else
    begin
      V := RatAdd(A[I].Val, RatMul(Fac, B[J].Val));
      if not RatIsZero(V) then
      begin
        Result[K].Col := A[I].Col;
        Result[K].Val := V;
        Inc(K);
      end;
      Inc(I);
      Inc(J);
    end;
  end;
  SetLength(Result, K);
end;

procedure AppendKroneckerSystem(var Rows: TSparseRows;
  const P, Q: TRatMatrix);
var
  N, I, J, K, L, Base: Integer;
  Row: TSparseRow;
begin
  N := Length(P);
  Base := Length(Rows);
  SetLength(Rows, Base + N * N);
  for I := 0 to N - 1 do
    for J := 0 to N - 1 do
    begin
      SetLength(Row, 0);
      { kron(P, I) term: (P x I)[(i,j),(k,l)] = P[i,k] * delta(j,l),
        so column k*N+j gets +P[i,k] for each nonzero P[i,k]. }
      for K := 0 to N - 1 do
        if not RatIsZero(P[I][K]) then
          SparseRowAddEntry(Row, K * N + J, P[I][K]);
      { -kron(I, Q) term: (I x Q)[(i,j),(k,l)] = delta(i,k) * Q[j,l],
        so column i*N+l gets -Q[j,l] for each nonzero Q[j,l]. }
      for L := 0 to N - 1 do
        if not RatIsZero(Q[J][L]) then
          SparseRowAddEntry(Row, I * N + L, RatNeg(Q[J][L]));
      Rows[Base + I * N + J] := Row;
    end;
end;

procedure SparseRREF(const ARows: TSparseRows; NumCols: Integer;
  out RREFRows: TSparseRows; out PivotCols: TIntArray);
var
  W: TSparseRows;
  I, C, PivRow, Cand, BestNNZ, E: Integer;
  Pv, Fac: TRational;
  TmpRow: TSparseRow;
begin
  SetLength(W, Length(ARows));
  for I := 0 to High(ARows) do
    W[I] := ARows[I];
  SetLength(PivotCols, 0);
  PivRow := 0;
  for C := 0 to NumCols - 1 do
  begin
    if PivRow > High(W) then
      Break;
    { pivot choice: fewest-nonzeros row first among the rows at or
      below PivRow that contain column C }
    Cand := -1;
    BestNNZ := MaxInt;
    for I := PivRow to High(W) do
      if SparseFind(W[I], C) >= 0 then
        if Length(W[I]) < BestNNZ then
        begin
          BestNNZ := Length(W[I]);
          Cand := I;
        end;
    if Cand < 0 then
      Continue;
    if Cand <> PivRow then
    begin
      TmpRow := W[PivRow];
      W[PivRow] := W[Cand];
      W[Cand] := TmpRow;
    end;
    E := SparseFind(W[PivRow], C);
    Pv := W[PivRow][E].Val;
    if not RatEq(Pv, Rat(1, 1)) then
      W[PivRow] := SparseRowScale(W[PivRow], RatInv(Pv));
    { eliminate column C from ALL other rows (RREF, not just REF) }
    for I := 0 to High(W) do
      if I <> PivRow then
      begin
        E := SparseFind(W[I], C);
        if E >= 0 then
        begin
          Fac := W[I][E].Val;
          W[I] := SparseRowAddScaled(W[I], W[PivRow], RatNeg(Fac));
        end;
      end;
    SetLength(PivotCols, Length(PivotCols) + 1);
    PivotCols[High(PivotCols)] := C;
    Inc(PivRow);
  end;
  RREFRows := W;
end;

function SparseKernelBasisRank(const ARows: TSparseRows;
  NumCols: Integer; out Rank: Integer): TRatVectors;
var
  RR: TSparseRows;
  Piv: TIntArray;
  IsPiv: array of Boolean;
  PivRowOf: array of Integer;
  I, J, F, E, NB: Integer;
  V: TRatVector;
begin
  SparseRREF(ARows, NumCols, RR, Piv);
  Rank := Length(Piv);
  SetLength(IsPiv, NumCols);
  SetLength(PivRowOf, NumCols);
  for J := 0 to NumCols - 1 do
  begin
    IsPiv[J] := False;
    PivRowOf[J] := -1;
  end;
  for I := 0 to High(Piv) do
  begin
    IsPiv[Piv[I]] := True;
    PivRowOf[Piv[I]] := I; { pivot rows are RR[0..Rank-1], in order }
  end;
  NB := 0;
  for J := 0 to NumCols - 1 do
    if not IsPiv[J] then
      Inc(NB);
  SetLength(Result, NB);
  NB := 0;
  for F := 0 to NumCols - 1 do
    if not IsPiv[F] then
    begin
      SetLength(V, NumCols);
      for J := 0 to NumCols - 1 do
        V[J] := Rat(0, 1);
      V[F] := Rat(1, 1);
      for I := 0 to High(Piv) do
      begin
        E := SparseFind(RR[PivRowOf[Piv[I]]], F);
        if E >= 0 then
          V[Piv[I]] := RatNeg(RR[PivRowOf[Piv[I]]][E].Val);
      end;
      Result[NB] := V;
      Inc(NB);
    end;
end;

function SparseKernelBasis(const ARows: TSparseRows;
  NumCols: Integer): TRatVectors;
var
  Rk: Integer;
begin
  Result := SparseKernelBasisRank(ARows, NumCols, Rk);
end;

{ ================= integer sqrt on limbs ================= }
{ Local limb arithmetic (base 10^9, little-endian, trimmed) used only
  for the exact perfect-square test on TBigInt magnitudes. Mirrors the
  representation inside RTheory.ExactField; the limbs of a TBigInt
  are read, never mutated. }

type
  TLimbs = array of UInt32;

procedure LTrim(var L: TLimbs);
var
  N: Integer;
begin
  N := Length(L);
  while (N > 0) and (L[N - 1] = 0) do
    Dec(N);
  SetLength(L, N);
end;

function LCmp(const A, B: TLimbs): Integer;
var
  I: Integer;
begin
  if Length(A) <> Length(B) then
    Exit(Length(A) - Length(B));
  for I := High(A) downto 0 do
    if A[I] <> B[I] then
    begin
      if A[I] > B[I] then
        Exit(1)
      else
        Exit(-1);
    end;
  Result := 0;
end;

function LAdd(const A, B: TLimbs): TLimbs;
var
  I, N: Integer;
  S, Carry: QWord;
begin
  N := Length(A);
  if Length(B) > N then
    N := Length(B);
  SetLength(Result, N + 1);
  Carry := 0;
  for I := 0 to N do
  begin
    S := Carry;
    if I < Length(A) then
      S := S + A[I];
    if I < Length(B) then
      S := S + B[I];
    Result[I] := S mod 1000000000;
    Carry := S div 1000000000;
  end;
  LTrim(Result);
end;

function LMul(const A, B: TLimbs): TLimbs;
var
  I, J: Integer;
  Carry, T: QWord;
begin
  if (Length(A) = 0) or (Length(B) = 0) then
  begin
    SetLength(Result, 0);
    Exit;
  end;
  SetLength(Result, Length(A) + Length(B));
  { the accumulation below reads Result cells, so zero them
    explicitly: SetLength alone preserves stale contents when the
    caller's result slot is reused }
  for I := 0 to High(Result) do
    Result[I] := 0;
  for I := 0 to High(A) do
  begin
    Carry := 0;
    for J := 0 to High(B) do
    begin
      T := Result[I + J] + QWord(A[I]) * B[J] + Carry;
      Result[I + J] := T mod 1000000000;
      Carry := T div 1000000000;
    end;
    Result[I + Length(B)] := Result[I + Length(B)] + Carry;
  end;
  LTrim(Result);
end;

function LHalve(const A: TLimbs): TLimbs;
var
  I: Integer;
  Rem, Cur: QWord;
begin
  SetLength(Result, Length(A));
  Rem := 0;
  for I := High(A) downto 0 do
  begin
    Cur := Rem * 1000000000 + A[I];
    Result[I] := Cur div 2;
    Rem := Cur mod 2;
  end;
  LTrim(Result);
end;

function LInc(const A: TLimbs): TLimbs;
var
  One: TLimbs;
begin
  SetLength(One, 1);
  One[0] := 1;
  Result := LAdd(A, One);
end;

{ Requires A >= 1. }
function LDec(const A: TLimbs): TLimbs;
var
  I: Integer;
  Borrow: Int64;
  D: Int64;
begin
  SetLength(Result, Length(A));
  Borrow := 1;
  for I := 0 to High(A) do
  begin
    D := Int64(A[I]) - Borrow;
    if D < 0 then
    begin
      D := D + 1000000000;
      Borrow := 1;
    end
    else
      Borrow := 0;
    Result[I] := UInt32(D);
  end;
  LTrim(Result);
end;

{ floor(sqrt(N)) by binary search over the limb representation:
  exact, and fast for the magnitudes that occur here. }
function LFloorSqrt(const N: TLimbs): TLimbs;
var
  Lo, Hi, Mid, Sq: TLimbs;
begin
  SetLength(Lo, 0);
  Hi := Copy(N, 0, Length(N));
  while LCmp(Lo, Hi) < 0 do
  begin
    Mid := LHalve(LInc(LAdd(Lo, Hi))); { upper midpoint }
    Sq := LMul(Mid, Mid);
    if LCmp(Sq, N) <= 0 then
      Lo := Mid
    else
      Hi := LDec(Mid);
  end;
  Result := Lo;
end;

function LimbsIsSquare(const N: TLimbs; out Root: TLimbs): Boolean;
begin
  Root := LFloorSqrt(N);
  Result := LCmp(LMul(Root, Root), N) = 0;
end;

function BigFromLimbs(const L: TLimbs; ASign: Integer): TBigInt;
var
  I: Integer;
begin
  if Length(L) = 0 then
  begin
    Result.Sign := 0;
    SetLength(Result.Mag, 0);
    Exit;
  end;
  Result.Sign := ASign;
  SetLength(Result.Mag, Length(L));
  for I := 0 to High(L) do
    Result.Mag[I] := L[I];
end;

function RatIsSquare(const R: TRational; out Root: TRational): Boolean;
var
  SN, SD: TLimbs;
  NumL, DenL: TLimbs;
begin
  if RatSign(R) < 0 then
    Exit(False);
  NumL := Copy(R.Num.Mag, 0, Length(R.Num.Mag));
  DenL := Copy(R.Den.Mag, 0, Length(R.Den.Mag));
  if not LimbsIsSquare(NumL, SN) then
    Exit(False);
  if not LimbsIsSquare(DenL, SD) then
    Exit(False);
  { R is normalized, so sqrt(num)/sqrt(den) is normalized too. }
  Root.Num := BigFromLimbs(SN, 1);
  if Length(SN) = 0 then
    Root.Num.Sign := 0;
  Root.Den := BigFromLimbs(SD, 1);
  Result := True;
end;

{ ================= dense helpers over Q ================= }

function RMatZero(Rows, Cols: Integer): TRatMatrix;
var
  I, J: Integer;
begin
  SetLength(Result, Rows, Cols);
  for I := 0 to Rows - 1 do
    for J := 0 to Cols - 1 do
      Result[I][J] := Rat(0, 1);
end;

function RMatIdent(N: Integer): TRatMatrix;
var
  I: Integer;
begin
  Result := RMatZero(N, N);
  for I := 0 to N - 1 do
    Result[I][I] := Rat(1, 1);
end;

function RMatAdd(const A, B: TRatMatrix): TRatMatrix;
var
  I, J: Integer;
begin
  Result := RMatZero(Length(A), Length(A[0]));
  for I := 0 to High(A) do
    for J := 0 to High(A[0]) do
      Result[I][J] := RatAdd(A[I][J], B[I][J]);
end;

function RMatSub(const A, B: TRatMatrix): TRatMatrix;
var
  I, J: Integer;
begin
  Result := RMatZero(Length(A), Length(A[0]));
  for I := 0 to High(A) do
    for J := 0 to High(A[0]) do
      Result[I][J] := RatSub(A[I][J], B[I][J]);
end;

function RMatNeg(const A: TRatMatrix): TRatMatrix;
var
  I, J: Integer;
begin
  Result := RMatZero(Length(A), Length(A[0]));
  for I := 0 to High(A) do
    for J := 0 to High(A[0]) do
      Result[I][J] := RatNeg(A[I][J]);
end;

function RMatScale(const A: TRatMatrix; const S: TRational): TRatMatrix;
var
  I, J: Integer;
begin
  Result := RMatZero(Length(A), Length(A[0]));
  for I := 0 to High(A) do
    for J := 0 to High(A[0]) do
      Result[I][J] := RatMul(A[I][J], S);
end;

function RMatMul(const A, B: TRatMatrix): TRatMatrix;
var
  R, K, C, I, J, L: Integer;
begin
  R := Length(A);
  K := Length(A[0]);
  C := Length(B[0]);
  Result := RMatZero(R, C);
  for I := 0 to R - 1 do
    for L := 0 to K - 1 do
      if not RatIsZero(A[I][L]) then
        for J := 0 to C - 1 do
          if not RatIsZero(B[L][J]) then
            Result[I][J] := RatAdd(Result[I][J],
              RatMul(A[I][L], B[L][J]));
end;

function RMatTranspose(const A: TRatMatrix): TRatMatrix;
var
  I, J: Integer;
begin
  Result := RMatZero(Length(A[0]), Length(A));
  for I := 0 to High(A) do
    for J := 0 to High(A[0]) do
      Result[J][I] := A[I][J];
end;

function RMatTrace(const A: TRatMatrix): TRational;
var
  I: Integer;
begin
  Result := Rat(0, 1);
  for I := 0 to Length(A) - 1 do
    Result := RatAdd(Result, A[I][I]);
end;

function RMatEq(const A, B: TRatMatrix): Boolean;
var
  I, J: Integer;
begin
  if (Length(A) <> Length(B)) or (Length(A[0]) <> Length(B[0])) then
    Exit(False);
  for I := 0 to High(A) do
    for J := 0 to High(A[0]) do
      if not RatEq(A[I][J], B[I][J]) then
        Exit(False);
  Result := True;
end;

function RMatIsZero(const A: TRatMatrix): Boolean;
var
  I, J: Integer;
begin
  for I := 0 to High(A) do
    for J := 0 to High(A[0]) do
      if not RatIsZero(A[I][J]) then
        Exit(False);
  Result := True;
end;

function RMatRank(const A: TRatMatrix): Integer;
var
  R, C, I, J, K, Rw, Piv: Integer;
  W: TRatMatrix;
  Fac: TRational;
  TmpRow: TRatVector;
begin
  R := Length(A);
  C := Length(A[0]);
  SetLength(W, R, C);
  for I := 0 to R - 1 do
    for J := 0 to C - 1 do
      W[I][J] := A[I][J];
  Rw := 0;
  for J := 0 to C - 1 do
  begin
    if Rw >= R then
      Break;
    Piv := -1;
    for I := Rw to R - 1 do
      if not RatIsZero(W[I][J]) then
      begin
        Piv := I;
        Break;
      end;
    if Piv < 0 then
      Continue;
    if Piv <> Rw then
    begin
      TmpRow := W[Rw];
      W[Rw] := W[Piv];
      W[Piv] := TmpRow;
    end;
    for I := Rw + 1 to R - 1 do
      if not RatIsZero(W[I][J]) then
      begin
        Fac := RatDiv(W[I][J], W[Rw][J]);
        for K := J to C - 1 do
          W[I][K] := RatSub(W[I][K], RatMul(Fac, W[Rw][K]));
      end;
    Inc(Rw);
  end;
  Result := Rw;
end;

end.
