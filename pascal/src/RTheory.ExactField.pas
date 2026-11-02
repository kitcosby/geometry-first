unit RTheory.ExactField;

{ Exact arithmetic in the biquadratic field Q(sqrt(2), sqrt(3)).

  Element: a + b*sqrt(2) + c*sqrt(3) + d*sqrt(6), with a..d in Q as
  normalized arbitrary-precision rationals. The rational type is
  built on a small base-10^9 bigint (this unit) because the e71/e72
  computations (structure constants, inverse Killing form) outgrow
  Int64 numerators -- that overflow was observed, loud, in the first
  Int64 version, and arbitrary precision is the honest fix. Every
  matrix entry in the e70/e71/e72 exact scripts lies in this field
  (the 15-degree values are (sqrt(6) +/- sqrt(2))/4), so the Pascal
  ports run with zero floating point anywhere.

  The field is ordered; FSign is exact by recursive reduction
  Q(sqrt(2),sqrt(3)) -> Q(sqrt(2)) -> Q. Sylvester inertia over this
  field is therefore exact -- it replaces the QQbar eigenvalue-sign
  counting in e71. }

{$mode objfpc}{$H+}

interface

uses
  SysUtils;

type
  ERationalOverflow = class(Exception);
  EFieldError = class(Exception);

  { Arbitrary-precision integer: sign + little-endian base-10^9
    limbs, no trailing zero limbs; zero is Sign 0, no limbs. }
  TBigInt = record
    Sign: Integer;
    Mag: array of UInt32;
  end;

  TRational = record
    Num, Den: TBigInt; { Den > 0; gcd(|Num|, Den) = 1; zero = 0/1 }
  end;

  { a + b*sqrt(2) + c*sqrt(3) + d*sqrt(6) }
  TField = record
    A, B, C, D: TRational;
  end;

  TIntArray = array of Integer;
  TFieldArray = array of TField;
  TFieldVectors = array of TFieldArray;
  TFMatrix = array of array of TField; { dynamic Rows x Cols }

{ --- rationals --- }

function Rat(ANum, ADen: Int64): TRational;
function RatIsZero(const R: TRational): Boolean;
function RatSign(const R: TRational): Integer;
function RatEq(const R, S: TRational): Boolean;
function RatAdd(const R, S: TRational): TRational;
function RatNeg(const R: TRational): TRational;
function RatSub(const R, S: TRational): TRational;
function RatMul(const R, S: TRational): TRational;
function RatInv(const R: TRational): TRational;
function RatDiv(const R, S: TRational): TRational;
function RatToString(const R: TRational): string;

{ --- the biquadratic field --- }

function Fld(const A, B, C, D: TRational): TField; overload;
function FldI(ANum, ADen: Int64): TField; overload; { rational a }
function FSqrt2: TField;
function FSqrt3: TField;
function FSqrt6: TField;
function FAdd(const X, Y: TField): TField;
function FNeg(const X: TField): TField;
function FSub(const X, Y: TField): TField;
function FMul(const X, Y: TField): TField;
function FInv(const X: TField): TField;
function FDiv(const X, Y: TField): TField;
function FIsZero(const X: TField): Boolean;
function FEq(const X, Y: TField): Boolean;
function FSign(const X: TField): Integer; { -1, 0, +1 -- exact }
function FToString(const X: TField): string;

{ --- exact linear algebra over the field --- }

function FMatZero(Rows, Cols: Integer): TFMatrix;
function FMatIdent(N: Integer): TFMatrix;
function FMatAdd(const A, B: TFMatrix): TFMatrix;
function FMatSub(const A, B: TFMatrix): TFMatrix;
function FMatMul(const A, B: TFMatrix): TFMatrix;
function FMatTrace(const A: TFMatrix): TField;
function FFlatten(const M: TFMatrix): TFieldArray;
{ Rank of a row list by exact Gauss-Jordan; PivotCols lists the
  pivot columns of the row-echelon form (like Sage's .pivots()). }
function FRowsRank(const Rows: TFieldVectors; out PivotCols: TIntArray): Integer;
function FRowsRankSimple(const Rows: TFieldVectors): Integer;
{ Exact inverse of a square matrix; raises EFieldError if singular. }
function FMatInvert(const N: TFMatrix): TFMatrix;
{ Exact nullspace basis of M (Rows x Cols): one vector per free column. }
function FMatNullBasis(const M: TFMatrix): TFieldVectors;
{ Sylvester inertia of a symmetric matrix: exact (Pos, Neg, Zero). }
procedure FSymmetricInertia(const S: TFMatrix; out Pos, Neg, Zero: Integer);
{ Faddeeva-LeVerrier char poly: l^N + C[1] l^(N-1) + ... + C[N]. }
function FMatCharPoly(const A: TFMatrix): TFieldArray;

implementation

const
  BIGBASE = 1000000000; { 10^9 }

type
  TLimbs = array of UInt32;

{ ================= bigint magnitudes ================= }

procedure MagTrim(var L: TLimbs);
var
  N: Integer;
begin
  N := Length(L);
  while (N > 0) and (L[N - 1] = 0) do
    Dec(N);
  SetLength(L, N);
end;

function MagCmp(const A, B: TLimbs): Integer;
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

function MagAdd(const A, B: TLimbs): TLimbs;
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
    Result[I] := S mod BIGBASE;
    Carry := S div BIGBASE;
  end;
  MagTrim(Result);
end;

{ Requires A >= B. }
function MagSub(const A, B: TLimbs): TLimbs;
var
  I: Integer;
  D, Borrow: Int64;
begin
  SetLength(Result, Length(A));
  Borrow := 0;
  for I := 0 to High(A) do
  begin
    D := Int64(A[I]) - Borrow;
    if I < Length(B) then
      D := D - Int64(B[I]);
    if D < 0 then
    begin
      D := D + BIGBASE;
      Borrow := 1;
    end
    else
      Borrow := 0;
    Result[I] := UInt32(D);
  end;
  MagTrim(Result);
end;

function MagMul(const A, B: TLimbs): TLimbs;
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
  for I := 0 to High(A) do
  begin
    Carry := 0;
    for J := 0 to High(B) do
    begin
      T := Result[I + J] + QWord(A[I]) * B[J] + Carry;
      Result[I + J] := T mod BIGBASE;
      Carry := T div BIGBASE;
    end;
    Result[I + Length(B)] := Result[I + Length(B)] + Carry;
  end;
  MagTrim(Result);
end;

function MagMulSmall(const A: TLimbs; M: UInt32): TLimbs;
var
  I: Integer;
  Carry, T: QWord;
begin
  if M = 0 then
  begin
    SetLength(Result, 0);
    Exit;
  end;
  SetLength(Result, Length(A) + 1);
  Carry := 0;
  for I := 0 to High(A) do
  begin
    T := QWord(A[I]) * M + Carry;
    Result[I] := T mod BIGBASE;
    Carry := T div BIGBASE;
  end;
  Result[High(Result)] := Carry;
  MagTrim(Result);
end;

{ Long division (Knuth D specialized): returns Q, R with A = B*Q + R,
  0 <= R < B. B must be nonzero. }
procedure MagDivMod(const A, B: TLimbs; out Q, R: TLimbs);
var
  N, M, J, I: Integer;
  D: UInt32;
  U, V: TLimbs;
  QHat, RHat, Prod, BorrowT: QWord;
  Borrow: Int64;
  T: Int64;
begin
  if MagCmp(A, B) < 0 then
  begin
    SetLength(Q, 0);
    R := Copy(A, 0, Length(A));
    Exit;
  end;
  N := Length(B);
  M := Length(A) - N;
  if N = 1 then
  begin
    { single-limb divisor: simple pass }
    SetLength(Q, Length(A));
    RHat := 0;
    for I := High(A) downto 0 do
    begin
      QHat := RHat * BIGBASE + A[I];
      Q[I] := QHat div B[0];
      RHat := QHat mod B[0];
    end;
    MagTrim(Q);
    if RHat = 0 then
      SetLength(R, 0)
    else
    begin
      SetLength(R, 1);
      R[0] := RHat;
    end;
    Exit;
  end;
  D := BIGBASE div (B[N - 1] + 1);
  U := MagMulSmall(A, D);
  SetLength(U, Length(A) + 1); { pad one zero limb on top }
  V := MagMulSmall(B, D);
  SetLength(Q, M + 1);
  for J := M downto 0 do
  begin
    QHat := (QWord(U[J + N]) * BIGBASE + U[J + N - 1]) div V[N - 1];
    RHat := (QWord(U[J + N]) * BIGBASE + U[J + N - 1]) mod V[N - 1];
    while (QHat >= BIGBASE) or
      (QHat * V[N - 2] > BIGBASE * RHat + U[J + N - 2]) do
    begin
      Dec(QHat);
      RHat := RHat + V[N - 1];
      if RHat >= BIGBASE then
        Break;
    end;
    { U[J..J+N] -= QHat * V }
    Borrow := 0;
    for I := 0 to N - 1 do
    begin
      Prod := QHat * V[I];
      T := Int64(U[J + I]) - Borrow - Int64(Prod mod BIGBASE);
      { Batch-2 disclosed infrastructure repair (2026-10-17): the
        single conditional add-back here was insufficient -- T can
        fall as low as -2*BIGBASE+2 (borrow-in up to BIGBASE-1 plus
        Prod mod BIGBASE up to BIGBASE-1 against a zero limb), so a
        negative limb survived and was wrapped by UInt32(T), silently
        corrupting the division (surfaced by the batch-2 m_refined
        closure data as ERationalOverflow in exact divisions). The
        while loop runs at most twice; Borrow stays <= BIGBASE,
        which the add-back / top-limb step below tolerates. }
      Borrow := Int64(Prod div BIGBASE);
      while T < 0 do
      begin
        T := T + BIGBASE;
        Inc(Borrow);
      end;
      U[J + I] := UInt32(T);
    end;
    T := Int64(U[J + N]) - Borrow;
    if T < 0 then
    begin
      { over-subtracted: add V back, decrement quotient digit }
      U[J + N] := UInt32(T + BIGBASE);
      Dec(QHat);
      BorrowT := 0;
      for I := 0 to N - 1 do
      begin
        BorrowT := QWord(U[J + I]) + V[I] + BorrowT;
        U[J + I] := BorrowT mod BIGBASE;
        BorrowT := BorrowT div BIGBASE;
      end;
      U[J + N] := (QWord(U[J + N]) + BorrowT) mod BIGBASE;
    end
    else
      U[J + N] := UInt32(T);
    Q[J] := QHat;
  end;
  MagTrim(Q);
  { remainder: U / D (exact small division) }
  SetLength(R, Length(U));
  RHat := 0;
  for I := High(U) downto 0 do
  begin
    QHat := RHat * BIGBASE + U[I];
    R[I] := QHat div D;
    RHat := QHat mod D;
  end;
  MagTrim(R);
end;

function MagDivExact(const A, B: TLimbs): TLimbs;
var
  Q, R: TLimbs;
begin
  MagDivMod(A, B, Q, R);
  if Length(R) <> 0 then
    raise ERationalOverflow.Create('bigint division not exact');
  Result := Q;
end;

function MagGCD(A, B: TLimbs): TLimbs;
var
  Q, R: TLimbs;
begin
  while Length(B) > 0 do
  begin
    MagDivMod(A, B, Q, R);
    A := B;
    B := R;
  end;
  Result := A;
end;

{ ================= bigint signed ================= }

function BigZero: TBigInt;
begin
  Result.Sign := 0;
  SetLength(Result.Mag, 0);
end;

function BigFromInt64(V: Int64): TBigInt;
var
  M: QWord;
begin
  if V = 0 then
    Exit(BigZero);
  if V < 0 then
  begin
    Result.Sign := -1;
    M := QWord(-(V + 1)) + 1;
  end
  else
  begin
    Result.Sign := 1;
    M := QWord(V);
  end;
  SetLength(Result.Mag, 0);
  while M > 0 do
  begin
    SetLength(Result.Mag, Length(Result.Mag) + 1);
    Result.Mag[High(Result.Mag)] := M mod BIGBASE;
    M := M div BIGBASE;
  end;
end;

function BigIsZero(const X: TBigInt): Boolean;
begin
  Result := X.Sign = 0;
end;

function BigNeg(const X: TBigInt): TBigInt;
begin
  Result := X;
  Result.Sign := -X.Sign;
end;

function BigAdd(const A, B: TBigInt): TBigInt;
var
  C: Integer;
begin
  if A.Sign = 0 then
    Exit(B);
  if B.Sign = 0 then
    Exit(A);
  if A.Sign = B.Sign then
  begin
    Result.Sign := A.Sign;
    Result.Mag := MagAdd(A.Mag, B.Mag);
    Exit;
  end;
  C := MagCmp(A.Mag, B.Mag);
  if C = 0 then
    Exit(BigZero);
  if C > 0 then
  begin
    Result.Sign := A.Sign;
    Result.Mag := MagSub(A.Mag, B.Mag);
  end
  else
  begin
    Result.Sign := B.Sign;
    Result.Mag := MagSub(B.Mag, A.Mag);
  end;
end;

function BigSub(const A, B: TBigInt): TBigInt;
begin
  Result := BigAdd(A, BigNeg(B));
end;

function BigMul(const A, B: TBigInt): TBigInt;
begin
  if (A.Sign = 0) or (B.Sign = 0) then
    Exit(BigZero);
  Result.Sign := A.Sign * B.Sign;
  Result.Mag := MagMul(A.Mag, B.Mag);
end;

function BigGCD(const A, B: TBigInt): TBigInt;
begin
  if A.Sign = 0 then
  begin
    Result := B;
    Result.Sign := Abs(B.Sign);
    Exit;
  end;
  if B.Sign = 0 then
  begin
    Result := A;
    Result.Sign := Abs(A.Sign);
    Exit;
  end;
  Result.Sign := 1;
  Result.Mag := MagGCD(A.Mag, B.Mag);
end;

{ Exact division: B must divide A. }
function BigDivExact(const A, B: TBigInt): TBigInt;
begin
  if A.Sign = 0 then
    Exit(BigZero);
  Result.Sign := A.Sign * B.Sign;
  Result.Mag := MagDivExact(A.Mag, B.Mag);
end;

function BigToString(const X: TBigInt): string;
var
  I: Integer;
begin
  if X.Sign = 0 then
    Exit('0');
  Result := IntToStr(X.Mag[High(X.Mag)]);
  for I := High(X.Mag) - 1 downto 0 do
    Result := Result + Format('%.9d', [X.Mag[I]]);
  if X.Sign < 0 then
    Result := '-' + Result;
end;

{ ================= rationals ================= }

function RatNorm(const ANum, ADen: TBigInt): TRational;
var
  G, ND: TBigInt;
begin
  if BigIsZero(ADen) then
    raise ERationalOverflow.Create('rational with zero denominator');
  ND := ADen;
  if ND.Sign < 0 then
    ND := BigNeg(ND);
  if BigIsZero(ANum) then
  begin
    Result.Num := BigZero;
    Result.Den := BigFromInt64(1);
    Exit;
  end;
  G.Sign := 1;
  G.Mag := MagGCD(ANum.Mag, ND.Mag);
  Result.Num := BigDivExact(ANum, G);
  if ADen.Sign < 0 then
    Result.Num := BigNeg(Result.Num);
  Result.Den := BigDivExact(ND, G);
end;

function Rat(ANum, ADen: Int64): TRational;
begin
  Result := RatNorm(BigFromInt64(ANum), BigFromInt64(ADen));
end;

function RatIsZero(const R: TRational): Boolean;
begin
  Result := BigIsZero(R.Num);
end;

function RatSign(const R: TRational): Integer;
begin
  Result := R.Num.Sign;
end;

function RatEq(const R, S: TRational): Boolean;
begin
  Result := (R.Num.Sign = S.Num.Sign) and (S.Num.Sign = R.Num.Sign) and
    (MagCmp(R.Num.Mag, S.Num.Mag) = 0) and
    (MagCmp(R.Den.Mag, S.Den.Mag) = 0);
end;

function RatAdd(const R, S: TRational): TRational;
var
  G, D1, D2, N1, N2: TBigInt;
begin
  { d = gcd(den1, den2); num = n1*(d2/d) + n2*(d1/d), reduced }
  G := BigGCD(R.Den, S.Den);
  D1 := BigDivExact(R.Den, G);
  D2 := BigDivExact(S.Den, G);
  N1 := BigMul(R.Num, D2);
  N2 := BigMul(S.Num, D1);
  Result := RatNorm(BigAdd(N1, N2), BigMul(R.Den, D2));
end;

function RatNeg(const R: TRational): TRational;
begin
  Result := R;
  Result.Num := BigNeg(R.Num);
end;

function RatSub(const R, S: TRational): TRational;
begin
  Result := RatAdd(R, RatNeg(S));
end;

function RatMul(const R, S: TRational): TRational;
var
  G1, G2: TBigInt;
  N, D: TBigInt;
begin
  { cross-cancel before multiplying to keep magnitudes small }
  G1 := BigGCD(R.Num, S.Den);
  G1.Sign := 1;
  G2 := BigGCD(S.Num, R.Den);
  G2.Sign := 1;
  N := BigMul(BigDivExact(R.Num, G1), BigDivExact(S.Num, G2));
  D := BigMul(BigDivExact(R.Den, G2), BigDivExact(S.Den, G1));
  Result := RatNorm(N, D);
end;

function RatInv(const R: TRational): TRational;
var
  N, D: TBigInt;
begin
  if BigIsZero(R.Num) then
    raise ERationalOverflow.Create('rational division by zero');
  N := R.Den;
  D := R.Num;
  if D.Sign < 0 then
  begin
    D := BigNeg(D);
    N := BigNeg(N);
  end;
  Result.Num := N;
  Result.Den := D;
end;

function RatDiv(const R, S: TRational): TRational;
begin
  Result := RatMul(R, RatInv(S));
end;

function RatToString(const R: TRational): string;
begin
  if (R.Den.Sign = 1) and (Length(R.Den.Mag) = 1) and (R.Den.Mag[0] = 1) then
    Result := BigToString(R.Num)
  else
    Result := BigToString(R.Num) + '/' + BigToString(R.Den);
end;

{ ================= the field ================= }
{ Basis products: s2^2=2, s3^2=3, s6^2=6, s2*s3=s6,
  s2*s6=2*s3, s3*s6=3*s2. }

function Fld(const A, B, C, D: TRational): TField;
begin
  Result.A := A;
  Result.B := B;
  Result.C := C;
  Result.D := D;
end;

function FldI(ANum, ADen: Int64): TField;
begin
  Result := Fld(Rat(ANum, ADen), Rat(0, 1), Rat(0, 1), Rat(0, 1));
end;

function FSqrt2: TField;
begin
  Result := Fld(Rat(0, 1), Rat(1, 1), Rat(0, 1), Rat(0, 1));
end;

function FSqrt3: TField;
begin
  Result := Fld(Rat(0, 1), Rat(0, 1), Rat(1, 1), Rat(0, 1));
end;

function FSqrt6: TField;
begin
  Result := Fld(Rat(0, 1), Rat(0, 1), Rat(0, 1), Rat(1, 1));
end;

function FAdd(const X, Y: TField): TField;
begin
  Result.A := RatAdd(X.A, Y.A);
  Result.B := RatAdd(X.B, Y.B);
  Result.C := RatAdd(X.C, Y.C);
  Result.D := RatAdd(X.D, Y.D);
end;

function FNeg(const X: TField): TField;
begin
  Result.A := RatNeg(X.A);
  Result.B := RatNeg(X.B);
  Result.C := RatNeg(X.C);
  Result.D := RatNeg(X.D);
end;

function FSub(const X, Y: TField): TField;
begin
  Result := FAdd(X, FNeg(Y));
end;

function FMul(const X, Y: TField): TField;
begin
  Result.A := RatAdd(RatAdd(RatMul(X.A, Y.A),
    RatMul(Rat(2, 1), RatMul(X.B, Y.B))),
    RatAdd(RatMul(Rat(3, 1), RatMul(X.C, Y.C)),
    RatMul(Rat(6, 1), RatMul(X.D, Y.D))));
  Result.B := RatAdd(RatAdd(RatMul(X.A, Y.B), RatMul(X.B, Y.A)),
    RatMul(Rat(3, 1), RatAdd(RatMul(X.C, Y.D), RatMul(X.D, Y.C))));
  Result.C := RatAdd(RatAdd(RatMul(X.A, Y.C), RatMul(X.C, Y.A)),
    RatMul(Rat(2, 1), RatAdd(RatMul(X.B, Y.D), RatMul(X.D, Y.B))));
  Result.D := RatAdd(RatAdd(RatMul(X.A, Y.D), RatMul(X.D, Y.A)),
    RatAdd(RatMul(X.B, Y.C), RatMul(X.C, Y.B)));
end;

type
  TQ2 = record
    U, V: TRational; { u + v*sqrt(2) }
  end;

function Q2Add(const P, Q: TQ2): TQ2;
begin
  Result.U := RatAdd(P.U, Q.U);
  Result.V := RatAdd(P.V, Q.V);
end;

function Q2Neg(const P: TQ2): TQ2;
begin
  Result.U := RatNeg(P.U);
  Result.V := RatNeg(P.V);
end;

function Q2Mul(const P, Q: TQ2): TQ2;
begin
  Result.U := RatAdd(RatMul(P.U, Q.U), RatMul(Rat(2, 1), RatMul(P.V, Q.V)));
  Result.V := RatAdd(RatMul(P.U, Q.V), RatMul(P.V, Q.U));
end;

function Q2Square(const P: TQ2): TQ2;
begin
  Result := Q2Mul(P, P);
end;

function Q2Scale(const P: TQ2; K: Integer): TQ2;
var
  R: TRational;
begin
  R := Rat(K, 1);
  Result.U := RatMul(P.U, R);
  Result.V := RatMul(P.V, R);
end;

function Q2IsZero(const P: TQ2): Boolean;
begin
  Result := RatIsZero(P.U) and RatIsZero(P.V);
end;

{ Exact sign of u + v*sqrt(2): the opposite-sign case is decided by
  comparing u^2 against 2*v^2 (squaring is valid on nonnegatives). }
function Q2Sign(const P: TQ2): Integer;
var
  Su, Sv, Sd: Integer;
  D: TRational;
begin
  if RatIsZero(P.V) then
    Exit(RatSign(P.U));
  if RatIsZero(P.U) then
    Exit(RatSign(P.V)); { sqrt(2) > 0 }
  Su := RatSign(P.U);
  Sv := RatSign(P.V);
  if Su = Sv then
    Exit(Su);
  D := RatSub(RatMul(P.U, P.U), RatMul(Rat(2, 1), RatMul(P.V, P.V)));
  Sd := RatSign(D);
  if Sd = 0 then
    Exit(0); { degenerate; impossible for (u,v)<>0, sqrt(2) irrational }
  { u>0,v<0: x>0 iff u^2 > 2v^2 iff D>0. u<0,v>0: x>0 iff D<0.
    In both cases sign(x) = sign(D)*sign(u). }
  Result := Sd * Su;
end;

function Q2Inv(const P: TQ2): TQ2;
var
  N2: TRational; { u^2 - 2 v^2, nonzero for P <> 0 }
  Inv: TRational;
begin
  N2 := RatSub(RatMul(P.U, P.U), RatMul(Rat(2, 1), RatMul(P.V, P.V)));
  if RatIsZero(N2) then
    raise EFieldError.Create('Q(sqrt(2)) inverse of degenerate element');
  Inv := RatInv(N2);
  Result.U := RatMul(P.U, Inv);
  Result.V := RatNeg(RatMul(P.V, Inv));
end;

{ Inverse via (p + q*sqrt(3))^-1 = (p - q*sqrt(3)) / (p^2 - 3 q^2),
  with p, q in Q(sqrt(2)) handled recursively. }
function FInv(const X: TField): TField;
var
  P, Q, N2, InvN2, RP, RQ: TQ2;
begin
  if FIsZero(X) then
    raise EFieldError.Create('field division by zero');
  P.U := X.A;
  P.V := X.B;
  Q.U := X.C;
  Q.V := X.D;
  N2 := Q2Add(Q2Square(P), Q2Neg(Q2Scale(Q2Square(Q), 3)));
  if Q2IsZero(N2) then
    raise EFieldError.Create('field inverse: zero norm');
  InvN2 := Q2Inv(N2);
  RP := Q2Mul(P, InvN2);
  RQ := Q2Mul(Q2Neg(Q), InvN2);
  Result.A := RP.U;
  Result.B := RP.V;
  Result.C := RQ.U;
  Result.D := RQ.V;
end;

function FDiv(const X, Y: TField): TField;
begin
  Result := FMul(X, FInv(Y));
end;

function FIsZero(const X: TField): Boolean;
begin
  Result := RatIsZero(X.A) and RatIsZero(X.B) and
    RatIsZero(X.C) and RatIsZero(X.D);
end;

function FEq(const X, Y: TField): Boolean;
begin
  Result := RatEq(X.A, Y.A) and RatEq(X.B, Y.B) and
    RatEq(X.C, Y.C) and RatEq(X.D, Y.D);
end;

{ Exact sign of p + q*sqrt(3), p,q in Q(sqrt(2)), by the same
  opposite-sign squaring argument (sqrt(3) > 0, irrational over
  Q(sqrt(2))). }
function FSign(const X: TField): Integer;
var
  P, Q, N2: TQ2;
  Sp, Sq, Sd: Integer;
begin
  P.U := X.A;
  P.V := X.B;
  Q.U := X.C;
  Q.V := X.D;
  if Q2IsZero(Q) then
    Exit(Q2Sign(P));
  if Q2IsZero(P) then
    Exit(Q2Sign(Q)); { sqrt(3) > 0 }
  Sp := Q2Sign(P);
  Sq := Q2Sign(Q);
  if Sp = Sq then
    Exit(Sp);
  N2 := Q2Add(Q2Square(P), Q2Neg(Q2Scale(Q2Square(Q), 3))); { p^2 - 3 q^2 }
  Sd := Q2Sign(N2);
  if Sd = 0 then
    Exit(0); { degenerate; impossible for (p,q) <> 0 }
  { p>0,q<0: x>0 iff p^2 > 3 q^2 iff N2>0. p<0,q>0: x>0 iff N2<0.
    In both cases sign(x) = sign(N2)*sign(p). }
  Result := Sd * Sp;
end;

function FToString(const X: TField): string;
begin
  Result := Format('(%s)+(%s)s2+(%s)s3+(%s)s6',
    [RatToString(X.A), RatToString(X.B), RatToString(X.C), RatToString(X.D)]);
end;

{ ================= exact linear algebra ================= }

function FMatZero(Rows, Cols: Integer): TFMatrix;
var
  I, J: Integer;
begin
  SetLength(Result, Rows, Cols);
  for I := 0 to Rows - 1 do
    for J := 0 to Cols - 1 do
      Result[I][J] := FldI(0, 1);
end;

function FMatIdent(N: Integer): TFMatrix;
var
  I: Integer;
begin
  Result := FMatZero(N, N);
  for I := 0 to N - 1 do
    Result[I][I] := FldI(1, 1);
end;

function FMatAdd(const A, B: TFMatrix): TFMatrix;
var
  R, C, I, J: Integer;
begin
  R := Length(A);
  C := Length(A[0]);
  Result := FMatZero(R, C);
  for I := 0 to R - 1 do
    for J := 0 to C - 1 do
      Result[I][J] := FAdd(A[I][J], B[I][J]);
end;

function FMatSub(const A, B: TFMatrix): TFMatrix;
var
  R, C, I, J: Integer;
begin
  R := Length(A);
  C := Length(A[0]);
  Result := FMatZero(R, C);
  for I := 0 to R - 1 do
    for J := 0 to C - 1 do
      Result[I][J] := FSub(A[I][J], B[I][J]);
end;

function FMatMul(const A, B: TFMatrix): TFMatrix;
var
  R, K, C, I, J, L: Integer;
  S: TField;
begin
  R := Length(A);
  K := Length(A[0]);
  C := Length(B[0]);
  Result := FMatZero(R, C);
  for I := 0 to R - 1 do
    for J := 0 to C - 1 do
    begin
      S := FldI(0, 1);
      for L := 0 to K - 1 do
        S := FAdd(S, FMul(A[I][L], B[L][J]));
      Result[I][J] := S;
    end;
end;

function FMatTrace(const A: TFMatrix): TField;
var
  I: Integer;
begin
  Result := FldI(0, 1);
  for I := 0 to Length(A) - 1 do
    Result := FAdd(Result, A[I][I]);
end;

function FFlatten(const M: TFMatrix): TFieldArray;
var
  R, C, I, J: Integer;
begin
  R := Length(M);
  C := Length(M[0]);
  SetLength(Result, R * C);
  for I := 0 to R - 1 do
    for J := 0 to C - 1 do
      Result[I * C + J] := M[I][J];
end;

function FRowsRank(const Rows: TFieldVectors; out PivotCols: TIntArray): Integer;
var
  R, C, I, J, K, PivRow, P: Integer;
  A: TFieldVectors;
  Inv, Factor, TmpF: TField;
begin
  R := Length(Rows);
  if R = 0 then
  begin
    SetLength(PivotCols, 0);
    Exit(0);
  end;
  C := Length(Rows[0]);
  SetLength(A, R);
  for I := 0 to R - 1 do
  begin
    SetLength(A[I], C);
    for J := 0 to C - 1 do
      A[I][J] := Rows[I][J];
  end;
  SetLength(PivotCols, 0);
  PivRow := 0;
  for J := 0 to C - 1 do
  begin
    P := -1;
    for I := PivRow to R - 1 do
      if not FIsZero(A[I][J]) then
      begin
        P := I;
        Break;
      end;
    if P < 0 then
      Continue;
    if P <> PivRow then
      for K := 0 to C - 1 do
      begin
        TmpF := A[PivRow][K];
        A[PivRow][K] := A[P][K];
        A[P][K] := TmpF;
      end;
    Inv := FInv(A[PivRow][J]);
    for K := 0 to C - 1 do
      A[PivRow][K] := FMul(A[PivRow][K], Inv);
    for I := 0 to R - 1 do
      if (I <> PivRow) and (not FIsZero(A[I][J])) then
      begin
        Factor := A[I][J];
        for K := 0 to C - 1 do
          A[I][K] := FSub(A[I][K], FMul(Factor, A[PivRow][K]));
      end;
    SetLength(PivotCols, Length(PivotCols) + 1);
    PivotCols[High(PivotCols)] := J;
    Inc(PivRow);
    if PivRow >= R then
      Break;
  end;
  Result := Length(PivotCols);
end;

function FRowsRankSimple(const Rows: TFieldVectors): Integer;
var
  Piv: TIntArray;
begin
  Result := FRowsRank(Rows, Piv);
end;

function FMatInvert(const N: TFMatrix): TFMatrix;
var
  Dim, I, J, K: Integer;
  A: TFMatrix;
  Inv, Factor, TmpF: TField;
  P: Integer;
begin
  Dim := Length(N);
  SetLength(A, Dim, 2 * Dim);
  for I := 0 to Dim - 1 do
    for J := 0 to Dim - 1 do
    begin
      A[I][J] := N[I][J];
      A[I][Dim + J] := FldI(0, 1);
    end;
  for I := 0 to Dim - 1 do
    A[I][Dim + I] := FldI(1, 1);
  for J := 0 to Dim - 1 do
  begin
    P := -1;
    for I := J to Dim - 1 do
      if not FIsZero(A[I][J]) then
      begin
        P := I;
        Break;
      end;
    if P < 0 then
      raise EFieldError.Create('FMatInvert: singular matrix');
    if P <> J then
      for K := 0 to 2 * Dim - 1 do
      begin
        TmpF := A[J][K];
        A[J][K] := A[P][K];
        A[P][K] := TmpF;
      end;
    Inv := FInv(A[J][J]);
    for K := 0 to 2 * Dim - 1 do
      A[J][K] := FMul(A[J][K], Inv);
    for I := 0 to Dim - 1 do
      if (I <> J) and (not FIsZero(A[I][J])) then
      begin
        Factor := A[I][J];
        for K := 0 to 2 * Dim - 1 do
          A[I][K] := FSub(A[I][K], FMul(Factor, A[J][K]));
      end;
  end;
  SetLength(Result, Dim, Dim);
  for I := 0 to Dim - 1 do
    for J := 0 to Dim - 1 do
      Result[I][J] := A[I][Dim + J];
end;

function FMatNullBasis(const M: TFMatrix): TFieldVectors;
var
  R, C, I, J, K: Integer;
  A: TFMatrix;
  Inv, Factor, TmpF: TField;
  IsPivot: array of Boolean;
  PivRowOf: array of Integer;
  P, PivRow, NF: Integer;
  V: TFieldArray;
begin
  R := Length(M);
  C := Length(M[0]);
  A := FMatZero(R, C);
  for I := 0 to R - 1 do
    for J := 0 to C - 1 do
      A[I][J] := M[I][J];
  SetLength(IsPivot, C);
  SetLength(PivRowOf, C);
  for J := 0 to C - 1 do
  begin
    IsPivot[J] := False;
    PivRowOf[J] := -1;
  end;
  PivRow := 0;
  for J := 0 to C - 1 do
  begin
    P := -1;
    for I := PivRow to R - 1 do
      if not FIsZero(A[I][J]) then
      begin
        P := I;
        Break;
      end;
    if P < 0 then
      Continue;
    if P <> PivRow then
      for K := 0 to C - 1 do
      begin
        TmpF := A[PivRow][K];
        A[PivRow][K] := A[P][K];
        A[P][K] := TmpF;
      end;
    Inv := FInv(A[PivRow][J]);
    for K := 0 to C - 1 do
      A[PivRow][K] := FMul(A[PivRow][K], Inv);
    for I := 0 to R - 1 do
      if (I <> PivRow) and (not FIsZero(A[I][J])) then
      begin
        Factor := A[I][J];
        for K := 0 to C - 1 do
          A[I][K] := FSub(A[I][K], FMul(Factor, A[PivRow][K]));
      end;
    IsPivot[J] := True;
    PivRowOf[J] := PivRow;
    Inc(PivRow);
    if PivRow >= R then
      Break;
  end;
  NF := 0;
  for J := 0 to C - 1 do
    if not IsPivot[J] then
      Inc(NF);
  SetLength(Result, NF);
  NF := 0;
  for J := 0 to C - 1 do
    if not IsPivot[J] then
    begin
      SetLength(V, C);
      for K := 0 to C - 1 do
        V[K] := FldI(0, 1);
      V[J] := FldI(1, 1);
      for K := 0 to C - 1 do
        if IsPivot[K] then
          V[K] := FNeg(A[PivRowOf[K]][J]);
      Result[NF] := V;
      Inc(NF);
    end;
end;

{ Exact Sylvester inertia: congruence-diagonalize S and count the
  signs of the diagonal pivots with the exact FSign. Zero-pivot
  handling: swap with a later nonzero diagonal; failing that, the
  standard add-row/add-column congruence (R_i += R_j, C_i += C_j),
  which makes the diagonal entry 2*S[i][j] <> 0 in that branch
  (all later diagonals are zero there). }

procedure FSymmetricInertia(const S: TFMatrix; out Pos, Neg, Zero: Integer);
var
  N, I, J, K: Integer;
  A: TFMatrix;
  Factor, TmpF: TField;
  Sg: Integer;
  Fixed: Boolean;
begin
  N := Length(S);
  A := FMatZero(N, N);
  for I := 0 to N - 1 do
    for J := 0 to N - 1 do
      A[I][J] := S[I][J];
  Pos := 0;
  Neg := 0;
  Zero := 0;
  for I := 0 to N - 1 do
  begin
    if FIsZero(A[I][I]) then
    begin
      Fixed := False;
      for J := I + 1 to N - 1 do
        if not FIsZero(A[J][J]) then
        begin
          for K := 0 to N - 1 do
          begin
            TmpF := A[I][K];
            A[I][K] := A[J][K];
            A[J][K] := TmpF;
          end;
          for K := 0 to N - 1 do
          begin
            TmpF := A[K][I];
            A[K][I] := A[K][J];
            A[K][J] := TmpF;
          end;
          Fixed := True;
          Break;
        end;
      if not Fixed then
        for J := I + 1 to N - 1 do
          if not FIsZero(A[I][J]) then
          begin
            for K := 0 to N - 1 do
              A[I][K] := FAdd(A[I][K], A[J][K]);
            for K := 0 to N - 1 do
              A[K][I] := FAdd(A[K][I], A[K][J]);
            Fixed := True;
            Break;
          end;
      if not Fixed then
      begin
        Zero := Zero + (N - I);
        Exit;
      end;
    end;
    Sg := FSign(A[I][I]);
    if Sg > 0 then
      Inc(Pos)
    else if Sg < 0 then
      Inc(Neg)
    else
    begin
      Inc(Zero);
      Continue;
    end;
    for J := I + 1 to N - 1 do
      if not FIsZero(A[J][I]) then
      begin
        Factor := FDiv(A[J][I], A[I][I]);
        for K := I to N - 1 do
          A[J][K] := FSub(A[J][K], FMul(Factor, A[I][K]));
        for K := I to N - 1 do
          A[K][J] := FSub(A[K][J], FMul(Factor, A[K][I]));
      end;
  end;
end;

{ Faddeeva-LeVerrier: M_0 = 0; M_k = A*M_{k-1} + c_{k-1}*I (c_0 = 1);
  c_k = -tr(A*M_k)/k. charpoly = l^N + c_1 l^{N-1} + ... + c_N. }

function FMatCharPoly(const A: TFMatrix): TFieldArray;
var
  N, K, I: Integer;
  B, T, AM: TFMatrix;
  CPrev, Ck: TField;
begin
  N := Length(A);
  SetLength(Result, N + 1);
  B := FMatZero(N, N);
  CPrev := FldI(1, 1);
  for K := 1 to N do
  begin
    T := FMatMul(A, B);
    for I := 0 to N - 1 do
      T[I][I] := FAdd(T[I][I], CPrev);
    AM := FMatMul(A, T);
    Ck := FNeg(FDiv(FMatTrace(AM), FldI(K, 1)));
    Result[K] := Ck;
    B := T;
    CPrev := Ck;
  end;
end;

end.
