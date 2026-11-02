unit RTheory.ExactClifford;

{ EXACT Pascal ports of the Q-linear Clifford family:
  q4_exact.py (T1), q5_exact.py (T2), q6a_exact.py (T3, the
  rational-search attempt), d2_structure.py (D2), d3b_forms.py (D3b).

  Every matrix entry is an exact rational (TRational, arbitrary
  precision, from RTheory.ExactField); every system is built and
  solved exactly with the sparse RREF engine in RTheory.RatSparse.
  The constructions replicate the Sage scripts literally, including
  their Kronecker-system expressions in Sage's tensor_product
  convention with the unknown vector read row-major -- see the
  TQ5ExactProof note on what that convention does to the basis
  matrices. No floating point anywhere; failures are reported as
  failures, with the discrepancy stated, never massaged. }

{$mode objfpc}{$H+}

interface

uses
  SysUtils, RTheory.Types, RTheory.Theorems, RTheory.ExactField,
  RTheory.RatSparse;

type
  TQ4ExactProof = class(TProofObject)
  public
    constructor Create; reintroduce;
    function VerifyNumeric(out ADetail: string): Boolean; override;
  end;

  TQ5ExactProof = class(TProofObject)
  public
    constructor Create; reintroduce;
    function VerifyNumeric(out ADetail: string): Boolean; override;
  end;

  TQ6aExactProof = class(TProofObject)
  public
    constructor Create; reintroduce;
    function VerifyNumeric(out ADetail: string): Boolean; override;
  end;

  TD2ExactProof = class(TProofObject)
  public
    constructor Create; reintroduce;
    function VerifyNumeric(out ADetail: string): Boolean; override;
  end;

  TD3bExactProof = class(TProofObject)
  public
    constructor Create; reintroduce;
    function VerifyNumeric(out ADetail: string): Boolean; override;
  end;

implementation

{ --- construction helpers, literal from the Sage scripts --- }

function Mat2(A, B, C, D: Int64): TRatMatrix;
begin
  Result := RMatZero(2, 2);
  Result[0][0] := Rat(A, 1);
  Result[0][1] := Rat(B, 1);
  Result[1][0] := Rat(C, 1);
  Result[1][1] := Rat(D, 1);
end;

{ Dense Kronecker product: (A x B)[(i,j),(k,l)] = A[i,k]*B[j,l],
  row index (i,j) -> i*rows(B)+j, col index (k,l) -> k*cols(B)+l. }
function KronDense(const A, B: TRatMatrix): TRatMatrix;
var
  RA, CA, RB, CB, I, J, K, L: Integer;
begin
  RA := Length(A);
  CA := Length(A[0]);
  RB := Length(B);
  CB := Length(B[0]);
  Result := RMatZero(RA * RB, CA * CB);
  for I := 0 to RA - 1 do
    for J := 0 to RB - 1 do
      for K := 0 to CA - 1 do
        if not RatIsZero(A[I][K]) then
          for L := 0 to CB - 1 do
            if not RatIsZero(B[J][L]) then
              Result[I * RB + J][K * CB + L] :=
                RatMul(A[I][K], B[J][L]);
end;

{ Builds exactly as q4_exact.py: 2x2 X, Y, Z (W is defined in the
  script but unused in the construction);
  e = [X x I2, Y x I2, Z x X, Z x Y] (4x4, signature +,-,+,-);
  Omega = e0*e1*e2*e3;
  E_mu = e_mu x I4 (mu=0..3), E[4+nu] = Omega x e_nu (16x16). }
procedure BuildClifford(out E: TRatMatrices; out SmallE: TRatMatrices;
  out Omega: TRatMatrix);
var
  I2, Mx, My, Mz: TRatMatrix;
  I4: TRatMatrix;
  Mu: Integer;
begin
  I2 := Mat2(1, 0, 0, 1);
  Mx := Mat2(0, 1, 1, 0);   { X^2 = +I }
  Mz := Mat2(1, 0, 0, -1);  { Z^2 = +I }
  My := Mat2(0, -1, 1, 0);  { Y^2 = -I }
  SetLength(SmallE, 4);
  SmallE[0] := KronDense(Mx, I2);
  SmallE[1] := KronDense(My, I2);
  SmallE[2] := KronDense(Mz, Mx);
  SmallE[3] := KronDense(Mz, My);
  Omega := RMatMul(RMatMul(SmallE[0], SmallE[1]),
    RMatMul(SmallE[2], SmallE[3]));
  I4 := RMatIdent(4);
  SetLength(E, 8);
  for Mu := 0 to 3 do
    E[Mu] := KronDense(SmallE[Mu], I4);
  for Mu := 0 to 3 do
    E[4 + Mu] := KronDense(Omega, SmallE[Mu]);
end;

procedure SplitSectors(const E: TRatMatrices; out R, J: TRatMatrices);
var
  Mu: Integer;
begin
  SetLength(R, 4);
  SetLength(J, 4);
  for Mu := 0 to 3 do
  begin
    R[Mu] := E[Mu];
    J[Mu] := E[4 + Mu];
  end;
end;

{ +1 if M = +I, -1 if M = -I, 0 otherwise (square matrices). }
function IdentSign(const M: TRatMatrix): Integer;
var
  I, J: Integer;
  Plus, Minus: Boolean;
begin
  Plus := True;
  Minus := True;
  for I := 0 to High(M) do
    for J := 0 to High(M[0]) do
      if I = J then
      begin
        if not RatEq(M[I][J], Rat(1, 1)) then
          Plus := False;
        if not RatEq(M[I][J], Rat(-1, 1)) then
          Minus := False;
      end
      else if not RatIsZero(M[I][J]) then
      begin
        Plus := False;
        Minus := False;
      end;
  if Plus then
    Exit(1);
  if Minus then
    Exit(-1);
  Result := 0;
end;

function AntiCommIsZero(const A, B: TRatMatrix): Boolean;
begin
  Result := RMatIsZero(RMatAdd(RMatMul(A, B), RMatMul(B, A)));
end;

function VecToMat16(const V: TRatVector): TRatMatrix;
var
  R, C: Integer;
begin
  Result := RMatZero(16, 16);
  for R := 0 to 15 do
    for C := 0 to 15 do
      Result[R][C] := V[R * 16 + C];
end;

function IsIntertwinerMat(const X: TRatMatrix;
  const R, J: TRatMatrices): Boolean;
var
  Mu: Integer;
begin
  for Mu := 0 to 3 do
    if not RMatEq(RMatMul(X, R[Mu]), RMatMul(J[Mu], X)) then
      Exit(False);
  Result := True;
end;

{ --- the q5 intertwiner basis, shared by T2 / T3 / D2 ---
  Built once per process (unit-level lazy singleton): the system is
  the Sage expression of q5_exact.py, literally:
  for each mu, the 256 rows of kron(R_mu^T, I16) - kron(I16, J_mu),
  1024 rows x 256 cols, solved exactly; basis vectors are reshaped
  row-major into 16x16 matrices, exactly as the scripts' vec_to_mat.

  CONVENTION NOTE (the math, stated plainly): with vec(X) read
  row-major, that system's equations are R_mu^T * B = B * J_mu^T --
  the transpose of the intertwiner equation B * R_mu = J_mu * B.
  The Sage scripts obtain the same kernel (their dimension claims
  are unaffected, and the Gram matrix of D2 comes out the same
  because the trace form is transpose-invariant), but the basis
  matrices themselves are the transposes of intertwiners, not
  intertwiners. TQ5ExactProof verifies both equations and reports
  the outcome honestly. }

var
  GCachedBasis: TRatMatrices;
  GCachedBasisDim: Integer;
  GCacheReady: Boolean = False;

function GetIntertwinerBasis(out B: TRatMatrices): Integer;
var
  E, SmallE, R, J: TRatMatrices;
  Omega: TRatMatrix;
  Rows: TSparseRows;
  Vecs: TRatVectors;
  Mu, I: Integer;
begin
  if not GCacheReady then
  begin
    BuildClifford(E, SmallE, Omega);
    SplitSectors(E, R, J);
    SetLength(Rows, 0);
    for Mu := 0 to 3 do
      AppendKroneckerSystem(Rows, RMatTranspose(R[Mu]), J[Mu]);
    Vecs := SparseKernelBasis(Rows, 256);
    GCachedBasisDim := Length(Vecs);
    SetLength(GCachedBasis, GCachedBasisDim);
    for I := 0 to GCachedBasisDim - 1 do
      GCachedBasis[I] := VecToMat16(Vecs[I]);
    GCacheReady := True;
  end;
  B := GCachedBasis;
  Result := GCachedBasisDim;
end;

{ M = C*I ? Returns True and sets C. }
function IsScalarMat(const M: TRatMatrix; out C: TRational): Boolean;
var
  I, J: Integer;
begin
  C := M[0][0];
  for I := 0 to High(M) do
    for J := 0 to High(M[0]) do
      if I = J then
      begin
        if not RatEq(M[I][J], C) then
          Exit(False);
      end
      else if not RatIsZero(M[I][J]) then
        Exit(False);
  Result := True;
end;

{ --- TQ4ExactProof (T1-EXACT, proofs/q4_exact.py) --- }

constructor TQ4ExactProof.Create;
begin
  inherited Create('T1-EXACT',
    'T1 exact port: Cl(2,2) and Cl(4,4) relations hold exactly over Q; Omega^2 = +I and {Omega, e_mu} = 0 exactly; commutant of the Cl(4,4) 16-dim representation has exact dimension 1 (scalars only)',
    psProved, 'proofs/q4_exact.py (+ Pascal exact port, RTheory.ExactClifford)');
end;

function TQ4ExactProof.VerifyNumeric(out ADetail: string): Boolean;
var
  E, SmallE: TRatMatrices;
  Omega, Sq: TRatMatrix;
  Rows: TSparseRows;
  Vecs: TRatVectors;
  I, J, Dim: Integer;
  SmallSqOk, SmallAntiOk, OmegaSqOk, OmegaAntiOk: Boolean;
  ESigOk, EAntiOk: Boolean;
  ExpSmall: array[0..3] of Integer;
  ExpE: array[0..7] of Integer;
  SigStr: string;
begin
  try
    ExpSmall[0] := 1;
    ExpSmall[1] := -1;
    ExpSmall[2] := 1;
    ExpSmall[3] := -1;
    for I := 0 to 7 do
      ExpE[I] := ExpSmall[I mod 4];
    BuildClifford(E, SmallE, Omega);
    SmallSqOk := True;
    for I := 0 to 3 do
    begin
      Sq := RMatMul(SmallE[I], SmallE[I]);
      if IdentSign(Sq) <> ExpSmall[I] then
        SmallSqOk := False;
    end;
    SmallAntiOk := True;
    for I := 0 to 3 do
      for J := I + 1 to 3 do
        if not AntiCommIsZero(SmallE[I], SmallE[J]) then
          SmallAntiOk := False;
    OmegaSqOk := IdentSign(RMatMul(Omega, Omega)) = 1;
    OmegaAntiOk := True;
    for I := 0 to 3 do
      if not AntiCommIsZero(Omega, SmallE[I]) then
        OmegaAntiOk := False;
    ESigOk := True;
    SigStr := '';
    for I := 0 to 7 do
    begin
      Sq := RMatMul(E[I], E[I]);
      if IdentSign(Sq) <> ExpE[I] then
        ESigOk := False;
      if IdentSign(Sq) = 1 then
        SigStr := SigStr + '+1,'
      else if IdentSign(Sq) = -1 then
        SigStr := SigStr + '-1,'
      else
        SigStr := SigStr + '??,';
    end;
    if SigStr <> '' then
      Delete(SigStr, Length(SigStr), 1);
    EAntiOk := True;
    for I := 0 to 7 do
      for J := I + 1 to 7 do
        if not AntiCommIsZero(E[I], E[J]) then
          EAntiOk := False;
    { Commutant: Sage system kron(I16, E_i) - kron(E_i^T, I16),
      appended blockwise via the RatSparse builder as
      kron(-E_i^T, I16) - kron(I16, -E_i), which is the same matrix:
      -E_i^T x I + I x E_i. 2048 rows x 256 cols. }
    SetLength(Rows, 0);
    for I := 0 to 7 do
      AppendKroneckerSystem(Rows,
        RMatNeg(RMatTranspose(E[I])), RMatNeg(E[I]));
    Vecs := SparseKernelBasis(Rows, 256);
    Dim := Length(Vecs);
    ADetail := Format(
      'EXACT (Q, zero-residual): Cl(2,2) squares (+,-,+,-) %s, anticommutators %s; Omega^2=+I %s, {Omega,e_mu}=0 %s; Cl(4,4) signature [%s] %s, all 28 anticommutators 0 %s; commutant dim %d (expect 1)',
      [BoolToStr(SmallSqOk, 'exact', 'FAILED'),
       BoolToStr(SmallAntiOk, 'exact', 'FAILED'),
       BoolToStr(OmegaSqOk, 'yes', 'NO'),
       BoolToStr(OmegaAntiOk, 'yes', 'NO'),
       SigStr, BoolToStr(ESigOk, 'as expected', 'FAILED'),
       BoolToStr(EAntiOk, 'yes', 'NO'), Dim]);
    Result := SmallSqOk and SmallAntiOk and OmegaSqOk and OmegaAntiOk
      and ESigOk and EAntiOk and (Dim = 1);
  except
    on Ex: Exception do
    begin
      ADetail := 'EXACT PORT FAILED: ' + Ex.ClassName + ': ' + Ex.Message;
      Result := False;
    end;
  end;
end;

{ --- TQ5ExactProof (T2-EXACT, proofs/q5_exact.py) --- }

constructor TQ5ExactProof.Create;
begin
  inherited Create('T2-EXACT',
    'T2 exact port: sector intertwiner space has exact dimension 16 over Q, and the RREF kernel basis (row-major reshape, as in the Sage script) is self-validated against B*R_mu = J_mu*B',
    psProved, 'proofs/q5_exact.py (+ Pascal exact port, RTheory.ExactClifford)');
end;

function TQ5ExactProof.VerifyNumeric(out ADetail: string): Boolean;
var
  E, SmallE, R, J, B: TRatMatrices;
  Omega: TRatMatrix;
  Dim, I, Mu, InterCount, SysCount: Integer;
  SysOk: Boolean;
begin
  try
    BuildClifford(E, SmallE, Omega);
    SplitSectors(E, R, J);
    Dim := GetIntertwinerBasis(B);
    InterCount := 0;
    SysCount := 0;
    for I := 0 to Dim - 1 do
    begin
      if IsIntertwinerMat(B[I], R, J) then
        Inc(InterCount);
      SysOk := True;
      for Mu := 0 to 3 do
        if not RMatEq(RMatMul(RMatTranspose(R[Mu]), B[I]),
          RMatMul(B[I], RMatTranspose(J[Mu]))) then
          SysOk := False;
      if SysOk then
        Inc(SysCount);
    end;
    ADetail := Format(
      'EXACT (Q, zero-residual): intertwiner-system kernel dim %d (expect 16); basis matrices satisfying B*R_mu = J_mu*B for all mu: %d/%d; satisfying the system equation R_mu^T*B = B*J_mu^T: %d/%d. The basis matrices are transposes of intertwiners (row-major vec vs the column-stacking Kronecker identity in the source construction); the dimension claim itself is exact.',
      [Dim, InterCount, Dim, SysCount, Dim]);
    Result := (Dim = 16) and (InterCount = Dim);
  except
    on Ex: Exception do
    begin
      ADetail := 'EXACT PORT FAILED: ' + Ex.ClassName + ': ' + Ex.Message;
      Result := False;
    end;
  end;
end;

{ --- TQ6aExactProof (T3-EXACT, proofs/q6a_exact.py) --- }

constructor TQ6aExactProof.Create;
begin
  inherited Create('T3-EXACT',
    'T3 exact port (q6a attempt 1): search the exact intertwiner basis (singles, then pairs B_i + B_j) for X with X^2 = +I over Q, verified exactly; the Groebner fallback of the script is class (c) and out of scope',
    psProved, 'proofs/q6a_exact.py (+ Pascal exact port, RTheory.ExactClifford)');
end;

function TQ6aExactProof.VerifyNumeric(out ADetail: string): Boolean;
var
  B: TRatMatrices;
  Dim, I, J: Integer;
  Found: Boolean;
  FoundDesc, FoundC, FoundRoot: string;

  function TryCandidate(const X0: TRatMatrix; const Desc: string): Boolean;
  var
    S, X: TRatMatrix;
    C, Root: TRational;
  begin
    Result := False;
    S := RMatMul(X0, X0);
    if not IsScalarMat(S, C) then
      Exit;
    if RatIsZero(C) then
      Exit;
    if RatSign(C) <= 0 then
      Exit; { negative c would give X^2 = -I, not the target }
    if not RatIsSquare(C, Root) then
      Exit;
    X := RMatScale(X0, RatInv(Root));
    if RMatEq(RMatMul(X, X), RMatIdent(16)) then
    begin
      FoundDesc := Desc;
      FoundC := RatToString(C);
      FoundRoot := RatToString(Root);
      Result := True;
    end;
  end;

begin
  try
    Dim := GetIntertwinerBasis(B);
    Found := False;
    FoundDesc := '';
    FoundC := '';
    FoundRoot := '';
    if Dim = 16 then
    begin
      for I := 0 to 15 do
        if TryCandidate(B[I], Format('B[%d]', [I])) then
        begin
          Found := True;
          Break;
        end;
      if not Found then
        for I := 0 to 14 do
        begin
          for J := I + 1 to 15 do
            if TryCandidate(RMatAdd(B[I], B[J]),
              Format('B[%d]+B[%d]', [I, J])) then
            begin
              Found := True;
              Break;
            end;
          if Found then
            Break;
        end;
    end;
    if Found then
      ADetail := Format(
        'EXACT (Q, zero-residual): explicit involution found: X = (%s)/%s with (%s)^2 square = %s*I, sqrt(c) = %s in Q; X^2 = I16 verified exactly',
        [FoundDesc, FoundRoot, FoundDesc, FoundC, FoundRoot])
    else
      ADetail := Format(
        'EXACT (Q, zero-residual): intertwiner basis dim %d; singles tried 16, pairs B_i+B_j tried 120 (i<j, lexicographic); no candidate has square a positive rational-square multiple of I16 -- search exhausted, FAIL honestly (the Groebner fallback is class (c), out of scope, not attempted)',
        [Dim]);
    Result := Found;
  except
    on Ex: Exception do
    begin
      ADetail := 'EXACT PORT FAILED: ' + Ex.ClassName + ': ' + Ex.Message;
      Result := False;
    end;
  end;
end;

{ --- TD2ExactProof (D2-EXACT, proofs/d2_structure.py) --- }

constructor TD2ExactProof.Create;
begin
  inherited Create('D2-EXACT',
    'D2 exact port: intertwiner space V (dim 16 over Q): I16 is not in V; V is not closed under 16x16 multiplication (first failing pair reported); Gram G[i,j] = tr(B_i^T B_j) is claimed in the registry to be exactly 4*I16',
    psProved, 'proofs/d2_structure.py (+ Pascal exact port, RTheory.ExactClifford)');
end;

function TD2ExactProof.VerifyNumeric(out ADetail: string): Boolean;
var
  E, SmallE, R, Js, B: TRatMatrices;
  Omega: TRatMatrix;
  G: TRatMatrix;
  Dim, I, J, Fi, Fj, OffNZ, GRank: Integer;
  IInV, Closed, GramOk: Boolean;
  TrList: string;
  Tr: TRational;
begin
  try
    BuildClifford(E, SmallE, Omega);
    SplitSectors(E, R, Js);
    Dim := GetIntertwinerBasis(B);
    IInV := IsIntertwinerMat(RMatIdent(16), R, Js);
    Closed := True;
    Fi := -1;
    Fj := -1;
    if Dim = 16 then
      for I := 0 to 15 do
      begin
        for J := 0 to 15 do
          if not IsIntertwinerMat(RMatMul(B[I], B[J]), R, Js) then
          begin
            Closed := False;
            Fi := I;
            Fj := J;
            Break;
          end;
        if not Closed then
          Break;
      end;
    { Gram matrix G[i,j] = tr(B_i^T * B_j) over Q }
    G := RMatZero(16, 16);
    if Dim = 16 then
      for I := 0 to 15 do
        for J := 0 to 15 do
          G[I][J] := RMatTrace(RMatMul(RMatTranspose(B[I]), B[J]));
    GramOk := Dim = 16;
    OffNZ := 0;
    for I := 0 to 15 do
      for J := 0 to 15 do
      begin
        if I = J then
        begin
          if not RatEq(G[I][J], Rat(4, 1)) then
            GramOk := False;
        end
        else
        begin
          if not RatIsZero(G[I][J]) then
          begin
            GramOk := False;
            Inc(OffNZ);
          end;
        end;
      end;
    GRank := RMatRank(G);
    TrList := '';
    if Dim = 16 then
      for I := 0 to 15 do
      begin
        Tr := RMatTrace(B[I]);
        if I > 0 then
          TrList := TrList + ',';
        TrList := TrList + RatToString(Tr);
      end;
    if GramOk then
      ADetail := Format(
        'EXACT (Q, zero-residual): dim V %d; I16 in V %s; V closed under multiplication %s (first failure B[%d]*B[%d] not an intertwiner); Gram = 4*I16 exact (rank %d, off-diagonal nonzero %d); tr(B_i) = [%s]',
        [Dim, BoolToStr(IInV, 'YES', 'no'),
         BoolToStr(Closed, 'YES', 'no'), Fi, Fj, GRank, OffNZ, TrList])
    else
      ADetail := Format(
        'EXACT (Q, zero-residual): dim V %d; I16 in V %s; V closed %s (first failure B[%d]*B[%d]); Gram is NOT 4*I16: rank %d, off-diagonal nonzero %d (diagonal values and full matrix in the run log); tr(B_i) = [%s]',
        [Dim, BoolToStr(IInV, 'YES', 'no'),
         BoolToStr(Closed, 'YES', 'no'), Fi, Fj, GRank, OffNZ, TrList]);
    Result := (Dim = 16) and (not IInV) and (not Closed) and GramOk;
  except
    on Ex: Exception do
    begin
      ADetail := 'EXACT PORT FAILED: ' + Ex.ClassName + ': ' + Ex.Message;
      Result := False;
    end;
  end;
end;

{ --- TD3bExactProof (D3B-EXACT, proofs/d3b_forms.py) --- }

constructor TD3bExactProof.Create;
begin
  inherited Create('D3B-EXACT',
    'D3b exact port: vol = E0*...*E7 has vol^2 = +I exactly and is not scalar; the Spin(4,4)-invariant form space (28 spin generators) has exact dimension 2, both basis forms symmetric of rank exactly 8, none skew',
    psProved, 'proofs/d3b_forms.py (+ Pascal exact port, RTheory.ExactClifford)');
end;

function TD3bExactProof.VerifyNumeric(out ADetail: string): Boolean;
var
  E, SmallE: TRatMatrices;
  Omega, Vol, ST: TRatMatrix;
  SList: TRatMatrices;
  Rows: TSparseRows;
  Vecs: TRatVectors;
  M: TRatMatrix;
  I, Mm, Nn, Dim, Cnt: Integer;
  VolSqOk, VolScalar, FormsOk: Boolean;
  Cv: TRational;
  FormInfo: string;
  Sym, Skew: Boolean;
  Rk: Integer;
begin
  try
    BuildClifford(E, SmallE, Omega);
    Vol := E[0];
    for I := 1 to 7 do
      Vol := RMatMul(Vol, E[I]);
    VolSqOk := RMatEq(RMatMul(Vol, Vol), RMatIdent(16));
    VolScalar := IsScalarMat(Vol, Cv);
    SetLength(SList, 28);
    Cnt := 0;
    for Mm := 0 to 7 do
      for Nn := Mm + 1 to 7 do
      begin
        SList[Cnt] := RMatSub(RMatMul(E[Mm], E[Nn]),
          RMatMul(E[Nn], E[Mm]));
        Inc(Cnt);
      end;
    { Invariant-form system of d3b_forms.py: for each S,
      kron(I16, S^T) + kron(S^T, I16); appended via the builder as
      kron(S^T, I16) - kron(I16, -S^T), the same sum. 7168 x 256. }
    SetLength(Rows, 0);
    for I := 0 to 27 do
    begin
      ST := RMatTranspose(SList[I]);
      AppendKroneckerSystem(Rows, ST, RMatNeg(ST));
    end;
    Vecs := SparseKernelBasis(Rows, 256);
    Dim := Length(Vecs);
    FormsOk := Dim = 2;
    FormInfo := '';
    for I := 0 to Dim - 1 do
    begin
      M := VecToMat16(Vecs[I]);
      Sym := RMatEq(RMatTranspose(M), M);
      Skew := RMatEq(RMatTranspose(M), RMatNeg(M));
      Rk := RMatRank(M);
      if not (Sym and (not Skew) and (Rk = 8)) then
        FormsOk := False;
      FormInfo := FormInfo + Format(' M%d: symmetric %s, skew %s, rank %d;',
        [I, BoolToStr(Sym, 'yes', 'NO'), BoolToStr(Skew, 'yes', 'no'), Rk]);
    end;
    ADetail := Format(
      'EXACT (Q, zero-residual): vol^2 = +I %s; vol scalar %s (vol[0,0] = %s); invariant-form dim %d (expect 2);%s',
      [BoolToStr(VolSqOk, 'yes', 'NO'),
       BoolToStr(VolScalar, 'YES', 'no'), RatToString(Vol[0][0]),
       Dim, FormInfo]);
    Result := VolSqOk and (not VolScalar) and (Dim = 2) and FormsOk;
  except
    on Ex: Exception do
    begin
      ADetail := 'EXACT PORT FAILED: ' + Ex.ClassName + ': ' + Ex.Message;
      Result := False;
    end;
  end;
end;

end.
