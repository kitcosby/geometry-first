unit RTheory.RatFunc;

{ Exact rational functions over Q, plus Laurent polynomials over Q(i).

  TPoly: polynomial over Q with TRational coefficients (arbitrary
  precision, from RTheory.ExactField), stored low-to-high, normalized
  (no trailing zero coefficients; the zero polynomial is empty).

  TRatFunc: Num/Den with TPoly parts, kept normalized: gcd cancelled
  and Den monic. Equality is therefore exact componentwise equality
  of the normalized forms. Mobius substitution t -> (a t + b)/(c t + d)
  is exact via reciprocal-style power clearing.

  TLaurent: Laurent polynomial over Q(i), coefficients indexed by
  integer powers (Offset is the exponent of Coeffs[0]). Only Add,
  Sub, Mul and ratio-equality by cross-multiplication are provided --
  that is all quadrant_flip_polar_a's w-identities need; no Laurent
  gcd is claimed or used.

  No floating point anywhere. }

{$mode objfpc}{$H+}

interface

uses
  SysUtils, RTheory.ExactField;

type
  TPoly = record
    Coeffs: array of TRational; { Coeffs[k] = coeff of t^k }
  end;

  TPolyArray = array of TPoly;

  TRatFunc = record
    Num, Den: TPoly;
  end;

  TRFMatrix = array of array of TRatFunc;

  TComplexRat = record
    Re, Im: TRational;
  end;

  TLaurent = record
    Offset: Integer; { exponent of Coeffs[0] }
    Coeffs: array of TComplexRat;
  end;

{ --- polynomials over Q --- }

function PolyZero: TPoly;
function PolyOne: TPoly;
function PolyConst(const R: TRational): TPoly;
function PolyVar: TPoly; { t }
function PolyMake(const Coeffs: array of TRational): TPoly;
function PolyNorm(const P: TPoly): TPoly;
function PolyDegree(const P: TPoly): Integer; { -1 for zero }
function PolyIsZero(const P: TPoly): Boolean;
function PolyEq(const A, B: TPoly): Boolean;
function PolyCoeff(const P: TPoly; K: Integer): TRational;
function PolyLead(const P: TPoly): TRational;
function PolyAdd(const A, B: TPoly): TPoly;
function PolySub(const A, B: TPoly): TPoly;
function PolyNeg(const P: TPoly): TPoly;
function PolyMul(const A, B: TPoly): TPoly;
function PolyScale(const P: TPoly; const R: TRational): TPoly;
procedure PolyDivMod(const A, B: TPoly; out Q, R: TPoly);
function PolyDivExact(const A, B: TPoly): TPoly;
function PolyGcd(const A, B: TPoly): TPoly; { monic, unless both zero }
function PolyDerivative(const P: TPoly): TPoly;
function PolyEvalField(const P: TPoly; const X: TField): TField;
function RatToField(const R: TRational): TField;
function PolyToString(const P: TPoly): string;

{ --- rational functions over Q --- }

function RFMake(const ANum, ADen: TPoly): TRatFunc;
function RFZero: TRatFunc;
function RFOne: TRatFunc;
function RFConst(const R: TRational): TRatFunc;
function RFInt(N, D: Int64): TRatFunc;
function RFVar: TRatFunc; { t }
function RFIsZero(const F: TRatFunc): Boolean;
function RFEq(const A, B: TRatFunc): Boolean;
function RFAdd(const A, B: TRatFunc): TRatFunc;
function RFSub(const A, B: TRatFunc): TRatFunc;
function RFNeg(const F: TRatFunc): TRatFunc;
function RFMul(const A, B: TRatFunc): TRatFunc;
function RFInv(const F: TRatFunc): TRatFunc;
function RFDiv(const A, B: TRatFunc): TRatFunc;
function RFDerivative(const F: TRatFunc): TRatFunc;
{ Substitute t -> (a t + b)/(c t + d), rational a,b,c,d, c t + d <> 0. }
function RFMobius(const F: TRatFunc; const A, B, C, D: TRational): TRatFunc;
function RFMatDet(const M: TRFMatrix): TRatFunc; { n <= 4, Leibniz }
function RFToString(const F: TRatFunc): string;

{ --- Laurent polynomials over Q(i) --- }

function CR(const ARe, AIm: TRational): TComplexRat;
function CRIsZero(const Z: TComplexRat): Boolean;
function CREq(const Z, W: TComplexRat): Boolean;
function LCZero: TLaurent;
function LCTerm(Power: Integer; const Coeff: TComplexRat): TLaurent;
function LCAdd(const A, B: TLaurent): TLaurent;
function LCSub(const A, B: TLaurent): TLaurent;
function LCMul(const A, B: TLaurent): TLaurent;
function LCScale(const A: TLaurent; const Z: TComplexRat): TLaurent;
function LCIsZero(const A: TLaurent): Boolean;
function LCEq(const A, B: TLaurent): Boolean;
{ num1/den1 == num2/den2, by cross-multiplication. }
function LCRatioEq(const N1, D1, N2, D2: TLaurent): Boolean;

implementation

function RZero: TRational;
begin
  Result := Rat(0, 1);
end;

function ROne: TRational;
begin
  Result := Rat(1, 1);
end;

{ ================= polynomials ================= }

function PolyCloneTrim(const P: TPoly; Len: Integer): TPoly;
var
  I: Integer;
begin
  SetLength(Result.Coeffs, Len);
  for I := 0 to Len - 1 do
    Result.Coeffs[I] := P.Coeffs[I];
end;

function PolyNorm(const P: TPoly): TPoly;
var
  N: Integer;
begin
  N := Length(P.Coeffs);
  while (N > 0) and RatIsZero(P.Coeffs[N - 1]) do
    Dec(N);
  Result := PolyCloneTrim(P, N);
end;

function PolyZero: TPoly;
begin
  SetLength(Result.Coeffs, 0);
end;

function PolyMake(const Coeffs: array of TRational): TPoly;
var
  I: Integer;
  Tmp: TPoly;
begin
  SetLength(Tmp.Coeffs, Length(Coeffs));
  for I := 0 to High(Coeffs) do
    Tmp.Coeffs[I] := Coeffs[I];
  Result := PolyNorm(Tmp);
end;

function PolyConst(const R: TRational): TPoly;
begin
  if RatIsZero(R) then
    Result := PolyZero
  else
  begin
    SetLength(Result.Coeffs, 1);
    Result.Coeffs[0] := R;
  end;
end;

function PolyOne: TPoly;
begin
  Result := PolyConst(ROne);
end;

function PolyVar: TPoly;
begin
  SetLength(Result.Coeffs, 2);
  Result.Coeffs[0] := RZero;
  Result.Coeffs[1] := ROne;
end;

function PolyDegree(const P: TPoly): Integer;
begin
  Result := Length(PolyNorm(P).Coeffs) - 1;
end;

function PolyIsZero(const P: TPoly): Boolean;
begin
  Result := PolyDegree(P) < 0;
end;

function PolyCoeff(const P: TPoly; K: Integer): TRational;
begin
  if (K < 0) or (K >= Length(P.Coeffs)) then
    Result := RZero
  else
    Result := P.Coeffs[K];
end;

function PolyLead(const P: TPoly): TRational;
var
  N: TPoly;
begin
  N := PolyNorm(P);
  if Length(N.Coeffs) = 0 then
    Result := RZero
  else
    Result := N.Coeffs[High(N.Coeffs)];
end;

function PolyEq(const A, B: TPoly): Boolean;
var
  NA, NB: TPoly;
  I: Integer;
begin
  NA := PolyNorm(A);
  NB := PolyNorm(B);
  if Length(NA.Coeffs) <> Length(NB.Coeffs) then
    Exit(False);
  for I := 0 to High(NA.Coeffs) do
    if not RatEq(NA.Coeffs[I], NB.Coeffs[I]) then
      Exit(False);
  Result := True;
end;

function PolyAdd(const A, B: TPoly): TPoly;
var
  I, N: Integer;
  Tmp: TPoly;
begin
  N := Length(A.Coeffs);
  if Length(B.Coeffs) > N then
    N := Length(B.Coeffs);
  SetLength(Tmp.Coeffs, N);
  for I := 0 to N - 1 do
    Tmp.Coeffs[I] := RatAdd(PolyCoeff(A, I), PolyCoeff(B, I));
  Result := PolyNorm(Tmp);
end;

function PolyNeg(const P: TPoly): TPoly;
var
  I: Integer;
begin
  SetLength(Result.Coeffs, Length(P.Coeffs));
  for I := 0 to High(P.Coeffs) do
    Result.Coeffs[I] := RatNeg(P.Coeffs[I]);
  Result := PolyNorm(Result);
end;

function PolySub(const A, B: TPoly): TPoly;
begin
  Result := PolyAdd(A, PolyNeg(B));
end;

function PolyScale(const P: TPoly; const R: TRational): TPoly;
var
  I: Integer;
begin
  if RatIsZero(R) then
    Exit(PolyZero);
  SetLength(Result.Coeffs, Length(P.Coeffs));
  for I := 0 to High(P.Coeffs) do
    Result.Coeffs[I] := RatMul(P.Coeffs[I], R);
  Result := PolyNorm(Result);
end;

function PolyMul(const A, B: TPoly): TPoly;
var
  I, J: Integer;
  Tmp: TPoly;
begin
  if PolyIsZero(A) or PolyIsZero(B) then
    Exit(PolyZero);
  SetLength(Tmp.Coeffs, Length(A.Coeffs) + Length(B.Coeffs) - 1);
  for I := 0 to High(Tmp.Coeffs) do
    Tmp.Coeffs[I] := RZero;
  for I := 0 to High(A.Coeffs) do
    for J := 0 to High(B.Coeffs) do
      Tmp.Coeffs[I + J] := RatAdd(Tmp.Coeffs[I + J],
        RatMul(A.Coeffs[I], B.Coeffs[J]));
  Result := PolyNorm(Tmp);
end;

procedure PolyDivMod(const A, B: TPoly; out Q, R: TPoly);
var
  NA, NB, RR: TPoly;
  DA, DB, K, J: Integer;
  InvLC, Coef: TRational;
  Tmp: TPoly;
begin
  NB := PolyNorm(B);
  if Length(NB.Coeffs) = 0 then
    raise EFieldError.Create('PolyDivMod: division by zero polynomial');
  NA := PolyNorm(A);
  DA := Length(NA.Coeffs) - 1;
  DB := Length(NB.Coeffs) - 1;
  if DA < DB then
  begin
    Q := PolyZero;
    R := NA;
    Exit;
  end;
  RR := PolyCloneTrim(NA, Length(NA.Coeffs));
  SetLength(RR.Coeffs, DA + 1);
  SetLength(Tmp.Coeffs, DA - DB + 1);
  for J := 0 to High(Tmp.Coeffs) do
    Tmp.Coeffs[J] := RZero;
  InvLC := RatInv(NB.Coeffs[DB]);
  for K := DA - DB downto 0 do
  begin
    Coef := RatMul(RR.Coeffs[DB + K], InvLC);
    Tmp.Coeffs[K] := Coef;
    if not RatIsZero(Coef) then
      for J := 0 to DB do
        RR.Coeffs[J + K] := RatSub(RR.Coeffs[J + K],
          RatMul(Coef, NB.Coeffs[J]));
  end;
  Q := PolyNorm(Tmp);
  R := PolyNorm(RR);
end;

function PolyDivExact(const A, B: TPoly): TPoly;
var
  Q, R: TPoly;
begin
  PolyDivMod(A, B, Q, R);
  if not PolyIsZero(R) then
    raise EFieldError.Create('PolyDivExact: nonzero remainder');
  Result := Q;
end;

function PolyGcd(const A, B: TPoly): TPoly;
var
  X, Y, Q, R: TPoly;
begin
  X := PolyNorm(A);
  Y := PolyNorm(B);
  if PolyIsZero(X) and PolyIsZero(Y) then
    Exit(PolyZero);
  if PolyIsZero(X) then
    Exit(PolyScale(Y, RatInv(PolyLead(Y))));
  if PolyIsZero(Y) then
    Exit(PolyScale(X, RatInv(PolyLead(X))));
  while not PolyIsZero(Y) do
  begin
    PolyDivMod(X, Y, Q, R);
    X := Y;
    Y := R;
  end;
  Result := PolyScale(X, RatInv(PolyLead(X)));
end;

function PolyDerivative(const P: TPoly): TPoly;
var
  I: Integer;
  Tmp: TPoly;
begin
  if Length(P.Coeffs) <= 1 then
    Exit(PolyZero);
  SetLength(Tmp.Coeffs, Length(P.Coeffs) - 1);
  for I := 1 to High(P.Coeffs) do
    Tmp.Coeffs[I - 1] := RatMul(P.Coeffs[I], Rat(I, 1));
  Result := PolyNorm(Tmp);
end;

function RatToField(const R: TRational): TField;
begin
  Result := Fld(R, RZero, RZero, RZero);
end;

function PolyEvalField(const P: TPoly; const X: TField): TField;
var
  I: Integer;
  N: TPoly;
begin
  N := PolyNorm(P);
  Result := FldI(0, 1);
  for I := High(N.Coeffs) downto 0 do
    Result := FAdd(FMul(Result, X), RatToField(N.Coeffs[I]));
end;

function PolyToString(const P: TPoly): string;
var
  N: TPoly;
  I: Integer;
begin
  N := PolyNorm(P);
  if Length(N.Coeffs) = 0 then
    Exit('0');
  Result := '[';
  for I := 0 to High(N.Coeffs) do
  begin
    if I > 0 then
      Result := Result + ',';
    Result := Result + RatToString(N.Coeffs[I]);
  end;
  Result := Result + '] (coeffs t^0..t^' + IntToStr(High(N.Coeffs)) + ')';
end;

{ ================= rational functions ================= }

function RFMake(const ANum, ADen: TPoly): TRatFunc;
var
  N, D, G: TPoly;
  LC: TRational;
begin
  D := PolyNorm(ADen);
  if PolyIsZero(D) then
    raise EFieldError.Create('RFMake: zero denominator');
  N := PolyNorm(ANum);
  if PolyIsZero(N) then
  begin
    Result.Num := PolyZero;
    Result.Den := PolyOne;
    Exit;
  end;
  G := PolyGcd(N, D);
  if PolyDegree(G) > 0 then
  begin
    { monic gcd of positive degree: cancel it. A constant (monic = 1)
      gcd is absorbed by the monic-denominator scaling below. }
    N := PolyDivExact(N, G);
    D := PolyDivExact(D, G);
  end;
  LC := PolyLead(D);
  if not RatEq(LC, ROne) then
  begin
    N := PolyScale(N, RatInv(LC));
    D := PolyScale(D, RatInv(LC));
  end;
  Result.Num := N;
  Result.Den := D;
end;

function RFZero: TRatFunc;
begin
  Result.Num := PolyZero;
  Result.Den := PolyOne;
end;

function RFOne: TRatFunc;
begin
  Result := RFConst(ROne);
end;

function RFConst(const R: TRational): TRatFunc;
begin
  Result.Num := PolyConst(R);
  Result.Den := PolyOne;
end;

function RFInt(N, D: Int64): TRatFunc;
begin
  Result := RFConst(Rat(N, D));
end;

function RFVar: TRatFunc;
begin
  Result.Num := PolyVar;
  Result.Den := PolyOne;
end;

function RFIsZero(const F: TRatFunc): Boolean;
begin
  Result := PolyIsZero(F.Num);
end;

function RFEq(const A, B: TRatFunc): Boolean;
begin
  Result := PolyEq(A.Num, B.Num) and PolyEq(A.Den, B.Den);
end;

function RFNeg(const F: TRatFunc): TRatFunc;
begin
  Result.Num := PolyNeg(F.Num);
  Result.Den := F.Den;
end;

function RFAdd(const A, B: TRatFunc): TRatFunc;
begin
  Result := RFMake(
    PolyAdd(PolyMul(A.Num, B.Den), PolyMul(B.Num, A.Den)),
    PolyMul(A.Den, B.Den));
end;

function RFSub(const A, B: TRatFunc): TRatFunc;
begin
  Result := RFAdd(A, RFNeg(B));
end;

function RFMul(const A, B: TRatFunc): TRatFunc;
begin
  Result := RFMake(PolyMul(A.Num, B.Num), PolyMul(A.Den, B.Den));
end;

function RFInv(const F: TRatFunc): TRatFunc;
begin
  if RFIsZero(F) then
    raise EFieldError.Create('RFInv: inverse of zero');
  Result := RFMake(F.Den, F.Num);
end;

function RFDiv(const A, B: TRatFunc): TRatFunc;
begin
  Result := RFMul(A, RFInv(B));
end;

function RFDerivative(const F: TRatFunc): TRatFunc;
begin
  Result := RFMake(
    PolySub(PolyMul(PolyDerivative(F.Num), F.Den),
      PolyMul(F.Num, PolyDerivative(F.Den))),
    PolyMul(F.Den, F.Den));
end;

function PolyPow(const P: TPoly; E: Integer): TPoly;
var
  I: Integer;
begin
  Result := PolyOne;
  for I := 1 to E do
    Result := PolyMul(Result, P);
end;

{ (c t + d)^deg(P) * P((a t + b)/(c t + d)) as a polynomial. }
function PolyMobiusTransform(const P: TPoly;
  const A, B, C, D: TRational): TPoly;
var
  N: TPoly;
  DegP, K: Integer;
  L1, L2: TPoly;
  Pow1, Pow2: TPolyArray;
  Term: TPoly;
begin
  N := PolyNorm(P);
  if PolyIsZero(N) then
    Exit(PolyZero);
  DegP := Length(N.Coeffs) - 1;
  L1 := PolyMake([B, A]); { a t + b }
  L2 := PolyMake([D, C]); { c t + d }
  SetLength(Pow1, DegP + 1);
  SetLength(Pow2, DegP + 1);
  for K := 0 to DegP do
  begin
    Pow1[K] := PolyPow(L1, K);
    Pow2[K] := PolyPow(L2, K);
  end;
  Result := PolyZero;
  for K := 0 to DegP do
    if not RatIsZero(N.Coeffs[K]) then
    begin
      Term := PolyMul(Pow1[K], Pow2[DegP - K]);
      Term := PolyScale(Term, N.Coeffs[K]);
      Result := PolyAdd(Result, Term);
    end;
end;

function RFMobius(const F: TRatFunc; const A, B, C, D: TRational): TRatFunc;
var
  NP, DP: TPoly;
  DN, DD, Diff: Integer;
  L2: TPoly;
begin
  NP := PolyMobiusTransform(F.Num, A, B, C, D);
  DP := PolyMobiusTransform(F.Den, A, B, C, D);
  DN := PolyDegree(F.Num);
  if DN < 0 then
    DN := 0;
  DD := PolyDegree(F.Den);
  if DD < 0 then
    DD := 0;
  Diff := DD - DN;
  L2 := PolyMake([D, C]);
  if Diff > 0 then
    NP := PolyMul(NP, PolyPow(L2, Diff))
  else if Diff < 0 then
    DP := PolyMul(DP, PolyPow(L2, -Diff));
  Result := RFMake(NP, DP);
end;

function RFMatDet(const M: TRFMatrix): TRatFunc;
var
  N, I, Counter, Total, X, Digit: Integer;
  Perm: array of Integer;
  Used: array of Boolean;
  Distinct: Boolean;
  Inversions: Integer;
  Term: TRatFunc;
  A2, B2: Integer;
begin
  N := Length(M);
  if N = 0 then
    Exit(RFOne);
  Total := 1;
  for I := 1 to N do
    Total := Total * N;
  SetLength(Perm, N);
  SetLength(Used, N);
  Result := RFZero;
  for Counter := 0 to Total - 1 do
  begin
    X := Counter;
    for I := 0 to N - 1 do
      Used[I] := False;
    Distinct := True;
    for I := 0 to N - 1 do
    begin
      Digit := X mod N;
      X := X div N;
      Perm[I] := Digit;
      if Used[Digit] then
        Distinct := False;
      Used[Digit] := True;
    end;
    if not Distinct then
      Continue;
    Inversions := 0;
    for A2 := 0 to N - 2 do
      for B2 := A2 + 1 to N - 1 do
        if Perm[A2] > Perm[B2] then
          Inc(Inversions);
    Term := RFOne;
    for I := 0 to N - 1 do
      Term := RFMul(Term, M[I][Perm[I]]);
    if Odd(Inversions) then
      Result := RFSub(Result, Term)
    else
      Result := RFAdd(Result, Term);
  end;
end;

function RFToString(const F: TRatFunc): string;
begin
  Result := '(' + PolyToString(F.Num) + ')/(' + PolyToString(F.Den) + ')';
end;

{ ================= Laurent over Q(i) ================= }

function CR(const ARe, AIm: TRational): TComplexRat;
begin
  Result.Re := ARe;
  Result.Im := AIm;
end;

function CRIsZero(const Z: TComplexRat): Boolean;
begin
  Result := RatIsZero(Z.Re) and RatIsZero(Z.Im);
end;

function CREq(const Z, W: TComplexRat): Boolean;
begin
  Result := RatEq(Z.Re, W.Re) and RatEq(Z.Im, W.Im);
end;

function CRAdd(const Z, W: TComplexRat): TComplexRat;
begin
  Result.Re := RatAdd(Z.Re, W.Re);
  Result.Im := RatAdd(Z.Im, W.Im);
end;

function CRNeg(const Z: TComplexRat): TComplexRat;
begin
  Result.Re := RatNeg(Z.Re);
  Result.Im := RatNeg(Z.Im);
end;

function CRMul(const Z, W: TComplexRat): TComplexRat;
begin
  Result.Re := RatSub(RatMul(Z.Re, W.Re), RatMul(Z.Im, W.Im));
  Result.Im := RatAdd(RatMul(Z.Re, W.Im), RatMul(Z.Im, W.Re));
end;

function LCNorm(const A: TLaurent): TLaurent;
var
  Lo, Hi, I: Integer;
begin
  Lo := 0;
  Hi := Length(A.Coeffs) - 1;
  while (Lo <= Hi) and CRIsZero(A.Coeffs[Lo]) do
    Inc(Lo);
  while (Hi >= Lo) and CRIsZero(A.Coeffs[Hi]) do
    Dec(Hi);
  if Lo > Hi then
  begin
    Result.Offset := 0;
    SetLength(Result.Coeffs, 0);
    Exit;
  end;
  Result.Offset := A.Offset + Lo;
  SetLength(Result.Coeffs, Hi - Lo + 1);
  for I := Lo to Hi do
    Result.Coeffs[I - Lo] := A.Coeffs[I];
end;

function LCZero: TLaurent;
begin
  Result.Offset := 0;
  SetLength(Result.Coeffs, 0);
end;

function LCTerm(Power: Integer; const Coeff: TComplexRat): TLaurent;
begin
  if CRIsZero(Coeff) then
    Exit(LCZero);
  Result.Offset := Power;
  SetLength(Result.Coeffs, 1);
  Result.Coeffs[0] := Coeff;
end;

function LCFetch(const A: TLaurent; Power: Integer): TComplexRat;
var
  Idx: Integer;
begin
  Idx := Power - A.Offset;
  if (Idx < 0) or (Idx >= Length(A.Coeffs)) then
    Result := CR(RZero, RZero)
  else
    Result := A.Coeffs[Idx];
end;

function LCAdd(const A, B: TLaurent): TLaurent;
var
  NA, NB: TLaurent;
  Lo, Hi, I: Integer;
  Tmp: TLaurent;
begin
  NA := LCNorm(A);
  NB := LCNorm(B);
  if Length(NA.Coeffs) = 0 then
    Exit(NB);
  if Length(NB.Coeffs) = 0 then
    Exit(NA);
  Lo := NA.Offset;
  if NB.Offset < Lo then
    Lo := NB.Offset;
  Hi := NA.Offset + Length(NA.Coeffs) - 1;
  if NB.Offset + Length(NB.Coeffs) - 1 > Hi then
    Hi := NB.Offset + Length(NB.Coeffs) - 1;
  Tmp.Offset := Lo;
  SetLength(Tmp.Coeffs, Hi - Lo + 1);
  for I := Lo to Hi do
    Tmp.Coeffs[I - Lo] := CRAdd(LCFetch(NA, I), LCFetch(NB, I));
  Result := LCNorm(Tmp);
end;

function LCSub(const A, B: TLaurent): TLaurent;
var
  NA, NB: TLaurent;
  Lo, Hi, I: Integer;
  Tmp: TLaurent;
  ZA, ZB: TComplexRat;
begin
  NA := LCNorm(A);
  NB := LCNorm(B);
  if Length(NB.Coeffs) = 0 then
    Exit(NA);
  if Length(NA.Coeffs) = 0 then
  begin
    Result.Offset := NB.Offset;
    SetLength(Result.Coeffs, Length(NB.Coeffs));
    for I := 0 to High(NB.Coeffs) do
      Result.Coeffs[I] := CRNeg(NB.Coeffs[I]);
    Exit;
  end;
  Lo := NA.Offset;
  if NB.Offset < Lo then
    Lo := NB.Offset;
  Hi := NA.Offset + Length(NA.Coeffs) - 1;
  if NB.Offset + Length(NB.Coeffs) - 1 > Hi then
    Hi := NB.Offset + Length(NB.Coeffs) - 1;
  Tmp.Offset := Lo;
  SetLength(Tmp.Coeffs, Hi - Lo + 1);
  for I := Lo to Hi do
  begin
    ZA := LCFetch(NA, I);
    ZB := LCFetch(NB, I);
    Tmp.Coeffs[I - Lo].Re := RatSub(ZA.Re, ZB.Re);
    Tmp.Coeffs[I - Lo].Im := RatSub(ZA.Im, ZB.Im);
  end;
  Result := LCNorm(Tmp);
end;

function LCMul(const A, B: TLaurent): TLaurent;
var
  NA, NB, Tmp: TLaurent;
  I, J: Integer;
begin
  NA := LCNorm(A);
  NB := LCNorm(B);
  if (Length(NA.Coeffs) = 0) or (Length(NB.Coeffs) = 0) then
    Exit(LCZero);
  Tmp.Offset := NA.Offset + NB.Offset;
  SetLength(Tmp.Coeffs, Length(NA.Coeffs) + Length(NB.Coeffs) - 1);
  for I := 0 to High(Tmp.Coeffs) do
    Tmp.Coeffs[I] := CR(RZero, RZero);
  for I := 0 to High(NA.Coeffs) do
    for J := 0 to High(NB.Coeffs) do
      Tmp.Coeffs[I + J] := CRAdd(Tmp.Coeffs[I + J],
        CRMul(NA.Coeffs[I], NB.Coeffs[J]));
  Result := LCNorm(Tmp);
end;

function LCScale(const A: TLaurent; const Z: TComplexRat): TLaurent;
var
  I: Integer;
begin
  if CRIsZero(Z) then
    Exit(LCZero);
  Result.Offset := A.Offset;
  SetLength(Result.Coeffs, Length(A.Coeffs));
  for I := 0 to High(A.Coeffs) do
    Result.Coeffs[I] := CRMul(A.Coeffs[I], Z);
  Result := LCNorm(Result);
end;

function LCIsZero(const A: TLaurent): Boolean;
begin
  Result := Length(LCNorm(A).Coeffs) = 0;
end;

function LCEq(const A, B: TLaurent): Boolean;
var
  NA, NB: TLaurent;
  I: Integer;
begin
  NA := LCNorm(A);
  NB := LCNorm(B);
  if Length(NA.Coeffs) <> Length(NB.Coeffs) then
    Exit(False);
  if Length(NA.Coeffs) = 0 then
    Exit(True);
  if NA.Offset <> NB.Offset then
    Exit(False);
  for I := 0 to High(NA.Coeffs) do
    if not CREq(NA.Coeffs[I], NB.Coeffs[I]) then
      Exit(False);
  Result := True;
end;

function LCRatioEq(const N1, D1, N2, D2: TLaurent): Boolean;
begin
  Result := LCIsZero(LCSub(LCMul(N1, D2), LCMul(N2, D1)));
end;

end.
