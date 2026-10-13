unit RTheory.Operators;

{ The core operators as a Pascal object.

  Definitions (Kit's canonical IO-word table; proofs/pythagorean_foundation.md):

    srx(x) = |csc x| + cot x          ("sininverx")
    sxp(x) = 1 / srx(x)               ("sinexpon";  = |csc x| - cot x, FLIP2)
    cxp(x) = |sec x| + tan x          ("cosinexpon")
    crx(x) = 1 / cxp(x)               ("cosinverx"; = |sec x| - tan x, FLIP1)
    urx(x) = srx(x) - crx(x)          ("unainverex")
    uxp(x) = cxp(x) - sxp(x)          ("unaexpon")
    sawup(x)   = 1 / urx(x)           (quasi-linear growth)
    sawdown(x) = 1 / uxp(x)           (quasi-linear decay)
    FlatWave(x) = sawup(x) + sawdown(x) = sgn(sin 2x)   (RECIP / FOUNDATION)

  A TOperatorPoint is an immutable snapshot of every operator at one x.
  The operators are singular where sin x = 0 or cos x = 0; constructing a
  snapshot there raises EOperatorDomain. }

{$mode objfpc}{$H+}

interface

uses
  SysUtils;

type
  EOperatorDomain = class(Exception);

  TOperatorPoint = class
  private
    FX: Double;
    FSrx, FSxp, FCxp, FCrx: Double;
    FUrx, FUxp: Double;
    FSawUp, FSawDown, FFlatWave: Double;
  public
    constructor Create(AX: Double);
    property X: Double read FX;
    { the primitive quartet }
    property Srx: Double read FSrx;
    property Sxp: Double read FSxp;
    property Cxp: Double read FCxp;
    property Crx: Double read FCrx;
    { the una pair }
    property Urx: Double read FUrx;
    property Uxp: Double read FUxp;
    { the saws and their sum }
    property SawUp: Double read FSawUp;
    property SawDown: Double read FSawDown;
    property FlatWave: Double read FFlatWave;
    { double-angle combinations (AB / COSRAT): A = 2 cot x, B = 2 tan x }
    function A: Double;   { srx - sxp }
    function B: Double;   { cxp - crx }
    function K2: Double;  { A - B }
    function D2: Double;  { A + B }
  end;

{ Raw class functions, for callers that want a single value. }
  TCoreOperators = class
  public
    class function Srx(AX: Double): Double;
    class function Sxp(AX: Double): Double;
    class function Cxp(AX: Double): Double;
    class function Crx(AX: Double): Double;
    class function Urx(AX: Double): Double;
    class function Uxp(AX: Double): Double;
    class function SawUp(AX: Double): Double;
    class function SawDown(AX: Double): Double;
    class function FlatWave(AX: Double): Double;
    class function FlatWaveSign(AX: Double): Integer; { sgn(sin 2x) }
  end;

implementation

const
  EPS_SINGULAR = 1e-12;

function SafeCsc(AX: Double): Double;
var
  S: Double;
begin
  S := Sin(AX);
  if Abs(S) < EPS_SINGULAR then
    raise EOperatorDomain.CreateFmt('csc singular at x=%g', [AX]);
  Result := 1.0 / S;
end;

function SafeSec(AX: Double): Double;
var
  C: Double;
begin
  C := Cos(AX);
  if Abs(C) < EPS_SINGULAR then
    raise EOperatorDomain.CreateFmt('sec singular at x=%g', [AX]);
  Result := 1.0 / C;
end;

function Cot(AX: Double): Double;
begin
  Result := Cos(AX) * SafeCsc(AX);
end;

function TanX(AX: Double): Double;
begin
  Result := Sin(AX) * SafeSec(AX);
end;

{ TOperatorPoint }

constructor TOperatorPoint.Create(AX: Double);
begin
  inherited Create;
  FX := AX;
  FSrx := Abs(SafeCsc(AX)) + Cot(AX);
  FSxp := 1.0 / FSrx;
  FCxp := Abs(SafeSec(AX)) + TanX(AX);
  FCrx := 1.0 / FCxp;
  FUrx := FSrx - FCrx;
  FUxp := FCxp - FSxp;
  if Abs(FUrx) < EPS_SINGULAR then
    raise EOperatorDomain.CreateFmt('urx = 0 at x=%g (sawup singular)', [AX]);
  if Abs(FUxp) < EPS_SINGULAR then
    raise EOperatorDomain.CreateFmt('uxp = 0 at x=%g (sawdown singular)', [AX]);
  FSawUp := 1.0 / FUrx;
  FSawDown := 1.0 / FUxp;
  FFlatWave := FSawUp + FSawDown;
end;

function TOperatorPoint.A: Double;
begin
  Result := FSrx - FSxp;
end;

function TOperatorPoint.B: Double;
begin
  Result := FCxp - FCrx;
end;

function TOperatorPoint.K2: Double;
begin
  Result := A - B;
end;

function TOperatorPoint.D2: Double;
begin
  Result := A + B;
end;

{ TCoreOperators }

class function TCoreOperators.Srx(AX: Double): Double;
begin
  Result := Abs(SafeCsc(AX)) + Cot(AX);
end;

class function TCoreOperators.Sxp(AX: Double): Double;
begin
  Result := 1.0 / Srx(AX);
end;

class function TCoreOperators.Cxp(AX: Double): Double;
begin
  Result := Abs(SafeSec(AX)) + TanX(AX);
end;

class function TCoreOperators.Crx(AX: Double): Double;
begin
  Result := 1.0 / Cxp(AX);
end;

class function TCoreOperators.Urx(AX: Double): Double;
begin
  Result := Srx(AX) - Crx(AX);
end;

class function TCoreOperators.Uxp(AX: Double): Double;
begin
  Result := Cxp(AX) - Sxp(AX);
end;

class function TCoreOperators.SawUp(AX: Double): Double;
begin
  Result := 1.0 / Urx(AX);
end;

class function TCoreOperators.SawDown(AX: Double): Double;
begin
  Result := 1.0 / Uxp(AX);
end;

class function TCoreOperators.FlatWave(AX: Double): Double;
begin
  Result := SawUp(AX) + SawDown(AX);
end;

class function TCoreOperators.FlatWaveSign(AX: Double): Integer;
var
  S: Double;
begin
  S := Sin(2.0 * AX);
  if S > 0 then
    Result := 1
  else if S < 0 then
    Result := -1
  else
    Result := 0;
end;

end.
