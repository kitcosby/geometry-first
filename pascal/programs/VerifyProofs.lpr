program VerifyProofs;

{ Geometry-first Pascal layer -- console verifier. (v2)

  1. Builds the proof registry (one Pascal object per theorems/INDEX.md
     entry, pre-Lie sections + the Lie campaign section + the wall
     section of docs/WALL.md) and prints it with its scope labels.
  2. Runs the numerically verifiable proof objects (RTheory.Verify,
     RTheory.LieVerify and RTheory.WallVerify) over the eight octant
     midpoints, including numeric Lie-closure checks (sl(8) dim 63,
     so(9) dim 36, sl(3) families, sl(2) block) and the wall probes
     (chain stall dim 10, Wick block closure dim 56, commutant dim 1).
  3. Exits 0 iff every verifiable object passes.

  Build (Free Pascal):
    cd pascal/programs
    fpc -Fu../src VerifyProofs.lpr
  or from pascal/:
    fpc -Fusrc programs/VerifyProofs.lpr }

{$mode objfpc}{$H+}

uses
  SysUtils, RTheory.Types, RTheory.Operators, RTheory.Theorems,
  RTheory.Verify, RTheory.Lie, RTheory.LieVerify,
  RTheory.Wall, RTheory.WallVerify;

var
  Registry: TProofRegistry;
  Checks: array of TProofObject;
  I, Failures: Integer;
  Detail: string;
  P: TOperatorPoint;

begin
  Registry := TProofRegistry.Create;
  try
    RegisterAllProofs(Registry);
    RegisterLieCampaign(Registry);
    RegisterWall(Registry);
    WriteLn('=== Proof registry (theorems/INDEX.md as Pascal objects) ===');
    for I := 0 to Registry.Count - 1 do
      WriteLn(Registry[I].Describe);
    WriteLn(Format('Total: %d objects; PROVED: %d; STRONG-EVIDENCE: %d; ASSERTED: %d; OPEN: %d',
      [Registry.Count, Registry.CountByStatus(psProved),
       Registry.CountByStatus(psStrongEvidence),
       Registry.CountByStatus(psAsserted),
       Registry.CountByStatus(psOpen)]));
    WriteLn('PROVED count includes mixed-scope entries whose exact split');
    WriteLn('(e.g. PROVED/CHECKED) is carried verbatim in their scope label.');
    WriteLn;

    WriteLn('=== Operator snapshot at x = pi/8 (octant 0 midpoint) ===');
    P := TOperatorPoint.Create(Pi / 8.0);
    try
      WriteLn(Format('srx=%.12f sxp=%.12f cxp=%.12f crx=%.12f',
        [P.Srx, P.Sxp, P.Cxp, P.Crx]));
      WriteLn(Format('urx=%.12f uxp=%.12f FlatWave=%.12f',
        [P.Urx, P.Uxp, P.FlatWave]));
    finally
      P.Free;
    end;
    WriteLn;

    WriteLn('=== Numeric verification (CHECKED-strength confirmation) ===');
    SetLength(Checks, 28);
    Checks[0] := NewRecipProof;
    Checks[1] := NewFlip1Proof;
    Checks[2] := NewFlip2Proof;
    Checks[3] := NewInvPairProof;
    Checks[4] := NewAbProof;
    Checks[5] := NewCosRatProof;
    Checks[6] := NewSinProdProof;
    Checks[7] := NewHyperProof;
    Checks[8] := NewCircleProof;
    Checks[9] := NewMobiusProof;
    Checks[10] := TE77Proof.Create;
    Checks[11] := TFlatWaveSignProof.Create;
    Checks[12] := TQ7Proof.Create('Q7',
      'dS/dx = (srx*|csc x| - crx*|sec x|) * S^2, S = sawup (Bernoulli, exact per-octant)',
      psProved, 'proofs/q7_dSdx.py');
    Checks[13] := TDurxProof.Create('DURX',
      'd(urx)/dx = crx*|sec x| - srx*|csc x| (exact)',
      psProved, 'Q7 + S = 1/urx');
    Checks[14] := THarmProof.Create('HARM',
      'F = (sin 2x / 4, cos 2x / 2), dF/dx = [[0,1],[-4,0]] * F',
      psProved, 'trig');
    { Lie campaign verifiable objects }
    Checks[15] := TCycleSl8Proof.Create;
    Checks[16] := TSo9PathProof.Create;
    Checks[17] := TTaFamilyProof.Create;
    Checks[18] := TE106Proof.Create;
    Checks[19] := TE107Proof.Create;
    Checks[20] := TCycProdProof.Create;
    Checks[21] := TM1M2ScalarProof.Create;
    Checks[22] := TDeltasProof.Create;
    Checks[23] := TXFamProof.Create;
    Checks[24] := TSo3ThetaProof.Create;
    { Wall probes (RTheory.WallVerify) }
    Checks[25] := TWallChain4Proof.Create;
    Checks[26] := TWallWickProof.Create;
    Checks[27] := TWallCommutantProof.Create;

    Failures := 0;
    for I := 0 to High(Checks) do
    begin
      if Checks[I].VerifyNumeric(Detail) then
        WriteLn(Format('PASS %-8s %s', [Checks[I].Id, Detail]))
      else
      begin
        WriteLn(Format('FAIL %-8s %s', [Checks[I].Id, Detail]));
        Inc(Failures);
      end;
    end;
    for I := 0 to High(Checks) do
      Checks[I].Free;
    WriteLn;
    if Failures = 0 then
    begin
      WriteLn('ALL CHECKS PASS');
      Halt(0);
    end
    else
    begin
      WriteLn(Format('%d CHECK(S) FAILED', [Failures]));
      Halt(1);
    end;
  finally
    Registry.Free;
  end;
end.
