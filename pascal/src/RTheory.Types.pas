unit RTheory.Types;

{ Geometry-first Pascal layer -- shared types.

  Every proof object in this codebase carries a scope label, exactly as
  in theorems/INDEX.md and docs/METHOD.md:

    psProved         complete mathematical proof (exact Sage / symbolic)
    psChecked        numerical verification, no proof yet
    psStrongEvidence strong evidence only -- NOT admitted as proved
                     (this is where T4* lives)
    psOpen           open item, no claim made
    psAsserted       declared input / dictionary, not derived

  The Pascal objects mirror the index; they do not upgrade any label. }

{$mode objfpc}{$H+}

interface

type
  TProofStatus = (
    psProved,
    psChecked,
    psStrongEvidence,
    psOpen,
    psAsserted
  );

function ProofStatusName(AStatus: TProofStatus): string;

implementation

function ProofStatusName(AStatus: TProofStatus): string;
begin
  case AStatus of
    psProved:         Result := 'PROVED';
    psChecked:        Result := 'CHECKED';
    psStrongEvidence: Result := 'STRONG-EVIDENCE';
    psOpen:           Result := 'OPEN';
    psAsserted:       Result := 'ASSERTED';
  else
    Result := 'UNKNOWN';
  end;
end;

end.
