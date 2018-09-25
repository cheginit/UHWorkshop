MRF.makeRelative(phi);
adjustPhi(phi, U, p);

tUEqn.clear();

// Non-orthogonal pressure corrector loop
while (simple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::ddt(p) == -c2*fvc::div(U)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve();
    adjustPhi(phi, U, p);
}


#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U.correctBoundaryConditions();
fvOptions.correct(U);
