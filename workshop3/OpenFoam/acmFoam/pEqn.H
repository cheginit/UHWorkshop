// Non-orthogonal pressure corrector loop
while (simple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::ddt(p) == -c2*fvc::div(U)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve();
}

// Explicitly relax pressure for momentum corrector
p.relax();
