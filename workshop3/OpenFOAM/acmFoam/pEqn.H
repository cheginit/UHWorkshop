    phi =  fvc::flux(U);
    MRF.makeRelative(phi);
    adjustPhi(phi, U, p);

    tUEqn.clear();

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(p) == -c2*fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            adjustPhi(phi, U, p);
        }
    }

    #include "continuityErrs.H"

    p.relax();
