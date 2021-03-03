{
    volVectorField HbyA(constrainHbyA(rAUc*UcEqn.H(), Uc, p));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
           fvc::flux(HbyA)
         + alphacf*rAUcf*fvc::ddtCorr(Uc, phic)
        )
    );


    adjustPhi(phiHbyA, Uc, p_rgh);
    

    phiHbyA += phicForces;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, Uc, phiHbyA, rAUcf);

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(alphacf*rAUcf, p_rgh)
         ==
            fvc::ddt(alphac) + fvc::div(alphacf*phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, pRefValue);

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phic = phiHbyA - p_rghEqn.flux()/alphacf;

            p_rghEqn.relax();

            Uc = HbyA
              + rAUc*fvc::reconstruct((phicForces - p_rghEqn.flux()/alphacf)/rAUcf);

            Uc.correctBoundaryConditions();

            fvOptions.correct(Uc);
        }
    }
    #include "continuityErrs.H"

    p == p_rgh + 1.0*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - 1.0*gh;
    }
}
