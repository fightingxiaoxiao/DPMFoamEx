volScalarField kCorrect(
    IOobject
    (
        "kCorrect",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "kCorrect",
        Uc.dimensions(),
        Zero
    )
    );

for(const auto& p : kinematicCloud)
{
    const label celli = p.cell();
    kCorrect[celli] += mag(p.UTurb());
}

kCorrect.correctBoundaryConditions();
