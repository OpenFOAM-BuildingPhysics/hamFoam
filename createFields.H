Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info<< "Reading field pc\n" << endl;
volScalarField pc
(
    IOobject
    (
        "pc",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField ws
(
    IOobject
    (
        "ws",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Ts
(
    IOobject
    (
        "Ts",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
); 

volScalarField Crel
(
    IOobject
    (
        "Crel",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Crel", dimensionSet(0, -2, 2, 0, 0, 0, 0), scalar(1E-3)),
    ("zeroGradient")
);	

volScalarField Krel
(
    IOobject
    (
        "Krel",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Krel", dimensionSet(0, 0, 1, 0, 0, 0, 0), scalar(1E-3)),
    ("zeroGradient")
);

volScalarField K_v
(
    IOobject
    (
        "K_v",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("K_v", dimensionSet(0, 0, 1, 0, 0, 0, 0), scalar(1E-3)),
    ("zeroGradient")
);

volScalarField K_pt
(
    IOobject
    (
        "K_pt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("K_pt", dimensionSet(1, -1, -1, -1, 0, 0, 0), scalar(1E-3)),
    ("zeroGradient")
);    

volScalarField rho_m
(
    IOobject
    (
        "rho_m",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rho_m", dimensionSet(1,-3,0,0,0,0,0), scalar(1E-3)),
    ("zeroGradient")
); 

volScalarField cap_m
(
    IOobject
    (
        "cap_m",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("cap_m", dimensionSet(0,2,-2,-1,0,0,0), scalar(1E-3)),
    ("zeroGradient")
);     

volScalarField lambda_m
(
    IOobject
    (
        "lambda_m",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("lambda_m", dimensionSet(1,1,-3,-1,0,0,0), scalar(1E-3)),
    ("zeroGradient")
);             	

uniformDimensionedVectorField g
(
    IOobject
    (
        "g",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
); 

