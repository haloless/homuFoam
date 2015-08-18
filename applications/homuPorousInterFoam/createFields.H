    
    // porosity
    Info << "Reading field varepsilon\n" << endl;
    // cell porosity
    volScalarField varepsilon
    (
        IOobject
        (
            "varepsilon",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    // porosity dictionary
    Info<< "Reading homuPorosityProperties\n" << endl;
    IOdictionary porosityProperties
    (
        IOobject
        (
            "homuPorosityProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    // MAX/MIN porosity
    dimensionedScalar varepsilonMin
    (
        "varepsilonMin",
        varepsilon.dimensions(),
        porosityProperties.lookup("varepsilonMin")
    );
    dimensionedScalar varepsilonMax
    (
        "varepsilonMax",
        varepsilon.dimensions(),
        porosityProperties.lookup("varepsilonMax")
    );
    // clamp the porosity
    forAll(varepsilon, icell)
    {
        if (varepsilon[icell]<varepsilonMin.value())
            varepsilon[icell] = varepsilonMin.value();
        if (varepsilon[icell]>varepsilonMax.value())
            varepsilon[icell] = varepsilonMax.value();
    }
    
    //
    dimensionedScalar dpEff
    (
        "dpEff",
        dimLength,
        porosityProperties.lookup("dpEff")
    );
    
    bool porousDragForceActive = true;
    porosityProperties.lookup("dragForceActive") >> porousDragForceActive;
    
    // face porosity
    surfaceScalarField varepsilonf
    (
        "varepsilonf",
        fvc::interpolate(varepsilon)
    );
    // solid volume fraction
    volScalarField varepsilonp = 1.0 - varepsilon;

    
    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"


    Info<< "Reading transportProperties\n" << endl;
    immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

    volScalarField& alpha1(mixture.alpha1());
    volScalarField& alpha2(mixture.alpha2());

    const dimensionedScalar& rho1 = mixture.rho1();
    const dimensionedScalar& rho2 = mixture.rho2();


    // Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        alpha1*rho1 + alpha2*rho2,
        alpha1.boundaryField().types()
    );
    rho.oldTime();


    // Mass flux
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rhoPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho)*phi
    );
    
    // intrinsic mass flux
    surfaceScalarField rhoPhi_varepsilon
    (
        IOobject
        (
            "rhoPhi_varepsilon",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoPhi / varepsilonf
    );


    // Construct incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, mixture)
    );

    #include "readGravitationalAcceleration.H"

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    Info<< "Setting full pressure p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*gh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        pimple.dict(),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }


    fv::IOoptionList fvOptions(mesh);


    // MULES flux from previous time-step
    surfaceScalarField phiAlpha
    (
        IOobject
        (
            "phiAlpha",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        phi*fvc::interpolate(alpha1)
    );

    // MULES Correction
    tmp<surfaceScalarField> tphiAlphaCorr0;
