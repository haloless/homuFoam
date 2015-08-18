    
    
    // The interFoam solver seems not to have laminar viscosity,
    // so the following code is added based on porousInterFoam
    
    // face effective viscosity (laminar + turbulent)
    surfaceScalarField muEff
    (
        "muEff",
        mixture.muf()
        + fvc::interpolate(rho * turbulence->nut())
    );
    // cell viscosity
    volScalarField mu(mixture.mu());
        
    
    fvVectorMatrix UEqn
    (
        /*
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::laplacian(muEff, U) //
      // + turbulence->divDevRhoReff(rho, U)
      - (fvc::grad(U) & fvc::grad(muEff))
        */
        
        1.0/varepsilon * fvm::ddt(rho, U)
        + 1.0/varepsilon * fvm::div(rhoPhi/varepsilonf, U)
        - fvm::laplacian(muEff/varepsilonf, U)
        - 1.0/varepsilon * (fvc::grad(U) & fvc::grad(muEff))
        
        // Ergun drag
        //+ ErgunCoef * U
     ==
        fvOptions(rho, U)
    );
    
    UEqn.relax();
    
    //
    if (porousDragForceActive) 
    {
        volScalarField ErgunCoef = 
        (
            150.0/(dpEff*dpEff) * (1.0-varepsilon)/varepsilon * mu
            + 1.75/dpEff * rho / varepsilon * mag(U)
        ) * (1.0-varepsilon) / pow(varepsilon,2);
        
        // explicit drag
        //UEqn -= fvm::SuSp(-ErgunCoef,U); // this also works
        //UEqn -= fvm::Su(-ErgunCoef*U, U);
        
        // implicit drag
        UEqn += fvm::Sp(ErgunCoef, U);
    }

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
