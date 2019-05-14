{
    volScalarField& he = thermo.he();
    const PtrList<volScalarField>& Hspecs = cTranModel.partialEnthalpies();
    const PtrList<volScalarField>& diffTCoeffs = cTranModel.diffTCoeffs();
    const PtrList<volScalarField>& diffCoeffs = cTranModel.diffCoeffs();
    const volScalarField& kappa = cTranModel.kappa();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
     ==
        Qdot
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);
    
    forAll(Hspecs, i)
    {
       EEqn -= fvc::laplacian(kappa/thermo.Cp()*Hspecs[i], Y[i]); // extra term from writing Fourier's law as a function of sensible enthalpy
       EEqn += fvc::laplacian(Hspecs[i]*rho*diffCoeffs[i], Y[i]); // energy flux due to mass diffusion (Fick's law) 
       EEqn += fvc::laplacian(Hspecs[i]*diffTCoeffs[i]/T, T); // energy flux due to mass diffusion (Soret effect)
    } 
    
    EEqn.relax();
    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
    cTranModel.updateState();
    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
