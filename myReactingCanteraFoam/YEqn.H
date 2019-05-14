tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    reaction->correct();
    Qdot = reaction->Qdot();
    volScalarField Yt(0.0*Y[0]);
    
    const PtrList<volScalarField>& diffCoeffs = cTranModel.diffCoeffs();
    const PtrList<volScalarField>& diffTCoeffs = cTranModel.diffTCoeffs();

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
//              - fvm::laplacian(turbulence->muEff(), Yi)
                - fvm::laplacian(rho*diffCoeffs[i], Yi) // mass diffusion due to mass fraction gradients
             ==
                reaction->R(Yi)
              + fvOptions(rho, Yi)
            );

            YiEqn += fvc::laplacian(diffTCoeffs[i]/T, T);   // mass diffusion due to temp gradients

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
                        //inert check        
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
