// Solve the Momentum equation
const volScalarField& mu = cTranModel.mu();

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U) + fvm::div(phi, U)
  + MRF.DDt(rho, U)
  //+ turbulence->divDevRhoReff(U)  // Check Removing turbulence
    - fvm::laplacian(mu, U)
    - fvc::div(mu*dev2(fvc::grad(U)().T()))
 ==
       
    fvOptions(rho, U)  // check Gravity
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}
