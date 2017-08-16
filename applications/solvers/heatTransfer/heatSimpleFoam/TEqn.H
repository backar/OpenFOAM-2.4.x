{
    volScalarField
        alphaEff("alphaEff", turbulence->nu()/Pr + turbulence->nut()/Prt);

    tmp<fvScalarMatrix> TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
     ==
        fvOptions(T)
    );

    TEqn().relax();

    fvOptions.constrain(TEqn());

    TEqn->solve();

    fvOptions.correct(T);
}
