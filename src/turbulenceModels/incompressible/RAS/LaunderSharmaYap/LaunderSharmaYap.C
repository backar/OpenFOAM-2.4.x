/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "LaunderSharmaYap.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDist.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LaunderSharmaYap, 0);
addToRunTimeSelectionTable(RASModel, LaunderSharmaYap, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> LaunderSharmaYap::fMu() const
{
    return exp(-3.4/sqr(scalar(1) + sqr(k_)/(nu()*epsilonTilda_)/50.0));
}


tmp<volScalarField> LaunderSharmaYap::f2() const
{
    return
        scalar(1)
      - 0.3*exp(-min(sqr(sqr(k_)/(nu()*epsilonTilda_)), scalar(50.0)));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

LaunderSharmaYap::LaunderSharmaYap
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            1.92
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.3
        )
    ),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    epsilonTilda_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateLowReNut("nut", mesh_)
    ),

    fvOptions_
    (
        mesh_,
        *this
    )
{
    bound(k_, kMin_);
    bound(epsilonTilda_, epsilonMin_);

    nut_ = Cmu_*fMu()*sqr(k_)/epsilonTilda_;
    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> LaunderSharmaYap::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> LaunderSharmaYap::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> LaunderSharmaYap::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}


tmp<fvVectorMatrix> LaunderSharmaYap::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}


bool LaunderSharmaYap::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void LaunderSharmaYap::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    tmp<volScalarField> S2 = 2*magSqr(symm(fvc::grad(U_)));
    volScalarField G(GName(), nut_*S2);

    const volScalarField E(2.0*nu()*nut_*fvc::magSqrGradGrad(U_));
    const volScalarField D(2.0*nu()*magSqr(fvc::grad(sqrt(k_))));

    // Implementation of Yap correction, Se
    // Dissipation rate equation
    const volScalarField
        F(pow(k_, 1.5) / \
          (pow(Cmu_, -0.75)*0.41*wallDist(mesh_).y()*epsilonTilda_)),
        Se(0.83*sqr(epsilonTilda_)/k_*(F - 1)*sqr(F));

    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilonTilda_)
      + fvm::div(phi_, epsilonTilda_)
      - fvm::laplacian(DepsilonEff(), epsilonTilda_)
     ==
        C1_*G*epsilonTilda_/k_
      - fvm::Sp(C2_*f2()*epsilonTilda_/k_, epsilonTilda_)
      + E
      + Se
      + fvOptions_(epsilonTilda_)
    );

    epsEqn().relax();
    fvOptions_.constrain(epsEqn());
    solve(epsEqn);
    bound(epsilonTilda_, epsilonMin_);
    fvOptions_.correct(epsilonTilda_);


    // Turbulent kinetic energy equation

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G - fvm::Sp((epsilonTilda_ + D)/k_, k_)
      + fvOptions_(k_)
    );

    kEqn().relax();
    fvOptions_.constrain(kEqn());
    solve(kEqn);
    bound(k_, kMin_);
    fvOptions_.correct(k_);


    // Re-calculate viscosity
    nut_ == Cmu_*fMu()*sqr(k_)/epsilonTilda_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
