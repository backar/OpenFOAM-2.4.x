/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "heatSolver.H"

#include "subCycle.H"
#include "singlePhaseTransportModel.H"
#include "incompressible/RAS/RASModel/RASModel.H"
#include "simpleControl.H"
#include "IOdictionary.H"


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatSolver::heatSolver
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    icoSsSolver(mesh, dict),
    Pr_(dict.lookupOrDefault<scalar>("Pr", 1)),
    Prt_(dict.lookupOrDefault<scalar>("Prt", 1)),
    Tpatch_(new volScalarField(fieldIOobject("T",
                                             *subMesh_,
                                             IOobject::READ_IF_PRESENT),
                               *subMesh_))
{
//    init(); // icoSolver::init() is executed twice by this
}


Foam::heatSolver::heatSolver
(
    const fvMesh& mesh
)
:
    icoSsSolver(mesh),
    Pr_(1),
    Prt_(1),
    Tpatch_(new volScalarField(fieldIOobject("T",
                                             *subMesh_,
                                             IOobject::READ_IF_PRESENT),
                               *subMesh_))
{
//    init(); // icoSolver::init() is executed twice by this
}


Foam::heatSolver::heatSolver
(
    const heatSolver& ptf
)
:
    icoSsSolver(ptf),
    Pr_(ptf.Prt_),
    Prt_(ptf.Prt_),
    Tpatch_(ptf.Tpatch_)
{}


Foam::heatSolver::~heatSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::heatSolver::solve()
{
    static unsigned iter_(0);
    if (iter_++%nSubCycles_) return;

    icoSsSolver::solve();

    volScalarField
        alphaEff("alphaEff", turbMdl_->nu()/Pr_ + turbMdl_->nut()/Prt_);

    tmp<fvScalarMatrix> TEqn
    (
        fvm::ddt(Tpatch_())
      + fvm::div(phiPatch_(), Tpatch_())
      - fvm::laplacian(alphaEff, Tpatch_())
      ==
        fvOptions_()(Tpatch_())
    );

    addParaFlux(TEqn());  // Better here than in fvOptions

    TEqn().relax();

    fvOptions_->constrain(TEqn());

    Foam::solve(TEqn());

    fvOptions_->correct(Tpatch_());
}


// ************************************************************************* //
