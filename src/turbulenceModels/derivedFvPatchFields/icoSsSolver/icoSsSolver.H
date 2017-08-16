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

Class
    Foam::icoSsSolver

Group
    grpWallBoundaryConditions

Description
    Calculate shear stress from sub grid solution.

SourceFiles
    icoSsSolver.C

\*---------------------------------------------------------------------------*/

#ifndef icoSsSolver_H
#define icoSsSolver_H

#include "subSolver.H"
#include "HashPtrTable.H"
#include "fvOptionList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class singlePhaseTransportModel;
namespace incompressible
{
class RASModel;
}

/*---------------------------------------------------------------------------*\
            Class icoSsSolver Declaration
\*---------------------------------------------------------------------------*/

class icoSsSolver
:
    public Foam::subSolver
{

protected:

    // Protected data

        //- Boolean to map or solve for pressure
        const bool mapP_, noAdv_;

        //- Relax gradient pressure mapping
        const scalar relaxGradP_;

        //- Number of sub cycles for sub solver
        const label nSubCycles_;

        //- Name of flux field (default = phi)
        const word phiName_;

        //- Sub fields
        autoPtr<volScalarField> pPatch_;
        autoPtr<volVectorField> Upatch_;
        autoPtr<surfaceScalarField> phiPatch_;
        // FIXME: Only for UMIST
        autoPtr<volScalarField> UxMain_;
        autoPtr<volScalarField> UxPatch_;

        //- Pointer to sub transport/thermo model
        autoPtr<singlePhaseTransportModel> transMdl_;

        //- Pointer to sub turbulence model
        autoPtr<incompressible::RASModel> turbMdl_;

        //- fvOption and needed IOdictionary
        autoPtr<IOdictionary> fvOptionsIO_;
        autoPtr<fv::optionList> fvOptions_;


    // Protected methods

        //- Initiate all fields not set in constructor
        virtual void init();

        //- Interpolate pressure from main to sub and set interface gradient
        virtual void update();

        //- Update and scale phi from phi main to ensure continuity
        void updatePhi();

        //- Map and interpolate main pressure on sub
        tmp<volVectorField> interpolateGradP();

        //- Solve streamwise component of velocity, return internal field
        tmp<vectorField> solveStreamwiseCmpt();

public:

//    //- Runtime type information
//    TypeName("icoSsSolver");


    // Constructors

        //- Construct from mesh and dictionary
        icoSsSolver
        (
            const fvMesh&,
            const dictionary&
        );

        //- Construct from mesh
        icoSsSolver
        (
            const fvMesh&
        );

        //- Construct as copy
        icoSsSolver
        (
            const icoSsSolver&
        );


        //- Destructor
        ~icoSsSolver();


    // Member functions

        //- Substep solver
        void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
