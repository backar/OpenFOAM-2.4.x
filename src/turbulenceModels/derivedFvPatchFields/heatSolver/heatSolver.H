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
    Foam::heatSolver

Group
    grpWallBoundaryConditions

Description
    Calculate shear stress a temperature from sub grid solution.

SourceFiles
    heatSolver.C

\*---------------------------------------------------------------------------*/

#ifndef heatSolver_H
#define heatSolver_H

#include "icoSsSolver.H"
#include "HashPtrTable.H"
#include "fvOptionList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class heatSolver Declaration
\*---------------------------------------------------------------------------*/

class heatSolver
:
    public Foam::icoSsSolver
{

protected:

    // Protected data

        //- Sub constants
        dimensionedScalar Pr_;
        dimensionedScalar Prt_;

        //- Sub fields
        autoPtr<volScalarField> Tpatch_;

public:

//    //- Runtime type information
//    TypeName("heatSolver");


    // Constructors

        //- Construct from mesh and dictionary
        heatSolver
        (
            const fvMesh&,
            const dictionary&
        );

        //- Construct from mesh
        heatSolver
        (
            const fvMesh&
        );

        //- Construct as copy
        heatSolver
        (
            const heatSolver&
        );


        //- Destructor
        ~heatSolver();


    // Member functions

        //- Substep solver
        void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
