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
    Foam::subSolver

Group
    grpWallBoundaryConditions

Description
    Calculate shear stress from sub grid solution.

SourceFiles
    subSolver.C

\*---------------------------------------------------------------------------*/

#ifndef subSolver_H
#define subSolver_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;
class IOobject;

/*---------------------------------------------------------------------------*\
            Class subSolver Declaration
\*---------------------------------------------------------------------------*/

class subSolver
{


protected:

    // Protected data

        //- Factor to be used in DoE
        const scalar fct_;

        //- Name of patch connected to the sub-grid
        const word numPatchName_;

        //- Name of interface patch
        const word intrfcPatchName_;

        //- Label of wall patch
        const label wallPatchID_;

        //- Reference to main mesh
        const fvMesh& mainMesh_;

        //- Reference to sub mesh (the same also if called from subMesh)
        fvMesh* subMesh_;
        static HashTable<fvMesh*, word, string::hash> subMeshes_;

    // Protected methods

        //- Initiate all fields not set in constructor
        virtual void init();

        //- An IOobject to initiate sub mesh connected to wall patch
        // mostly a conveniance method to overcome private constructor
        IOobject meshIOobject(const fvMesh&);

        //- An IOobject to initiate sub fields connected to wall patch
        IOobject fieldIOobject
        (
            const word,
            const fvMesh&,
            const IOobject::readOption rOpt=IOobject::MUST_READ
        );

        label interPatchID() const;

        label wallPatchID() const;

        //- Labels of para sub boundary faces
        const labelList paraFaces() const;

        //- Labels of orto sub boundary faces
        const labelList ortoFaces() const;

        //- Labels of corner faces
        const labelPairList cornerFaces() const;

        //- Labels of main faces, owners and neighbours at and around interface
        const labelListList interFaceInfo() const;

        //- Labels of main faces along wall patch
        const labelList wallFacesToMainFaces() const;

        //- Labels of sub cells overlapping main wall cell
        const labelListList mainWallFaceToSubCells() const;

        //- Labels of internal faces as pair; (face, cell) per sub wall faces
        const List<labelPairList> wallFaceToSubCellFaces() const;

        //- Labels of overlapping faces as pair; (sub, main) per
        //    boundary / internal index and per para face
        const List<List<labelPairList> > mainFaceToSubFaces() const;

        // FIXME: Transfer to mainOrtoToSubCells
        //- Labels of sub faces partly overlapping main wall faces
        const GeometricField<scalar, fvsPatchField, surfaceMesh>
            faceLabelField();

        //- Cell center ratio at para patch, used for linear interpolation
        const scalarField weights() const;

        //- advection and diffusion face flux over interPatch
        tmp<vectorField> faceFlux(const word) const;

        //- Set value on interface of sub field from main and sub field
        template<typename Type>
        void interpolFaceValue
        (
            GeometricField<Type, fvPatchField, volMesh>&,
            const bool onlyOrto = false,
            const bool averageCorner = false
        );

        //- Set value on interface of sub field from main
        template<typename Type>
        void setInterfaceValue
        (
            GeometricField<Type, fvPatchField, volMesh>&,
            const bool onlyOrto = false,
            const bool averageCorner = false
        );

        //- Reset value on interface of sub field
        template<typename Type>
        void resetInterfaceValue(GeometricField<Type, fvPatchField,volMesh>&);

        //- Zero value on interface of sub field
        template<typename Type>
        void zeroInterfaceValue(GeometricField<Type, fvPatchField, volMesh>&);

        //- Add flux source at patch assuming orthogonal mesh
        //  around the patch
        template<typename Type>
        void addParaFlux(fvMatrix<Type>& eqn);

        //- Volume average sub field onto main (only wall cells)
        template<class Type>
        tmp<Field<Type> > averageField(const Field<Type>& subField);

        //- Plot source of p-equation at interface
        void plotSrcP(Field<scalar>&);
        template<typename Type>
        void plotSrcP(Field<Type>&);

        //- Remove normal component of field towards patch with label
        //- all types must be supported through the templated version
        //- Sub field version looping over face cells
        tmp<Field<vector> > projectOnPatch(const Field<vector>&, const label);
        template<typename Type>
        tmp<Field<Type> > projectOnPatch(const Field<Type>&, const label);
        //-
        //- Main field version looping over wall faces
        tmp<Field<vector> > projectOnWall(const Field<vector>&, const label);
        template<typename Type>
        tmp<Field<Type> > projectOnWall(const Field<Type>&, const label);

        //- Set wall values in main to get correct face values at interface
        //-  needed by production of k (magSqr(symm(fvc::grad(U))))
        template<class Type>
        tmp<Field<Type> > correctFaceValue
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );
        template<typename Type>
        void fixInterpolation
        (
            GeometricField<Type, fvPatchField, volMesh>&,
            const bool onlyPositive = false
        );

public:

//    //- Runtime type information
//    TypeName("subSolver");


    // Constructors

        //- Construct from mesh and dictionary
        subSolver
        (
            const fvMesh&,
            const dictionary&
        );

        //- Construct from mesh
        subSolver
        (
            const fvMesh&
        );

        //- Construct as copy
        subSolver
        (
            const subSolver&
        );

        //- Destruct - declaration must be in same file which declare autoPtr
        virtual ~subSolver();


    // Member functions

        //- Wall source difference between main and sub
        template<class Type>
        tmp<Field<Type> > wallSource(const word, const Type, const bool=false);

        //- FIXME: make this pure virtual =0 when merged fvO:numShear / numSrc
        //virtual void update() = 0;
        virtual void update() {};

        //- FIXME: make this pure virtual =0 when merged fvO:numShear / numSrc
        //- Substep solver
        //virtual void solve() = 0;
        virtual void solve() {};
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "subSolverTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
