/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Jon-Anders Baeckar
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

Class
    Foam::fv::NumericFaceFluxSource

Description
    Numeric face flux source, calculated by a wall sub solver, and described
    using an input dictionary.

    Example of the source specification:

    \verbatim
    <Type>NumericFaceFluxSourceCoeffs
    {
        fieldNames_     (U);                // names of fields to apply source (replace with two below!)
        wallSrcFields   (U, T);             // apply wallSrc to fields from sub field
        fluxSrcFields   (k, epsion);        // apply fluxSrc to fields from sub field
        wallName        'wall';             // name of wall patch
        interfaceName   'oldInternalFaces'; // name of interface patch
        subCycles       1;                  // number of sub solver steps
        pName           'p';                // name of p field
        phiName         'phi';              // name of phi field
        rhoName         'rho';              // name of rho field (not used yet)
        UName           'U';                // name of U field
        mapPressure     yes;
        faceSource      no;
        wallSource      no;
        turbField       no;
        subField        no;
    }
    \endverbatim

SeeAlso
    Foam::fvOption
    Foam::turbulenceModels::subSolver

SourceFiles
    NumericFaceFluxSource.C

\*---------------------------------------------------------------------------*/

#ifndef NumericFaceFluxSource_H
#define NumericFaceFluxSource_H

#include "fvOption.H"
// FIXME: Choose type of subSolver at run time.
//#include "subSolver.H"
#include "heatSolver.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

// Forward declaration of classes

template<class Type>
class NumericFaceFluxSource;


// Forward declaration of friend functions

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const NumericFaceFluxSource<Type>&
);


/*---------------------------------------------------------------------------*\
                     Class NumericFaceFluxSource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class NumericFaceFluxSource
:
    // FIXME: This is NOT a good design ...
    //        NumericFaceFluxSource is templated and several instances exist
    //        at the same time: U(vector), T(scalar), [k, e](list of scalars).
    //        Each of these instances is both an option and a subSolver, the
    //        U(vector) does also contain an icoSolver which is a subSolver.
    //        On top of this, there may be several NumericFaceFluxSource,
    //        and thus, several ico/subSolver connected to separate walls.
    //        To support different subSolver for different walls, many methods
    //        use a static hash container where the has depends on the wall.
    //        Thus, I have implemented singleton methods.
    //
    //        The subSolver should be the same for all template types but
    //        unique per each wall instance, a singleton per wall.
    //        All fields, templates, do only need the base class, subSolver,
    //        except for U(vector) which needs the youngest class.
    //
    //        A.
    //        Another implementation might be a registry of singleton classes.
    //        The registry key would then be the name of the wall, or in the
    //        future a name of a cell set.
    //        U(vector) could use youngest heat/icoSolver and the other fields
    //        could type_cast to subSolver.
    //        However, it is important that U(vector) instantiates it.
    //        Potential problems is what types of constructors is needed. A
    //        copy constructor might give a reference of an existing, but
    //        what to do if only main mesh is given without coeffs?
    //
    //        B.
    //        Neither globals nor singleton should be needed. To communicate
    //        between main and sub an IOobject could be put on the sub mesh...
    //        ... or could it? It is main who should write and sub to read!
    //        The problem I did solve was that sub could NOT read main fields.
    //        Maybe go back and solve this problem for turbulent entities.
    //        Why cannot I reach main k from sub->mainMesh->turb->k()?
    //
    public option, public subSolver
//    public option
{
private:

    // Private data

        //_ Conveniance typedefs
        typedef GeometricField<Type, Foam::fvPatchField, Foam::surfaceMesh>
            surfTypeField;
        typedef GeometricField<Type, Foam::fvPatchField, Foam::volMesh>
            volTypeField;

    // Private Member Functions

        //- Check data
        void checkData() const;

        //- Disallow default bitwise copy construct
        NumericFaceFluxSource(const NumericFaceFluxSource&);

        //- Disallow default bitwise assignment
        void operator=(const NumericFaceFluxSource&);

protected:

    // Protected data

        //- some booleans
        const bool
            mapPressure_, subField_, turbField_, faceSource_, wallSource_;

        //- field names
        const word pName_, phiName_, UName_;

        //- reference to saved equation coefficients container
        HashPtrTable<fvMatrix<Type>, word, Hash<word> >& savedEqnCoeff_;

        //- pointer to sub solver
        heatSolver* sub_;

    // Protected Member Functions

        //- Give address to equation container based on type
        HashPtrTable<fvMatrix<Type>, word, Hash<word> >& eqnContainerAddress();

public:

    //- Runtime type information
    TypeName("NumericFaceFluxSource");


    // Constructors

        //- Construct from components
        NumericFaceFluxSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        ~NumericFaceFluxSource();

    // Member Functions

        // Evaluation

            //- Source face flux term to field equation
            virtual void addSup
            (
                fvMatrix<Type>& eqn,
                const label fieldI
            );

            //- Source face flux term to compressible field equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<Type>& eqn,
                const label fieldI
            );

            //- Remove wall cells from main matrix and add flux from sub
            virtual void setValue(fvMatrix<Type>& eqn, const label fieldI);

            //- Calculate flux from main at para patch and save for later
            virtual void correct(volTypeField& field);

        // I-O

            //- Write the source properties
            virtual void writeData(Ostream&) const;

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NumericFaceFluxSource.C"
#   include "NumericFaceFluxSourceIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
