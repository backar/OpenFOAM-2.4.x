/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include <assert.h>

// * * * * * * * * * * * * * Protected Static Data * * * * * * * * * * * * * //

//- Containers of saved equation coefficients (from last iteration)
static Foam::HashPtrTable
<
    Foam::fvMatrix<scalar>,
    Foam::word,
    Foam::Hash<Foam::word>
>
savedScalarEqnCoeff_;
static Foam::HashPtrTable
<
    Foam::fvMatrix<vector>,
    Foam::word,
    Foam::Hash<Foam::word>
>
savedVectorEqnCoeff_;
static Foam::HashPtrTable
<
    Foam::fvMatrix<sphericalTensor>,
    Foam::word,
    Foam::Hash<Foam::word>
>
savedSphericalTensorEqnCoeff_;
static Foam::HashPtrTable
<
    Foam::fvMatrix<symmTensor>,
    Foam::word,
    Foam::Hash<Foam::word>
>
savedSymmTensorEqnCoeff_;
static Foam::HashPtrTable
<
    Foam::fvMatrix<tensor>,
    Foam::word,
    Foam::Hash<Foam::word>
>
savedTensorEqnCoeff_;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
HashPtrTable<fvMatrix<Type>, word, Hash<word> >& eqnContainerRef()
{
    return NULL;
}

template<>
HashPtrTable<fvMatrix<scalar>, word, Hash<word> >& eqnContainerRef()
{
    return savedScalarEqnCoeff_;
}

template<>
HashPtrTable<fvMatrix<vector>, word, Hash<word> >& eqnContainerRef()
{
    return savedVectorEqnCoeff_;
}

template<>
HashPtrTable<fvMatrix<sphericalTensor>, word, Hash<word> >& eqnContainerRef()
{
    return savedSphericalTensorEqnCoeff_;
}

template<>
HashPtrTable<fvMatrix<symmTensor>, word, Hash<word> >& eqnContainerRef()
{
    return savedSymmTensorEqnCoeff_;
}

template<>
HashPtrTable<fvMatrix<tensor>, word, Hash<word> >& eqnContainerRef()
{
    return savedTensorEqnCoeff_;
}

static unsigned iter_(0);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::NumericFaceFluxSource<Type>::NumericFaceFluxSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    subSolver(mesh, coeffs_),
    mapPressure_(coeffs_.lookupOrDefault<bool>("mapPressure", true)),
    subField_(coeffs_.lookupOrDefault<bool>("subField", false)),
    turbField_(coeffs_.lookupOrDefault<bool>("turbField", false)),
    faceSource_(coeffs_.lookupOrDefault<bool>("faceSource", true)),
    wallSource_(coeffs_.lookupOrDefault<bool>("wallSource", false)),
    pName_(coeffs_.lookupOrDefault<word>("pName", "p")),
    phiName_(coeffs_.lookupOrDefault<word>("phiName", "phi")),
    UName_(coeffs_.lookupOrDefault<word>("UName", "U")),
    savedEqnCoeff_(eqnContainerRef<Type>())
{
    if (subField_ or turbField_)
    {
        // These are not used for sub but assert anyway
        assert(faceSource_);
        assert(not wallSource_);
    }
    coeffs_.lookup("fieldNames") >> fieldNames_;
    applied_.setSize(fieldNames_.size(), false);

    bool UinMainFieldNames = false;
    forAll(fieldNames_, i)
        if (not fieldNames_.empty() and fieldNames_[i] == UName_ and not subField_)
            UinMainFieldNames = true;
    if (UinMainFieldNames)
        sub_ = new heatSolver(mesh, coeffs_);
    else
        sub_ = NULL;

    Info<< "    - creating numeric source: " << this->name()
        << ": F(" << faceSource_ << "), W(" << wallSource_
        << "), S(" << subField_ << "), T(" << turbField_
        << ")" << endl;
}


template<class Type>
Foam::fv::NumericFaceFluxSource<Type>::~NumericFaceFluxSource()
{
    if (sub_) delete sub_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fv::NumericFaceFluxSource<Type>::addSup
(
    fvMatrix<Type>& eqn,
    const label fieldI
)
{
    if (debug)
    {
        Info<< "NumericFaceFluxSource<" << pTraits<Type>::typeName
            << ">::addSup for source " << name_ << endl;
    }

    // Called from fvOptions() before equation is solved
    const word key(numPatchName_ + fieldNames_[fieldI]);

    if (not savedEqnCoeff_.found(key)) return;

    if (sub_ and fieldNames_[fieldI] == UName_)
    {
        sub_->solve();
        iter_++;
    }

    if (subField_ and turbField_)
    {
        // Add flux source from main to sub field at interface to main.
        // Main fields are not reachable from here for sub. Thus, set in
        // solver if possible, and for turbulent fields the source has to be
        // set when called from main, see below in correct
        fvMatrix<Type>* mainEqn(savedEqnCoeff_[key]);
        eqn.diag() = mainEqn->diag();
        eqn.source() = mainEqn->source();
    }
    else
    if (iter_>1 and not subField_)
    {
        // Switch contribution what goes into the next cells,
        // from main wall cells to sub top cells
        const volTypeField& subField
            (
               subMesh_->thisDb().lookupObject<volTypeField>
               (
                   fieldNames_[fieldI]
               )
            );
        const vectorField flux(faceFlux(eqn.psi().name()));

        fvMatrix<Type>* mainEqn(savedEqnCoeff_[key]);

        const labelListList& interFaceInfo_(interFaceInfo());
        const fvPatch& interPatch(subMesh_->boundary()[interPatchID()]);
        labelList faceLabels(ortoFaces());
        if (faceSource_) faceLabels.append(paraFaces());
        forAll(faceLabels, indx)
        {
            const label
                sFacI(faceLabels[indx]),
                nextI(interFaceInfo_[3][sFacI]),
                wOwnr(interFaceInfo_[4][sFacI]),
                direc(interFaceInfo_[5][sFacI]),
                sCelI(interPatch.faceCells()[sFacI]);

            // Add contribution from sub top cell
            const scalar
                Diff(flux[sFacI].x()),
                uAdv(flux[sFacI].y()),
                lAdv(flux[sFacI].z());

            // Remove original diag contribution from wall cell,
            // using the one calculated in setValue from PREVIOUS
            // time step. The corresponding source contribution is
            // zeroed in setValue.
            /*
             * TEST PROCEDURE:
             *    1. Small mesh, one sub layer, but big enough to contain
             *        neg/pos phi for own/nei cells.
             *    2. Use uncorrected laplacian and comment out the div term
             *        in laminar.C.
             *    3. Hardcode mainEqn->lower/upper for a few iterations
             *        to NOT lag an iteration. This is useful to test
             *        diag, A, and H, and maybe dP and phi.
             *        OBSERVE: To find these correct do not run addSup.
             *    4. U can only be tested at convergence as solving
             *        U explicitly or implicitly will differ even with inner
             *        iterations set to only one. Best is to check Cf and Nu.
             */
            scalar dLU(0), sLU(0);
            // UD: INVESTIGATE FURTHER, DOES NOT GIVE SAME GOOD RESULTS AS CD
            /*
            if (wOwnr)
            {
            dLU = Diff + lAdv;
            sLU = Diff + uAdv;
            eqn.diag()[nextI] = -dLU;
                eqn.diag()[nextI] -= mainEqn->upper()[sFacI];
            }
            else
            {
            dLU = Diff + uAdv;
            sLU = Diff - lAdv;
            eqn.diag()[nextI] = -dLU;
                eqn.diag()[nextI] -= mainEqn->lower()[sFacI];
            }
            */
            // BEST for IMP and UPP for CD (using neg/pos and own/nei)
            dLU = Diff + direc*lAdv;
            sLU = Diff + direc*uAdv;
            eqn.diag()[nextI] = -dLU;
            if (wOwnr)
            {
                eqn.diag()[nextI] -= mainEqn->upper()[sFacI];
            }
            else
            {
                eqn.diag()[nextI] -= mainEqn->lower()[sFacI];
            }
            eqn.source()[nextI] = -sLU * subField[sCelI];
        }

        // FIXME: Assert not Dirichlet condition
        if (wallSource_)
        {
          // Change homogeneous Dirichlet to slip condition
          // This gives better results than applying OF slip condition?! WHY?
          // Dirichlet (U) internalCoeffs looks like (nu*A/y, nu*A/y, nu*A/y)
          // while slip looks like (nu*A/y, 0, 0) for a wall orthogonal with x
          // which makes it possible to set instead of project on wall
            /*
            fvMatrix<Type>* mainEqn(savedEqnCoeff_[key]);
            const Field<Type>&
                internalCoeffs(mainEqn->internalCoeffs()[wallPatchID()]);
            eqn.internalCoeffs()[wallPatchID_] =
                projectOnWall(internalCoeffs, wallPatchID_);
            */

            const fvPatch& wallPatch(mainMesh_.boundary()[wallPatchID_]);

            // Add source calculated in sub solver as wall shear stress
            // and compensate Dirichlet condition with an explicit tau
            const Field<Type>
                tauA(wallSource(fieldNames_[fieldI], pTraits<Type>::zero));
            forAll(wallPatch, mFacI)
                eqn.source()[wallPatch.faceCells()[mFacI]] = tauA[mFacI];
        }
    }
}


template<class Type>
void Foam::fv::NumericFaceFluxSource<Type>::addSup
(
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const label fieldI
)
{
    if (debug)
    {
        Info<< "NumericFaceFluxSource<" << pTraits<Type>::typeName
            << ">::addSup for source " << name_ << endl;
    }

    return this->addSup(eqn, fieldI);
}


template<class Type>
void Foam::fv::NumericFaceFluxSource<Type>::setValue
(
    fvMatrix<Type>& eqn,
    const label fieldI
)
{
    if (debug)
    {
        Info<< "NumericFaceFluxSource<" << pTraits<Type>::typeName
            << ">::setValue for source " << name_ << endl;
    }

    // Called from fvOptions.constrain before equation is solved
    const word key(numPatchName_ + fieldNames_[fieldI]);

    const volTypeField&
        subField
        (
            subMesh_->thisDb().lookupObject<volTypeField>(fieldNames_[fieldI])
        );

    if (not subField_)
    {
        const fvPatch &interPatch(subMesh_->boundary()[interPatchID()]);

        if (fieldNames_[fieldI] == pName_)
        {
            plotSrcP(eqn.source());
            return;
        }
        // Init container for flux source main->sub
        if (not savedEqnCoeff_.found(key))
        {
            // Creating a new matrix triggers updateCoeffs, for omega and eps
            // wall functions this gives a FOAM FATAL ERROR as it fails to
            // find GName in the database of the sub-region.
            // FIXME: For now I check in the wall function if G exists, better
            //        to fix something from here if possible.
            savedEqnCoeff_.insert
                (
                    key,
                    new fvMatrix<Type>(subField, eqn.dimensions())
                );
            if (eqn.asymmetric() and eqn.hasLower())
                savedEqnCoeff_[key]->lower(interPatch.size());
            if (eqn.hasUpper())
                savedEqnCoeff_[key]->upper(interPatch.size());

            // Save away Dirichlet boundary condition
            savedEqnCoeff_[key]->internalCoeffs()[wallPatchID()] =
                eqn.internalCoeffs()[wallPatchID_];
            savedEqnCoeff_[key]->boundaryCoeffs()[wallPatchID()] =
                eqn.boundaryCoeffs()[wallPatchID_];
        }
        fvMatrix<Type>* mainEqn(savedEqnCoeff_[key]);

        const labelListList& interFaceInfo_(interFaceInfo());
        if (faceSource_)
        {
            // 1. Save wall coefficient contributing to next from wall
            // 2. Zero out coeffficient contributing to next from wall
            //    which is used to add to source in setValues
            //    (actually src -= coeff*val and diag=D, src=D*val)

            forAll(interPatch, sFacI)
            {
                const label
                    wOwnr(interFaceInfo_[4][sFacI]),
                    mFacI(interFaceInfo_[0][sFacI]);
                if (eqn.symmetric())
                {
                    if (eqn.hasUpper())
                    {
                        mainEqn->upper()[sFacI] = eqn.upper()[mFacI];
                        if (iter_>1)
                        {
                            eqn.upper()[mFacI] = 0;
                        }
                    }
                    if (eqn.hasLower())
                    {
                        mainEqn->upper()[sFacI] = eqn.lower()[mFacI];
                        if (iter_>1)
                        {
                            eqn.lower()[mFacI] = 0;
                        }
                    }
                }
                else
                {
            // FIXME: Test to only reduce part of contribution from wall cell
            //        and add the rest from sub (average between sub and wall)
            //        Motive: get better pressure prediction
                    mainEqn->upper()[sFacI] = eqn.upper()[mFacI];
                    mainEqn->lower()[sFacI] = eqn.lower()[mFacI];
                    if (iter_>1)
                    {
                        if (wOwnr) eqn.lower()[mFacI] = 0;
                        else eqn.upper()[mFacI] = 0;
                    }
                }
            }

            // Remove wall cells from matrix and set values from sub
            if (fieldNames_[fieldI] != UName_)
            {
                // FIXME: Possibly make mainCells a static hash of lists
                labelList mainCells(interPatch.size(), -1);
                forAll(interPatch, sFacI)
                {
                    mainCells[sFacI] = interFaceInfo_[2][sFacI];
                }
                eqn.setValues(mainCells, averageField(subField));
            }
        }
    }
}


template<class Type>
void Foam::fv::NumericFaceFluxSource<Type>::correct(volTypeField& field)
{
    if (debug)
    {
        Info<< "NumericFaceFluxSource<" << pTraits<Type>::typeName
            << ">::correct for source " << name_ << endl;
    }

    // Called from fvOptions.correct after equation is solved
    const word key(numPatchName_ + field.name());
    // FIXME: As today this method is called twice for U with the motive of
    //        reseting the wall values before U is solved for main.
    //        A less intrusive way would be to reset the wall values before
    //        main U matrix is assembled ... this will be equally intrusive,
    //        as fvOptions(U) must be called in the beginning of the assemble.
    if (not subField_ and field.name() == UName_)
    {
        // Set cell and wall face value in wall main cell to evaluate
        // interpolation correctly, needed by e.g. fvc::grad()
        // No improvement seen if doing the same for turbulent fields
        // not even with Spalart-Allmaras or SST which really need it
        static HashTable<unsigned, word, string::hash> hashIter;
        static HashTable<Field<Type>, word, string::hash> hashU, hashUw;
        fvPatchField<Type>& wField = field.boundaryField()[wallPatchID_];
        if (not hashU.found(numPatchName_))
        {
            hashIter.insert(numPatchName_, 0);
            Field<Type> cacheU(wField);
            hashUw.insert(numPatchName_, cacheU);
            hashU.insert(numPatchName_, cacheU);
        }
        if (hashIter[numPatchName_] < iter_)
        {
            hashIter[numPatchName_] = iter_;
            forAll(wField, wFacI)
            {
                hashUw[numPatchName_][wFacI] = wField[wFacI];
                hashU[numPatchName_][wFacI] =
                    field[wField.patch().faceCells()[wFacI]];
            }
            fixInterpolation(field);
        }
        else
        {
            Info << "Reset U in " << numPatchName_
                << " wall cells from cache." << endl;
            forAll(wField, wFacI)
            {
                wField[wFacI] = hashUw[numPatchName_][wFacI];
                field[wField.patch().faceCells()[wFacI]] =
                    hashU[numPatchName_][wFacI];
            }
        }
    }
    else
    if (not subField_ and turbField_ and savedEqnCoeff_.found(key))
    {
        // Add flux source from main solution and save in container, these
        // could equally be set above in setValue as only field has changed
        fvMatrix<Type>* mainEqn(savedEqnCoeff_[key]);

        const polyPatch&
            interPatch(subMesh_->boundary()[interPatchID()].patch());
        forAll(interPatch, sFacI)
        {
            // Zero first as corner cells get multiple contributions
            const label sCelI(interPatch.faceCells()[sFacI]);
            mainEqn->diag()[sCelI] = 0;
            mainEqn->source()[sCelI] = pTraits<Type>::zero;
        }

        // Add diffusion source from sub solution
        const vectorField flux(faceFlux(field.name()));
        const labelListList& interFaceInfo_(interFaceInfo());
        forAll(interPatch, sFacI)
        {
            const label
                nextI(interFaceInfo_[3][sFacI]),
                direc(interFaceInfo_[5][sFacI]),
                sCelI(interPatch.faceCells()[sFacI]);

            const scalar
                Diff(flux[sFacI].x()),
                uAdv(flux[sFacI].y()),
                lAdv(flux[sFacI].z());

            const scalar
                dLU(Diff + direc*uAdv),
                sLU(Diff + direc*lAdv);

            // FIXME: Mind corners get multiple input, test below
            mainEqn->diag()[sCelI] -= dLU;
            mainEqn->source()[sCelI] -= sLU * field[nextI];
        }
    }
}


// ************************************************************************* //
