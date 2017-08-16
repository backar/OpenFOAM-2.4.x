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

#include <assert.h>
#include "fvPatchField.H"
#include "incompressible/RAS/RASModel/RASModel.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


template<typename Type>
void Foam::subSolver::interpolFaceValue
(
    GeometricField<Type, fvPatchField, volMesh>& field,
    const bool onlyOrto,
    const bool averageCorner
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> volTypeField;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> surfTypeField;

    const volTypeField&
        mainField(mainMesh_.thisDb().lookupObject<volTypeField>(field.name()));

    fvPatchField<Type>& subInterface =
        refCast<fvPatchField<Type> >
        (
            field.boundaryField()[interPatchID()]
        );
    const vectorField
        mC(mainMesh_.C()),
        sC(subMesh_->C()),
        Cf(subInterface.patch().Cf());
    const labelListList& interFaceInfo_(interFaceInfo());

    const labelList& faceCells(subInterface.patch().faceCells());
    labelList faceLabels(ortoFaces());
    if (not onlyOrto) faceLabels.append(paraFaces());
    forAll(faceLabels, indx)
    {
        const label
            sFacI(faceLabels[indx]),
            sCelI(faceCells[sFacI]),
            nextI(interFaceInfo_[3][sFacI]);
        const scalar
            dNf(mag(Cf[sFacI] - mC[nextI])),
            wS(dNf/mag(sC[sCelI] - mC[nextI]));

        // Calculate face value from main and sub
        subInterface[sFacI] = (1-wS)*mainField[nextI] + wS*field[sCelI];
    }
}


template<typename Type>
void Foam::subSolver::setInterfaceValue
(
    GeometricField<Type, fvPatchField, volMesh>& field,
    const bool onlyOrto,
    const bool averageCorner
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> volTypeField;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> surfTypeField;

    const volTypeField&
        mainField(mainMesh_.thisDb().lookupObject<volTypeField>(field.name()));
    const surfTypeField mainValue(fvc::interpolate(mainField));

    fvPatchField<Type>& subInterface =
        refCast<fvPatchField<Type> >
        (
            field.boundaryField()[interPatchID()]
        );
    const labelListList& interFaceInfo_(interFaceInfo());
    labelList faceLabels(ortoFaces());
    if (not onlyOrto) faceLabels.append(paraFaces());
    forAll(faceLabels, indx)
    {
        const label
            sFacI(faceLabels[indx]),
            mFacI(interFaceInfo_[0][sFacI]);
        subInterface[sFacI] = mainValue[mFacI];
    }

    // Update corner faces with an area averaged mean value (when solving p)
    if (averageCorner)
    {
        const scalarField faceArea(subInterface.patch().magSf());
        forAll(cornerFaces(), indx)
        {
            const scalar
                A0(faceArea[cornerFaces()[indx][0]]),
                A1(faceArea[cornerFaces()[indx][1]]);
            const Type
                v0(subInterface[cornerFaces()[indx][0]]),
                v1(subInterface[cornerFaces()[indx][1]]),
                vMean((v0*A0 + v1*A1)/(A0 + A1));
            subInterface[cornerFaces()[indx][0]] = vMean;
            subInterface[cornerFaces()[indx][1]] = vMean;
        }
    }
}


template<typename Type>
void Foam::subSolver::resetInterfaceValue
(
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    fvPatchField<Type>& subInterface =
        refCast<fvPatchField<Type> >
        (
            field.boundaryField()[interPatchID()]
        );
    const polyPatch &interPatch(subMesh_->boundary()[interPatchID()].patch());
    const labelList &faceCells(interPatch.faceCells());
    forAll(interPatch, sFacI)
    {
        subInterface[sFacI] = field[faceCells[sFacI]];
    }
}


template<typename Type>
void Foam::subSolver::zeroInterfaceValue
(
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    fvPatchField<Type>& subInterface =
        refCast<fvPatchField<Type> >
        (
            field.boundaryField()[interPatchID()]
        );
    forAll(subInterface, sFacI)
    {
        subInterface[sFacI] = pTraits<Type>::zero;
    }
}


template<typename Type>
void Foam::subSolver::addParaFlux(Foam::fvMatrix<Type>& eqn)
{
    const GeometricField<Type, Foam::fvPatchField, Foam::volMesh>&
        main
        (
            mainMesh_.thisDb().lookupObject
            <
                GeometricField<Type, Foam::fvPatchField, Foam::volMesh>
            >(eqn.psi().name())
        );
    Foam::Field<scalar> &diag(eqn.diag());
    Foam::Field<Type> &source(eqn.source());
    const vectorField flux(faceFlux(eqn.psi().name()));
    const labelListList& interFaceInfo_(interFaceInfo());

    const fvPatch& interPatch(subMesh_->boundary()[interPatchID()]);
    const labelList &faceCells(interPatch.faceCells());
    forAll(interPatch, sFacI)
    {
        const label
            nextI(interFaceInfo_[3][sFacI]),
//            wOwnr(interFaceInfo_[4][sFacI]),
            direc(interFaceInfo_[5][sFacI]),
            sCelI(faceCells[sFacI]);

        const scalar  // Scale source with volume ratio
            Diff(flux[sFacI].x()),
            uAdv(flux[sFacI].y()),
            lAdv(flux[sFacI].z());

        scalar dLU(0), sLU(0);
//        WORKS FOR CD both imp and upp considering pos/neg phi for own/nei
        dLU = Diff + direc*uAdv;
        sLU = Diff + direc*lAdv;
//        WORKS FOR UD both imp and upp considering pos/neg phi for own/nei
/*
        if (wOwnr)
        {
              // INVESTIGATE FURTHER, DOES NOT GIVE SAME GOOD RESULTS AS CD
            dLU = Diff + uAdv;
            sLU = Diff + lAdv;
        }
        else
        {
            dLU = Diff - lAdv;
            sLU = Diff + uAdv;
        }
*/

        // FIXME: Mind corners get multiple input, check
        diag[sCelI] += dLU;
        source[sCelI] += sLU * main[nextI];
    }
}


template<typename Type>
Foam::tmp<Foam::Field<Type> > Foam::subSolver::wallSource
(
    const word fieldName,
    const Type zero,
    const bool compressible
)
{
    typedef GeometricField<Type,Foam::fvPatchField,Foam::volMesh> volTypeField;
    const volTypeField
        &subField(subMesh_->thisDb().lookupObject<volTypeField>(fieldName));

    const incompressible::RASModel&
        turbulence
        (
            subMesh_->thisDb().lookupObject<incompressible::RASModel>
            (
                "turbulenceModel"
            )
        );

    const fvPatch
        &numPatch(mainMesh_.boundary()[wallPatchID_]),
        &patch(subField.boundaryField()[wallPatchID()].patch());
    const labelList
        &wallFaces(wallFacesToMainFaces()),
        &faceCells(patch.faceCells());
    const scalarField
        &subRY(patch.deltaCoeffs()),
        &subFA(patch.magSf()),
        subNu(turbulence.nu()->boundaryField()[wallPatchID()]);

    Foam::tmp<Foam::Field<Type> > tTauA
        (
            new Field<Type>(numPatch.size(), pTraits<Type>::zero)
        );
    Foam::Field<Type> &tauA(tTauA());

    forAll(patch, sFacI)
    {
        const label
            mFacI(wallFaces[sFacI]),
            sCelI(faceCells[sFacI]);
        tauA[mFacI] = subNu[sFacI]*subFA[sFacI]
            * subField[sCelI]*subRY[sFacI];
    }
    return projectOnPatch(tauA, wallPatchID());
//    return tTauA;
}


template<typename Type>
Foam::tmp<Foam::Field<Type> > Foam::subSolver::averageField
(
    const Foam::Field<Type>& subField
)
{
    const fvPatch& interPatch(subMesh_->boundary()[interPatchID()]);
    Foam::tmp<Foam::Field<Type> > tMainWallField
        (
            new Field<Type>(interPatch.size(), pTraits<Type>::zero)
        );
    Foam::Field<Type>& mainWallField(tMainWallField());
    bool multiMainOrtoCells(false);

    const labelList
        &faceCells(interPatch.faceCells()),
        &paraFaces_(paraFaces()),
        &ortoFaces_(ortoFaces());
    const labelListList
        &interFaceInfo_(interFaceInfo()),
        &mainWallFacesToSubCells_(mainWallFaceToSubCells());
    const scalarField
        &mainV(mainMesh_.V()),
        &subV(subMesh_->V());
    const vectorField
        &mainC(mainMesh_.C()),
        &subC(subMesh_->C());
    forAll(paraFaces_, indx)
    {
        const label
            sFacI(paraFaces_[indx]),
            wFacI(interFaceInfo_[1][sFacI]);
        if (wFacI<0) // multiple main orto cells (overlapping sub grid)
        {
            multiMainOrtoCells = true;
            continue;
        }
        assert(not mainWallFacesToSubCells_[wFacI].empty());
        forAll(mainWallFacesToSubCells_[wFacI], jndx)
        {
            const label sCelI(mainWallFacesToSubCells_[wFacI][jndx]);
            mainWallField[sFacI] += subField[sCelI]*subV[sCelI];
        }
        mainWallField[sFacI] /= mainV[interFaceInfo_[2][sFacI]];
    }
    // Map rather than average
    if (multiMainOrtoCells)
    {
      Info << "MAP FIELD ON MAIN FOR ORTO" << endl;
        const scalar SML(1e-9);
        forAll(ortoFaces_, indx)
        {
            const label
                sFacI(ortoFaces_[indx]),
                mCelI(interFaceInfo_[2][sFacI]),
                sCelI(faceCells[sFacI]);
            // Assume perfect one-to-one overlap between main and sub cells
            const scalar
                sV(subV[sCelI]),
                sCbRV(pow(sV, 1./3.));
            assert(mag(mainC[mCelI]-subC[sCelI])/sCbRV < SML);
            assert((mainV[mCelI] - sV)/sV < SML);
            mainWallField[sFacI] = subField[sCelI];
        }
    }
    return tMainWallField;
}


template<typename Type>
void Foam::subSolver::plotSrcP(Field<Type>&)
{
    // Add no correction, a scalar version overrides this one
}


template<class Type>
Foam::tmp<Field<Type> > Foam::subSolver::projectOnPatch
(
    const Foam::Field<Type>& field,
    const label
)
{
    // Add no correction, a vector version overrides this one
    Foam::tmp<Foam::Field<Type> > tField
        (
            new Field<Type>(field)
        );
    return tField;
}


template<class Type>
Foam::tmp<Field<Type> > Foam::subSolver::projectOnWall
(
    const Foam::Field<Type>& field,
    const label
)
{
    // Add no correction, a vector version overrides this one
    Foam::tmp<Foam::Field<Type> > tField
        (
            new Field<Type>(field)
        );
    return tField;
}


template<typename Type>
void Foam::subSolver::fixInterpolation
(
    Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& mainField,
    const bool onlyPositive
)
{
    // FIXME: For each orto face a main cell is needed to store the velocity
    //        to be able to calculate a correct gradU in turb::S2.
    //        To achieve this some meshing manipulations need to be done:
    //        - refined main, create sub mesh, coarsen main except orto cells
    //        - split main, copy sub, split one sub in orto and rest, refine,
    //          merge and stitch (alternative use Arbitrary Mesh Interface)
    //          in OF this is named AMI and OFext it is called GGI
    // FIXME: Arrange in python to use MBdP with refinements above.
    //        Govern from cfd what patch uses MBdP and NWT, respectively.
    //
    // FIXME: To stitch in certain conditions does not work, see below
    //        http://www.cfd-online.com/Forums/openfoam-meshing-utilities/113863-when-can-stitchmesh-used.html
    //        which is updated Jan 2016 with a solution about using cyclicAMI
    //        instead of stitchMesh. The important thing is to use
    //        "coincidentFullMatch" for the "transform" value, see
    //        http://www.cfdsupport.com/water-turbine-cfd-manual/node42.html
    //        There is also a user tool
    //        https://openfoamwiki.net/index.php/Contrib/SplitMeshWithSets
    //        My beleif is that a while loop in stitchMesh fails when it does
    //        not come to an end, as the cylinder in a donut, or a square in
    //        an L-shape or similar.
    //        stitchMesh uses polyTopoChanger in dynamicMesh where I beleive
    //        the loop reside.

    Info << "Update " << mainField.name() << " in " << numPatchName_
        << " wall cells to evaluate interpolation correctly at interface."
        << endl;
    typedef GeometricField<Type, fvPatchField, volMesh> volTypeField;
    const volTypeField&
        subField
        (
            subMesh_->thisDb().lookupObject<volTypeField>(mainField.name())
        );

//- Setting main wall face value improves E term in Launder-Sharma but breaks
//  something else resulting in overall worse results
//- Keep away from the wall and the E does not come in at all for main.
    const polyPatch&
        interPatch(subMesh_->boundary()[interPatchID()].patch()),
        wallPatch(mainMesh_.boundary()[wallPatchID_].patch());
    fvPatchField<Type>& mainWallBnd(mainField.boundaryField()[wallPatchID_]);

    const vectorField
        &mC(mainMesh_.C()),
        &sC(subMesh_->C()),
        &Cf(interPatch.faceCentres()),
        &Sf(interPatch.faceCentres());
    const labelListList& interFaceInfo_(interFaceInfo());
    const labelList
        &faceCells(interPatch.faceCells()),
        &wallCells(wallPatch.faceCells());

    // Fix main cell value
    forAll(interPatch, sFacI)
    {
        const label
            sCelI(faceCells[sFacI]),
            wallI(interFaceInfo_[2][sFacI]),
            nextI(interFaceInfo_[3][sFacI]);
        const scalar
            dNf = mag(mC[nextI] - Cf[sFacI]),
            wM = dNf/mag(mC[nextI] - mC[wallI]),
            wS = dNf/mag(mC[nextI] - sC[sCelI]);

        // Calculate face value from main and sub
        const Type faceValue((1-wS)*mainField[nextI] + wS*subField[sCelI]);

        // Set main wall cell so interpolation gives correct face value
        mainField[wallI] = (faceValue - (1-wM)*mainField[nextI]) / wM;

        if (onlyPositive) // Needed for turbulent entities
        {
            const Type small(1e-12*pTraits<Type>::one);
            mainField[wallI] = max(small, mainField[wallI]);
        }
    }

    // Average main corner cells
    const labelPairList &cornerFaces_(cornerFaces());
    forAll(cornerFaces_, indx)
    {
        label wallI(-1);
        scalar sumArea(0);
        Type aveValue(pTraits<Type>::zero);
        forAll(cornerFaces_[indx], jndx)
        {
            const label
                sFacI(cornerFaces_[indx][jndx]),
                sCelI(faceCells[sFacI]),
                nextI(interFaceInfo_[3][sFacI]);
            wallI = interFaceInfo_[2][sFacI];
            const scalar
                faceArea(mag(Sf[sFacI])),
                dNf(mag(mC[nextI] - Cf[sFacI])),
                wM(dNf/mag(mC[nextI] - mC[wallI])),
                wS(dNf/mag(mC[nextI] - sC[sCelI]));
            const Type
                faceValue((1-wS)*mainField[nextI] + wS*subField[sCelI]);
            aveValue += faceArea * (faceValue - (1-wM)*mainField[nextI]) / wM;
            sumArea += faceArea;
        }
        // Interpolate main cell from each face value with face area as weight
        mainField[wallI] = aveValue / sumArea;
    }

    // Fix main boundary value
    forAll(wallPatch, wFacI)
    {
        const label wCelI(wallCells[wFacI]);

        // Simulate a homogeneous Neumann, needed for Dirichlet conditions
        mainWallBnd[wFacI] = mainField[wCelI];

        // FIXME: The result on turbulent channel gets worse by this
        //            May this be suitable for WSS not using Dirichlet?
        // Set bnd value so interpolation gives correct wall value
        //mainWallBnd[wFacI] = 2*mainField[wCelI] - faceValue;
    }
}



// ************************************************************************* //
