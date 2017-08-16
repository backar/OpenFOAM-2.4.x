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
#include "IOobjectList.H"
#include "subSolver.H"


// * * * * * * * * * * * * * * Static Data Member  * * * * * * * * * * * * * //

// Cannot use HashPtrTable as fvMesh cannot be copied, but fvMesh* can
HashTable<fvMesh*> subSolver::subMeshes_ = HashTable<fvMesh*>();


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::subSolver::init()
{
    if (not subMeshes_.found(numPatchName_))
    {
        subMeshes_.insert(numPatchName_, new fvMesh(meshIOobject(mainMesh_)));
    }
    subMesh_ = subMeshes_[numPatchName_];
    //    CAUTION to place stuff here that is NOT meant for the sub
    //            e.g. assert problem in yPlus has been observed
//    cornerFaces();
}


IOobject Foam::subSolver::meshIOobject(const fvMesh& mesh)
{
    return IOobject
    (
        numPatchName_,
        mesh.time().timeName(),
        mesh.time(),
        IOobject::MUST_READ
    );
}


IOobject Foam::subSolver::fieldIOobject
(
    const word fieldName,
    const fvMesh& mesh,
    const IOobject::readOption rOpt
)
{
    return IOobject
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        rOpt,
        IOobject::AUTO_WRITE
    );
}


label Foam::subSolver::interPatchID() const
{
    return subMesh_->boundaryMesh().findPatchID(intrfcPatchName_);
}


label Foam::subSolver::wallPatchID() const
{
    return subMesh_->boundaryMesh().findPatchID(numPatchName_);
}


const labelList Foam::subSolver::paraFaces() const
{
    static HashTable<labelList, word, string::hash> paraFaces;
    if (not paraFaces.found(numPatchName_))
    {
        const polyPatch
            &interPatch(subMesh_->boundary()[interPatchID()].patch()),
            &wallPatch(subMesh_->boundary()[wallPatchID()].patch());
        labelList labels(0);
        forAll(interPatch, sFacI)
        {
            if (interFaceInfo()[6][sFacI]) labels.append(sFacI);
        }
        assert(labels.size() == wallPatch.size());
        paraFaces.insert(numPatchName_, labels);
    }
    return paraFaces[numPatchName_];
}


const labelList Foam::subSolver::ortoFaces() const
{
    // Require one-to-one main and sub for orto faces
    static HashTable<labelList, word, string::hash> ortoFaces;
    if (not ortoFaces.found(numPatchName_))
    {
        const polyPatch
            &interPatch(subMesh_->boundary()[interPatchID()].patch()),
            &wallPatch(subMesh_->boundary()[wallPatchID()].patch());
        labelList labels(0);
        forAll(interPatch, sFacI)
        {
            if (not interFaceInfo()[6][sFacI]) labels.append(sFacI);
        }
        assert(labels.size() == (interPatch.size() - wallPatch.size()));
        HashTable<label, label, string::hash> mainCells;
        forAll(labels, indx)
        {
            const label
                sFacI(labels[indx]),
                mCelI(interFaceInfo()[2][sFacI]);
            if (not mainCells.found(mCelI))
            {
                mainCells.insert(mCelI, sFacI);
            }
            else
            {
                assert(false);  //  orto requires one-to-one mapping
            }
        }
        ortoFaces.insert(numPatchName_, labels);
    }
    return ortoFaces[numPatchName_];
}


const labelPairList Foam::subSolver::cornerFaces() const
{
    // ortoFaces are stored in first and paraFaces in second
    static HashTable<labelPairList, word, string::hash> cornerFaces;
    if (not cornerFaces.found(numPatchName_))
    {
        const fvPatch &interPatch(subMesh_->boundary()[interPatchID()]);
        labelPairList labels(0);
        forAll(ortoFaces(), indx)
        {
            forAll(paraFaces(), jndx)
            {
                if (interPatch.faceCells()[ortoFaces()[indx]] ==
                    interPatch.faceCells()[paraFaces()[jndx]])
                {
                    const labelPair
                        corner(ortoFaces()[indx], paraFaces()[jndx]);
                    labels.append(corner);
                }
            }
        }
        cornerFaces.insert(numPatchName_, labels);
        label nDir(0);
        forAll(paraFaces(), i)
        {
            if (interFaceInfo()[5][paraFaces()[i]]>0) nDir++;
        }
        Info
            << numPatchName_ << ": para; "
            << paraFaces().size() << ", " << nDir << endl;
        nDir = 0;
        forAll(ortoFaces(), i)
        {
            if (interFaceInfo()[5][ortoFaces()[i]]>0) nDir++;
        }
        Info
            << numPatchName_ << ": orto; "
            << ortoFaces().size() << ", " << nDir << endl;
    }
    return cornerFaces[numPatchName_];
}


const labelListList Foam::subSolver::interFaceInfo() const
{
    static HashTable<labelListList, word, string::hash> interfaceInfo;
    if (not interfaceInfo.found(numPatchName_))
    {
        // FIXME: Add opposingFaceLabel(faceLabel, meshFaces) placed in 1
        // Place internal faces in 0, wall faces in 1,
        // wall cells in 2, inner cells in 3, mark in 4 if wall is owner,
        // parallel direction in 5, and mark in 6 for para
        const scalar SML(1e-12);
        const unsigned size = 7;
        const fvPatch
            &interPatch(subMesh_->boundary()[interPatchID()]),
            &wallPatch(mainMesh_.boundary()[wallPatchID_]);
        const vectorField
            pCf(interPatch.Cf()),
            wCf(wallPatch.Cf()),
            mCf(mainMesh_.Cf()),
            pSf(interPatch.Sf()),
            wSf(wallPatch.Sf()),
            mSf(mainMesh_.Sf()),
            pNf(interPatch.nf()),
            wNf(wallPatch.nf()),
            mNf(mainMesh_.Sf()/mainMesh_.magSf());
        const labelList wallCellLbls(wallPatch.patch().faceCells());
        labelListList infoLbls(size);
        for (unsigned i=0; i<size; i++)
        {
            // Initiate to correct size
            infoLbls[i].setSize(interPatch.size());
            label dummy(-1);
            if (i==5) dummy = 0;
            for (label j=0; j<interPatch.size(); j++) infoLbls[i][j] = dummy;
        }
        forAll(interPatch, sFacI)
        {
            label foundMainFaceI(-1);
            forAll(mCf, mFI)
            {
                // Do face centres coincide, and face area the same?
                // Finds all para (and possibly some orto) faces
                // Same center and size - perfect overlap
                if ((mag(pCf[sFacI] - mCf[mFI]) < SML) and
                    ((mag(pSf[sFacI] - mSf[mFI]) < SML) or
                     (mag(pSf[sFacI] + mSf[mFI]) < SML)))
                {
                    foundMainFaceI = mFI;
                    break;
                }
            }
            infoLbls[0][sFacI] = foundMainFaceI;
            if (infoLbls[0][sFacI]>-1)
            {
                const scalar mFacI(foundMainFaceI);
                bool isWallCellOwner = false;
                // bnd face is neither set for corner nor orto faces
                // except the one closest to bnd
                forAll(wallCellLbls, wFacI)
                {
                    if (wallCellLbls[wFacI]==mainMesh_.owner()[mFacI])
                    {
                        isWallCellOwner = true;
                        infoLbls[1][sFacI] = wFacI;
                        break;
                    }
                    else
                    if (wallCellLbls[wFacI]==mainMesh_.neighbour()[mFacI])
                    {
                        isWallCellOwner = false;
                        infoLbls[1][sFacI] = wFacI;
                        break;
                    }
                }
                if (isWallCellOwner)
                {
                    infoLbls[2][sFacI] = mainMesh_.owner()[mFacI];
                    infoLbls[3][sFacI] = mainMesh_.neighbour()[mFacI];
                    infoLbls[4][sFacI] = 1;
                } else
                {
                    infoLbls[2][sFacI] = mainMesh_.neighbour()[mFacI];
                    infoLbls[3][sFacI] = mainMesh_.owner()[mFacI];
                    infoLbls[4][sFacI] = 0;
                }
                const scalar direction(pNf[sFacI] & mNf[mFacI]);
                // sub and main share the same face, dir must be -1 or 1
                assert(fabs(direction) - 1.0 < SML);
                infoLbls[5][sFacI] = static_cast<label>(direction);
                infoLbls[6][sFacI] = 1;  // May be a para face
                forAll(wallPatch, wFacI)
                {
                    // Mark orto faces with 0
                    // Is face normal not parallel with main wall patch?
                    if (fabs(pNf[sFacI] & wNf[wFacI]) < 0.5) // Skewed cells?
                    {
                        infoLbls[6][sFacI] = 0;  // Must be an orto face
                    }
                }
            }
            // FIXME: Test again if still valid
            //        Cannot have overlapping sub-grid e.g. out and wall
            //        at least not for non-uniform main mesh.
            //        As the problem has been most obvious for n-grid, the
            //        problem might be in the update of the main mesh. Other
            //        reasons might be that in a corner the next main cell is
            //        also a wall cell, mind how T and k treat the next cell.
            //        The corner cell gets (correctly?) wall shear stress
            //        contribution from two walls in the U-equation.

            // Check for orto face coinciding with main inner face
            if ((infoLbls[0][sFacI] == -1) or not infoLbls[6][sFacI])
            {
                scalar
                    minDiffCfs(1e9),
                    mFacI(-1);
                forAll(mCf, mFI)
                {
                    // Same normal and orthogonal dist gives same plane,
                    // loop to find min dist
                    const vector Diff(pCf[sFacI]-mCf[mFI]);
                    if ((mag(pSf[sFacI] ^ mSf[mFI]) < SML) and
                        (mag(pSf[sFacI] & Diff) < SML) and
                        (mag(Diff) < minDiffCfs))
                    {
                        minDiffCfs = mag(Diff);
                        mFacI = mFI;
                    }
                }
                assert(mFacI>-1);
                infoLbls[0][sFacI] = mFacI;
                infoLbls[1][sFacI] = -1;  // Marker of multiple main orto cells
                const label
                    nei(mainMesh_.neighbour()[mFacI]),
                    own(mainMesh_.owner()[mFacI]),
                    sub(interPatch.faceCells()[sFacI]);
                const vector
                    nDist(mainMesh_.C()[nei] - pCf[sFacI]),
                    oDist(mainMesh_.C()[own] - pCf[sFacI]),
                    sDist(subMesh_->C()[sub] - pCf[sFacI]);
                // Is own-face or nei-face aligned with sub-face
                if ((oDist & sDist) > 0)
                {
                    infoLbls[2][sFacI] = mainMesh_.owner()[mFacI];
                    infoLbls[3][sFacI] = mainMesh_.neighbour()[mFacI];
                    infoLbls[4][sFacI] = 1;
                } else
                {
                    assert((nDist & sDist) > 0);
                    infoLbls[2][sFacI] = mainMesh_.neighbour()[mFacI];
                    infoLbls[3][sFacI] = mainMesh_.owner()[mFacI];
                    infoLbls[4][sFacI] = 0;
                }
                const scalar direction(pNf[sFacI] & mNf[mFacI]);
                // sub and main are in same plane, dir must be -1 or 1
                assert(fabs(direction) - 1.0 < SML);
                infoLbls[5][sFacI] = static_cast<label>(direction);
                infoLbls[6][sFacI] = 0;
            }
            for (unsigned i=0; i<size; i++)
            {
                if (i==5) assert(abs(infoLbls[i][sFacI]) == 1);
                else if (i!=1) assert(infoLbls[i][sFacI] > -1);
            }
        }
        interfaceInfo.insert(numPatchName_, infoLbls);
    }
    // FIXME: Returning interFaceInfo structure may take considerable time
    return interfaceInfo[numPatchName_];
}


const labelList Foam::subSolver::wallFacesToMainFaces() const
{
    // Reach the main faces at the wall from the overlapping sub faces,
    // used in wallSource
    static HashTable<labelList, word, string::hash> wallFace2MainFaces;
    if (not wallFace2MainFaces.found(numPatchName_))
    {
        const scalar SML(1e-12);
        const polyPatch
            &subPatch(subMesh_->boundary()[wallPatchID()].patch()),
            &mainPatch(mainMesh_.boundary()[wallPatchID_].patch());
        // Should be equal unless the sub subgrid is shortened
        assert(subPatch.size() <= mainPatch.size());
        labelList wallFaceLbls(subPatch.size(), -1);
        forAll(subPatch, sFacI)
        {
            forAll(mainPatch, mFacI)
            {
                // Do face centres coincide?
                if (mag(subPatch.faceCentres()[sFacI]-
                        mainPatch.faceCentres()[mFacI]) < SML)
                {
                    wallFaceLbls[sFacI] = mFacI;
                }
            }
            assert(wallFaceLbls[sFacI]>-1);
        }
        wallFace2MainFaces.insert(numPatchName_, wallFaceLbls);
    }
    return wallFace2MainFaces[numPatchName_];
}


const labelListList Foam::subSolver::mainWallFaceToSubCells() const
{
    // To be used to map pressure from main next or next wall to sub
    static HashTable<labelListList, word, string::hash> wallFace2SubCells;
    if (not wallFace2SubCells.found(numPatchName_))
    {
        const scalar SML(1e-10);
        const volVectorField& subCC(subMesh_->C());
        const polyPatch
            &interPatch(subMesh_->boundary()[interPatchID()].patch()),
            &wallPatch(mainMesh_.boundary()[wallPatchID_].patch());

        const label nSubLayers(subCC.size()/paraFaces().size());
        labelListList subCellLbls(wallPatch.size());

        forAll(wallPatch, wFacI)
        {
            // Connect para face with wall face (main)
            label pFacI(-1);
            forAll(paraFaces(), indx)
            {
                const label
                    sFacI(paraFaces()[indx]),
                    wCelI(interFaceInfo()[2][sFacI]);
                if (wCelI == wallPatch.faceCells()[wFacI])
                {
                    pFacI = sFacI;
                    break;
                }
            }

            // In case multiple main corner cells sub cells should overlap
            // one to one, thus, only wall cell can be mapped here.
            // Use ortoFaces to map the rest of the cells
            if (pFacI>-1)
            {
                // interface and wall position vectors
                const vector
                    &paraCf(interPatch.faceCentres()[pFacI]),
                    &wallCf(wallPatch.faceCentres()[wFacI]),
                    wallSf(wallPatch.faceAreas()[wFacI]),
                    wallNf(wallSf/mag(wallSf)),
                    CfCf(paraCf - wallCf),
                    uCfCf(CfCf/mag(CfCf));  // Unit vector
                assert(mag(uCfCf ^ wallNf) < SML);

                // Places list of sub cells overlapping main cell
                Tuple2<label, scalar> bigErr(-1, 1e30);
                label bigErrI(0);
                List<Tuple2<label, scalar> > subSubCells(nSubLayers, bigErr);
                forAll(subCC, sCelI)
                {
                    const vector
                        CcCf(subCC[sCelI] - wallCf),
                        orthErr(CcCf - uCfCf * (uCfCf & CcCf));
                    const scalar rOrthErr(mag(orthErr) / mag(CcCf));

                    if ((subSubCells[bigErrI].first() < 0) or
                        (subSubCells[bigErrI].second() > rOrthErr))
                    {
                        subSubCells[bigErrI] =
                            Tuple2<label, scalar>(sCelI, rOrthErr);

                        // Find largest error currently in subSubCells
                        bigErrI = 0;
                        bigErr = subSubCells[bigErrI];
                        forAll(subSubCells, indx)
                        {
                            if ((subSubCells[indx].first() < 0) or
                                (subSubCells[indx].second() > bigErr.second()))
                            {
                                bigErrI = indx;
                                bigErr = subSubCells[bigErrI];
                            }
                        }
                    }
                }
                // Collapse subCellLbls to zero size
                subCellLbls[wFacI].setSize(0);
                // Append the index having the smallest rOrthErr
                forAll(subSubCells, indx)
                {
                    assert(subSubCells[indx].first() > -1);
                    assert(subSubCells[indx].second() < mag(CfCf)/nSubLayers);
                    subCellLbls[wFacI].append(subSubCells[indx].first());
                }
                assert(subCellLbls[wFacI].size() == nSubLayers);
            }
            else
            {
                // Assert empty list of sub cells if multiple main cells found
                assert(subCellLbls[wFacI].empty());
            }
        }
        wallFace2SubCells.insert(numPatchName_, subCellLbls);
    }
    return wallFace2SubCells[numPatchName_];
}


const List<labelPairList> Foam::subSolver::wallFaceToSubCellFaces() const
{
    // To be used to scale phi sub from main
    // store as wallFace2SubFaces[indx][subIndx](sCelI, sFacI)
    // OBSERVE:
    //    indx:     sub wall faces
    //    subIndx:  index similar to mainWallToSubCells()
    //    sCelI:    cell connected with sFacI with wallDistCc < wallDistCf
    //    sFacI:    internal faces in colon of sub cells
    static HashTable<List<labelPairList>, word, string::hash>
        wallFace2SubFaces;
    if (not wallFace2SubFaces.found(numPatchName_))
    {
        const polyPatch
            &wallPatch(subMesh_->boundary()[wallPatchID()].patch());
        List<labelPairList> subFaceCellLbls(wallPatch.size());
        forAll(wallPatch, wFacI)
        {
            label   // Initiate current face and cell labels
                cFacI(wFacI + wallPatch.start()), // local to global index
                cCelI(wallPatch.faceCells()[wFacI]);
            const cell wCell(subMesh_->cells()[cCelI]);
            cFacI = wCell.opposingFaceLabel(cFacI, subMesh_->faces());
            subFaceCellLbls[wFacI].append(labelPair(cCelI, cFacI));
            while(subMesh_->boundaryMesh().whichPatch(cFacI) < 0) // Internal
            {
                if (mag(cCelI) == subMesh_->neighbour()[cFacI])
                {
                    cCelI = subMesh_->owner()[cFacI];
                }
                else
                if (mag(cCelI) == subMesh_->owner()[cFacI])
                {
                    cCelI = subMesh_->neighbour()[cFacI];
                }
                else
                {
                    assert(false);
                }
                const cell cCell(subMesh_->cells()[cCelI]);
                cFacI = cCell.opposingFaceLabel(cFacI, subMesh_->faces());
                subFaceCellLbls[wFacI].append(labelPair(cCelI, cFacI));
            }
        }
        wallFace2SubFaces.insert(numPatchName_, subFaceCellLbls);
    }
    return wallFace2SubFaces[numPatchName_];
}


// FIXME: Move to icoSolver, why it concerns the mesh more than the solver?
const List<List<labelPairList> > Foam::subSolver::mainFaceToSubFaces() const
{
    /*
      To be used to scale phi sub from main
      store as paraFace2SubFaces[bndPatch][Indx][subIndx](mFacI, sFacI)
      with internal faces at bndPatch=nPatch as given by
      mesh.boundaryMesh().whichPatch(c[i]), where c is a cell
      OBSERVE: Indx equals:
         1. sFacI for "the wall" boundary, wallFacesToMainFaces (one-to-one)
         2. index of orto/paraFaces() for interFace patch       (one-to-one)
         3. mFacI for over-lapping internal sub faces          (one-to-many)
         4. sFacI according to orto/paraFaces() for other bnds (one-to-many)
    */
    static HashTable<List<List<labelPairList> >, word, string::hash>
        paraFace2SubFaces;
    if (not paraFace2SubFaces.found(numPatchName_))
    {
        const scalar SML(1e-12);
        const label nSubLayers(subMesh_->cells().size()/paraFaces().size());

        const vectorField
            mCf(mainMesh_.Cf()),
            mSf(mainMesh_.Sf()),
            mNf(mainMesh_.Sf()/mainMesh_.magSf()),
            sCf(subMesh_->Cf()),
            sSf(subMesh_->Sf()),
            sNf(subMesh_->Sf()/subMesh_->magSf());
        const polyBoundaryMesh
            &mBndMesh(mainMesh_.boundaryMesh()),
            &sBndMesh(subMesh_->boundaryMesh());

        const label sIntI(subMesh_->boundaryMesh().size());
        List<List<labelPairList> > subFaceLbls(sIntI + 1);
        forAll(subFaceLbls, ptchI)
        {
            // Collapse inner list to zero size
            subFaceLbls[ptchI].setSize(0);
        }

        // Append all wall faces
        // OBSERVE: 1. Ordered by sub wall faces (indx of wallFacesToMainFaces)
        forAll(wallFacesToMainFaces(), sFacI)
        {
            const label wFacI(wallFacesToMainFaces()[sFacI]);
            const labelPairList wallPair(1, labelPair(wFacI, sFacI));
            subFaceLbls[wallPatchID()].append(wallPair);
        }

        forAll(paraFaces(), indx)
        {
            const label
                sIFcI(paraFaces()[indx]),
                mIFcI(interFaceInfo()[0][sIFcI]),
                mWFcI(interFaceInfo()[1][sIFcI]), // -1 corner if multiple orto
                mCelI(interFaceInfo()[2][sIFcI]);

            const labelPairList interPair(1, labelPair(mIFcI, sIFcI));
            subFaceLbls[interPatchID()].append(interPair);

            bool partOfCorner(false);
            const labelPairList& cornerFaces_(cornerFaces());
            forAll(cornerFaces_, jndx)
            {
                if (sIFcI == cornerFaces_[jndx].second())
                {
                    partOfCorner = true;  // para is part of corner
                    break;
                }
            }
            if (partOfCorner)
            {
                continue;   // All faces on corner cells is added in orto loop
            }
            assert(mWFcI > -1);   // From now on only "para only" faces

            const cell mFaces(mainMesh_.cells()[mCelI]);
            forAll(mFaces, mI)  // Loop through all faces of the main wall cell
            {
                const label
                    mGlbFI(mFaces[mI]),
                    mPtcI(mBndMesh.whichPatch(mGlbFI));
                const word
                    patchName((mPtcI==-1) ? "" : mBndMesh.names()[mPtcI]);
                const label                     // Set mFacI to global face:
                    mFacI((mPtcI==-1) ? mGlbFI :                  // internal
                          mBndMesh[patchName].whichFace(mGlbFI)); // bnd

                // wall or para interFace - no need to loop over sub cells
                if ((mPtcI == wallPatchID_) or
                    ((mPtcI == -1) and (mFacI == mIFcI)))
                {
                    continue;   // para or wall face already appended above
                }
                else
                if ((mPtcI == -1) and (subFaceLbls[sIntI].size() > mFacI) and
                    (subFaceLbls[sIntI][mFacI].size() > 0))
                {
                  // FIXME: Find condition separating normal case from when
                  //        current cell is next to the orto cells (one-to-one)
                  //        then size is equal to unity instead of nLayers
                    assert((subFaceLbls[sIntI][mFacI].size() == 1) or
                           (subFaceLbls[sIntI][mFacI].size() == nSubLayers));
                    continue;   // internal main face, mFacI, appended earlier
                }

                const word intPtchName(intrfcPatchName_);   // Short name
                scalar subFacesArea(0);
                forAll(mainWallFaceToSubCells()[mWFcI], jndx)
                {
                    const label sCelI(mainWallFaceToSubCells()[mWFcI][jndx]);
                    const cell sFaces(subMesh_->cells()[sCelI]);
                    forAll(sFaces, sI)
                    {
                        const label
                            sGlbFI(sFaces[sI]),
                            sPtcI(sBndMesh.whichPatch(sGlbFI));
                        // orto face
                        if ((mPtcI == -1) and (sPtcI==interPatchID()) and
                            (sIFcI!=sBndMesh[intPtchName].whichFace(sGlbFI)))
                        {
                            assert(false); // Corner is omitted, see above
                            // FIXME: Check if only para needs it, maybe
                            assert(not ortoFaces().size());
                            const label
                                sFacI(sBndMesh[intPtchName].whichFace(sGlbFI));
                            const polyPatch
                                &sPatch(sBndMesh[intPtchName]);
                            const vector
                                spCf(sPatch.faceCentres()[sFacI]),
                                spSf(sPatch.faceAreas()[sFacI]),
                                spNf(spSf/mag(sPatch.faceAreas()[sFacI]));
                            // Same normal and orthog dist gives same plane
                            if ((mag(mNf[mFacI] ^ spNf) < SML) and
                                (mag(mNf[mFacI] & (mCf[mFacI]-spCf)) < SML))
                            {
                                assert(mag(mNf[mFacI] & spNf) - 1 < SML);
                                const labelPair bndPair(mFacI, sFacI);
                                // para at 0, all ortos appended after
                                // OBSERVE: 2. Ordered by indx of paraFaces()
                                subFaceLbls[sPtcI][indx].append(bndPair);
                                subFacesArea += mag(spSf);
                            }
                        }
                        else
                        if ((sPtcI==interPatchID()) or (sPtcI==wallPatchID()))
                        {
                            continue;   // already appended above
                        }
                        else
                        if ((mPtcI == -1) and (sPtcI == -1))  // internal face
                        {
                            const label sFacI(sGlbFI);
                            // Same normal and orthogonal dist for same plane
                            if ((mag(mNf[mFacI] ^ sNf[sFacI]) < SML) and
                                (mag(mNf[mFacI] & (mCf[mFacI]-sCf[sFacI]))<SML)
                               )
                            {
                                assert(mag(mNf[mFacI] & sNf[sFacI]) - 1 <SML);
                                const labelPair internPair(mFacI, sFacI);
                                if (subFaceLbls[sIntI].size() < 1+mFacI)
                                {
                                    subFaceLbls[sIntI].setSize(1+mFacI);
                                }
                                // OBSERVE: 3. Ordered by mFacI for internals
                                subFaceLbls[sIntI][mFacI].append(internPair);
                                subFacesArea += subMesh_->magSf()[sFacI];
                            }
                        }
                        else
                        if ((mPtcI == -1) or (sPtcI == -1))
                        {
                            continue;   // both need to be internal faces
                        }
                        else
                        if (patchName == sBndMesh.names()[sPtcI]) //  boundary
                        {
                            const label
                                sFacI(sBndMesh[patchName].whichFace(sGlbFI));
                            assert(mBndMesh.names()[mPtcI] == patchName);
                            const polyPatch
                                &mPatch(mBndMesh[patchName]),
                                &sPatch(sBndMesh[patchName]);
                            const vector
                                mpCf(mPatch.faceCentres()[mFacI]),
                                mpSf(mPatch.faceAreas()[mFacI]),
                                mpNf(mpSf/mag(mPatch.faceAreas()[mFacI])),
                                spCf(sPatch.faceCentres()[sFacI]),
                                spSf(sPatch.faceAreas()[sFacI]),
                                spNf(spSf/mag(sPatch.faceAreas()[sFacI]));
                            // Same normal and orthog dist gives same plane
                            if ((mag(mpNf ^ spNf) < 1000*SML) and
                                (mag(mpNf & (mpCf-spCf)) < SML)
                               )
                            {
                                assert(mag(mpNf & spNf) - 1 < SML);
                                const labelPair bndPair(mFacI, sFacI);
                                subFaceLbls[sPtcI].setSize(1+sIFcI);
                                // OBSERVE: 4. Ordered by para sFacI for bnds
                                subFaceLbls[sPtcI][sIFcI].append(bndPair);
                                subFacesArea += mag(spSf);
                            }
                        }
                    }
                }
                if (mPtcI == -1)        // main internal face
                {
                    assert((mag(mainMesh_.magSf()[mFacI])-subFacesArea) < SML);
                }
                else
                if (subFacesArea > 0)   // patch is part of cell mCelI
                {
                    const vector
                        mainFaceArea(mBndMesh[patchName].faceAreas()[mFacI]);
                    assert((mag(mainFaceArea) - subFacesArea) < SML);
                }
            }
        }

        /*
        const fvPatch& interPatch(subMesh_->boundary()[interPatchID()]);
        // FIXME: only needed if ortoFaces().size()>1
        // FIXME: orto not tested lately, merge para and orto
        forAll(ortoFaces(), indx)
        {
            const label
                sIFcI(ortoFaces()[indx]),
                sCelI(interPatch.faceCells()[sIFcI]),
                mIFcI(interFaceInfo()[0][sIFcI]),
              //mWFcI(interFaceInfo()[1][sIFcI]), // -1 as orto, except bottom
                mCelI(interFaceInfo()[2][sIFcI]);

            const labelPairList interPair(1, labelPair(mIFcI, sIFcI));
            // FIXME: at appending para, it says para at 0 orto after
            //        but here orto is first?!?
            // OBSERVE: 2. Ordered by indx of ortoFaces()
            subFaceLbls[interPatchID()].append(interPair);

            const cell mFaces(mainMesh_.cells()[mCelI]);
            forAll(mFaces, mI)  // Loop through all faces of the main wall cell
            {
                const label
                    mGlbFI(mFaces[mI]),
                    mPtcI(mBndMesh.whichPatch(mGlbFI));
                const word
                    patchName((mPtcI==-1) ? "" : mBndMesh.names()[mPtcI]);
                const label                     // Set mFacI to global face:
                    mFacI((mPtcI==-1) ? mGlbFI :                  // internal
                          mBndMesh[patchName].whichFace(mGlbFI)); // bnd

                // orto interFace - no need to loop over sub cells
                if ((mPtcI == wallPatchID_) or
                    ((mPtcI == -1) and (mFacI == mIFcI)))
                {
                    continue;   // orto or wall face, already appended above
                }
                else
                if ((mPtcI == -1) and (subFaceLbls[sIntI].size() > mFacI) and
                    (subFaceLbls[sIntI][mFacI].size() > 0))
                {
                    assert(subFaceLbls[sIntI][mFacI].size() == 1);
                    continue;   // internal main face, mFacI, appended earlier
                }

                const word intPtchName(intrfcPatchName_);   // Short name
                scalar subFacesArea(0);
                const cell sFaces(subMesh_->cells()[sCelI]);
                forAll(sFaces, sI)
                {
                    const label
                        sGlbFI(sFaces[sI]),
                        sPtcI(sBndMesh.whichPatch(sGlbFI));
                    // para face
                    if ((mPtcI == -1) and (sPtcI==interPatchID()) and
                        (sIFcI!=sBndMesh[intPtchName].whichFace(sGlbFI)))
                    {
                        continue;   // This will be added below in para loop
                    }
                    else
                    if ((sPtcI==interPatchID()) or (sPtcI==wallPatchID()))
                    {
                        continue;   // already appended above, or will be
                    }
                    else
                    if ((mPtcI == -1) and (sPtcI == -1))  // internal face
                    {
                        const label sFacI(sGlbFI);
                        // Same normal and orthogonal dist for same plane
                        if ((mag(mNf[mFacI] ^ sNf[sFacI]) < SML) and
                            (mag(mNf[mFacI] & (mCf[mFacI]-sCf[sFacI]))<SML)
                           )
                        {
                            assert(mag(mNf[mFacI] & sNf[sFacI]) - 1 <SML);
                            const labelPair internPair(mFacI, sFacI);
                            if (subFaceLbls[sIntI].size() < 1+mFacI)
                            {
                                subFaceLbls[sIntI].setSize(1+mFacI);
                            }
                            // OBSERVE: 3. Ordered by mFacI for internals
                            subFaceLbls[sIntI][mFacI].append(internPair);
                            subFacesArea += subMesh_->magSf()[sFacI];
                        }
                    }
                    else
                    if ((mPtcI == -1) or (sPtcI == -1))
                    {
                        continue;   // both need to be internal faces
                    }
                    else
                    if (patchName == sBndMesh.names()[sPtcI]) //  boundary
                    {
                        const label
                            sFacI(sBndMesh[patchName].whichFace(sGlbFI));
                        assert(mBndMesh.names()[mPtcI] == patchName);
                        const polyPatch
                            &mPatch(mBndMesh[patchName]),
                            &sPatch(sBndMesh[patchName]);
                        const vector
                            mpCf(mPatch.faceCentres()[mFacI]),
                            mpSf(mPatch.faceAreas()[mFacI]),
                            mpNf(mpSf/mag(mPatch.faceAreas()[mFacI])),
                            spCf(sPatch.faceCentres()[sFacI]),
                            spSf(sPatch.faceAreas()[sFacI]),
                            spNf(spSf/mag(sPatch.faceAreas()[sFacI]));
                        // Same normal and orthog dist gives same plane
                        if ((mag(mpNf ^ spNf) < 1000*SML) and
                            (mag(mpNf & (mpCf-spCf)) < SML)
                           )
                        {
                            assert(mag(mpNf & spNf) - 1 < SML);
                            const labelPair bndPair(mFacI, sFacI);
                            subFaceLbls[sPtcI].setSize(1+sIFcI);
                            // OBSERVE: 4. Ordered by orto sFacI for bnds
                            subFaceLbls[sPtcI][sIFcI].append(bndPair);
                            subFacesArea += mag(spSf);
                        }
                    }
                }
            }
        }
        */

        // Check that map seems OK
        forAll(subFaceLbls, sPtcI)
        {
            const word patchName((sPtcI<sBndMesh.size()) ?
                                 sBndMesh.names()[sPtcI] : "Internal");
            long bndSize(0);  // Boundary size per boundary
            vectorField mCf, mNf, sCf, sNf;
            if (sPtcI == subMesh_->boundaryMesh().size()) // internal faces
            {
                mCf = mainMesh_.Cf();
                sCf = subMesh_->Cf();
                mNf = mainMesh_.Sf()/mag(mainMesh_.Sf());
                sNf = subMesh_->Sf()/mag(subMesh_->Sf());
                sNf = subMesh_->Sf()/mag(subMesh_->Sf());
            }
            else
            if (sPtcI == interPatchID())                  // inter patch faces
            {
                mCf = mainMesh_.Cf();
                sCf = subMesh_->boundaryMesh()[interPatchID()].faceCentres();
                mNf = mainMesh_.Sf()/mag(mainMesh_.Sf());
                const vectorField
                    sSf(subMesh_->boundaryMesh()[interPatchID()].faceAreas());
                sNf = sSf/mag(sSf);
            }
            else                                          // boundary faces
            {
                const word patchName(subMesh_->boundaryMesh().names()[sPtcI]);
                const label
                    mPtcI(mainMesh_.boundaryMesh().findPatchID(patchName));
                mCf = mainMesh_.boundaryMesh()[mPtcI].faceCentres();
                sCf = subMesh_->boundaryMesh()[sPtcI].faceCentres();
                const vectorField
                    mSf(mainMesh_.boundaryMesh()[mPtcI].faceAreas()),
                    sSf(subMesh_->boundaryMesh()[sPtcI].faceAreas());
                mNf = mSf/mag(mSf);
                sNf = sSf/mag(sSf);
            }
            forAll(subFaceLbls[sPtcI], indx)
            {
                scalar errCum(0), errMx(1e-30);
                forAll(subFaceLbls[sPtcI][indx], jndx)
                {
                    const label
                        mFacI(subFaceLbls[sPtcI][indx][jndx].first()),
                        sFacI(subFaceLbls[sPtcI][indx][jndx].second());
                    assert((mag(mNf[mFacI] - sNf[sFacI]) < 1e-11) or
                           (mag(mNf[mFacI] + sNf[sFacI]) < 1e-11));
                    errCum += (vector(1,1,1)&(mCf[mFacI]-sCf[sFacI]));
                    errMx = max(errMx, mag(mCf[mFacI]-sCf[sFacI]));
                }
                assert(mag(errCum)<errMx*nSubLayers);
                if (subFaceLbls[sPtcI][indx].size()>0)
                {
                    bndSize++;
                }
                if (sPtcI == wallPatchID())
                {
                    assert(subFaceLbls[sPtcI][indx].size()==1);
                }
                else
                // FIXME: Identify indx for corner cell, check orto separatly
                if (sPtcI == interPatchID())
                {
                    assert
                    (
                        (subFaceLbls[sPtcI][indx].size()==1) or
                        (subFaceLbls[sPtcI][indx].size()==1+ortoFaces().size())
                    );
                }
                else
                // Int and bnds except wall and inter is a multiple of nLayers
                {
                    // OBSERVE: size is not necessary equal to nSubLayers
                    //          if bnd contains faces from two sides of a cell
                    //          And if it is also close to orto cells
                    //          (one-to-one) size could be equal to 2.
                    if (not(subFaceLbls[sPtcI][indx].size() % nSubLayers == 0))
                    {
                        if (sPtcI == sIntI)
                        {
                            // Internal face of one-to-one orto cell
                            assert(subFaceLbls[sPtcI][indx].size() == 1);
                            continue;
                        }
                    }
                    assert((subFaceLbls[sPtcI][indx].size() % nSubLayers == 0)
                           or (subFaceLbls[sPtcI][indx].size() == 2));
                }
            }
            if (bndSize and (sPtcI == sIntI))
            {
                // Estimate number of internal faces
                const scalar root(sqrt((paraFaces().size()) + 0.5));
                const long n(root);
                // FIXME: Asserts for 3D or 2D with extChannel
                if (false and bndSize != paraFaces().size() - 1)
                {
                    assert(bndSize >= paraFaces().size() - 1);
                    assert(bndSize <= 2*n*(n - 1));  // 3D
                }
            }
            else    // FIXME: What to assert for orto faces?
            if (bndSize and not (sPtcI == interPatchID())
                                 and (ortoFaces().size() < 2))
            {
                assert(subMesh_->boundaryMesh()[sPtcI].size() % bndSize==0);
            }
        }
        paraFace2SubFaces.insert(numPatchName_, subFaceLbls);
    }
    return paraFace2SubFaces[numPatchName_];
}


// FIXME: Transfer to mainOrtoToSubCells
// FIXME: change from surfLabelField to Hash(bndName) of labelList
const GeometricField<scalar, fvsPatchField, surfaceMesh>
Foam::subSolver::faceLabelField()
{
    typedef GeometricField<scalar, fvsPatchField, surfaceMesh> surfLabelField;
    static HashTable<surfLabelField, word, string::hash> subFace2MainFace;
    if (not subFace2MainFace.found(numPatchName_))
    {
        const scalar
            GREAT(1e30),
            SMALL(1e-12);
        surfLabelField labelField
            (
                IOobject
                (
                    "faceLabels",
                    subMesh_->magSf().instance(),
                    *subMesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE),
                *subMesh_,
                dimensionedScalar("label", dimTime/dimTime, GREAT)
            );
        labelField.internalField().resize(subMesh_->nInternalFaces());

        const surfaceVectorField
            sCf(subMesh_->Cf()),
            sSf(subMesh_->Sf()),
            mCf(mainMesh_.Cf()),
            mSf(mainMesh_.Sf()),
            sNf(sSf / subMesh_->magSf()),
            mNf(mSf / mainMesh_.magSf());

        const polyPatch &wallPatch(mainMesh_.boundary()[wallPatchID_].patch());
        forAll(mainWallFaceToSubCells(), wFacI)
        {
            const label mCelI(wallPatch.faceCells()[wFacI]);
            const cell mFaces(mainMesh_.cells()[mCelI]);

            const labelList subCells(mainWallFaceToSubCells()[wFacI]);
            forAll(subCells, indx)
            {
                const label sCelI(subCells[indx]);
                const cell sFaces(subMesh_->cells()[sCelI]);
                forAll(sFaces, jndx)
                {
                    const label sFacI(sFaces[jndx]);
                    if (not subMesh_->isInternalFace(sFacI)) continue;
                    if (labelField[sFacI] < GREAT) continue;
                    forAll(mFaces, kndx)
                    {
                        const label mFacI(mFaces[kndx]);
                        if (mainMesh_.isInternalFace(mFacI))
                        {
                            // Same center and size - perfect overlap
                            if ((mag(mCf[mFacI] - sCf[sFacI]) < SMALL) and
                                ((mag(mSf[mFacI] - sSf[sFacI]) < SMALL) or
                                 (mag(mSf[mFacI] + sSf[sFacI]) < SMALL)))
                            {
                                const scalar direc(mNf[mFacI] & sNf[sFacI]);
                                assert(fabs(direc) - 1.0 < SMALL);
                                labelField[sFacI] = direc * mFacI;
                                break;
                            }
                        }
                    }
                    if (labelField[sFacI] == GREAT)
                    {
                        // Same normal and orthogonal dist gives same plane,
                        // loop to find min dist (recall, only one main cell)
                        scalar minDist(1e30);
                        forAll(mFaces, kndx)
                        {
                            const label mFacI(mFaces[kndx]);
                            if (mainMesh_.isInternalFace(mFacI))
                            {
                                const vector Dist(mCf[mFacI] - sCf[sFacI]);
                                if ((mag(mSf[mFacI] ^ sSf[sFacI])<SMALL) and
                                    mag(Dist & mSf[mFacI]) < SMALL)
                                {
                                    if (mag(Dist) < minDist)
                                    {
                                        minDist = mag(Dist);
                                        const scalar
                                            direc(mNf[mFacI]&sNf[sFacI]);
                                        assert(fabs(direc) - 1.0 < SMALL);
                                        labelField[sFacI] = direc*mFacI;
                                    }
                                }
                            }
                        }
                    }
                    //assert(labelField[sFacI] < GREAT);
                    // Assume left-overs are non-overlapping within a main cell
                    // Leave those to have GREAT as label
                }
            }
        }
        forAll(subMesh_->boundaryMesh(), sPatI)
        {
            // The surface vector from a patch is always aimed outwards
            const word patchName(subMesh_->boundaryMesh()[sPatI].name());
            const label mPatI(mainMesh_.boundaryMesh().findPatchID(patchName));
            const polyPatch
                &subPatch(subMesh_->boundaryMesh()[sPatI]),
                &mainPatch(mainMesh_.boundaryMesh()[mPatI]);

            const vectorField
                spCf(subPatch.faceCentres()),
                spSf(subPatch.faceAreas()),
                spNf(subPatch.faceAreas()/mag(subPatch.faceAreas())),
                mpCf(mainPatch.faceCentres()),
                mpSf(mainPatch.faceAreas()),
                mpNf(mainPatch.faceAreas()/mag(mainPatch.faceAreas()));

            scalarField& labelPatch(labelField.boundaryField()[sPatI]);
            labelPatch.resize(subPatch.size());
            forAll(subPatch, sFacI)
            {
                if (sPatI == interPatchID())
                {
                    const label direc(interFaceInfo()[5][sFacI]);
                    assert(fabs(direc) - 1.0 < SMALL);
                    labelPatch[sFacI] = // Add SMALL to support label = 0
                        direc * (interFaceInfo()[0][sFacI] + SMALL);
                }
                else
                {
                    labelPatch[sFacI] = GREAT;
                    forAll(mainPatch, mFacI)
                    {
                        // Same center and size - perfect overlap
                        if ((mag(spCf[sFacI] - mpCf[mFacI]) < SMALL) and
                            ((mag(spSf[sFacI] - mpSf[mFacI]) < SMALL) or
                             (mag(spSf[sFacI] + mpSf[mFacI]) < SMALL)))
                        {
                            labelPatch[sFacI] = mFacI;
                            break;
                        }
                    }
                    if (labelPatch[sFacI] == GREAT)
                    {
                        scalar minDist(1e30);
                        forAll(mainPatch, mFacI)
                        {
                            // Same normal and orthogonal dist gives same plane
                            // Loop to find min dist, may find neighbour!!!
                            const vector Dist(mpCf[mFacI] - spCf[sFacI]);
                            if ((mag(mpSf[mFacI] ^ spSf[sFacI]) < SMALL) and
                                mag(Dist & mpSf[mFacI]) < SMALL)
                            {
                                if (mag(Dist) < minDist)
                                {
                                    minDist = mag(Dist);
                                    labelPatch[sFacI] = mFacI;
                                }
                            }
                        }
                    }
                }
                assert(labelPatch[sFacI] < GREAT);
            }
        }
        subFace2MainFace.insert(numPatchName_, labelField);
    }
    return subFace2MainFace[numPatchName_];
}


const scalarField Foam::subSolver::weights() const
{
    // Cf OpenFOAM surfaceInterpolation.C
    static HashTable<scalarField, word, string::hash> yRatio;
    if (not yRatio.found(numPatchName_))
    {
        const polyPatch
            &interPatch(subMesh_->boundary()[interPatchID()].patch());
        const vectorField Cf(interPatch.faceCentres());
        const vectorField Sf(interPatch.faceAreas());
        scalarField yRat(interPatch.size());
        forAll(interPatch, faceI)
        {
            const label
                sub(interPatch.faceCells()[faceI]),
                nxt(interFaceInfo()[3][faceI]);
            const scalar
                Y(mag(Sf[faceI] & (mainMesh_.C()[nxt]-Cf[faceI]))),
                y(mag(Sf[faceI] & (subMesh_->C()[sub]-Cf[faceI])));
            yRat[faceI] = y / (y + Y);
        }
        yRatio.insert(numPatchName_, yRat);
    }
    return yRatio[numPatchName_];
}


tmp<vectorField> Foam::subSolver::faceFlux(const word fieldName) const
{
    // OBSERVE: turbulent properties in general subSolver, OK?
    const incompressible::RASModel& mainTurb =
        mainMesh_.thisDb().lookupObject<incompressible::RASModel>
        (
            "turbulenceModel"
        );
    const incompressible::RASModel& subTurb =
        subMesh_->thisDb().lookupObject<incompressible::RASModel>
        (
            "turbulenceModel"
        );

    const scalarField
        w(weights()),
        mainDeff(mainTurb.DphiEff(fieldName)),
        subDeff(subTurb.DphiEff(fieldName));

    // Place diffusion in x, advection in y (upper/sub) and z (lower/main)
    tmp<vectorField> tFlux(new vectorField(w.size(), pTraits<vector>::zero));
    vectorField &flux(tFlux());

    const surfaceScalarField
//        &sPhi(subMesh_->thisDb().lookupObject<surfaceScalarField>("phi")),
        &mPhi(mainMesh_.thisDb().lookupObject<surfaceScalarField>("phi"));
//    const scalarField &sPhiP(sPhi.boundaryField()[interPatchID()]);

    const fvPatch &interPatch(subMesh_->boundary()[interPatchID()]);
    const labelList& faceCells(interPatch.faceCells());
    const labelListList& interFaceInfo_(interFaceInfo());
    const scalarField faceArea(interPatch.magSf());
    const vectorField
        &mainC(mainMesh_.C()),
        &subC(subMesh_->C());
    forAll(interPatch, sFacI)
    {
        const label
            sCelI(faceCells[sFacI]),
            mFacI(interFaceInfo_[0][sFacI]),
//            direc(interFaceInfo_[5][sFacI]),
            nextI(interFaceInfo_[3][sFacI]);
        flux[sFacI].x() =
            // FIXME: Add cross-terms and non-orthogonal part to diffusion,
            // see Jasak's PhD p84, and divDevReff in any turbulence model
              ((w[sFacI])*mainDeff[nextI] + (1-w[sFacI])*subDeff[sCelI])
            * faceArea[sFacI] / mag(mainC[nextI] - subC[sCelI]);
        //
        // Sf:  a positive surface-normal vector points away from the owner
        //      cell and towards the neighbour cell,
        //      it always points out of a region.
        // phi: That means for positive flux owner is upstream cell.
        //
        // In case of no perfect overlap phi needs to be scaled with surface.
        // For CD we use w, i.e. w = Cs-Cf / (Cs-Cf - Cm-Cf)

        // FIXME: Remove
        /*
        scalar phi(0);
        // Use Sf from main to be consistent with signs, ie direc*Sf_sub
        switch(int(fct_))
        {
          case 0:
            phi = mPhi[mFacI];
          case 1:
            phi = direc*sPhiP[sFacI];
          case 2:
            phi = w[sFacI] * mPhi[mFacI] + (1-w[sFacI]) * direc*sPhiP[sFacI];
        }
        */
        const scalar
            phi(mPhi[mFacI]), // Make sure phi = mPhi at assembling of main
            lAdvCD(-w[sFacI] * phi),      // lower
            uAdvCD((1 - w[sFacI]) * phi), // upper
            lAdvUD(min(0, phi)),
            uAdvUD(max(0, phi)),
            //f(0.5*(abs(fct_)-1));
            //f(0.5*fct_);
            f(1);
        flux[sFacI].y() = (f*uAdvCD + (1-f)*uAdvUD);
        flux[sFacI].z() = (f*lAdvCD + (1-f)*lAdvUD);
    }
    return tFlux;
}


void Foam::subSolver::plotSrcP(Foam::Field<scalar>& src)
{
    const label N(5);
    scalarField
        maxFct(N-2, 0),
        minFct(N-2, 0);
    vectorField
        maxCf(N-2, pTraits<vector>::zero),
        minCf(N-2, pTraits<vector>::zero);
    static List<scalarField> cache(N);
    cache[0] = src;
    const labelListList& interFaceInfo_(interFaceInfo());
    const fvPatch &interPatch(subMesh_->boundary()[interPatchID()]);
    forAll(interPatch, sFacI)
    //forAll(ortoFaces(), indx)
    {
        const label
    //        sFacI(ortoFaces()[indx]),
            mCelI(interFaceInfo_[2][sFacI]);
        for(label i=2; i<N; i++)
        {
            if (not cache[i].empty())
            {
                const scalar diff(cache[0][mCelI] - cache[i][mCelI]);
                if (maxFct[i-2] < diff)
                {
                    maxCf[i-2] = interPatch.Cf()[sFacI];
                    maxFct[i-2] = diff;
                }
                else
                if (minFct[i-2] > diff)
                {
                    minCf[i-2] = interPatch.Cf()[sFacI];
                    minFct[i-2] = diff;
                }
            }
        }
    }
    Info << "Cf, mxDifSrcP: " << maxCf << ", " << maxFct << endl;
    Info << "Cf, mnDifSrcP: " << minCf << ", " << minFct << endl;
    for(label i=N-1; i>0; i--)
    {
        cache[i] = cache[i-1];
    }
}


Foam::tmp<vectorField>
Foam::subSolver::projectOnPatch(const vectorField& field, const label patchID)
{
    Foam::tmp<vectorField> tParallalField
        (
            new vectorField(field)
        );
    vectorField& paraField(tParallalField());

    const labelList& mainFaces(wallFacesToMainFaces());
    const vectorField nf(subMesh_->boundary()[patchID].nf());
    forAll(nf, sFacI)
    {
        const label mFacI(mainFaces[sFacI]);
        paraField[mFacI] -= nf[sFacI] * (nf[sFacI] & field[mFacI]);
    }
    return tParallalField;
}


Foam::tmp<vectorField>
Foam::subSolver::projectOnWall(const vectorField& field, const label patchID)
{
    Foam::tmp<vectorField> tParallalField
        (
            new vectorField(field)
        );
    vectorField& paraField(tParallalField());

    const vectorField nf(mainMesh_.boundary()[patchID].nf());

    forAll(nf, faceI)
    {
        paraField[faceI] -= nf[faceI] * (nf[faceI] & field[faceI]);
    }
    return tParallalField;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::subSolver::subSolver
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fct_(dict.lookupOrDefault<scalar>("factor", 1.)),
    numPatchName_(dict.lookup("wallPatch")),
    intrfcPatchName_(dict.lookupOrDefault<word>("interfacePatch",
                                                "oldInternalFaces")),
    wallPatchID_(mesh.boundaryMesh().findPatchID(numPatchName_)),
    mainMesh_(mesh)
{
    subSolver::init();
}


Foam::subSolver::subSolver
(
    const fvMesh& mesh
)
:
    fct_(1.),
    numPatchName_(""),
    intrfcPatchName_(""),
    wallPatchID_(0),
    mainMesh_(mesh)
{
    subSolver::init();
}


Foam::subSolver::subSolver
(
    const subSolver& ptf
)
:
    fct_(ptf.fct_),
    numPatchName_(ptf.numPatchName_),
    intrfcPatchName_(ptf.intrfcPatchName_),
    wallPatchID_(ptf.wallPatchID_),
    mainMesh_(ptf.mainMesh_)
{}


Foam::subSolver::~subSolver()
{
    subMesh_ = NULL;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
