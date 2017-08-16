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
#include <sys/stat.h>
#include "icoSsSolver.H"

#include "subCycle.H"
#include "fixedGradientFvPatchField.H"
#include "singlePhaseTransportModel.H"
#include "incompressible/RAS/RASModel/RASModel.H"
#include "simpleControl.H"
#include "IOdictionary.H"


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


void Foam::icoSsSolver::init()
{
    // FIXME: Only needed for UMIST - if Ux exists run in UMIST mode
    struct stat buffer;
    fileName
        tPath(mainMesh_.time().path()+"/"+mainMesh_.time().timeName()),
        mPath(tPath + "/Ux"),
        sPath(tPath + "/" + numPatchName_ + "/Ux");
    if (stat (mPath.toAbsolute().c_str(), &buffer) == 0)
    {
        UxMain_.clear();
        UxMain_.set(new volScalarField(fieldIOobject("Ux", mainMesh_),
                                       mainMesh_));
    }
    if (stat (sPath.toAbsolute().c_str(), &buffer) == 0)
    {
        UxPatch_.clear();
        UxPatch_.set(new volScalarField(fieldIOobject("Ux",*subMesh_),
                                        *subMesh_));
    }
}


void Foam::icoSsSolver::update()
{
    if (mapP_)
    {
        // Scale sub phi from main phi before assembling matrice.
        updatePhi();
    }
    else
    {
        // Set Dirichlet condition at interface for pressure.
        setInterfaceValue(pPatch_());
    }
}


void Foam::icoSsSolver::updatePhi()
{
    // Scaling of phi from phiMain to ensure continuity:
    // - First update sub phi from sub U to acknowledge internal changes
    phiPatch_() = (fvc::interpolate(Upatch_()) & subMesh_->Sf());

    // - Go through all overlaping internal and boundary faces and update phi
    const surfaceScalarField&
        phiMain(mainMesh_.thisDb().lookupObject<surfaceScalarField>("phi"));

    const List<List<labelPairList> >&
        mainFaceToSubFaces_(mainFaceToSubFaces());
    forAll(mainFaceToSubFaces_, sPtcI)
    {
        const scalarField *mPhi(NULL);
        scalarField *sPhi(NULL);
        vectorField mNf, sNf;
        if (sPtcI == subMesh_->boundaryMesh().size())   // internal faces
        {
            mPhi = &phiMain.internalField();
            sPhi = &phiPatch_->internalField();
            mNf = mainMesh_.Sf()/mag(mainMesh_.Sf());
            sNf = subMesh_->Sf()/mag(subMesh_->Sf());
        }
        else
        if (sPtcI == interPatchID())                    // inter patch faces
        {
            mPhi = &phiMain.internalField();
            sPhi = &phiPatch_->boundaryField()[interPatchID()];
            mNf = mainMesh_.Sf()/mag(mainMesh_.Sf());
            const vectorField
                sSf(subMesh_->boundaryMesh()[interPatchID()].faceAreas());
            sNf = sSf/mag(sSf);
        }
        else                                            // boundary faces
        {
            const word patchName(subMesh_->boundaryMesh().names()[sPtcI]);
            const label mPtcI(mainMesh_.boundaryMesh().findPatchID(patchName));
            mPhi = &phiMain.boundaryField()[mPtcI];
            sPhi = &phiPatch_->boundaryField()[sPtcI];
            const vectorField
                mSf(mainMesh_.boundaryMesh()[mPtcI].faceAreas()),
                sSf(subMesh_->boundaryMesh()[sPtcI].faceAreas());
            mNf = mSf/mag(mSf);
            sNf = sSf/mag(sSf);
        }

        /* Scale all overlapping sub faces from main
          OBSERVE: Indx equals (copy and paste from subSolver.C):
          1. sFacI for "the wall" boundary, wallFacesToMainFaces (one-to-one)
          2. index of orto/paraFaces() for interFace patch       (one-to-one)
          3. mFacI for over-lapping internal sub faces          (one-to-many)
          4. sFacI according to orto/paraFaces() for other bnds (one-to-many)
          FIXME: This loop goes through at least some faces twice
        */
        forAll(mainFaceToSubFaces_[sPtcI], indx)
        {
            if (sPhi->empty()) break;   // Boundaries not relevant for sub

            // Calculate total mass flux from sub cells per main face, mFacI
            scalar sPhiTot(0);
            forAll(mainFaceToSubFaces_[sPtcI][indx], jndx)
            {
                const label
                    mFacI(mainFaceToSubFaces_[sPtcI][indx][jndx].first()),
                    sFacI(mainFaceToSubFaces_[sPtcI][indx][jndx].second());
                const vector UxMain(mNf[mFacI] * (*mPhi)[mFacI]);
                if (sPtcI==interPatchID() and jndx==1) // first of ortoFaces()
                {
                    // Skip para (jndx==0), only sum up phi from orto faces
                    sPhiTot = 0;
                }
                // Assert wall parallel mass flow has the same sign as main.
                // This is important for UMIST but is valid in general.
                if ((UxMain & (sNf[sFacI]*(*sPhi)[sFacI])) < 0)
                {
                    // Reset flux in wrong direction, needed for UMIST (NWS)
                    (*sPhi)[sFacI] = 0;   // UMIST with NFF does not converge
                }
                sPhiTot += (mNf[mFacI] & sNf[sFacI]) * (*sPhi)[sFacI];

            }
            if (sPhiTot != 0.)
            {
                forAll(mainFaceToSubFaces_[sPtcI][indx], jndx)
                {
                    const label
                        mFacI(mainFaceToSubFaces_[sPtcI][indx][jndx].first()),
                        sFacI(mainFaceToSubFaces_[sPtcI][indx][jndx].second());
                    if (noAdv_)
                    {
                        (*sPhi)[sFacI] = 0;
                    }
                    else
                    if (sPtcI==interPatchID())  // inter face
                    {
                        (*sPhi)[sFacI] =  // For inter phiTot is not needed
                            (mNf[mFacI] & sNf[sFacI]) * (*mPhi)[mFacI];
                    }
                    else
                    {
                        (*sPhi)[sFacI] *= (*mPhi)[mFacI] / sPhiTot;
                    }
                }
            }
        }
    }

    // - Update all non overlapping sub faces for divergence free sub cells
    // FIXME: For compression acknowledge non-zero divPhi for main.
    const polyPatch
        &interPatch(subMesh_->boundary()[interPatchID()].patch()),
        &wallPatch(subMesh_->boundary()[wallPatchID()].patch());
    const labelList
        &faceCells(interPatch.faceCells()),
        &wallCells(wallPatch.faceCells()),
        &ortoFaces_(ortoFaces());

    volScalarField curDivPhi(fvc::surfaceIntegrate(phiPatch_()));
    const scalarField
        &curDivPhiInt(curDivPhi.internalField()),
        &subV(subMesh_->V());

    const List<labelPairList>&
        wallFaceToSubCellFaces_(wallFaceToSubCellFaces());
    forAll(wallFaceToSubCellFaces_, wFacI)
    {
        // Omitting orto colon has only a marginal effect on phi
        const label wCelI(wallCells[wFacI]);
        bool ortoWallCell(false);
        forAll(ortoFaces_, indx)
        {
            const label oCelI(faceCells[ortoFaces_[indx]]);
            if (wCelI == oCelI)
            {
                ortoWallCell = true;
                break;
            }
        }
        if (ortoWallCell)
        {
            continue;   // Orto cells is mapped one-to-one and divergence free
        }
        scalar curNetFlow(0);
        forAll(wallFaceToSubCellFaces_[wFacI], indx)
        {
            const label
                sCelI(wallFaceToSubCellFaces_[wFacI][indx].first()),
                sFacI(wallFaceToSubCellFaces_[wFacI][indx].second());
            // Update divPhi for current cell from former phi update
            curDivPhi[sCelI] += curNetFlow / subV[sCelI];
            curNetFlow = curDivPhiInt[sCelI] * subV[sCelI];
            if (subMesh_->boundaryMesh().whichPatch(sFacI) < 0) // Internal
            {
                const label sgn((sCelI == subMesh_->owner()[sFacI]) ? 1 : -1);
                phiPatch_()[sFacI] -= sgn * curNetFlow;
            }
            else  // Must be face at interface, find local from global index
            {
                const label pFcI(interPatch.whichFace(sFacI));
                phiPatch_->boundaryField()[interPatchID()][pFcI] -= curNetFlow;
            }
        }
    }
}


Foam::tmp<volVectorField> Foam::icoSsSolver::interpolateGradP()
{
    // Mapping grad(pMain) directly avoids trouble with bnd of subP
    const volScalarField&
        pMain(mainMesh_.thisDb().lookupObject<volScalarField>("p"));
    const volVectorField gradP(fvc::grad(pMain));
    const surfaceVectorField surfGradP(fvc::interpolate(gradP));
    const scalarField
        &mainV(mainMesh_.V()),
        &subV(subMesh_->V());
    const vectorField
        &wallGradP(gradP.boundaryField()[wallPatchID_]),
        &mainWallCf(mainMesh_.boundary()[wallPatchID_].Cf()),
        &mainCf(mainMesh_.Cf()),
        &mainC(mainMesh_.C()),
        &subC(subMesh_->C());

    Foam::tmp<volVectorField> tSubGradP(new volVectorField(Upatch_));
    tSubGradP->dimensions().reset(tSubGradP->dimensions()/dimTime);
    Foam::Field<vector>& subGradP = tSubGradP();
    volScalarField& subP(pPatch_());
    bool multiMainOrtoCells(false);

    // Use main para and wall face values to interpolate to sub field
    const labelListList& interFaceInfo_(interFaceInfo());
    const labelListList& mainWallFacesToSubCells_(mainWallFaceToSubCells());
    const labelList
        &paraFaces_(paraFaces()),
        &ortoFaces_(ortoFaces());
    forAll(paraFaces_, indx)
    {
        const label
            sFacI(paraFaces_[indx]),
            pFacI(interFaceInfo_[0][sFacI]),
            wFacI(interFaceInfo_[1][sFacI]),
            mCelI(interFaceInfo_[2][sFacI]);
        if (wFacI<0) // multiple main orto cells (overlapping sub grid)
        {
            multiMainOrtoCells = true;
            continue;
        }
        assert(not mainWallFacesToSubCells_[wFacI].empty());
        const scalar
            distCP(mag(mainC[mCelI] - mainCf[pFacI])),
            distCW(mag(mainC[mCelI] - mainWallCf[wFacI]));
        forAll(mainWallFacesToSubCells_[wFacI], jndx)
        {
            const label sCelI(mainWallFacesToSubCells_[wFacI][jndx]);
            const vector cC(subC[sCelI] - mainC[mCelI]);

            if (true) // Constant pressure gradient
            {
                subGradP[sCelI] = gradP[mCelI];
            }
            else      // Piecewise linear pressure gradient
            {
                const vector
                    cP(subC[sCelI] - mainCf[pFacI]),
                    cW(subC[sCelI] - mainWallCf[wFacI]);
                const scalar
                    wP(mag(cP)/distCP),
                    wW(mag(cW)/distCW);
                subGradP[sCelI] = (wP < wW) ?
                    (wP*gradP[mCelI] + (1 - wP)*surfGradP[pFacI]) :
                    (wW*gradP[mCelI] + (1 - wW)*wallGradP[wFacI]);
            }
            subP[sCelI] = pMain[mCelI] + relaxGradP_ * (cC & subGradP[sCelI]);
        }
    }

    const fvPatch& interPatch(subMesh_->boundary()[interPatchID()]);
    if (multiMainOrtoCells)
    {
      Info << "MAP GRAD P FOR ORTO" << endl;
        const scalar SML(1e-9);
        forAll(ortoFaces_, indx)
        {
            const label
                sFacI(ortoFaces_[indx]),
                mCelI(interFaceInfo_[2][sFacI]),
                sCelI(interPatch.faceCells()[sFacI]);
            // Assume perfect one-to-one overlap between main and sub cells
            assert(mag(mainC[mCelI]-subC[sCelI])/pow(subV[sCelI], 1./3.)<SML);
            assert((mainV[mCelI] - subV[sCelI])/subV[sCelI] < SML);
            subGradP[sCelI] = gradP[mCelI];
            subP[sCelI] = pMain[mCelI];
        }
    }
    return tSubGradP;
}


Foam::tmp<vectorField> Foam::icoSsSolver::solveStreamwiseCmpt()
{
    const polyPatch &wallPatch(mainMesh_.boundary()[wallPatchID_].patch());
    const vectorField &wallSf(wallPatch.faceAreas());
    const vector  // Assume 2D with empty boundaries in y-direction
        unitX(1,0,0),
        unitY(0,1,0),
        unitZ(0,0,1),
        streamCmpt(cmptMag(1./mag(wallSf[0])*wallSf[0] ^ unitY));
    assert((mag(streamCmpt-unitX) < 1e-12) or (mag(streamCmpt-unitZ) < 1e-12));

    // Transfer Umain to UxMain, needed to transfer face flux from main
    const volVectorField&
        mainU(mainMesh_.thisDb().lookupObject<volVectorField>("U"));
    // FIXME: Maybe use component instead for scalar product
    const vectorField &mU(mainU.internalField());
    scalarField &mUx(UxMain_->internalField());
    mUx = (streamCmpt & mU);

    tmp<fvScalarMatrix> UEqn
    (
        fvm::ddt(UxPatch_())
      + fvm::div(phiPatch_(), UxPatch_())
      - fvm::laplacian(turbMdl_->nuEff(), UxPatch_())
      ==
        fvOptions_()(UxPatch_())
    );
    addParaFlux(UEqn());  // Better here than in fvOptions

    UEqn->relax();

    fvOptions_->constrain(UEqn());

    // FIXME: add assumeZeroGradP as internal variable
//    bool assumeZeroGradP(false);
    bool assumeZeroGradP(fct_==0);
    if (assumeZeroGradP)
    {
        Foam::solve(UEqn());
    }
    else
    {
        Foam::solve(UEqn() == - (streamCmpt & interpolateGradP()) );
    }
    return streamCmpt * UxPatch_->internalField();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::icoSsSolver::icoSsSolver
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    subSolver(mesh, dict),
    mapP_(dict.lookupOrDefault<bool>("mapPressure", true)),
    noAdv_(dict.lookupOrDefault<bool>("noAdvection", false)),
    relaxGradP_(dict.lookupOrDefault<scalar>("relaxGradP", 0.9)),
    nSubCycles_(dict.lookupOrDefault<label>("nSubCycles", 1)),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    pPatch_(new volScalarField(fieldIOobject("p", *subMesh_), *subMesh_)),
    Upatch_(new volVectorField(fieldIOobject("U", *subMesh_), *subMesh_)),
    phiPatch_
    (
        new surfaceScalarField
        (
            fieldIOobject(phiName_, *subMesh_, IOobject::NO_READ),
            linearInterpolate(Upatch_()) & subMesh_->Sf()
        )
    ),
    transMdl_(new singlePhaseTransportModel(Upatch_, phiPatch_)),
    turbMdl_
    (
        incompressible::RASModel::New
        (
            Upatch_(),
            phiPatch_(),
            transMdl_()
        ).ptr()
    ),
    fvOptionsIO_
    (
        new IOdictionary
        (
            IOobject
            (
                "fvOptions",
                subMesh_->time().system(),
                *subMesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        )
    ),
    fvOptions_(new fv::optionList(*subMesh_, fvOptionsIO_()))
{
    init();
}


Foam::icoSsSolver::icoSsSolver
(
    const fvMesh& mesh
)
:
    subSolver(mesh),
    mapP_(true),
    noAdv_(false),
    relaxGradP_(0.9),
    nSubCycles_(1),
    phiName_("phi"),
    pPatch_(new volScalarField(fieldIOobject("p", *subMesh_), *subMesh_)),
    Upatch_(new volVectorField(fieldIOobject("U", *subMesh_), *subMesh_)),
    phiPatch_
    (
        new surfaceScalarField
        (
            fieldIOobject(phiName_, *subMesh_, IOobject::READ_IF_PRESENT),
            linearInterpolate(Upatch_()) & subMesh_->Sf()
        )
    ),
    transMdl_(new singlePhaseTransportModel(Upatch_, phiPatch_)),
    turbMdl_
    (
        incompressible::RASModel::New
        (
            Upatch_,
            phiPatch_,
            transMdl_()
        ).ptr()
    ),
    fvOptionsIO_
    (
        new IOdictionary
        (
            IOobject
            (
                "fvOptions",
                subMesh_->time().system(),
                *subMesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        )
    ),
    fvOptions_(new fv::optionList(*subMesh_, fvOptionsIO_()))
{
    init();
}


Foam::icoSsSolver::icoSsSolver
(
    const icoSsSolver& ptf
)
:
    subSolver(ptf),
    mapP_(ptf.mapP_),
    noAdv_(ptf.noAdv_),
    relaxGradP_(ptf.relaxGradP_),
    nSubCycles_(ptf.nSubCycles_),
    phiName_(ptf.phiName_),
    pPatch_(ptf.pPatch_),
    Upatch_(ptf.Upatch_),
    phiPatch_(ptf.phiPatch_),
    transMdl_(ptf.transMdl_),
    turbMdl_(ptf.turbMdl_),
    fvOptionsIO_(ptf.fvOptionsIO_),
    fvOptions_(ptf.fvOptions_)
{}


Foam::icoSsSolver::~icoSsSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::icoSsSolver::solve()
{
    static HashTable<unsigned, word, string::hash> hashIter;
    if (not hashIter.found(numPatchName_))
    {
        hashIter.insert(numPatchName_, 0);
    }
    // Ensure main has runned as default a couple of time steps
    if (hashIter[numPatchName_] < 3)
    {
        const volVectorField
            &Umain(mainMesh_.thisDb().lookupObject<volVectorField>("U"));
        const polyPatch
            &wallPatch(mainMesh_.boundary()[wallPatchID_].patch());
        const labelListList &subCells(mainWallFaceToSubCells());
        const labelList &wallCells(wallPatch.faceCells());
        forAll(wallPatch, wFacI)
        {
            const label mCelI(wallCells[wFacI]);
            forAll(subCells[wFacI], indx)
            {
                const label sCelI(subCells[wFacI][indx]);
                Upatch_()[sCelI] = Umain[mCelI];
            }
        }
        Upatch_->correctBoundaryConditions();
        phiPatch_() = (fvc::interpolate(Upatch_()) & subMesh_->Sf());
    }
    if (hashIter[numPatchName_]++ % nSubCycles_) return;

    // FIXME: add solveOnlyStreamwise as internal variable
    bool solveOnlyStreamwise(UxPatch_.valid() and (hashIter[numPatchName_]>1));
    if (solveOnlyStreamwise)
    {
        Upatch_().internalField() = solveStreamwiseCmpt();
    }

    update();

    if (solveOnlyStreamwise)
    {
        Upatch_().internalField() =
            fvc::reconstruct(phiPatch_())->internalField();
    }
    else
    {
        tmp<fvVectorMatrix> UEqn
        (
            fvm::ddt(Upatch_())
          + fvm::div(phiPatch_(), Upatch_())
          + turbMdl_->divDevReff(Upatch_())
          ==
            fvOptions_()(Upatch_())
        );
        addParaFlux(UEqn());  // Better here than in fvOptions

        UEqn->relax();

        fvOptions_->constrain(UEqn());

        // FIXME: add assumeZeroGradP as internal variable
//        bool assumeZeroGradP(false);
        bool assumeZeroGradP(fct_==0);
        if (mapP_ and assumeZeroGradP)
        {
            Foam::solve(UEqn());
        }
        else
        if (mapP_)
        {
            Foam::solve(UEqn() == - interpolateGradP());
        }
        else  // solve P
        {
            if (fct_) interpolFaceValue(pPatch_(), true);
            else interpolFaceValue(pPatch_(), false, true);

            Foam::solve(UEqn() == - fvc::grad(pPatch_()));

            volScalarField rAU(1.0/UEqn().A());
            volVectorField HbyA("HbyA", Upatch_());
            HbyA = rAU*UEqn().H();
            UEqn.clear();

            surfaceScalarField
                phiHbyA("phiHbyA", fvc::interpolate(HbyA) & subMesh_->Sf());
            adjustPhi(phiHbyA, Upatch_(), pPatch_());
            fvOptions_->makeRelative(phiHbyA);

            Info<< "Assembling p and correcting phi with pressure flux."
                << endl;
            fvScalarMatrix pEqn
            (
                fvm::laplacian(rAU, pPatch_()) == fvc::div(phiHbyA)
            );
            pPatch_->storePrevIter();
            pEqn.setReference(0, 0.0);

            pEqn.solve();

            // pEqn.flux() non-zero even if not solving but hardly affect p,T,U
            phiPatch_() = phiHbyA - pEqn.flux();

            // Momentum corrector
            pPatch_->relax(); // Only meaningful after p is solved
            Upatch_() = HbyA
                - rAU*fvc::reconstruct(fvc::snGrad(pPatch_())*subMesh_->magSf());
        }
    }

    Upatch_->correctBoundaryConditions();
    fvOptions_->correct(Upatch_());   // Does nothing

    // Need to get correct gradU when calculating G=nut*magSqr(fvc::grad(U))
    interpolFaceValue(Upatch_());
    turbMdl_->correct();
    Upatch_().correctBoundaryConditions();  // Homogeneous Neumann

//    runTime.write();
//    #include "initContinuityErrs.H"
}


// ************************************************************************* //
