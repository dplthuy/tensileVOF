/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
              Modified work | 2019 Andy Spitzenberger
              Modified work | 2022 Dennis Thuy
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
#include "addToRunTimeSelectionTable.H"
#include "binned.H"
#include "coupledPolyPatch.H"
#include "dynamicMultiDimRefineFvMesh.H"
#include "dynamicFvMesh.H"
#include "emptyPolyPatch.H"
#include "extractEulerianParticlesLBAMR.H"
#include "mathematicalConstants.H"
#include "objectMap.H"
#include "pairPatchAgglomeration.H"
#include "regionSplit2Dv1906.H"
#include "surfaceFields.H"
#include "surfaceInterpolate.H"
#include "volFields.H"

#include <queue>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(extractEulerianParticlesLBAMR, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        extractEulerianParticlesLBAMR,
        dictionary
    );
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
void Foam::functionObjects::extractEulerianParticlesLBAMR::createFiles()
{
    if (writeToFile() && !extractEulerianFilePtr_)
    {
    	extractEulerianFilePtr_ = createFile("eulerianParticleData_" + faceZoneName_);
    	writeDataFileHeader(extractEulerianFilePtr_());
    }
}

void Foam::functionObjects::extractEulerianParticlesLBAMR::writeDataFileHeader
(
    Ostream& os
) const
{
    writeHeader(os, "extracted Eulerian Particles");
    writeHeaderValue(os, "FaceZone", faceZoneName_);
    writeHeaderValue(os, "Extracted field", alphaName_);
    writeHeaderValue(os, "Alpha threshold", alphaThreshold_);
    writeHeaderValue(os, "Minimum diameter", minDiameter_);
    writeHeaderValue(os, "Maximum diameter", maxDiameter_);
    writeHeader(os, "");
    writeCommented(os, "ID");
    writeTabbed(os, "x");
    writeTabbed(os, "y");
    writeTabbed(os, "z");
    writeTabbed(os, "time");
    writeTabbed(os, "Ux");
    writeTabbed(os, "Uy");
    writeTabbed(os, "Uz");
    writeTabbed(os, "D");
    writeTabbed(os, "V");
    
    os << endl;
}

void Foam::functionObjects::extractEulerianParticlesLBAMR::checkFaceZone()
{
    DebugInFunction << endl;
    
    if(!isA<dynamicFvMesh>(mesh_))
    {
        FatalErrorInFunction
            << "Mesh is  not of type dynamicFvMesh, functionObject extractEulerianParticlesLBAMR cannot be used"
            << exit(FatalError);
    }

    zoneID_ = mesh_.faceZones().findZoneID(faceZoneName_);
    if (zoneID_ == -1)
    {
        FatalErrorInFunction
            << "Unable to find faceZone " << faceZoneName_
            << ".  Available faceZones are: " << mesh_.faceZones().names()
            << exit(FatalError);
    }

    const faceZone& fz = mesh_.faceZones()[zoneID_];
    const label nFaces = fz.size();
    const label allFaces = returnReduce(nFaces, sumOp<label>());

    if (allFaces < nInjectorLocations_)
    {
        FatalErrorInFunction
            << "faceZone " << faceZoneName_
            << ": Number of faceZone faces (" << allFaces
            << ") is less than the number of requested locations ("
            << nInjectorLocations_ << ")."
            << exit(FatalError);
    }

    Info<< type() << " " << name() << " output:" << nl
        << "    faceZone : " << faceZoneName_ << nl
        << "    faces    : " << allFaces << nl
        << endl;

    // Initialise old iteration 
    // Note: for restart, this info needs to be written/read
    regions0_.setSize(fz.size(), -1);
    faceZone0_.setSize(fz.size(), -1);
    
    forAll(fz, localFacei)
    {
    	faceZone0_[localFacei] = fz[localFacei];
    }
}


void Foam::functionObjects::extractEulerianParticlesLBAMR::initialiseBins()
{
    DebugInFunction << endl;

    if (!nInjectorLocations_)
    {
        return;
    }

    const faceZone& fz = mesh_.faceZones()[zoneID_];

    // Agglomerate faceZone faces into nInjectorLocations_ global locations
    const indirectPrimitivePatch patch
    (
        IndirectList<face>(mesh_.faces(), fz),
        mesh_.points()
    );

    const label nFaces = fz.size();
    label nLocations = nInjectorLocations_;

    if (Pstream::parRun())
    {
        label nGlobalFaces = returnReduce(nFaces, sumOp<label>());
        scalar fraction = scalar(nFaces)/scalar(nGlobalFaces);
        nLocations = ceil(fraction*nInjectorLocations_);
        if (debug)
        {
            Pout<< "nFaces:" << nFaces
                << ", nGlobalFaces:" << nGlobalFaces
                << ", fraction:" << fraction
                << ", nLocations:" << nLocations
                << endl;
        }
    }

    pairPatchAgglomeration ppa
    (
        patch.localFaces(),
        patch.localPoints(),
        10,
        50,
        nLocations,
        labelMax,
        180
    );

    ppa.agglomerate();

    label nCoarseFaces = 0;
    if (nFaces != 0)
    {
        fineToCoarseAddr_ = ppa.restrictTopBottomAddressing();
        nCoarseFaces = max(fineToCoarseAddr_) + 1;
    }

    globalCoarseFaces_ = globalIndex(nCoarseFaces);

    Info<< "Created " << returnReduce(nCoarseFaces, sumOp<label>())
        << " coarse faces" << endl;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::functionObjects::extractEulerianParticlesLBAMR::phiU() const
{
    DebugInFunction << endl;

    const surfaceScalarField& phi
    (
        mesh_.lookupObject<surfaceScalarField>(phiName_)
    );

    if (phi.dimensions() == dimMass/dimTime)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        return phi/fvc::interpolate(rho);
    }

    return phi;
}


void Foam::functionObjects::extractEulerianParticlesLBAMR::setBlockedFaces
(
    const surfaceScalarField& alphaf,
    const faceZone& fz,
    boolList& blockedFaces
)
{
    DebugInFunction << endl;

    // Initialise storage for patch and patch-face indices where faceZone
    // intersects mesh patch(es)
    patchIDs_.setSize(fz.size(), -1);
    patchFaceIDs_.setSize(fz.size(), -1);

    label nBlockedFaces = 0;
    forAll(fz, localFacei)
    {
        const label facei = fz[localFacei];

        if (mesh_.isInternalFace(facei))
        {
            if (alphaf[facei] > alphaThreshold_)
            {
                blockedFaces[localFacei] = true;
            }
        }
        else
        {
            label patchi = mesh_.boundaryMesh().whichPatch(facei);
            label patchFacei = -1;

            const polyPatch& pp = mesh_.boundaryMesh()[patchi];
            const scalarField& alphafp = alphaf.boundaryField()[patchi];

            if (isA<coupledPolyPatch>(pp))
            {
                const coupledPolyPatch& cpp =
                    refCast<const coupledPolyPatch>(pp);

                if (cpp.owner())
                {
                    patchFacei = cpp.whichFace(facei);
                }
            }
            else if (!isA<emptyPolyPatch>(pp))
            {
                patchFacei = pp.whichFace(facei);
            }

            if (patchFacei == -1)
            {
                patchi = -1;
            }
            else if (alphafp[patchFacei] > alphaThreshold_)
            {
                blockedFaces[localFacei] = true;
            }

            patchIDs_[localFacei] = patchi;
            patchFaceIDs_[localFacei] = patchFacei;
        }
    }

    DebugInFunction << "Number of blocked faces: " << nBlockedFaces << endl;
}

void Foam::functionObjects::extractEulerianParticlesLBAMR::collectParticles
(
    const scalar time,
    const boolList& collectParticle
)
{
    DebugInFunction << "collectParticle: " << collectParticle << endl;

    // Collect particles on local processors that have passed through faceZone
    forAll(collectParticle, particlei)
    {
        if (!collectParticle[particlei])
        {
            continue;
        }

        eulerianParticleMod p = particles_[particlei];

        if (p.faceIHit != -1 && nInjectorLocations_)
        {
            // Use coarse face index for tag output
            label coarseFacei = fineToCoarseAddr_[p.faceIHit];
            p.faceIHit = globalCoarseFaces_.toGlobal(coarseFacei);
        }

        reduce(p, sumParticleOp<eulerianParticleMod>());

        scalar pDiameter = 0.0;
        
        if (mesh_.nSolutionD() == 3)
        {
        	pDiameter = cbrt(6.0*p.V/constant::mathematical::pi);
        }
        else if (mesh_.nSolutionD() == 2)
        {
        	const boundBox& bBox = mesh_.bounds();
        	const point& minBox = bBox.min();
        	const point& maxBox = bBox.max();
        	scalar dimEmpty = 0.0; //Size of the domain in the empty direction
        	
        	forAll(mesh_.solutionD(), dimI)
        	{
        	   if (mesh_.solutionD()[dimI] == -1)
        	   {
        	   	dimEmpty = maxBox[dimI] - minBox[dimI];
        	   }
        	}
        	
        	pDiameter = sqrt(4.0*p.V/(constant::mathematical::pi*dimEmpty));
        }

        if ((pDiameter > minDiameter_) && (pDiameter < maxDiameter_))
        {
            if (Pstream::master())
            {
                const point position = p.VC/(p.V + ROOTVSMALL);
                const vector U = p.VU/(p.V + ROOTVSMALL);
                
                if (writeToFile())
                {
                    extractEulerianFilePtr_() << nCollectedParticles_ << tab << position[0] << tab
                    	<< position[1] << tab << position[2] << tab << time << tab << U.x() << tab
                    	<< U.y() << tab << U.z() << tab << pDiameter << tab << p.V << endl;
                }

                collectedVolume_ += p.V;
            }

            nCollectedParticles_++;
        }
        else
        {
            // Discard particles over/under diameter threshold
            nDiscardedParticles_++;
            discardedVolume_ += p.V;
        }
    }
}

void Foam::functionObjects::extractEulerianParticlesLBAMR::applyFaceZoneFaceMap
(
    const faceZone& fz,
    const boolList& blockedFaces,
    labelList& faceZoneFaceMap
)
{
    labelList faceMap = mesh_.dynamicMultiDimRefineFvMeshFaceMap();
    
    forAll(fz, localFacei)
    {
    	const label globalFacei = fz[localFacei];
    	
    	//Get corresponding face from previous time-step
    	label mappedFaceOld = faceMap[globalFacei];
    	 	
    	//TODO: MAKE THIS LOOKUP PROCEDURE MORE EFFICIENT
    	forAll(faceZone0_, oldFacei)
    	{
    	   if (faceZone0_[oldFacei] == mappedFaceOld)
    	   {
    	   	faceZoneFaceMap[localFacei] = oldFacei;
    	   	break;
    	   }
    	}
    }
}

void Foam::functionObjects::extractEulerianParticlesLBAMR::calculateAddressing
(
    const label nRegionsOld,
    const label nRegionsNew,
    const scalar time,
    const bool meshChanged,
    const labelList& faceZoneFaceMap,
    labelList& regionFaceIDs
)
{
    DebugInFunction << endl;        

    // Determine mapping between old and new regions so that we can
    // accumulate particle info
    boolListList oldToNewRegions(max(nRegionsOld,1));
    forAll(oldToNewRegions, oldRegioni)
    {
        boolList newRegions(max(nRegionsNew,1), false);
        oldToNewRegions[oldRegioni] =  newRegions;
    }
    
    forAll(regionFaceIDs, facei)
    {
        label newRegioni = regionFaceIDs[facei];

        label oldRegioni = -1;
        if (meshChanged)
        {
            label oldFacei = faceZoneFaceMap[facei];
            
            oldRegioni = regions0_[oldFacei];
        }
        else
        {
            oldRegioni = regions0_[facei];
        }
        
        // if face is part of an old and new region, those regions are connected and belong to the same particle
        if (newRegioni != -1 && oldRegioni != -1)
        {
            oldToNewRegions[oldRegioni][newRegioni] = true;
        }
    }
    
    // Ensure all old regions point to the same new regions over all processors
    Pstream::listCombineGather(oldToNewRegions, maxEqOp<boolList>());
    Pstream::listCombineScatter(oldToNewRegions);
    

    // CCL of connected Regions to get Particles
    Map<label> newRegionToParticleMap(nRegionsNew);
    List<eulerianParticleMod> newParticles(nRegionsNew);
    label particlei = 0;
    
    boolList visitedOld(nRegionsOld, false);
    boolList visitedNew(nRegionsNew, false);
    

    if (nRegionsNew) { forAll(oldToNewRegions[0], newRegioni) { if (!visitedNew[newRegioni])
    {
        std::queue<label> Q;
        std::queue<bool> newFlag;
        std::queue<label> oldRegionsToParticle;

        Q.push(newRegioni);
        newFlag.push(true);
        visitedNew[newRegioni] = true;
        newRegionToParticleMap.insert(newRegioni, particlei);

        while(!Q.empty())
        {
            label top = Q.front();
            Q.pop();

            if (newFlag.front())
            {
                // find connections from new to old Region
                newFlag.pop();
                forAll(oldToNewRegions, oldRegioni)
                {
                    if (oldToNewRegions[oldRegioni][top] && !visitedOld[oldRegioni])
                    {
                        Q.push(oldRegioni);
                        oldRegionsToParticle.push(oldRegioni);
                        newFlag.push(false);
                        visitedOld[oldRegioni]=true;
                    }
                }
            }

            else
            {
                // find connections from old to new region
                newFlag.pop();
                forAll(oldToNewRegions[top], newRegionii)
                {
                    if (oldToNewRegions[top][newRegionii] && !visitedNew[newRegionii])
                    {
                        Q.push(newRegionii);
                        newFlag.push(true);
                        visitedNew[newRegionii]=true;
                        newRegionToParticleMap.insert(newRegionii, particlei);
                    }
                }

                // find connections between oldRegions from regionToParticleMap_
                label oldParticlei = regionToParticleMap_[top];
                forAll(oldToNewRegions, oldRegioni)
                {
                    if (regionToParticleMap_[oldRegioni] == oldParticlei && !visitedOld[oldRegioni])
                    {
                        Q.push(oldRegioni);
                        oldRegionsToParticle.push(oldRegioni);
                        newFlag.push(false);
                        visitedOld[oldRegioni]=true;
                    }
                }
            }

        }

        label nConnectedOldRegions = oldRegionsToParticle.size();

        if (nConnectedOldRegions == 1)
        {
            // n:1 matching between new and old regions
            label oldParticlei = regionToParticleMap_[oldRegionsToParticle.front()];
            newParticles[particlei] = particles_[oldParticlei];
        }
        else if (nConnectedOldRegions > 1)
        {
            // n:m matching between new and old regions
            Map<label> oldParticleToRegionMap(nConnectedOldRegions);
            for (label i=0; i<nConnectedOldRegions; i++)
            {
                label oldRegioniConnected = oldRegionsToParticle.front();
                oldRegionsToParticle.pop();
                label oldParticlei = regionToParticleMap_[oldRegioniConnected];
                if (oldParticleToRegionMap.insert(oldParticlei, oldRegioniConnected))
                {
                    newParticles[particlei] =
                    sumParticleOp<eulerianParticleMod>()(newParticles[particlei], particles_[oldParticlei]);
                }
            }
        }


        particlei++;
    }}}
    newParticles.resize(particlei);
    
    //Collect Particles
    if (nRegionsOld)
    {
        boolList collectParticleFlag(particles_.size(), true);
        forAll(oldToNewRegions, oldRegioni)
        {
            label oldParticlei = regionToParticleMap_[oldRegioni];
            if (visitedOld[oldRegioni]) collectParticleFlag[oldParticlei]=false;
        }
        collectParticles(time, collectParticleFlag);
    }

    // Reset the particles list and addressing for latest available info
    particles_.transfer(newParticles);
    regionToParticleMap_ = newRegionToParticleMap;
}


void Foam::functionObjects::extractEulerianParticlesLBAMR::accumulateParticleInfo
(
    const surfaceScalarField& alphaf,
    const surfaceScalarField& phi,
    const labelList& regionFaceIDs,
    const faceZone& fz
)
{
    DebugInFunction << endl;

    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);
    const surfaceVectorField Uf(fvc::interpolate(U));

    const scalar deltaT = mesh_.time().deltaTValue();
    const pointField& faceCentres = mesh_.faceCentres();

    scalar dVtotal = 0;
    scalar dVavailable = 0;
    
    forAll(regionFaceIDs, localFacei)
    {
        const label newRegioni = regionFaceIDs[localFacei];
        const label meshFacei = fz[localFacei];
        
        scalar magPhii = mag(faceValue(phi, localFacei, meshFacei));

        scalar alphai = faceValue(alphaf, localFacei, meshFacei);
        
        scalar dV = magPhii*deltaT*alphai;

        dVtotal += dV;
        if (newRegioni != -1)
        {
            dVavailable += dV;
            const label particlei = regionToParticleMap_[newRegioni];
            eulerianParticleMod& p = particles_[particlei];

            if (p.faceIHit < 0)
            {
                // New particle - does not exist in particles_ list
                p.faceIHit = localFacei;
                p.V = 0;
                p.VC = vector::zero;
                p.VU = vector::zero;
            }

            // Accumulate particle properties
            vector Ufi = faceValue(Uf, localFacei, meshFacei);
            p.V += dV;
            p.VC += dV*faceCentres[meshFacei];
            p.VU += dV*Ufi;
        }
    }
    
    reduce(dVtotal, sumOp<scalar>());
    reduce(dVavailable, sumOp<scalar>());
    totalVolume_ += dVtotal;
    availableVolume_ += dVavailable;

}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::functionObjects::extractEulerianParticlesLBAMR::extractEulerianParticlesLBAMR
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    dynamicMultiDimRefineFvMeshFunctionObject(name, runTime, dict),
    writeFile(runTime, name),
    faceZoneName_(word::null),
    zoneID_(-1),
    patchIDs_(),
    patchFaceIDs_(),
    faceZone0_(),
    alphaName_("alpha"),
    alphaThreshold_(1e-8),
    UName_("U"),
    rhoName_("rho"),
    phiName_("phi"),
    nInjectorLocations_(0),
    fineToCoarseAddr_(),
    globalCoarseFaces_(),
    regions0_(),
    nRegions0_(0),
    particles_(),
    regionToParticleMap_(),
    minDiameter_(ROOTVSMALL),
    maxDiameter_(GREAT),
    nCollectedParticles_(getProperty<label>("nCollectedParticles", 0)),
    collectedVolume_(getProperty<scalar>("collectedVolume", 0)),
    nDiscardedParticles_(getProperty<label>("nDiscardedParticles", 0)),
    discardedVolume_(getProperty<scalar>("discardedVolume", 0)),
    availableVolume_(getProperty<scalar>("availableVolume", 0)),
    totalVolume_(getProperty<scalar>("totalVolume", 0))
{
    read(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::functionObjects::extractEulerianParticlesLBAMR::read
(
    const dictionary& dict
)
{
    DebugInFunction << endl;

    if (dynamicMultiDimRefineFvMeshFunctionObject::read(dict) && writeFile::read(dict))
    {
        dict.readEntry("faceZone", faceZoneName_);
        dict.readEntry("alpha", alphaName_);

        dict.readIfPresent("alphaThreshold", alphaThreshold_);
        dict.readIfPresent("U", UName_);
        dict.readIfPresent("rho", rhoName_);
        dict.readIfPresent("phi", phiName_);
        dict.readIfPresent("nLocations", nInjectorLocations_);
        dict.readIfPresent("minDiameter", minDiameter_);
        dict.readIfPresent("maxDiameter", maxDiameter_);

        checkFaceZone();

        if (nInjectorLocations_)
        {
            initialiseBins();
        }

        return true;
    }

    return false;
}


bool Foam::functionObjects::extractEulerianParticlesLBAMR::execute()
{
    DebugInFunction << endl;

    Log << type() << " " << name() << " output:" << nl;
    
    createFiles();

    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>(alphaName_);

    const surfaceScalarField alphaf
    (
        typeName + ":alphaf",
        fvc::interpolate(alpha)
    );

    const faceZone& fz = mesh_.faceZones()[zoneID_];
    const indirectPrimitivePatch patch
    (
        IndirectList<face>(mesh_.faces(), fz),
        mesh_.points()
    );
    
    // Set the blocked faces, i.e. where alpha > alpha threshold value
    boolList blockedFaces(fz.size(), false);
    setBlockedFaces(alphaf, fz, blockedFaces);
    
    // Split the  faceZone according to the blockedFaces
    // - Returns a list of (disconnected) region index per face zone face
    regionSplit2Dv1906 regionFaceIDs(mesh_, patch, blockedFaces);
    
    // Global number of regions
    const label nRegionsNew = regionFaceIDs.nRegions();
    
    // Map of the face zone local faces
    labelList faceZoneFaceMap(fz.size(), -1);
    
    // If the faceZone has changed, map the old regions
    bool meshChanged = false;
    if (mesh_.changing())
    {
    	applyFaceZoneFaceMap(fz, blockedFaces, faceZoneFaceMap);
    	meshChanged = true;
    } 
    
    // Calculate the addressing between the old and new region information
    // Also collects particles that have traversed the faceZone
    // - Note: may also update regionFaceIDs
    calculateAddressing
    (
        nRegions0_,
        nRegionsNew,
        mesh_.time().value(),
        meshChanged,
        faceZoneFaceMap,
        regionFaceIDs
    );
    
    // Process latest region information
    tmp<surfaceScalarField> tphi = phiU();
    accumulateParticleInfo(alphaf, tphi(), regionFaceIDs, fz);

    // Reset the blocked faces for the next integration step
    nRegions0_ = nRegionsNew;
    regions0_ = regionFaceIDs;
    // Save the faceZone faces for the next step
    faceZone0_.setSize(fz.size(), -1);
    
    forAll(fz, localFacei)
    {
    	faceZone0_[localFacei] = fz[localFacei];
    }


    Log << "    Collected particles   : " << nCollectedParticles_ << nl
        << "    Discarded particles   : " << nDiscardedParticles_ << nl
        << "    Particles in progress : " << particles_.size() << nl
        << "    Collected volume      : " << collectedVolume_ << nl
        << "    Available volume      : " << availableVolume_ << nl
        << "    Total volume          : " << totalVolume_ << nl
        << "    Discarded volume      : " << discardedVolume_ << nl
        << endl;

    
    return true;
}


bool Foam::functionObjects::extractEulerianParticlesLBAMR::write()
{
    DebugInFunction << endl;

    setProperty("nCollectedParticles", nCollectedParticles_);
    setProperty("nDiscardedParticles", nDiscardedParticles_);
    setProperty("collectedVolume", collectedVolume_);
    setProperty("availableVolume", availableVolume_);
    setProperty("totalVolume", totalVolume_);
    setProperty("discardedVolume", discardedVolume_);

    return true;
}


// ************************************************************************* //
