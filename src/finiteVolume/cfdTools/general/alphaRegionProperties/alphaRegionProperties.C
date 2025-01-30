/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
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

#include "alphaRegionProperties.H"
#include "boundBox.H"
#include "OFstream.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
void Foam::alphaRegionProperties::setBlockedFaces
(
    boolList& liquidBlockedFaces
)
{
    forAll(mesh_.cells(), meshCelli)
    {
    
        const labelList& cellFaces = mesh_.cells()[meshCelli];
        forAll(cellFaces, meshFacei)
        {
            label facei = cellFaces[meshFacei];
            if (mesh_.isInternalFace(facei))
            {
                const label own = mesh_.faceOwner()[facei];
                const label nb = mesh_.faceNeighbour()[facei];
                if ((alpha_[own] > alphaThreshold_) != (alpha_[nb] > alphaThreshold_))
                {
                    liquidBlockedFaces[facei] = true;
                }
            }
            else
            {
                label patchi = mesh_.boundaryMesh().whichPatch(facei);
                const polyPatch& pp = mesh_.boundaryMesh()[patchi];
                if (pp.coupled())
                {
                    label patchFacei = facei - pp.start();
                    label own = mesh_.faceOwner()[facei];
                    const scalarField alphaPNb = alpha_.boundaryField()[patchi].patchNeighbourField();

                    if ((alpha_[own] > alphaThreshold_) != (alphaPNb[patchFacei] > alphaThreshold_))
                    {
                        liquidBlockedFaces[facei] = true;
                    }
                }
                else
                {
                    liquidBlockedFaces[facei] = true;
                }
            }
        }
    }
}

void Foam::alphaRegionProperties::removeAlphaRegions
(
    DynamicList<DynamicList<label>>& alphaRegionLocalLabels
)
{
    scalar volumeRemoved = 0.0;
    
    forAll(alphaRegionLocalLabels, regioni)
    {
        label localRegionSize = alphaRegionLocalLabels[regioni].size();
        label globalRegionSize = localRegionSize;
        reduce(globalRegionSize, sumOp<label>());

        if (globalRegionSize < NCellsMin_ && localRegionSize > 0)
        {           
            forAll (alphaRegionLocalLabels[regioni], regionCelli)
            {
                label localCelli = alphaRegionLocalLabels[regioni][regionCelli];

                scalar dV = mesh_.V()[localCelli] * alpha_[localCelli];
                volumeRemoved += dV;
                alpha_[localCelli] = 0.0;
            }

            alphaRegionLocalLabels[regioni].remove();
        }
    }
    reduce(volumeRemoved, sumOp<scalar>());
    removedVolumeTotal_ += volumeRemoved;
    Info << "Total volume removed: " << removedVolumeTotal_ << endl;
}

void Foam::alphaRegionProperties::getInterfaceLabelsFromRegions
(
    const DynamicList<DynamicList<label>>& regionLabels,
    DynamicList<DynamicList<label>>& interfaceLabels
)
{
    interfaceLabels.clear();
    interfaceLabels.setSize(regionLabels.size());

    labelList interfaceSize(regionLabels.size(), 0);

    const boolList& interfaceCell = surf_.interfaceCell();

    forAll(regionLabels, regioni)
    {
        forAll(regionLabels[regioni], cellidx)
        {
            label localCelli = regionLabels[regioni][cellidx];

            if (interfaceCell[localCelli])
            {
                interfaceSize[regioni]++;
                interfaceLabels[regioni].append(localCelli);
            }
        }
    }
}

void Foam::alphaRegionProperties::calculateAlphaRegionProperties
(
    const DynamicList<DynamicList<label>>& alphaRegionLabels
)
{
    List<Foam::eulerianParticleMod> alphaRegions(alphaRegionLabels.size());
    
    accumulateRegionProperties(alphaRegionLabels, alphaRegions);

    scalarList regionD(alphaRegions.size());
    volumeEquivalentDiameter(alphaRegions, regionD);

    scalarList regionRMin(alphaRegions.size());
    scalarList regionRMax(alphaRegions.size());

    minMaxRadii(alphaRegionLabels, alphaRegions, regionRMin, regionRMax);

    List<vector> contPhaseU(alphaRegions.size());
    continuousPhaseVelocity(alphaRegions, regionD, contPhaseU);

    writeRegionProperties(alphaRegions, regionD, regionRMin, regionRMax, contPhaseU);
}

void Foam::alphaRegionProperties::accumulateRegionProperties
(
    const DynamicList<DynamicList<label>>& alphaRegionLabels,
    List<Foam::eulerianParticleMod>& alphaRegions
)
{    
    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);
    const volScalarField& T = mesh_.lookupObject<volScalarField>(TName_);

    forAll(alphaRegionLabels, regioni)
    {
        DynamicList<label> regionILabels = alphaRegionLabels[regioni];
        Foam::eulerianParticleMod& p = alphaRegions[regioni];

        if (regionILabels.size())
        {
            const pointField& cellCentres = mesh_.cellCentres();

            forAll(regionILabels, regionCelli)
            {
                label localCelli = alphaRegionLabels[regioni][regionCelli];
                        
                scalar dV = mesh_.V()[localCelli] * alpha_[localCelli];
                p.V += dV;
                p.VC += dV * cellCentres[localCelli];
                p.VU += dV * U[localCelli];
                p.VT += dV * T[localCelli];
                
            }
        }

        reduce(p, Foam::sumParticleOp<Foam::eulerianParticleMod>());
        p.VC = p.VC /(p.V + ROOTVSMALL);
        p.VU = p.VU /(p.V + ROOTVSMALL);
        p.VT = p.VT /(p.V + ROOTVSMALL);
    }
}

void Foam::alphaRegionProperties::volumeEquivalentDiameter
(
    const List<Foam::eulerianParticleMod>& alphaRegions,
    scalarList& regionD
)
{
    forAll(alphaRegions, regioni)
    {
        const Foam::eulerianParticleMod& p = alphaRegions[regioni];
        
        if (mesh_.nSolutionD() == 3)
        {
        	regionD[regioni] = cbrt(6.0*p.V/constant::mathematical::pi);
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
        	regionD[regioni] = sqrt(4.0*p.V/(constant::mathematical::pi*dimEmpty));
        }
    }            
}

void Foam::alphaRegionProperties::minMaxRadii
(
    const DynamicList<DynamicList<label>>& alphaRegionLabels,
    const List<Foam::eulerianParticleMod>& alphaRegions,
    scalarList& regionRMin,
    scalarList& regionRMax
)
{
    const boolList& interfaceCell = surf_.interfaceCell();

    forAll(alphaRegions,regioni)
    {
        DynamicList<label> regionILabels = alphaRegionLabels[regioni];

        scalar rMin2 = GREAT;
        scalar rMax2 = 0;

        if (regionILabels.size())
        {
            const Foam::eulerianParticleMod& p = alphaRegions[regioni];
            point COMi(p.VC);

            forAll(regionILabels, regionCelli)
            {
                label localCelli = alphaRegionLabels[regioni][regionCelli];

                if (interfaceCell[localCelli])
                {
                    point surfCenteri(surf_.centre()[localCelli]);

                    vector ctrToCOM = surfCenteri - COMi;
                    scalar r2 = ctrToCOM.x()*ctrToCOM.x() + ctrToCOM.y()*ctrToCOM.y() + ctrToCOM.z()*ctrToCOM.z();
                    if (r2 < rMin2)
                    {
                        rMin2 = r2;
                    }
                    if (r2 > rMax2)
                    {
                        rMax2 = r2;
                    }
                }
            }
        }
        
        regionRMin[regioni] = sqrt(rMin2);
        regionRMax[regioni] = sqrt(rMax2);
    }

    reduce(regionRMin, minOp<scalarList>());
    reduce(regionRMax, maxOp<scalarList>());
}

void Foam::alphaRegionProperties::continuousPhaseVelocity
(
    const List<Foam::eulerianParticleMod>& alphaRegions,
    const scalarList& regionD,
    List<vector>& contPhaseU
)
{    
    const pointField& ctrs = mesh_.cellCentres();

    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);

    vector zeroVec(0,0,0);
    List<vector> alphaUdV(alphaRegions.size(), zeroVec);
    scalarList alphadV(alphaRegions.size(),SMALL);

    forAll(ctrs, elemi)
    {
        forAll(alphaRegions, regioni)
        {
            const Foam::eulerianParticleMod& p = alphaRegions[regioni];
            scalar dist2 = magSqr(ctrs[elemi] - p.VC);
            scalar maxDist = 1.5 * regionD[regioni];
            scalar maxDist2 = maxDist * maxDist;

            if ((maxDist > 0.0) && (dist2 < maxDist2))
            {
                scalar dV = mesh_.V()[elemi];
                alphadV[regioni] += (1 - alpha_[elemi]) * dV;
                alphaUdV[regioni] += (1 - alpha_[elemi]) * dV * U[elemi];
            }
        }

    }
    reduce(alphadV, sumOp<scalarList>());
    reduce(alphaUdV, sumOp<List<vector>>());

    contPhaseU = alphaUdV / alphadV;
}


void Foam::alphaRegionProperties::writeRegionProperties
(
    const List<Foam::eulerianParticleMod>& alphaRegions,
    const scalarList& regionD,
    const scalarList& regionRMin,
    const scalarList& regionRMax,
    const List<vector>& contPhaseU
)
{    

    if (Pstream::master())
    {
        fileName VOFCloudFileName = "VOFCloudData/" + mesh_.time().timeName()  ;

        OFstream OS(VOFCloudFileName);
        if (OS.opened())
        {
            OS << "Volume" << tab << "COM" << tab << "D" << tab << "RMin" << tab << "RMax" << tab
                << "UDrop" << tab << "UCont" << tab << "T" << endl;
        }

        forAll(alphaRegions, regioni)
        {
            const Foam::eulerianParticleMod& p = alphaRegions[regioni];

            if (OS.opened() && p.V > 0.0)
            {
                OS << p.V << tab << p.VC << tab << regionD[regioni]
                    << tab << regionRMin[regioni] << tab << regionRMax[regioni] << tab << p.VU
                    << tab << contPhaseU[regioni] << tab << p.VT << endl;
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::alphaRegionProperties::alphaRegionProperties
(
    const dictionary& dict,
    isoAdvection& advector,
    volScalarField& alpha
)
:
    NCellsMin_(0),
    removeRegions_(false),
    calculateRegionProperties_(false),
    alphaThreshold_(1e-8),
    removedVolumeTotal_(0.0),
    UName_("U"),
    rhoName_("rho"),
    TName_("T"),
    mesh_(alpha.mesh()),
    cloudWriteInterval_(GREAT),
    cloudWriteTimeIndex_(-1),
    alpha_(alpha),
    surf_(advector.surf())
{
    
    if (dict.found("NCellsMin"))
    {
        removeRegions_ = true;
        dict.readEntry("NCellsMin", NCellsMin_);
    }

    if (dict.found("alphaThreshold"))
    {
        dict.readEntry("alphaThreshold", alphaThreshold_);
    }

    if (dict.found("VOFCloudWriteInterval"))
    {
        calculateRegionProperties_ = true;
        dict.readEntry("VOFCloudWriteInterval",cloudWriteInterval_);
        mkDir("VOFCloudData");
    }
    
    /* Print info */
    Info<< "alphaRegionProperties" << nl;
    if (removeRegions_)
    {
        Info<< "Removing all alpha regions smaller than " << NCellsMin_ << " cells" << endl;
    }
    if (calculateRegionProperties_)
    {
        Info << "Calculating properties of the alpha regions " << endl;
    }
    Info<< "Threshold for alpha " << alphaThreshold_ << endl;

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::alphaRegionProperties::execute
(
    DynamicList<DynamicList<label>>& interfaceLabels,
    bool firstIter
)
{     
    interfaceLabels.clearStorage();

    // Block the faces of the regions
    boolList liquidBlockedFaces(mesh_.faces().size(), false);
    setBlockedFaces(liquidBlockedFaces);

    // INTERFACE FROM THE LIQUID REGIONS
    modRegionSplit liquidRegions(mesh_, liquidBlockedFaces);
    DynamicList<DynamicList<label>>& liquidRegionLocalLabels = liquidRegions.getLocalRegionCellLabels();

    //Remove liquid regions if they are small
    if (removeRegions_)
    {
        removeAlphaRegions(liquidRegionLocalLabels);
    }
    
    getInterfaceLabelsFromRegions(liquidRegionLocalLabels, interfaceLabels);
    
    /*********************************************************/
    
    if (calculateRegionProperties_ && firstIter)
    {
        scalar timeValue = mesh_.time().value();
        scalar startTime = mesh_.time().startTime().value();
        scalar dt = mesh_.time().deltaT().value();

        if (timeValue == dt)
        {
            cloudWriteTimeIndex_ = 0;
        }
        label writeIndex = ((timeValue - startTime) + 0.5*dt)/cloudWriteInterval_;
        if (writeIndex > cloudWriteTimeIndex_)
        {
            calculateAlphaRegionProperties(liquidRegionLocalLabels);
            cloudWriteTimeIndex_ = writeIndex;
        }
        
    }

}

// ************************************************************************* //
