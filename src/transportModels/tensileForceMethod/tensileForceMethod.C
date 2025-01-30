/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "tensileForceMethod.H"
#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "syncTools.H"
#include "slicedSurfaceFields.H"
#include "unitConversion.H"
#include "symmetryPolyPatch.H"
#include "reconstructedDistanceFunction.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::tensileForceMethod::correctContactAngle() const
{
    // Cell gradient of alpha
    const volVectorField gradAlpha(fvc::grad(alpha1_, "nHat"));

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));

    surfaceVectorField::Boundary& nHatb = nHatfv.boundaryFieldRef();
    const surfaceVectorField::Boundary& gradAlphabf = gradAlphaf.boundaryField();

    const volScalarField::Boundary& abf = alpha1_.boundaryField();

    const fvBoundaryMesh& boundary = mesh_.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleTwoPhaseFvPatchScalarField& acap =
                const_cast<alphaContactAngleTwoPhaseFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleTwoPhaseFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                degToRad() * acap.theta(U_.boundaryField()[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            acap.gradient() = (nf & nHatp)*mag(gradAlphabf[patchi]);
            acap.evaluate();
        }
    }
}

void Foam::tensileForceMethod::calculateTensileForce
(
    const DynamicList<DynamicList<label>>& connectedInterfaceLabels,
    List<scalar>& interfacePressureJump
)
{    
    Fsigma_ = dimensionedVector(dimMass/(dimArea*sqr(dimTime)), Zero);
    sigma_ = sigmaPtr_->sigma().ref();

    // Initialize CPC stencil for nb normal search. Has to be initialized here because of parallel communication (all procs must reach this point)
    zoneDistribute& CPCDistribute = zoneDistribute::New(mesh_);
    CPCDistribute.setUpCommforZone(nextToInterface_, false);
    Map<vector> mapSurfNormal(CPCDistribute.getDatafromOtherProc(nextToInterface_, surfNormal_));
    Map<vector> mapCellCentre(CPCDistribute.getDatafromOtherProc(nextToInterface_, mesh_.C()));

    // Interfacial summations required for pressure jump calculation
    List<scalar> interfaceFsigmaSum(interfacePressureJump.size(), Zero);
    List<scalar> interfaceAreaSum(interfacePressureJump.size(), Zero);

    // Loop over all cells that contain an interface
    // Calculates first Fsigma, then calculate required components to distribute Fpj
    forAll(connectedInterfaceLabels, interfacei)
    {

        List<label> interfaceILabels = connectedInterfaceLabels[interfacei];

        if (interfaceILabels.size())
        {
            forAll(interfaceILabels, interfaceCelli)
            {
                label localCelli = interfaceILabels[interfaceCelli];

                if (mag(surfNormal_[localCelli]) != 0)
                {
                    vector FsigmaCelli = plicTensileForce(localCelli, mapSurfNormal, mapCellCentre, CPCDistribute);
                    Fsigma_[localCelli] = FsigmaCelli / mesh_.V()[localCelli];

                    vector n = surfNormal_[localCelli]/mag(surfNormal_[localCelli]);
                    interfaceFsigmaSum[interfacei] += FsigmaCelli & n;
                    interfaceAreaSum[interfacei] += mag(surfNormal_[localCelli]);
                }
            }
        }
    }

    reduce(interfaceFsigmaSum, sumOp<List<scalar>>());
    reduce(interfaceAreaSum, sumOp<List<scalar>>());

    forAll(connectedInterfaceLabels, interfacei)
    {
        if (interfaceAreaSum[interfacei] != 0)
        {
            interfacePressureJump[interfacei] = interfaceFsigmaSum[interfacei] / interfaceAreaSum[interfacei];

            /*Info << "For interface " << interfacei << " with surface area " << interfaceAreaSum[interfacei]
                << ", the pressure jump equals " << interfacePressureJump[interfacei] << endl;*/
        }
    }
}

Foam::vector Foam::tensileForceMethod::plicTensileForce
(
    const label celli,
    const Map<vector>& mapSurfNormal,
    const Map<vector>& mapCellCentre,
    zoneDistribute& CPCDistribute
)
{
    const labelList& cellFaces = mesh_.cells()[celli];
    cutFacePLIC faceCut(mesh_);
    DynamicList<point> facePlicPts;

    label nCutFaces = 0;
    label nTensileForces = 0;
    List<vector> tensileForce(cellFaces.size(), Zero);
    List<scalar> edgeLength(cellFaces.size(), Zero);

    vector n = surfNormal_[celli]/mag(surfNormal_[celli]);
    scalar surfTensionCoeff = sigma_[celli];

    forAll(cellFaces, cellFacei)
    {
        label meshFacei = cellFaces[cellFacei];

        // Calculate the tensile force for faces that are cut by interface element edge
        if (faceCut.calcSubFace(meshFacei, n, surfCentre_[celli]) == 0)
        {
            // Get the cut points and calculate the tangent. Check for tangent direction
            facePlicPts.clear();
            facePlicPts.append(faceCut.surfacePoints());
            
            if (facePlicPts.size() == 2) //Always the case, but check just in case.
            {
                vector nbNormal(Zero);   
                bool emptyFace = neighbourNormal(celli, meshFacei, nbNormal);

                if(!emptyFace)
                {   
                    if(mag(nbNormal) == 0.0)
                    {
                        nbNormal = neighbourNormalFromCPCStencil(celli, meshFacei, facePlicPts, mapSurfNormal, mapCellCentre, CPCDistribute);
                    }
                    nCutFaces++;

                    edgeLength[cellFacei] = mag(facePlicPts[0]-facePlicPts[1]);
                }

                if (mag(nbNormal) != 0)
                {
                    //Tangent vector to the interface element needs to be oriented anticlockwise
                    vector plicTangent = antiClockWisePlicTangent(facePlicPts, celli);

                    tensileForce[cellFacei] = plicTangent ^ nbNormal;
                    nTensileForces++;
                }
            }
        }
    }

    if ((nCutFaces != nTensileForces) && nTensileForces > 0)
    {
        return balanceTensileForce(celli, tensileForce, edgeLength, surfTensionCoeff);
    }
    return 0.5 * surfTensionCoeff * sum(tensileForce);

}

Foam::vector Foam::tensileForceMethod::antiClockWisePlicTangent
(
    const DynamicList<point>& facePlicPts,
    const label celli
)
{
    vector p0 = facePlicPts[0];
    vector p1 = facePlicPts[1];

    vector v0 = surfCentre_[celli] - p0;
    vector v1 = surfCentre_[celli] - p1;

    vector crossProduct = v0 ^ v1;

    scalar orientation = surfNormal_[celli] & crossProduct;

    vector plicTangent(Zero);
    if (orientation > 0)
    {
        plicTangent = p1 - p0;
    }
    else
    {
        plicTangent = p0 - p1;
    }

    return plicTangent;
}

bool Foam::tensileForceMethod::neighbourNormal
(
    const label celli,
    const label facei,
    vector& nbNormal
)
{

    if (mesh_.isInternalFace(facei))
    {
        const label own = mesh_.faceOwner()[facei];
        const label nb = mesh_.faceNeighbour()[facei];

        if (celli == own)
        {
            nbNormal = surfNormal_[nb];
        }
        else if (celli == nb)
        {
            nbNormal = surfNormal_[own];
        }
    }
    else    // Face on patch, different treatment based on patch type
    {
        label patchi = mesh_.boundaryMesh().whichPatch(facei);
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        label patchFacei = facei - pp.start();

        if (isA<emptyPolyPatch>(pp))
        {// No force calculation on empty patch
            return true;
        }
        else if (isA<symmetryPolyPatch>(pp))
        { // Mirror the normal vector in the face
            const surfaceVectorField::Boundary& SfHatBf = SfHat_.boundaryField();
            vector SfHati = SfHatBf[patchi][patchFacei];

            vector n = surfNormal_[celli];

            nbNormal = n - (2 * n & SfHati) * SfHati / (mag(SfHati)*mag(SfHati));
        }
        else if (pp.coupled())
        {// Get the neigbour from the patch.
            const vectorField normalNbField = surfNormal_.boundaryField()[patchi].patchNeighbourField();

            nbNormal = normalNbField[patchFacei];
        }
        else
        {
            nbNormal = surfNormal_[celli];
        }
    }

    if (mag(nbNormal) !=0)
    {
        nbNormal = nbNormal/mag(nbNormal);
    }

    return false;
}

Foam::vector Foam::tensileForceMethod::neighbourNormalFromCPCStencil
(
    const label celli,
    const label facei,
    const DynamicList<point>& facePlicPts,
    const Map<vector>& mapSurfNormal,
    const Map<vector>& mapCellCentre,
    zoneDistribute& CPCDistribute
)
{
    vector nbNormal(Zero);
    vector SfHati(Zero);
    const surfaceVectorField::Boundary& SfHatBf = SfHat_.boundaryField();

    if(mesh_.isInternalFace(facei))
    {
        SfHati = SfHat_[facei];
    }
    else
    {
        label patchi = mesh_.boundaryMesh().whichPatch(facei);
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];

        label patchFacei = facei - pp.start();
        SfHati = SfHatBf[patchi][patchFacei];
    }

    const labelListList& CPCStencil = CPCDistribute.getStencil();
    
    point facePlicEdgeCentre = 0.5*sum(facePlicPts);
    vector cellCenterToPlicEdge = facePlicEdgeCentre - mesh_.C()[celli];
    scalar distCellCentreToPlicEdge = SfHati & cellCenterToPlicEdge;

    scalar maxCosTheta = -1.0;

    for (const label gblIdx : CPCStencil[celli])
    {
        vector stencilCellCentrei = CPCDistribute.getValue(mesh_.C(), mapCellCentre, gblIdx);
        vector stencilCentreToPlicEdge = facePlicEdgeCentre - stencilCellCentrei;
        scalar distStencilCentreToPlicEdge = SfHati & stencilCentreToPlicEdge;

        if (neg(distCellCentreToPlicEdge*distStencilCentreToPlicEdge)) // Cell centers must be on opposite sides of the face
        {
            vector stencilNormali = CPCDistribute.getValue(surfNormal_, mapSurfNormal, gblIdx);

            if (mag(stencilNormali) != 0)
            {
                scalar cosNormalAngle = (surfNormal_[celli] & stencilNormali) / (mag(surfNormal_[celli])*mag(stencilNormali)); // Angle between normal of celli and normal from stencil

                if (cosNormalAngle > 0 && cosNormalAngle < 1) // Angle between normals is between 0 and 90 degrees
                {
                    // Cell with normal closest to celli gets assigned
                    if(cosNormalAngle > maxCosTheta)
                    {
                        nbNormal = stencilNormali;
                        maxCosTheta = cosNormalAngle;
                    }
                }
            }
        }
    }

    if (mag(nbNormal) !=0)
    {
        nbNormal = nbNormal/mag(nbNormal);
    }

    return nbNormal;
}

Foam::vector Foam::tensileForceMethod::balanceTensileForce
(
    const label celli,
    const List<vector>& tensileForce,
    const List<scalar>& edgeLength,
    const scalar surfTensionCoeff
) const
{
    vector forceSum = sum(tensileForce);
    vector n = surfNormal_[celli]/mag(surfNormal_[celli]);

    if (mag(forceSum) > 0.0)
    {
        if(mesh_.nSolutionD() == 3)
        {
            scalar totalEdgeLength = sum(edgeLength);
            scalar forceEdgeLength = 0.0;

            forAll(tensileForce, cellFacei)
            {
                if (mag(tensileForce[cellFacei]) != 0)
                {
                    forceEdgeLength += edgeLength[cellFacei];
                }
            }

            vector forceProjection = (forceSum & n) * n;

            if (mag(forceProjection) != 0.0)
            {
                scalar sign = (n & forceProjection) / (mag(forceProjection)*mag(n));
                scalar magnitude = mag(forceProjection) * totalEdgeLength / forceEdgeLength;

                return 0.5 * surfTensionCoeff * sign * n * magnitude;
            }
        }
        else
        {
            return (forceSum & n) * n;
        }
    }
    
    return forceSum;

}

void Foam::tensileForceMethod::pressureJumpCorrection
(
    const DynamicList<DynamicList<label>>& connectedInterfaceLabels,
    const List<scalar>& interfacePressureJump
)
{
    Fpj_ = dimensionedVector(dimMass/(dimArea*sqr(dimTime)), Zero);

    forAll(connectedInterfaceLabels, interfacei)
    {
        if (interfacePressureJump[interfacei] !=0 )
        {
            List<label> interfaceILabels = connectedInterfaceLabels[interfacei];

            forAll(interfaceILabels, interfaceCelli)
            {
                label localCelli = interfaceILabels[interfaceCelli];

                if (mag(surfNormal_[localCelli]) != 0)
                {
                    Fpj_[localCelli] = - interfacePressureJump[interfacei] * surfNormal_[localCelli] / mesh_.V()[localCelli];
                }
            }
        }
    }
}

void Foam::tensileForceMethod::mapInterfaceForceDensity()
{
    zoneDistribute& CPCDistribute = zoneDistribute::New(mesh_);
    CPCDistribute.setUpCommforZone(nextToInterface_, false);

    Map<vector> mapSurfCentre(CPCDistribute.getDatafromOtherProc(nextToInterface_, surfCentre_));
    Map<vector> mapCellCentre(CPCDistribute.getDatafromOtherProc(nextToInterface_, mesh_.C()));
    Map<scalar> mapRho(CPCDistribute.getDatafromOtherProc(nextToInterface_, rho_));

    const labelListList& CPCStencil = CPCDistribute.getStencil();
    
    volVectorField weightField
    (
        IOobject
        (
            "weightField",
            alpha1_.time().timeName(),
            mesh_
        ),
        mesh_,
        vector(Zero)
    );
    interfaceForce_ = dimensionedVector(dimMass/(sqr(dimTime)*dimArea), Zero);
    const boolList& interfaceCell = surf_.interfaceCell();

    // Set weight totals for interface cells
    forAll(mesh_.cells(), celli)
    {
        scalar sumWeight = 0.0;

        if(interfaceCell[celli])
        {
            // Loop over all cells in the stencil of celli
            // CAREFUL: gblIdx is not the global cell index, only position in the Map!
            for (const label gblIdx : CPCStencil[celli])
            {
                scalar rho = CPCDistribute.getValue(rho_, mapRho, gblIdx);
                vector cellCentre = CPCDistribute.getValue(mesh_.C(), mapCellCentre, gblIdx);
                
                scalar d = mag(cellCentre - surfCentre_[celli]);
                if (d == 0.0)
                {
                    d = SMALL;
                }

                sumWeight += rho/d;
                               
            }

            vector forceTotal = Fsigma_[celli] + Fpj_[celli];
            weightField[celli] = forceTotal / sumWeight;
        }
    }

    Map<vector> mapWeightField(CPCDistribute.getDatafromOtherProc(nextToInterface_, weightField));

    forAll(mesh_.cells(), celli)
    {
        if(nextToInterface_[celli] || interfaceCell[celli])
        {
            // Loop over all cells in the stencil of celli
            // CAREFUL: gblIdx is not the global cell index, only position in the Map!
            for (const label gblIdx : CPCStencil[celli])
            {
                vector weightForce = CPCDistribute.getValue(weightField, mapWeightField, gblIdx);

                if (mag(weightForce) != 0)
                {
                    scalar rho = rho_[celli];
                    scalar d = mag(mesh_.C()[celli] - CPCDistribute.getValue(surfCentre_, mapSurfCentre, gblIdx));

                    if (d == 0.0)
                    {
                        d = SMALL;
                    }
                    interfaceForce_[celli] += (rho/d) * weightForce;
                }
            }
        }
    }
}

void Foam::tensileForceMethod::decomposeForceToCellFaces()
{
    interfaceForcef_ = dimensionedVector(dimMass/(dimArea*sqr(dimTime)), Zero);

    // Create the objects required for syncronizing the surface field over the different processors
    vectorField allSurfForcef(mesh_.nFaces(), Zero);
    const surfaceVectorField::Boundary& SfBf = mesh_.Sf().boundaryField();
    const surfaceVectorField::Boundary& SfHatBf = SfHat_.boundaryField();

    surfaceScalarField rhof = fvc::interpolate(rho_);
    const surfaceScalarField::Boundary& rhofBf = rhof.boundaryField();

    slicedSurfaceVectorField intForce
    (
        IOobject
        (
            "intForcef",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimMass/(sqr(dimTime)*dimArea),
        allSurfForcef,
        false
    );

    // Split in internal and boundary field
    vectorField& intForceIf = intForce;
    surfaceVectorField::Boundary& intForceBf = intForce.boundaryFieldRef();
    
    forAll(mesh_.cells(), celli)
    {
        if (mag(interfaceForce_[celli]) != 0)
        {

            const labelList& cellFaces = mesh_.cells()[celli];

            scalarList forceDecompositionX(cellFaces.size(), Zero);
            scalarList forceDecompositionY(cellFaces.size(), Zero);
            scalarList forceDecompositionZ(cellFaces.size(), Zero);

            forAll(cellFaces, cellFacei)
            {
                label meshFacei = cellFaces[cellFacei];
                scalar rhofi = 0.0;
                vector Sfi(Zero);
                vector SfHati(Zero);

                if(mesh_.isInternalFace(meshFacei))
                {
                    rhofi = rhof[meshFacei];
                    Sfi = mesh_.Sf()[meshFacei];
                    SfHati = SfHat_[meshFacei];

                }
                else
                {
                    label patchi = mesh_.boundaryMesh().whichPatch(meshFacei);
                    const polyPatch& pp = mesh_.boundaryMesh()[patchi];

                    if (!isA<emptyPolyPatch>(pp))
                    {
                        label patchFacei = meshFacei - pp.start();

                        rhofi = rhofBf[patchi][patchFacei];
                        Sfi = SfBf[patchi][patchFacei];
                        SfHati = SfHatBf[patchi][patchFacei];
                    }
                }

                scalar faceNormalComponent = SfHati & interfaceForce_[celli];
                forceDecompositionX[cellFacei] = faceNormalComponent * rhofi * SfHati.x() * mag(Sfi);
                forceDecompositionY[cellFacei] = faceNormalComponent * rhofi * SfHati.y() * mag(Sfi);
                forceDecompositionZ[cellFacei] = faceNormalComponent * rhofi * SfHati.z() * mag(Sfi);
            }

            scalar forceDecompositionSumX = sum(forceDecompositionX);
            scalar forceDecompositionSumY = sum(forceDecompositionY);
            scalar forceDecompositionSumZ = sum(forceDecompositionZ);

            forAll(cellFaces, cellFacei)
            {
                label meshFacei = cellFaces[cellFacei];

                vector weighedForce = interfaceForce_[celli];
                if (forceDecompositionSumX != 0.0)
                {
                    weighedForce.x() *= mag(forceDecompositionX[cellFacei] / forceDecompositionSumX);
                }
                if (forceDecompositionSumY != 0.0)
                {
                    weighedForce.y() *= mag(forceDecompositionY[cellFacei] / forceDecompositionSumY);
                }
                if (forceDecompositionSumZ != 0.0)
                {
                    weighedForce.z() *= mag(forceDecompositionZ[cellFacei] / forceDecompositionSumZ);
                }

                if (mesh_.isInternalFace(meshFacei))
                {
                    intForceIf[meshFacei] += weighedForce;   
                }
                else
                {
                    label patchi = mesh_.boundaryMesh().whichPatch(meshFacei);
                    vectorField& intForcePf = intForceBf[patchi];
                    const polyPatch& pp = mesh_.boundaryMesh()[patchi];
                    if (!isA<emptyPolyPatch>(pp))
                    {
                        label patchFacei = meshFacei - pp.start();

                        intForcePf[patchFacei] += weighedForce;
                    }
                }
            }
        }
    }
    
    //Synchronize across processors
    syncTools::syncFaceList(mesh_, allSurfForcef,plusEqOp<vector>());

    interfaceForcef_ += intForce;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tensileForceMethod::tensileForceMethod
(
    const volScalarField& alpha1,
    const volVectorField& U,
    const volScalarField& rho,
    const reconstructionSchemes& surf,
    const IOdictionary& dict
)
:
    mesh_(alpha1.mesh()),
    transportPropertiesDict_(dict),
    sigmaPtr_(surfaceTensionModel::New(dict, mesh_)),
    deltaN_
    (
        "deltaN",
        1e-8/cbrt(average(mesh_.V()))
    ),
    alpha1_(alpha1),
    U_(U),
    rho_(rho),
    sigma_
    (
        IOobject
        (
            "sigma",
            alpha1_.time().timeName(),
            mesh_
        ),
        sigmaPtr_->sigma().ref()
    ),
    SfHat_
    (
        IOobject
        (
            "SfHat",
            alpha1_.time().timeName(),
            mesh_
        ),
        mesh_.Sf()/mesh_.magSf()
    ),
    Fsigma_
    (
        IOobject
        (
            "Fsigma",
            alpha1_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimMass/(dimArea*sqr(dimTime)), Zero)
    ),
    Fpj_
    (
        IOobject
        (
            "Fpj",
            alpha1_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimMass/(dimArea*sqr(dimTime)), Zero)
    ),
    interfaceForce_
    (
        IOobject
        (
            "interfaceForce",
            alpha1_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimMass/(dimArea*sqr(dimTime)), Zero)
    ),
    interfaceForcef_
    (
        IOobject
        (
            "interfaceForcef",
            alpha1_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimMass/(dimArea*sqr(dimTime)), Zero)
    ),
    surf_(surf),
    surfNormal_(surf.normal()),
    surfCentre_(surf.centre())
{
    reconstructedDistanceFunction RDF(mesh_);
    RDF.markCellsNearSurf(surf_.interfaceCell(), 1);
    nextToInterface_ = RDF.nextToInterface();

    correctContactAngle();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField>
Foam::tensileForceMethod::nearInterface() const
{
    return pos0(alpha1_ - 0.01)*pos0(0.99 - alpha1_);
}

void Foam::tensileForceMethod::interfaceForce
(
    const DynamicList<DynamicList<label>>& connectedInterfaceLabels
)
{
    correctContactAngle();
    
    reconstructedDistanceFunction RDF(mesh_);
    RDF.markCellsNearSurf(surf_.interfaceCell(), 1);
    nextToInterface_ = RDF.nextToInterface();
    SfHat_ = mesh_.Sf()/mesh_.magSf();

    List<scalar> interfacePressureJump(connectedInterfaceLabels.size(), Zero);

    calculateTensileForce(connectedInterfaceLabels, interfacePressureJump);
    pressureJumpCorrection(connectedInterfaceLabels, interfacePressureJump);
    mapInterfaceForceDensity();
    decomposeForceToCellFaces();
}

bool Foam::tensileForceMethod::read()
{
    sigmaPtr_->readDict(transportPropertiesDict_);

    return true;
}


// ************************************************************************* //
