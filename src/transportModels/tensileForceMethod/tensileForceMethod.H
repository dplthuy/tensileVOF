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

Class
    Foam::tensileForceMethod

Description
    Contains the interface properties.

    Properties to aid interFoam:
    -# Correct the alpha boundary condition for dynamic contact angle.
    -# Calculate interface curvature.

SourceFiles
    tensileForceMethod.C

\*---------------------------------------------------------------------------*/

#ifndef tensileForceMethod_H
#define tensileForceMethod_H

#include "IOdictionary.H"
#include "surfaceTensionModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "reconstructionSchemes.H"
#include "cutFacePLIC.H"
#include "zoneDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class tensileForceMethod Declaration
\*---------------------------------------------------------------------------*/

class tensileForceMethod
{
    // Private data
        //- Reference to mesh
        const fvMesh& mesh_;

        //- Keep a reference to the transportProperties dictionary
        const dictionary& transportPropertiesDict_;

        //- Surface tension
        autoPtr<surfaceTensionModel> sigmaPtr_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

        const volScalarField& alpha1_;  
        const volVectorField& U_;       
        const volScalarField& rho_;     
        volScalarField sigma_;          
        surfaceVectorField SfHat_;      

        // Vector field of the surface tension force 
        volVectorField Fsigma_;

        // Vector field of the pressure jump force 
        volVectorField Fpj_;

        // Vector fields of resulting force density at the interface
        volVectorField interfaceForceDensity_;
        surfaceVectorField interfaceForceDensityf_;
        
        //geometric vof surface access
        const reconstructionSchemes& surf_;
        const volVectorField& surfNormal_;
        const volVectorField& surfCentre_;

        // Is cell next to interface cell?
        boolList nextToInterface_;

    // Private Member Functions

        //- No copy construct
        tensileForceMethod(const tensileForceMethod&) = delete;

        //- No copy assignment
        void operator=(const tensileForceMethod&) = delete;

        //- Correction for the boundary condition on the unit normal nHat on
        //  walls to produce the correct contact dynamic angle
        //  calculated from the component of U parallel to the wall
        void correctContactAngle() const;

        void calculateTensileForce
        (
            const DynamicList<DynamicList<label>>& connectedInterfaceLabels,
            List<scalar>& interfacePressureJump
        );

        vector plicTensileForce
        (
            const label celli,
            const Map<vector>& mapSurfNormal,
            const Map<vector>& mapCellCentre,
            zoneDistribute& CPCDistribute
        );

        vector antiClockWisePlicTangent
        (
            const DynamicList<point>& facePlicPts,
            const label celli
        );

        bool neighbourNormal
        (
            const label celli,
            const label facei,
            vector& nbNormal
        );
        
        // Get the normal vector of neighbor element from CPC stencil. 
        //  - Considers the position of the nb element with respect to the considered face
        //  - Only selects a nb normal if the angle between the normals is 0 < angle < 90 degree
        vector neighbourNormalFromCPCStencil
        (
            const label celli,
            const label facei,
            const DynamicList<point>& facePlicPts,
            const Map<vector>& mapSurfNormal,
            const Map<vector>& mapCellCentre,
            zoneDistribute& CPCDistribute
        );
        
        // For unbalanced elements, have the surface tension force point in the same direction as the normal
        vector balanceTensileForce
        (
            const label celli,
            const List<vector>& tensileForce,
            const List<scalar>& edgeLength,
            const scalar surfTensionCoeff
        ) const;

        void pressureJumpCorrection
        (
            const DynamicList<DynamicList<label>>& connectedInterfaceLabels,
            const List<scalar>& interfacePressureJump
        );

        void mapInterfaceForce();

        void decomposeForceToCellFaces();

        
public:

    // Constructors

        //- Construct from volume fraction field gamma and IOdictionary
        tensileForceMethod
        (
            const volScalarField& alpha1,
            const volVectorField& U,
            const volScalarField& rho,
            const reconstructionSchemes& surf,
            const IOdictionary&
        );


    // Member Functions

        const dimensionedScalar& deltaN() const
        {
            return deltaN_;
        }

        tmp<volVectorField> surfaceTensionForceDensity() const
        {
            return interfaceForceDensity_;
        }

        tmp<surfaceScalarField> cellFaceSurfaceTensionForceDensity() const
        {
            return interfaceForceDensityf_ & SfHat_;
        }

        
      
        //- Indicator of the proximity of the interface
        //  Field values are 1 near and 0 away for the interface.
        tmp<volScalarField> nearInterface() const;

        void interfaceForce
        (
            const DynamicList<DynamicList<label>>& connectedInterfaceLabels
        );

        //- Read transportProperties dictionary
        bool read();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
