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

Class
    Foam::alphaRegionProperties

Description
    Provides controls for the pressure reference is closed-volume simulations
    and a general method for limiting the pressure during the startup of
    steady-state simulations.

SourceFiles
    alphaRegionPropertiesI.H
    alphaRegionProperties.C

\*---------------------------------------------------------------------------*/

#ifndef alphaRegionProperties_H
#define alphaRegionProperties_H

#include "isoAdvection.H"
#include "reconstructionSchemes.H"
#include "volFields.H"
#include "modRegionSplit.H"
#include "eulerianParticleMod.H"
#include "writeFile.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                         Class alphaRegionProperties Declaration
\*---------------------------------------------------------------------------*/

class alphaRegionProperties
{
    // Private data

        //- Number of cells lower limit
        label NCellsMin_;

        //- Switch to remove regions with size below NCellsMin
        bool removeRegions_;

        //- Switch for calculating region properties
        bool calculateRegionProperties_;

        //- Threshold for marking cell as blocked
        scalar alphaThreshold_;

        //- Total removed volume
        scalar removedVolumeTotal_;

        // Data field properties
        word UName_;

        word rhoName_;

        word TName_;

        // Pointer to the mesh
        const fvMesh& mesh_;

        //- Write interval for the cloud data (runTime)
        scalar cloudWriteInterval_;

        //- Write time index for the cloud data
        label cloudWriteTimeIndex_;

        //- Volume fraction field
        volScalarField& alpha_;

        //- Geometric VoF surface access
        reconstructionSchemes& surf_;

        //- Block a face if its owner and neighbor are not both in the same alpha-phase. Separately for gas and liquid
        void setBlockedFaces
        (
            boolList& liquidBlockedFaces
        );

        //- Remove alpha regions if they have size below NCellsMin
        void removeAlphaRegions
        (
            DynamicList<DynamicList<label>>& alphaRegionLocalLabels
        );

        void getInterfaceLabelsFromRegions
        (
            const DynamicList<DynamicList<label>>& regionLabels,
            DynamicList<DynamicList<label>>& interfaceLabels
        );

        //- Calculate the properties of the alpha regions
        void calculateAlphaRegionProperties
        (
            const DynamicList<DynamicList<label>>& alphaRegionLabels
        );

        //- Accumulate the properties for a region based on the cell labels
        void accumulateRegionProperties
        (
            const DynamicList<DynamicList<label>>& alphaRegionLabels,
            List<Foam::eulerianParticleMod>& alphaRegions
        );

        //- Calculate the diameter of the region for a sphere with equal volume
        void volumeEquivalentDiameter
        (
            const List<Foam::eulerianParticleMod>& alphaRegions,
            scalarList& regionD
        );

         //- Find the minimum and maximum radii of the alphaRegion
        void minMaxRadii
        (
            const DynamicList<DynamicList<label>>& alphaRegionLabels,
            const List<Foam::eulerianParticleMod>& alphaRegions,
            scalarList& regionRMin,
            scalarList& regionRMax
        );

        //- Calculate the velocity of the continuous phase surrounding the alphaRegion
        void continuousPhaseVelocity
        (
            const List<Foam::eulerianParticleMod>& alphaRegions,
            const scalarList& regionD,
            List<vector>& contPhaseU
        );

        //- write the region properties
        void writeRegionProperties
        (
            const List<Foam::eulerianParticleMod>& alphaRegions,
            const scalarList& regionD,
            const scalarList& regionRMin,
            const scalarList& regionRMax,
            const List<vector>& contPhaseU
        );

public:

    // Constructors

        //- Construct from the SIMPLE/PIMPLE sub-dictionary
        alphaRegionProperties
        (
            const dictionary& dict,
            isoAdvection& advector,
            volScalarField& alpha
        );


    // Member Functions

        //- Remove structures and collect their properties based on dictionary input
        void execute
        (
            DynamicList<DynamicList<label>>& interfaceLabels,
            bool firstIter
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
