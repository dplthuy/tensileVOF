/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::regionSplit2Dv1906

Description
    Splits a patch into regions based on a mask field.  Result is a globally
    consistent label list of region index per patch face.

SourceFiles
    regionSplit2Dv1906.C

\*---------------------------------------------------------------------------*/

#ifndef regionSplit2Dv1906_H
#define regionSplit2Dv1906_H

#include "DynamicList.H"
#include "boolList.H"
#include "labelList.H"
#include "indirectPrimitivePatch.H"
#include "patchEdgeFaceRegion.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class faceZone;
class polyMesh;

/*---------------------------------------------------------------------------*\
                        Class regionSplit2Dv1906 Declaration
\*---------------------------------------------------------------------------*/

class regionSplit2Dv1906
:
    public labelList
{
    // Private data

        //- Number of regions
        label nRegions_;


    // Private Member Functions

        //- No copy construct
        regionSplit2Dv1906(const regionSplit2Dv1906&) = delete;

        //- No copy assignment
        void operator=(const regionSplit2Dv1906&) = delete;


public:

    // Constructors

        //- Construct from mesh and list of blocked faces
        regionSplit2Dv1906
        (
            const polyMesh& mesh,
            const indirectPrimitivePatch& patch,
            const boolList& blockedFaces,
            const label offset = 0
        );


    //- Destructor
    ~regionSplit2Dv1906();


    // Member Functions

        //- Return the global number of regions
        label nRegions() const
        {
            return nRegions_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
