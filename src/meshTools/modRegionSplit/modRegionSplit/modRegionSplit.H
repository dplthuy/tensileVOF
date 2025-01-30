/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::modRegionSplit

Description
    This class separates the mesh into distinct unconnected regions,
    each of which is then given a label according to globalNumbering().


    Say 6 cells, 3 processors, with single baffle on proc1.

    \verbatim
              baffle
                |
    +---+---+---+---+---+---+
    |   |   |   |   |   |   |
    +---+---+---+---+---+---+
      proc0 | proc1 | proc2

    \endverbatim


    \verbatim

    1: determine local regions (uncoupled)

    +---+---+---+---+---+---+
    | 0 | 0 | 0 | 1 | 0 | 0 |
    +---+---+---+---+---+---+
      proc0 | proc1 | proc2

    \endverbatim


    2: make global

    \verbatim

    +---+---+---+---+---+---+
    | 0 | 0 | 1 | 2 | 3 | 3 |
    +---+---+---+---+---+---+
      proc0 | proc1 | proc2

    \endverbatim


    3: merge connected across procs

    \verbatim

    +---+---+---+---+---+---+
    | 0 | 0 | 0 | 2 | 2 | 2 |
    +---+---+---+---+---+---+
      proc0 | proc1 | proc2

    \endverbatim


    4. determine locally owner regions.

    Determine compact numbering for the local regions and send these to
    all processors that need them:

    - proc0 uses regions:
      - 0 which is local to it.
    - proc1 uses regions
      - 0 which originates from proc0
      - 2 which is local to it
    - proc2 uses regions
      - 2 which originates from proc1

    So proc1 needs to get the compact number for region 0 from proc0 and proc2
    needs to get the compact number for region 2 from proc1:

    \verbatim

    +---+---+---+---+---+---+
    | 0 | 0 | 0 | 1 | 1 | 1 |
    +---+---+---+---+---+---+
      proc0 | proc1 | proc2

    \endverbatim

    Can optionally keep all regions local to the processor.

Note
    does not walk across cyclicAMI/cyclicACMI - since these are not
    \c coupled() at the patch level.

SourceFiles
    modRegionSplit.C

\*---------------------------------------------------------------------------*/

#ifndef modRegionSplit_H
#define modRegionSplit_H

#include "globalIndex.H"
#include "labelPair.H"
#include "bitSet.H"
#include "boolList.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;

/*---------------------------------------------------------------------------*\
                         Class modRegionSplit Declaration
\*---------------------------------------------------------------------------*/

class modRegionSplit
:
    public MeshObject<polyMesh, TopologicalMeshObject, modRegionSplit>,
    public labelList
{
    // Private Data

        //- Indexing into the regions
        globalIndex globalNumbering_;

        //- Temporary list of cells that have changed
        mutable DynamicList<label> changedCells_;

        //- Temporary list of faces that have changed
        mutable DynamicList<label> changedFaces_;

        //- Dynamic list of the global cells that make up a region
        DynamicList<DynamicList<label>> regionGlobalCellLabels_;

        //- Dynamic list of the local cells that make up a region
        DynamicList<DynamicList<label>> regionLocalCellLabels_;

    // Private Member Functions

        //- Check that boundary faces are synchronised, or Fatal
        void checkBoundaryFaceSync(const boolList& blockedFace) const;

        //- Update faceRegion data between (non-processor) coupled faces.
        void updateFacePair
        (
            const label face0,
            const label face1,
            labelList& faceRegion,
            DynamicList<label>& facesChanged
        ) const;

        //- Given a seed cell label, fill cellRegion/faceRegion with markValue
        //- for contiguous region around it
        void fillSeedMask
        (
            const UList<labelPair>& explicitConnections,
            const label seedCellID,
            const label markValue,
            labelList& cellRegion,
            labelList& faceRegion
        ) const;


        //- Calculate the local region split.
        //  \return number of processor-local regions,
        //      without consolidation between procesors
        label localRegionSplit
        (
            const UList<labelPair>& explicitConnections,
            labelList& cellRegion,
            labelList& faceRegion
        ) const;

        //- Manually consolidate regions globally by swapping information
        //  between processor domains and reducing the regions accordingly.
        //
        //  \return globalIndex into the local regions after reduction.
        globalIndex reduceRegionsImpl
        (
            const label numLocalRegions,
            const boolList& blockedFace,
            labelList& cellRegion
        ) const;

        //- Calculate the size (number of cells) of the different regions locally
        void cellRegionLabels
        (
            const polyMesh& mesh,
            const labelList& cellRegion
        );

        // - Calculate the size (number of cells) of the different regions globally from local information
        void cellRegionLabelsToGlobal();


public:

    //- Runtime type information
    ClassName("modRegionSplit");


    // Constructors

        //- Construct from mesh and whether face is blocked, optionally
        //- with additional explicit connections between
        //- normal boundary faces.
        //
        //  \note blockedFace must be consistent across coupled faces!
        modRegionSplit
        (
            const polyMesh& mesh,
            const boolList& blockedFace,
            const List<labelPair>& explicitConnections = List<labelPair>(),
            const bool doGlobalRegions = Pstream::parRun()
        );


    // Member Functions

        //- Return global region numbering
        const globalIndex& globalNumbering() const noexcept
        {
            return globalNumbering_;
        }

        //- Return local number of regions
        label nLocalRegions() const
        {
            return globalNumbering().localSize(Pstream::myProcNo());
        }

        //- Return total number of regions
        label nRegions() const
        {
            return globalNumbering().size();
        }

        //- Return global cell labels per region
        DynamicList<DynamicList<label>>& getGlobalRegionCellLabels()
        {
            cellRegionLabelsToGlobal();
            return regionGlobalCellLabels_;
        }

        //- Return local cell labels per region
        DynamicList<DynamicList<label>>& getLocalRegionCellLabels()
        {
            return regionLocalCellLabels_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
