/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::functionObjects::dynamicMultiDimRefineFvMeshFunctionObject

Description
    Specialization of Foam::functionObject for an Foam::fvMesh, providing a
    reference to the Foam::fvMesh.

    If the selected region is not an Foam::fvMesh a Foam::FatalError will be
    generated.

See also
    Foam::regionFunctionObject
    Foam::functionObject

SourceFiles
    dynamicMultiDimRefineFvMeshFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_dynamicMultiDimRefineFvMeshFunctionObject_H
#define functionObjects_dynamicMultiDimRefineFvMeshFunctionObject_H

#include "regionFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class dynamicMultiDimRefineFvMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class dynamicMultiDimRefineFvMeshFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class dynamicMultiDimRefineFvMeshFunctionObject
:
    public regionFunctionObject
{
protected:

    // Protected Member Data

        //- Reference to the fvMesh
        const dynamicMultiDimRefineFvMesh& mesh_;


    // Protected Member Functions

        //- No copy construct
        dynamicMultiDimRefineFvMeshFunctionObject(const dynamicMultiDimRefineFvMeshFunctionObject&) = delete;

        //- No copy assignment
        void operator=(const dynamicMultiDimRefineFvMeshFunctionObject&) = delete;


public:

    //- Runtime type information
    TypeName("dynamicMultiDimRefineFvMeshFunctionObject");


    // Constructors

        //- Construct from Time and dictionary
        dynamicMultiDimRefineFvMeshFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from the region objectRegistry and dictionary
        dynamicMultiDimRefineFvMeshFunctionObject
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );


    //- Destructor
    virtual ~dynamicMultiDimRefineFvMeshFunctionObject() = default;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
