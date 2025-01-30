/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::LESModels::Vreman

Group
    grpLESTurbulence

Description
    The Smagorinsky SGS model.

    Reference:
    \verbatim
        Vreman, A. (2004)
        An eddy-viscosity subgrid-scale model for turbulent shear flow:
        Algebraic theory and applications
        Physics of Fluids, 16(10), 3670-3681.
    \endverbatim

    The form of the Vreman model implemented is obtained from the
    k-equation model assuming local equilibrium which provides estimates of both
    k and epsilon:
    For k calculations, see Vreman(2004) pp.3677 left column.

    \verbatim
        k = 2*nut*|S|

    where

        S = symm(grad(U));
        nut = c * sqrt(B_beta / (alpha_ij * alpha_ij))
    \endverbatim

    The default model coefficients are
    \verbatim
        SmagorinskyCoeffs
        {
            CVreman             0.07225;
            Ce                  1.048;
        }
    \endverbatim

See also
    Foam::LESModels::kEqn

SourceFiles
    Vreman.C

\*---------------------------------------------------------------------------*/

#ifndef Vreman_H
#define Vreman_H

#include "LESModel.H"
#include "LESeddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class Vreman Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class Vreman
:
    public LESeddyViscosity<BasicTurbulenceModel>
{
    // Private Member Functions

        //- No copy construct
        Vreman(const Vreman&) = delete;

        //- No copy assignment
        void operator=(const Vreman&) = delete;


protected:

    // Protected data
	
	//Fields
	volScalarField k_;
	
	// Cvreman = 2.5 * Csmagorinsky^2; Csmagorinsky = 0.17
        dimensionedScalar CVreman_;

    // Protected Member Functions

        //- Return aBeta
        //  calculated from the given velocity gradient
        volScalarField aBeta(const tmp<volTensorField>& tmpGradU) const;
        
        //- Return bBeta
        //  calculated from the given velocity gradient
        volScalarField bBeta(const tmp<volTensorField>& tmpGradU) const;
        
        //- Return SGS kinetic energy
        //  calculated from the given velocity gradient
        tmp<volScalarField> k(const tmp<volTensorField>& gradU, const volScalarField& nut) const;
        

        //- Update the SGS eddy viscosity
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("Vreman");


    // Constructors

        //- Construct from components
        Vreman
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~Vreman() = default;


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
            //return k(fvc::grad(this->U_), this->nut_);
        }

        //- Return sub-grid dissipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Return sub-grid specific dissipation rate
        virtual tmp<volScalarField> omega() const;

        //- Correct Eddy-Viscosity and related properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Vreman.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
