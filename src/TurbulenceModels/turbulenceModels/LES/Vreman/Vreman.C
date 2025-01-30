/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "Vreman.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
volScalarField Vreman<BasicTurbulenceModel>::aBeta
(
    const tmp<volTensorField>& tmpGradU
) const
{
    volTensorField gradU = tmpGradU();
    
    volScalarField d1u1(gradU.component(tensor::XX));
    volScalarField d2u1(gradU.component(tensor::YX));
    volScalarField d3u1(gradU.component(tensor::ZX));
    
    volScalarField d1u2(gradU.component(tensor::XY));
    volScalarField d2u2(gradU.component(tensor::YY));
    volScalarField d3u2(gradU.component(tensor::ZY));
    
    volScalarField d1u3(gradU.component(tensor::ZX));
    volScalarField d2u3(gradU.component(tensor::ZY));
    volScalarField d3u3(gradU.component(tensor::ZZ));
    
    volScalarField aBeta(d1u1*d1u1 + d1u2*d1u2 + d1u3*d1u3 + d2u1*d2u1 + d2u2*d2u2 + d2u3*d2u3 + d3u1*d3u1 + d3u2*d3u2 + d3u3*d3u3);

    return aBeta;
}

template<class BasicTurbulenceModel>
volScalarField Vreman<BasicTurbulenceModel>::bBeta
(
    const tmp<volTensorField>& tmpGradU
) const
{
    volTensorField gradU = tmpGradU();
    
    volScalarField d1u1(gradU.component(tensor::XX));
    volScalarField d2u1(gradU.component(tensor::YX));
    volScalarField d3u1(gradU.component(tensor::ZX));
    
    volScalarField d1u2(gradU.component(tensor::XY));
    volScalarField d2u2(gradU.component(tensor::YY));
    volScalarField d3u2(gradU.component(tensor::ZY));
    
    volScalarField d1u3(gradU.component(tensor::ZX));
    volScalarField d2u3(gradU.component(tensor::ZY));
    volScalarField d3u3(gradU.component(tensor::ZZ));
    
    volScalarField delta2(this->delta() * this->delta());
    
    volScalarField beta11(delta2*(d1u1*d1u1 + d2u1*d2u1 + d3u1*d3u1));
    volScalarField beta12(delta2*(d1u1*d1u2 + d2u1*d2u2 + d3u1*d3u2));
    volScalarField beta13(delta2*(d1u1*d1u3 + d2u1*d2u3 + d3u1*d3u3));
    volScalarField beta22(delta2*(d1u2*d1u2 + d2u2*d2u2 + d3u2*d3u2));
    volScalarField beta23(delta2*(d1u2*d1u3 + d2u2*d2u3 + d3u2*d3u3));
    volScalarField beta33(delta2*(d1u3*d1u3 + d2u3*d2u3 + d3u3*d3u3));

    volScalarField bBeta(beta11*beta22 - beta12*beta12 + beta11*beta33 - beta13*beta13 + beta22*beta33 - beta23*beta23);

    return bBeta;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> Vreman<BasicTurbulenceModel>::k
(
    const tmp<volTensorField>& gradU,
    const volScalarField& nut
) const
{
    volSymmTensorField S(symm(gradU));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("k", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->mesh_
            ),
            2*nut*sqrt(S&&S)
        )
    );
}

template<class BasicTurbulenceModel>
void Vreman<BasicTurbulenceModel>::correctNut()
{
    
    volScalarField aBeta(this->aBeta(fvc::grad(this->U_)));
    volScalarField bBeta(this->bBeta(fvc::grad(this->U_)));
    
    forAll(this->nut_, cellI)
    {
    	if (aBeta[cellI] == 0 || bBeta[cellI] < 1e-8) 
    	{
    		this->nut_[cellI] = 0;
    	} else {
    		this->nut_[cellI] = 0.07225 * sqrt(bBeta[cellI]/aBeta[cellI]);
    	}
    }
    
    //volScalarField k(this->k(fvc::grad(this->U_), this->nut_));
    k_ = this->k(fvc::grad(this->U_), this->nut_);

    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
    
    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Vreman<BasicTurbulenceModel>::Vreman
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    LESeddyViscosity<BasicTurbulenceModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    
    k_
    (
    	IOobject
    	(
    	    IOobject::groupName("k", this->alphaRhoPhi_.group()),
    	    this->runTime_.timeName(),
    	    this->mesh_,
    	    IOobject::NO_READ,
    	    IOobject::AUTO_WRITE
    	),
    	this->mesh_
    ),

    CVreman_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "CVreman",
            this->coeffDict_,
            0.07225				
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool Vreman<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        CVreman_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> Vreman<BasicTurbulenceModel>::epsilon() const
{
    volScalarField k(this->k(fvc::grad(this->U_), this->nut_));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->mesh_
            ),
            this->Ce_*k*sqrt(k)/this->delta()
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> Vreman<BasicTurbulenceModel>::omega() const
{
    volScalarField k(this->k(fvc::grad(this->U_), this->nut_));
    volScalarField epsilon(this->Ce_*k*sqrt(k)/this->delta());

    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("omega", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_
        ),
        epsilon/(0.09*k)
    );
}


template<class BasicTurbulenceModel>
void Vreman<BasicTurbulenceModel>::correct()
{
    LESeddyViscosity<BasicTurbulenceModel>::correct();
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
