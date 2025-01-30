/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenCFD Ltd.
     \\/     M anipulation  |
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

namespace Foam
{
//namespace functionObjects
//{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class VOFParticle>
class sumParticleOp
{
    public:
    eulerianParticleMod operator()
    (
        const eulerianParticleMod& p0,
        const eulerianParticleMod& p1
    ) const
    {
        if ((p0.faceIHit != -1) && (p1.faceIHit == -1))
        {
            return p0;
        }
        else if ((p0.faceIHit == -1) && (p1.faceIHit != -1))
        {
            return p1;
        }
        else if ((p0.faceIHit != -1) && (p1.faceIHit != -1))
        {
            // Choose particle with the largest collected volume and
            // accumulate total volume
            if (p0.V > p1.V)
            {
                eulerianParticleMod p = p0;
                p.V = p0.V + p1.V;
                p.VC = p0.VC + p1.VC;
                p.VU = p0.VU + p1.VU;
                p.VT = p0.VT + p1.VT;
                return p;
            }
            else
            {
                eulerianParticleMod p = p1;
                p.V = p0.V + p1.V;
                p.VC = p0.VC + p1.VC;
                p.VU = p0.VU + p1.VU;
                p.VT = p0.VT + p1.VT;
                return p;
            }
        }
        else
        {
            ///eulerianParticleMod p;
            //return p;
            if (p0.V > p1.V)
            {
                eulerianParticleMod p = p0;
                p.V = p0.V + p1.V;
                p.VC = p0.VC + p1.VC;
                p.VU = p0.VU + p1.VU;
                p.VT = p0.VT + p1.VT;
                return p;
            }
            else
            {
                eulerianParticleMod p = p1;
                p.V = p0.V + p1.V;
                p.VC = p0.VC + p1.VC;
                p.VU = p0.VU + p1.VU;
                p.VT = p0.VT + p1.VT;
                return p;
            }
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
