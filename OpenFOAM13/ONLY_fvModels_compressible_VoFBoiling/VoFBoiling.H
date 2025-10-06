/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

Class
    Foam::fv::compressible::VoFboiling

Description
    boiling fvModel

Usage
    Example usage:
    \verbatim
    VoFboiling
    {
        type    VoFboiling;

        model   SchnerrSauer;

        KunzCoeffs
        {
            pSat    2300;   // Saturation pressure

            UInf    20.0;
            tInf    0.005; // L = 0.1 m
            Cc      1000;
            Cv      1000;
        }

        MerkleCoeffs
        {
            pSat    2300;   // Saturation pressure

            UInf    20.0;
            tInf    0.005;  // L = 0.1 m
            Cc      80;
            Cv      1e-03;
        }

        SchnerrSauerCoeffs
        {
            pSat    2300;   // Saturation pressure

            n       1.6e+13;
            dNuc    2.0e-06;
            Cc      1;
            Cv      1;
        }
    }
    \endverbatim

SourceFiles
    VoFboiling.C

\*---------------------------------------------------------------------------*/

#ifndef VoFboiling_H
#define VoFboiling_H

#include "fvModel.H"
#include "boilingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class compressibleTwoPhaseVoFMixture;

namespace fv
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
                      Class VoFboiling Declaration
\*---------------------------------------------------------------------------*/

class VoFBoiling
:
    public fvModel
{
    // Private Data

        //- Reference to the mixture properties
        const compressibleTwoPhaseVoFMixture& mixture_;

        //- The boiling model
        autoPtr<Foam::compressible::boilingModel> boiling_;


public:

    //- Runtime type information
    TypeName("compressible::VoFBoiling");


    // Constructors

        //- Construct from explicit source name and mesh
        VoFBoiling
        (
            const word& sourceName,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        VoFBoiling
        (
            const VoFBoiling&
        ) = delete;


    // Member Functions

        //- Return the list of fields for which the option adds source term
        //  to the transport equation
        virtual wordList addSupFields() const;

        //- Add a source to the phase continuity equation
        virtual void addSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvMatrix<scalar>& eqn
        ) const;


        // Mesh changes

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);

            //- Update for mesh motion
            virtual bool movePoints();


        //- Correct the boiling model
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const VoFBoiling&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
