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

\*---------------------------------------------------------------------------*/

#ifndef VoFBoilingIncompressible_H
#define VoFBoilingIncompressible_H

#include "fvModel.H"
#include "boilingIncompressibleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class incompressibleTwoPhaseVoFMixture;

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class VoFBoilingIncompressible Declaration
\*---------------------------------------------------------------------------*/

class VoFBoilingIncompressible
:
    public fvModel
{
    // Private Data

        //- Reference to the mixture properties
        const incompressibleTwoPhaseVoFMixture& mixture_;

        //- The BoilingIncompressible model
        autoPtr<boilingIncompressibleModel> boilingIncompressible_;


public:

    //- Runtime type information
    TypeName("VoFBoilingIncompressible");


    // Constructors

        //- Construct from explicit source name and mesh
        VoFBoilingIncompressible
        (
            const word& sourceName,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        VoFBoilingIncompressible
        (
            const VoFBoilingIncompressible&
        ) = delete;


    // Member Functions

        //- Return the list of fields for which the option adds source term
        //  to the transport equation
        virtual wordList addSupFields() const;

        //- Add a source to the phase continuity equation
        virtual void addSup
        (
            const volScalarField& alpha,
            fvMatrix<scalar>& eqn
        ) const;

        //- Add a source to the mixture momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            const volVectorField& U,
            fvMatrix<vector>& eqn
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


        //- Correct the BoilingIncompressible model
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const VoFBoilingIncompressible&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
