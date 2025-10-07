/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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

#include "incompressibleBoilingVoF.H"
#include "localEulerDdtScheme.H"
#include "fvCorrectPhi.H"
#include "geometricZeroField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(incompressibleBoilingVoF, 0);
    addToRunTimeSelectionTable(solver, incompressibleBoilingVoF, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::incompressibleBoilingVoF::incompressibleBoilingVoF(fvMesh& mesh)
:
    twoPhaseVoFSolver
    (
        mesh,
        autoPtr<twoPhaseVoFMixture>(new incompressibleTwoPhaseVoFMixture(mesh))
    ),

    mixture
    (
        refCast<incompressibleTwoPhaseVoFMixture>(twoPhaseVoFSolver::mixture)
    ),

    p
    (
        IOobject
        (
            "p",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*buoyancy.gh
    ),

    pressureReference_
    (
        p,
        p_rgh,
        pimple.dict()
    ),

    momentumTransport
    (
        U,
        phi,
        alphaPhi1,
        alphaPhi2,
        mixture
    )
{
    if (correctPhi || mesh.topoChanging())
    {
        rAU = new volScalarField
        (
            IOobject
            (
                "rAU",
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimTime/dimDensity, 1)
        );
    }

    if (!runTime.restart() || !divergent())
    {
        correctUphiBCs(U_, phi_, true);

        fv::correctPhi
        (
            phi_,
            U,
            p_rgh,
            rAU,
            autoPtr<volScalarField>(),
            pressureReference(),
            pimple
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::incompressibleBoilingVoF::~incompressibleBoilingVoF()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::incompressibleBoilingVoF::prePredictor()
{
    twoPhaseVoFSolver::prePredictor();

    const dimensionedScalar& rho1 = mixture.rho1();
    const dimensionedScalar& rho2 = mixture.rho2();

    // Calculate the mass-flux
    rhoPhi = alphaPhi1*rho1 + alphaPhi2*rho2;
}


void Foam::solvers::incompressibleBoilingVoF::momentumTransportPredictor()
{
    momentumTransport.predict();
}


void Foam::solvers::incompressibleBoilingVoF::thermophysicalTransportPredictor()
{}


void Foam::solvers::incompressibleBoilingVoF::pressureCorrector()
{
    incompressiblePressureCorrector(p);
}


void Foam::solvers::incompressibleBoilingVoF::thermophysicalPredictor()
{}


void Foam::solvers::incompressibleBoilingVoF::momentumTransportCorrector()
{
    momentumTransport.correct();
}


void Foam::solvers::incompressibleBoilingVoF::thermophysicalTransportCorrector()
{}


// ************************************************************************* //
