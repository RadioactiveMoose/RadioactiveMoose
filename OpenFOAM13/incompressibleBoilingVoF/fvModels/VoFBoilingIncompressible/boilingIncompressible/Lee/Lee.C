/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "Lee.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace boilingIncompressibleModels
{
    defineTypeNameAndDebug(Lee, 0);
    addToRunTimeSelectionTable(boilingIncompressibleModel, Lee, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boilingIncompressibleModels::Lee::Lee
(
    const dictionary& dict,
    const incompressibleTwoPhases& phases
)
:
    boilingIncompressibleModel(dict, phases),
    tInf_("tInf", dimTime, dict),
    Cv_("Cv", dimless, dict), //changer la dimension pour la passer en 1/t ? pour le modèle simplifie, oui
    Cc_("Cc", dimless, dict), //changer la dimension pour la passer en 1/t ?
    T0_("T0", dimTemperature, 0.01),               // Small temperature difference

    // Calculate coefficients once (like Kunz does)
    mcCoeff_(Cc_ * rhov() / tInf_),  // dimensionedScalar
    mvCoeff_(Cv_ * rhol() / tInf_)   // dimensionedScalar
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::boilingIncompressibleModels::Lee::mDotcvAlpha() const
{
    // Get current temperature field
    const volScalarField& TField = phases_.mesh().lookupObject<volScalarField>("T");
    const volScalarField::Internal& T = TField();

    const volScalarField::Internal limitedAlphal
    (
        min(max(alphal(), scalar(0)), scalar(1))
    );

    const volScalarField::Internal limitedAlphav
    (
        min(max(alphav(), scalar(0)), scalar(1))
    );

    // Get saturation temperatures
    
    const volScalarField::Internal Tsat
    (
	volScalarField::Internal::New
	(
		"Tsat",
		phases_.mesh(),
		Tsat_  // This is the dimensionedScalar from base class
    	)
    );
    
    // Lee model implementation
    // Condensation: when T < T_sat (vapor -> liquid)
    // Vaporisation: when T > T_sat (liquid -> vapor)

    return Pair<tmp<volScalarField::Internal>>
    (
    	// Condensation rate: rho_v * Cc * alpha_v * max(T_sat - T, 0) / T_sat
    	mcCoeff_ * limitedAlphav * max(Tsat - T, T0_) / max(Tsat, T0_),  // Use Tsat
    	// Vaporisation rate: rho_l * Cv * alpha_l * max(T - T_sat, 0) / T_sat  
    	mvCoeff_ * limitedAlphal * max(T - Tsat, T0_) / max(Tsat, T0_)   // Use Tsat
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::boilingIncompressibleModels::Lee::mDotcvT() const
{
    // Get current temperature field
    const volScalarField::Internal& T =
        phases_.mesh().lookupObject<volScalarField>("T")();

    const volScalarField::Internal limitedAlphal
    (
        min(max(alphal(), scalar(0)), scalar(1))
    );

    const volScalarField::Internal limitedAlphav
    (
        min(max(alphav(), scalar(0)), scalar(1))
    );


    // Get saturation temperatures
    const volScalarField::Internal Tsat
    (
	volScalarField::Internal::New("Tsat", phases_.mesh(), Tsat_)
    );

    // Return coefficients to multiply (T - T_sat) terms
    return Pair<tmp<volScalarField::Internal>>
    (
        // Condensation rate coefficient: active when T < T_sat
        mcCoeff_ * limitedAlphav * neg(T - Tsat) / max(Tsat, T0_),
        
        // Vaporisation rate coefficient: active when T > T_sat
        mvCoeff_ * limitedAlphal * pos(T - Tsat) / max(Tsat, T0_)
    );
}


void Foam::boilingIncompressibleModels::Lee::correct()
{
    // Update any time-dependent parameters if needed
    // Currently no dynamic updates required for Lee model
}


bool Foam::boilingIncompressibleModels::Lee::read
(
    const dictionary& dict
)
{
    if (boilingIncompressibleModel::read(dict))
    {
/*
        dict.lookup("diffL") >> diffL_;
        dict.lookup("diffG") >> diffG_;
*/
        dict.lookup("tInf") >> tInf_;
        dict.lookup("Cv") >> Cv_;
        dict.lookup("Cc") >> Cc_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
