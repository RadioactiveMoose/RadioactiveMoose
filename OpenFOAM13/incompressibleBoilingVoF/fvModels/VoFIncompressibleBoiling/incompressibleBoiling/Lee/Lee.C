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
namespace boilingModels
{
    defineTypeNameAndDebug(Lee, 0);
    addToRunTimeSelectionTable(boilingModel, Lee, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boilingModels::Lee::Lee
(
    const dictionary& dict,
    const incompressibleTwoPhases& phases
)
:
    boilingModel(dict, phases),
/* NOT NEEDED in simplified Lee Model - i.e. Cv_ = diffL/A , A being area of contact
    diffL_("diffL", dimArea/dimTime, dict),
    diffG_("diffG", dimArea/dimTime, dict),
*/
    tInf_("tInf", dimTime, dict),
    Cv_("Cv", dimless, dict), //changer la dimension pour la passer en 1/t ? pour le modèle simplifie, oui
    Cc_("Cc", dimless, dict), //changer la dimension pour la passer en 1/t ?
    T0_("T0", dimTemperature, 0.01)                // Small temperature difference
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::boilingModels::Lee::mvCoeff() const
{
    // Vaporisation coefficient: rho_l * Cv / tInf
    // Using thermal diffusivity scaling
    return Cv_ * rhol() / tInf_;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::boilingModels::Lee::mcCoeff() const
{
    // Condensation coefficient: rho_v * Cc / tInf  
    // Using thermal diffusivity scaling
    return Cc_ * rhov() / tInf_;
}


Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::boilingModels::Lee::mDotcvAlphal() const
{
    // Get current temperature field
    const volScalarField& TField = phases_.mesh().lookupObject<volScalarField>("T");
    const volScalarField::Internal& T = TField();

    const volScalarField::Internal alphal
    (
        min(max(this->alphal(), scalar(0)), scalar(1))
    );

    const volScalarField::Internal alphav
    (
        min(max(this->alphav(), scalar(0)), scalar(1))
    );

    // Get saturation temperatures
    const tmp<volScalarField::Internal> tTsatl = this->Tsatl();
    const volScalarField::Internal& Tsatl = tTsatl();
    const tmp<volScalarField::Internal> tTsatv = this->Tsatv();
    const volScalarField::Internal& Tsatv = tTsatv();

    // Lee model implementation
    // Condensation: when T < T_sat (vapor -> liquid)
    // Vaporisation: when T > T_sat (liquid -> vapor)

    return Pair<tmp<volScalarField::Internal>>
    (
        // Condensation rate: rho_v * Cc * alpha_v * max(T_sat - T, 0) / T_sat
        mcCoeff() * alphav * max(Tsatv - T, T0_) / max(Tsatv, T0_),
        
        // Vaporisation rate: rho_l * Cv * alpha_l * max(T - T_sat, 0) / T_sat  
        mvCoeff() * alphal * max(T - Tsatl, T0_) / max(Tsatl, T0_)
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::boilingModels::Lee::mDotcvT() const
{
    // Get current temperature field
    const volScalarField::Internal& T =
        phases_.mesh().lookupObject<volScalarField>("T")();

    const volScalarField::Internal alphal
    (
        min(max(this->alphal(), scalar(0)), scalar(1))
    );

    const volScalarField::Internal alphav
    (
        min(max(this->alphav(), scalar(0)), scalar(1))
    );

    // Get saturation temperatures
    const volScalarField::Internal Tsatl(this->Tsatl());
    const volScalarField::Internal Tsatv(this->Tsatv());

    // Return coefficients to multiply (T - T_sat) terms
    return Pair<tmp<volScalarField::Internal>>
    (
        // Condensation rate coefficient: active when T < T_sat
        mcCoeff() * alphav * neg(T - Tsatv) / max(Tsatv, T0_),
        
        // Vaporisation rate coefficient: active when T > T_sat
        mvCoeff() * alphal * pos(T - Tsatl) / max(Tsatl, T0_)
    );
}


void Foam::boilingModels::Lee::correct()
{
    // Update any time-dependent parameters if needed
    // Currently no dynamic updates required for Lee model
}


bool Foam::boilingModels::Lee::read
(
    const dictionary& dict
)
{
    if (boilingModel::read(dict))
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
