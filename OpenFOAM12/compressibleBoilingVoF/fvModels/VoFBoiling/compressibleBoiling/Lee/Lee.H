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

Class
    Foam::compressible::boilingModels::Lee

Description
    Lee boiling model for temperature-driven phase change.

    The Lee model calculates mass transfer rates based on temperature difference
    from saturation temperature, incorporating thermal diffusivity effects.
    
    Mass transfer rates:
    - Evaporation: rho_l * Cv * alpha_l * max(T - T_sat, 0) / T_sat
    - Condensation: rho_v * Cc * alpha_v * max(T_sat - T, 0) / T_sat

    Reference:
    \verbatim
    Lee, W.H. (1979). A Pressure Iteration Scheme for Two-Phase Flow Modeling.
    Multiphase Transport Fundamentals, Reactor Safety, Applications
    \endverbatim

Usage:
    \table
        Property | Description                     | Required | Default value
        liquid   | Name of the liquid phase             | yes      |
        TSat     | Saturation vapor temperature         | yes      |
        Uinf     | Free-stream velocity                 | yes      |
        tInf     | Vaporisation/condensation time scale | yes      |
        Cv       | Vapourisation rate coefficient       | yes      |
        Cc       | Condensation rate coefficient        | yes      |
    \endtable

    Example:
    \verbatim
        model       Lee;

        liquid      liquid;

        TSat
        {
            type        constant;
            value       373.15;
        }

        LeeCoeffs
        {
            diffL       1.4e-7;     // m^2/s
            diffG       2.0e-5;     // m^2/s
            tInf        0.001;      // s
            Cv          100;        // 1/s
            Cc          100;        // 1/s
        }
    \endverbatim

SourceFiles
    Lee.C

\*---------------------------------------------------------------------------*/

#ifndef Lee_H
#define Lee_H

#include "boilingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace boilingModels
{

/*---------------------------------------------------------------------------*\
                                 Class Lee
\*---------------------------------------------------------------------------*/

class Lee
:
    public boilingModel
{
    // Private Data

/* Pas besoin de la vitesse
        //- Free-stream velocity
        dimensionedScalar UInf_;
       
	//- Thermal difusivity term for liquid
	dimensionedScalar diffL_;
	
	//- Thermal difusivity term for gas
	dimensionedScalar diffG_;
*/ 	
        //- Vaporisation time scale
        dimensionedScalar tInf_;

// Cv_ and Cc_ need to be divided by tInf to have the same dimension as x_i in the Lee model provided my Mudawar et al.
        
        //- Vaporisation rate coefficient [1]
        dimensionedScalar Cv_;

        //- Condensation rate coefficient [1]
        dimensionedScalar Cc_;

        //- Small temperature difference for numerical stability [K]
        dimensionedScalar T0_;

    // Private Member functions

        tmp<volScalarField::Internal> mvCoeff() const;

        tmp<volScalarField::Internal> mcCoeff() const;


public:

    //- Runtime type information
    TypeName("Lee");


    // Constructors

        //- Construct for phases
        Lee
        (
            const dictionary& dict,
            const compressibleTwoPhases& phases
        );


    //- Destructor
    virtual ~Lee()
    {}

    // Member Functions

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply alphav for the condensation rate and a
        //  coefficient to multiply alphal for the vaporisation rate
        virtual Pair<tmp<volScalarField::Internal>> mDotcvAlphal() const;

        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply (T - TSat)
        virtual Pair<tmp<volScalarField::Internal>> mDotcvT() const;

        //- Correct the Lee phaseChange model
        virtual void correct();

        //- Read the dictionary and update
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace boilingModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
