/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::HAMexternalHeatFluxFvPatchScalarField

Description
    This boundary condition supplies a heat flux condition based on:
    - convective heat flux
    - latent and sensible heat transfer due to vapor exchange
    - radiative heat flux
    - sensible heat flux due to rain

    Example usage:
        patchName
        {
            type            HAMexternalHeatFlux;
            value           uniform 300.0;
        }
    
    Input files are provided under subfolder 0/patchName/
    - Tambient: ambient air temperature
    - alpha: convective heat transfer coefficient
    - rad: net raditive flux
    - beta: convective moisture transfer coefficient
    - pv_o: exterior vapor pressure
    - gl: rain flux
    - rainTemp: rain water temperature

SourceFiles
    HAMexternalHeatFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef HAMexternalHeatFluxFvPatchScalarField_H
#define HAMexternalHeatFluxFvPatchScalarField_H

#include "mixedFvPatchFields.H"

#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class HAMexternalHeatFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class HAMexternalHeatFluxFvPatchScalarField
:
    public mixedFvPatchScalarField
{
public:

private:

    // Private data
    
        //- table time vs ambient temperature [K]
        interpolationTable<scalar> Tambient;

        //- table time vs convective heat transfer coefficient [W/m2K]
        interpolationTable<scalar> alpha;
        
        //- table time vs radiative heat flux [W/m2]
        interpolationTable<scalar> rad;   
        
        //- table time vs convective moisture transfer coefficient [s/m]
        interpolationTable<scalar> beta;
        
        //- table time vs exterior vapor pressure [Pa]
        interpolationTable<scalar> pv_o;
        
        //- table time vs rain flux [kg/m2s]
        interpolationTable<scalar> gl;
        
        //- table time vs rain temperature [K]
        interpolationTable<scalar> rainTemp;                                     

public:

    //- Runtime type information
    TypeName("HAMexternalHeatFlux");


    // Constructors

        //- Construct from patch and internal field
        HAMexternalHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        HAMexternalHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // HAMexternalHeatFluxFvPatchScalarField
        // onto a new patch
        HAMexternalHeatFluxFvPatchScalarField
        (
            const HAMexternalHeatFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        HAMexternalHeatFluxFvPatchScalarField
        (
            const HAMexternalHeatFluxFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new HAMexternalHeatFluxFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        HAMexternalHeatFluxFvPatchScalarField
        (
            const HAMexternalHeatFluxFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new HAMexternalHeatFluxFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Mapping functions


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
