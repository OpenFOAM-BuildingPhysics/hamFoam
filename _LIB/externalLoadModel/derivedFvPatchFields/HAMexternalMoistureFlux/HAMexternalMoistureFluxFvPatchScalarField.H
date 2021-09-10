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
    Foam::HAMexternalMoistureFluxFvPatchScalarField

Description
    This boundary condition supplies a moisture flux condition based on:
    - convective moisture flux
    - rain flux

    Example usage:
        patchName
        {
            type            HAMexternalMoistureFlux;
            value           uniform -1E7;
        }
    
    Input files are provided under subfolder 0/patchName/
    - Tambient: ambient air temperature
    - beta: convective moisture transfer coefficient
    - pv_o: exterior vapor pressure
    - gl: rain flux

SourceFiles
    HAMexternalMoistureFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef HAMexternalMoistureFluxFvPatchScalarField_H
#define HAMexternalMoistureFluxFvPatchScalarField_H

#include "mixedFvPatchFields.H"

#include "TableFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class HAMexternalMoistureFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class HAMexternalMoistureFluxFvPatchScalarField
:
    public mixedFvPatchScalarField
{
public:

private:

    // Private data  
        
        //- table time vs convective moisture transfer coefficient [s/m]
        const autoPtr<Function1s::TableFile<scalar>> beta;
        
        //- table time vs exterior vapor pressure [Pa]
        const autoPtr<Function1s::TableFile<scalar>> pv_o;
        
        //- table time vs rain flux [kg/m2s]
        const autoPtr<Function1s::TableFile<scalar>> gl; 

public:

    //- Runtime type information
    TypeName("HAMexternalMoistureFlux");


    // Constructors

        //- Construct from patch and internal field
        HAMexternalMoistureFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        HAMexternalMoistureFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // HAMexternalMoistureFluxFvPatchScalarField
        // onto a new patch
        HAMexternalMoistureFluxFvPatchScalarField
        (
            const HAMexternalMoistureFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        HAMexternalMoistureFluxFvPatchScalarField
        (
            const HAMexternalMoistureFluxFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new HAMexternalMoistureFluxFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        HAMexternalMoistureFluxFvPatchScalarField
        (
            const HAMexternalMoistureFluxFvPatchScalarField&,
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
                new HAMexternalMoistureFluxFvPatchScalarField(*this, iF)
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
