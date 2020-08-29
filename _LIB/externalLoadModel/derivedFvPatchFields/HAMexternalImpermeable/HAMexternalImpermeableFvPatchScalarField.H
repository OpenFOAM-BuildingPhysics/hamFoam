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
    Foam::HAMexternalImpermeableFvPatchScalarField

Description
    This boundary condition supplies a moisture flux condition based on:
    - convective moisture flux
    - rain flux

    Example usage:
        patchName
        {
            type            HAMexternalImpermeable;
            value           uniform -1E7;
        }
    
    Input files are provided under subfolder 0/patchName/
    - Tambient: ambient air temperature
    - beta: convective moisture transfer coefficient
    - pv_o: exterior vapor pressure
    - gl: rain flux

SourceFiles
    HAMexternalImpermeableFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef HAMexternalImpermeableFvPatchScalarField_H
#define HAMexternalImpermeableFvPatchScalarField_H

#include "mixedFvPatchFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class HAMexternalImpermeableFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class HAMexternalImpermeableFvPatchScalarField
:
    public mixedFvPatchScalarField
{
public:

private:

    // Private data


public:

    //- Runtime type information
    TypeName("HAMexternalImpermeable");


    // Constructors

        //- Construct from patch and internal field
        HAMexternalImpermeableFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        HAMexternalImpermeableFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // HAMexternalImpermeableFvPatchScalarField
        // onto a new patch
        HAMexternalImpermeableFvPatchScalarField
        (
            const HAMexternalImpermeableFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        HAMexternalImpermeableFvPatchScalarField
        (
            const HAMexternalImpermeableFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new HAMexternalImpermeableFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        HAMexternalImpermeableFvPatchScalarField
        (
            const HAMexternalImpermeableFvPatchScalarField&,
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
                new HAMexternalImpermeableFvPatchScalarField(*this, iF)
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
