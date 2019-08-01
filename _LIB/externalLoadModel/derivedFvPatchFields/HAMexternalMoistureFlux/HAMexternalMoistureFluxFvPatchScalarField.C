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

\*---------------------------------------------------------------------------*/

#include "HAMexternalMoistureFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "uniformDimensionedFields.H"

#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HAMexternalMoistureFluxFvPatchScalarField::
HAMexternalMoistureFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


Foam::HAMexternalMoistureFluxFvPatchScalarField::
HAMexternalMoistureFluxFvPatchScalarField
(
    const HAMexternalMoistureFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::HAMexternalMoistureFluxFvPatchScalarField::
HAMexternalMoistureFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF)
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


Foam::HAMexternalMoistureFluxFvPatchScalarField::
HAMexternalMoistureFluxFvPatchScalarField
(
    const HAMexternalMoistureFluxFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf)
{}


Foam::HAMexternalMoistureFluxFvPatchScalarField::
HAMexternalMoistureFluxFvPatchScalarField
(
    const HAMexternalMoistureFluxFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::HAMexternalMoistureFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalar rhol=1.0e3; scalar Rv=8.31451*1000/(18.01534);   

    scalarField& pcp = *this;

    scalarField Ts(pcp.size(), 0.0);
        Ts = patch().lookupPatchField<volScalarField, scalar>("Ts");
    scalarField Krel(pcp.size(), 0.0);
        Krel = patch().lookupPatchField<volScalarField, scalar>("Krel");
    scalarField K_v(pcp.size(), 0.0);
        K_v = patch().lookupPatchField<volScalarField, scalar>("K_v");
    const fvPatchScalarField&
        fieldTs = refCast
            <const fvPatchScalarField>
            (
                patch().lookupPatchField<volScalarField, scalar>("Ts")
            );
            
    const fvPatchScalarField&
        fieldpc = refCast
            <const fvPatchScalarField>
            (
                patch().lookupPatchField<volScalarField, scalar>("pc")
            );            

    const polyPatch& p = this->patch().patch();
    const polyMesh& mesh = p.boundaryMesh().mesh();
    Time& time = const_cast<Time&>(mesh.time());
    interpolationTable<scalar> Tambient
    (
        "$FOAM_CASE/0/Tambient"
    ); 
    interpolationTable<scalar> beta
    (
        "$FOAM_CASE/0/beta"
    ); 
    interpolationTable<scalar> pv_o
    (
        "$FOAM_CASE/0/pv_o"
    ); 
    interpolationTable<scalar> gl
    (
        "$FOAM_CASE/0/gl"
    ); 

    scalarField pvsat_s = exp(6.58094e1-7.06627e3/Ts-5.976*log(Ts));
    scalarField pv_s = pvsat_s*exp((pcp)/(rhol*Rv*Ts));

    scalarField g_conv = beta(time.value())*(pv_o(time.value())-pv_s);
    scalarField g_cond = (Krel+K_v)*fieldpc.snGrad();
    scalar gl_ = gl(time.value());

    // term with temperature gradient:
    scalarField K_pt(pcp.size(), 0.0);
        K_pt = patch().lookupPatchField<volScalarField, scalar>("K_pt");                 
    scalarField X = K_pt*fieldTs.snGrad();
    //////////////////////////////////

    //-- Gravity flux --//
    //lookup gravity vector
    uniformDimensionedVectorField g = db().lookupObject<uniformDimensionedVectorField>("g");
    scalarField gn = g.value() & patch().nf();

    scalarField phiG = Krel*rhol*gn;
    //////////////////////////////////

    refGrad() = (g_conv + gl_ + phiG - X)/(Krel+K_v);
    refValue() =  -500.0 + 1.0;
    
    forAll(valueFraction(),faceI)
    {
        if(pcp[faceI] > -500.0 && gl_ > 0.0 && (gl_ > g_cond[faceI] - g_conv[faceI] - phiG[faceI] + X[faceI]) )
        {
            valueFraction()[faceI] = 1.0;
        }
        else{valueFraction()[faceI] = 0.0;}
    }
    //valueFraction() = 0.0;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::HAMexternalMoistureFluxFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    //os.writeKeyword("pc")<< pcName_ << token::END_STATEMENT << nl;
    //os.writeKeyword("T")<< TName_ << token::END_STATEMENT << nl;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        HAMexternalMoistureFluxFvPatchScalarField
    );
}

// ************************************************************************* //
