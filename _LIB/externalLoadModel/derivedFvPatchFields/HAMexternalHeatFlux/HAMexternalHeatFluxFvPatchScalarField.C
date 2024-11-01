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

#include "HAMexternalHeatFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HAMexternalHeatFluxFvPatchScalarField::
HAMexternalHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    Tambient(),
    alpha(),
    rad(),
    beta(),
    pv_o(),
    gl(),
    rainTemp()
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


Foam::HAMexternalHeatFluxFvPatchScalarField::
HAMexternalHeatFluxFvPatchScalarField
(
    const HAMexternalHeatFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::HAMexternalHeatFluxFvPatchScalarField::
HAMexternalHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    Tambient(new Function1s::Table<scalar>("Tambient", dict.subDict("Tambient"))),
    alpha(new Function1s::Table<scalar>("alpha", dict.subDict("alpha"))),
    rad(new Function1s::Table<scalar>("rad", dict.subDict("rad"))),
    beta(new Function1s::Table<scalar>("beta", dict.subDict("beta"))),        
    pv_o(new Function1s::Table<scalar>("pv_o", dict.subDict("pv_o"))),        
    gl(new Function1s::Table<scalar>("gl", dict.subDict("gl"))),        
    rainTemp(new Function1s::Table<scalar>("rainTemp", dict.subDict("rainTemp")))                
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


Foam::HAMexternalHeatFluxFvPatchScalarField::
HAMexternalHeatFluxFvPatchScalarField
(
    const HAMexternalHeatFluxFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::HAMexternalHeatFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalar rhol=1.0e3; scalar Rv=8.31451*1000/(18.01534);   
    scalar cap_v = 1880; scalar Tref = 273.15; scalar L_v = 2.5e6; scalar cap_l = 4182;

    scalarField& Tp = *this;
    
    scalarField lambda_m(Tp.size(), 0.0);
        lambda_m = patch().lookupPatchField<volScalarField, scalar>("lambda_m");  
    scalarField K_pt(Tp.size(), 0.0);
        K_pt = patch().lookupPatchField<volScalarField, scalar>("K_pt"); 
    scalarField Krel(Tp.size(), 0.0);
        Krel = patch().lookupPatchField<volScalarField, scalar>("Krel");
    scalarField K_v(Tp.size(), 0.0);
        K_v = patch().lookupPatchField<volScalarField, scalar>("K_v");                    
    scalarField pc(Tp.size(), 0.0);
        pc = patch().lookupPatchField<volScalarField, scalar>("pc");
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

    scalarField q_conv = alpha->value(time.value())*(Tambient->value(time.value())-Tp);

    scalarField pvsat_s = exp(6.58094e1-7.06627e3/Tp-5.976*log(Tp));
    scalarField pv_s = pvsat_s*exp((pc)/(rhol*Rv*Tp));

    scalarField g_conv = beta->value(time.value())*(pv_o->value(time.value())-pv_s);    
    scalarField LE = (cap_v*(Tp-Tref)+L_v)*g_conv;//Latent and sensible heat transfer due to vapor exchange    
    
    //-- Gravity-enthalpy flux --//
    //lookup gravity vector
    uniformDimensionedVectorField g = db().lookupObject<uniformDimensionedVectorField>("g");
    scalarField gn = g.value() & patch().nf();

    scalarField phiG = Krel*rhol*gn;
    scalarField phiGT = (cap_l*(Tp-Tref))*phiG;

    // term with capillary moisture gradient:                          
    scalarField X = ((cap_l*(Tp-Tref)*Krel)+(cap_v*(Tp-Tref)+L_v)*K_v)*fieldpc.snGrad();
    // moisture flux term with temperature gradient:               
    scalarField Xmoist = K_pt*fieldTs.snGrad();
    //////////////////////////////////      

    scalarField CR(Tp.size(), 0.0);    
    scalar gl_ = gl->value(time.value());
    if(gl_ > 0)
    {
        //scalarField g_cond = (Krel+K_v)*fieldpc.snGrad();
        scalarField g_cond = (Krel+K_v)*(-10.0-fieldpc.patchInternalField())*patch().deltaCoeffs();       
        forAll(CR,faceI)
        {
            scalar rainFlux = 0;
            //if(pc[faceI] > -100.0 && (gl_ > g_cond[faceI] - g_conv[faceI] - phiG[faceI] + Xmoist[faceI]) )
            if( (gl_ > g_cond[faceI] - g_conv[faceI] - phiG[faceI] + Xmoist[faceI]) )
            {
                rainFlux = g_cond[faceI] - g_conv[faceI] - phiG[faceI] + Xmoist[faceI];
            }
            else
            {
                rainFlux = gl_;
            }
            CR[faceI] = rainFlux * cap_l*(rainTemp->value(time.value()) - Tref);
        }
    }    

    if(fieldpc.type() == "HAMexternalImpermeable")
    {
        valueFraction() = 0;
        refValue() = 0;
        refGrad() = (q_conv + rad->value(time.value()))/(lambda_m);
    }
    else
    {    
        refGrad() = (q_conv + LE + rad->value(time.value()) + CR + phiGT - X)/(lambda_m+(cap_v*(Tp-Tref)+L_v)*K_pt);
        refValue() =  0;
        valueFraction() = 0.0;
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::HAMexternalHeatFluxFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    writeEntry(os, Tambient());
    writeEntry(os, alpha());
    writeEntry(os, rad());    
    writeEntry(os, beta());
    writeEntry(os, pv_o());
    writeEntry(os, gl());
    writeEntry(os, rainTemp());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        HAMexternalHeatFluxFvPatchScalarField
    );
}

// ************************************************************************* //
