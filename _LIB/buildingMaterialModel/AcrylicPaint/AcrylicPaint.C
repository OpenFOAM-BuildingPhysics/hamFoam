/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "AcrylicPaint.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace buildingMaterialModels
{
    defineTypeNameAndDebug(AcrylicPaint, 0);

    addToRunTimeSelectionTable
    (
        buildingMaterialModel,
        AcrylicPaint,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::buildingMaterialModels::AcrylicPaint::AcrylicPaint
(
    const word& name,
    const dictionary& buildingMaterialDict,
    const word& cellZoneModel
)
:
    buildingMaterialModel(name, buildingMaterialDict, cellZoneModel)
{
    
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//- Correct the buildingMaterial moisture content (cell)
void Foam::buildingMaterialModels::AcrylicPaint::update_w_C_cell(const volScalarField& pc, volScalarField& w, volScalarField& Crel, label& celli)
{
    w.ref()[celli] = 0;   
    Crel.ref()[celli] = SMALL;   
}

//- Correct the buildingMaterial liquid permeability (cell)
void Foam::buildingMaterialModels::AcrylicPaint::update_Krel_cell(const volScalarField& pc, const volScalarField& w, volScalarField& Krel, label& celli)
{

    Krel.ref()[celli] = 0;

}

//- Correct the buildingMaterial vapor permeability (cell)
void Foam::buildingMaterialModels::AcrylicPaint::update_Kvap_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_v, volScalarField& K_pt, label& celli)

{
    scalar rho_l = 1.0e3;
    scalar R_v = 8.31451 * 1000 / (18.01534);
    scalar L_v = 2.5e6;
    scalar a = 4.68e-14;
    scalar b = 2.52e-15;
    scalar c = 5.65;

    scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3 / T.internalField()[celli] - 5.976*Foam::log(T.internalField()[celli])); // saturation vapour pressure [Pa]
    scalar relhum = Foam::exp(pc.internalField()[celli] / (rho_l*R_v*T.internalField()[celli])); // relative humidity [-]

    scalar delta = a + b*Foam::exp(c*relhum); // Water vapour diffusion coefficient "for acrylic paint" [s]

    K_v.ref()[celli] = (delta*p_vsat*relhum) / (rho_l*R_v*T.internalField()[celli]);
    K_pt.ref()[celli] = ( (delta*p_vsat*relhum)/(rho_l*R_v*pow(T.internalField()[celli],2)) ) * (rho_l*L_v - pc.internalField()[celli]);

}

//*********************************************************** //
