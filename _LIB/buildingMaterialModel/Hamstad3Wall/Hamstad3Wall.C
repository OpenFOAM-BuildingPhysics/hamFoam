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

#include "Hamstad3Wall.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace buildingMaterialModels
{
    defineTypeNameAndDebug(Hamstad3Wall, 0);

    addToRunTimeSelectionTable
    (
        buildingMaterialModel,
        Hamstad3Wall,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::buildingMaterialModels::Hamstad3Wall::Hamstad3Wall
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
void Foam::buildingMaterialModels::Hamstad3Wall::update_w_C_cell(const volScalarField& pc, volScalarField& w, volScalarField& Crel, label& celli)
{
    List<scalar> reta; reta.setSize(2); reta[0] = 6.0e-3; reta[1] = 1.2e-2;
    List<scalar> retn; retn.setSize(2); retn[0] = 2.4883; retn[1] = 2.3898;
    List<scalar> retm; retm.setSize(2); retm[0] = -0.59812; retm[1] = -0.58155;
    List<scalar> retw; retw.setSize(2); retw[0] = 0.41; retw[1] = 0.59;
    scalar w_tmp = 0; scalar tmp = 0; scalar C_tmp = 0; scalar tmp2 = 0;
    for (int i = 0; i <= 1; i++)
    {
        tmp = pow( (reta[i]*pc.internalField()[celli]) , retn[i] );
        w_tmp = w_tmp + retw[i] / ( pow( (1 + tmp) , retm[i] ));
        tmp2 = pow( (1 + tmp) , retm[i] );
        C_tmp = C_tmp - retw[i]/tmp2 * retm[i]*retn[i]*tmp/((1 + tmp)*pc.internalField()[celli]); 
    } 
    w.ref()[celli] = w_tmp*871;   
    Crel.ref()[celli] = mag( C_tmp*871 );   
}

//- Correct the buildingMaterial liquid permeability (cell)
void Foam::buildingMaterialModels::Hamstad3Wall::update_Krel_cell(const volScalarField& pc, const volScalarField& w, volScalarField& Krel, label& celli)

{
    List<scalar> a; a.setSize(6); a[0] = -46.245; a[1] = 294.506; a[2] = -1439; a[3] = 3249; a[4] = -3370; a[5] = 1305;
    scalar rho_w = 1.0e3;
    scalar tmp = 0;
    for (int i = 0; i <= 6; i++)
    {
        tmp = tmp + (a[i] * pow((w.internalField()[celli] / rho_w), i));
    }
    Krel.ref()[celli] = exp(tmp);
}

//- Correct the buildingMaterial vapor permeability (cell)
void Foam::buildingMaterialModels::Hamstad3Wall::update_Kvap_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_v, volScalarField& K_pt, label& celli)
{
    scalar rho_l = 1.0e3; 
    scalar R_v = 8.31451 * 1000 / (18.01534);
    scalar L_v = 2.5e6;
    scalar mu_dry = 5.6;
    scalar p = 0.2;
    scalar w_sat = 871;

    scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3/T.internalField()[celli] - 5.976*Foam::log(T.internalField()[celli])); // saturation vapour pressure [Pa]
    scalar relhum = Foam::exp(pc.internalField()[celli]/(rho_l*R_v*T.internalField()[celli])); // relative humidity [-]

    scalar tmp = 1 - (w.internalField()[celli] / w_sat);
    scalar delta = 2.61e-5 * tmp / (R_v*T.internalField()[celli] * mu_dry * ((1-p)*tmp*tmp + p)); // Water vapour diffusion coefficient "for Wall3" [s]

    K_v.ref()[celli] = (delta*p_vsat*relhum) / (rho_l*R_v*T.internalField()[celli]);
    K_pt.ref()[celli] = ( (delta*p_vsat*relhum)/(rho_l*R_v*pow(T.internalField()[celli],2)) ) * (rho_l*L_v - pc.internalField()[celli]);

}

//*********************************************************** //
