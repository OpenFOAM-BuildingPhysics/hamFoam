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

#include "MichelaBrick.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace buildingMaterialModels
{
    defineTypeNameAndDebug(MichelaBrick, 0);

    addToRunTimeSelectionTable
    (
        buildingMaterialModel,
        MichelaBrick,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::buildingMaterialModels::MichelaBrick::MichelaBrick
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
void Foam::buildingMaterialModels::MichelaBrick::update_w_C_cell(const volScalarField& pc, volScalarField& w, volScalarField& Crel, label& celli)
{
    List<scalar> reta; reta.setSize(3); reta[0]=-8e-6; reta[1]=-1e-6; reta[2]=-5e-7; //reta[3]=-6.5e-4; reta[4]=-6.5e-4; 
    List<scalar> retn; retn.setSize(3); retn[0]=2.6; retn[1]=1.55; retn[2]=3; //retn[3]=4.00; retn[4]=4.00;
    List<scalar> retm; retm.setSize(3); retm[0]=0.615384615; retm[1]=0.35483871; retm[2]=0.66666667; //retm[3]=0.75; retm[4]=0.75;
    List<scalar> retw; retw.setSize(3); retw[0]=0.9; retw[1]=0.05; retw[2]=0.05; //retw[3]=0; retw[4]=0; //the last 2 are zero for wetting retention curve
    scalar w_tmp = 0; scalar tmp = 0; scalar C_tmp = 0; scalar tmp2 = 0;    
    for (int i=0; i<=2; i++)
    {
        tmp = pow( (reta[i]*pc.internalField()[celli]) , retn[i] );
        w_tmp = w_tmp + retw[i] / ( pow( (1 + tmp) , retm[i] ));
        tmp2 = pow( (1 + tmp) , retm[i] );
        C_tmp = C_tmp - retw[i]/tmp2 * retm[i]*retn[i]*tmp/((1 + tmp)*pc.internalField()[celli]); 
    } 
    w.ref()[celli] = w_tmp*250;   
    Crel.ref()[celli] = mag( C_tmp*250);   
}

//- Correct the buildingMaterial liquid permeability (cell)
void Foam::buildingMaterialModels::MichelaBrick::update_Krel_cell(const volScalarField& pc, const volScalarField& w, volScalarField& Krel, label& celli)
{
    scalar logpc = log10(-pc.internalField()[celli]);
    scalar logKl = 0;
    int i;
    double logpc_M[]={2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
        3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9,
        4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9,
        5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9,
        6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9,
        7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9,
        8.0}; 
    double logKl_M[]={-9.41411, -9.41411, -9.41411, -9.41411, -9.41411, -9.41411, -9.41411, -9.41411, -9.41411, -9.41411,
        -9.41411, -9.41411, -9.41411, -9.41411, -9.41411, -9.41411, -9.41411, -9.41411, -9.41411, -9.41411,
        -9.41411, -9.41411, -9.41411, -9.41411, -9.41411, -9.41411, -9.41411, -9.41411, -9.41411, -9.52628,
        -9.75476, -10.10900, -10.55747, -11.04712, -11.52984, -11.97805, -12.38269, -12.74402, -13.06338, -13.33836,
        -13.5635, -13.74067, -13.89334, -14.05679, -14.25228, -14.47765, -14.71888, -14.96252, -15.20012, -15.42765,
        -15.6439, -15.84929, -16.04488, -16.23211, -16.41243, -16.58719, -16.75754, -16.92446, -17.08872, -17.25093,
        -17.4116};

    if (logpc < scalar(2.0))
    {
        i = 0;
        logKl = logKl_M[i] + (((logKl_M[i+1] - logKl_M[i])/(logpc_M[i+1] - logpc_M[i]))*(logpc - logpc_M[i])) ;
    }
    else if (logpc >= scalar(8.0))
    {
        i = 59;
        logKl = logKl_M[i] + (((logKl_M[i+1] - logKl_M[i])/(logpc_M[i+1] - logpc_M[i]))*(logpc - logpc_M[i])) ;
    }
    else
    {
        for (i=0; i<=59; ++i)
        {
            if ( (logpc_M[i] <= logpc) && (logpc < logpc_M[i+1]) )
            {
                logKl = logKl_M[i] + (((logKl_M[i+1] - logKl_M[i])/(logpc_M[i+1] - logpc_M[i]))*(logpc - logpc_M[i])) ;
                break;
            }
        }
    }
    Krel.ref()[celli] = pow(10,logKl);
}

//- Correct the buildingMaterial vapor permeability (cell)
void Foam::buildingMaterialModels::MichelaBrick::update_Kv_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_v, label& celli)
{
    scalar rho_l = 1.0e3; 
    scalar R_v = 8.31451*1000/(18.01534); 

    scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3/T.internalField()[celli] - 5.976*Foam::log(T.internalField()[celli])); // saturation vapour pressure [Pa]
    scalar relhum = Foam::exp(pc.internalField()[celli]/(rho_l*R_v*T.internalField()[celli])); // relative humidity [-]
    
    scalar tmp = 1 - (w.internalField()[celli]/250); 
    scalar delta = 2.61e-5 * tmp/(R_v*T.internalField()[celli]*16*(0.503*tmp*tmp + 0.497)); // Water vapour diffusion coefficient [s]
    
    K_v.ref()[celli] = (delta*p_vsat*relhum)/(rho_l*R_v*T.internalField()[celli]);
}

//- Correct the buildingMaterial K_pt (cell)
void Foam::buildingMaterialModels::MichelaBrick::update_Kpt_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_pt, label& celli)
{
    scalar rho_l = 1.0e3; 
    scalar R_v = 8.31451*1000/(18.01534); 
    scalar L_v = 2.5e6;

    scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3/T.internalField()[celli] - 5.976*Foam::log(T.internalField()[celli])); // saturation vapour pressure [Pa]
    //scalar dpsatdt = (7.06627e3/(T.internalField()[celli]*T.internalField()[celli]) - 5.976/T.internalField()[celli]) * p_vsat; // saturation vapour pressure [Pa]
        
    scalar relhum = Foam::exp(pc.internalField()[celli]/(rho_l*R_v*T.internalField()[celli])); // relative humidity [-]
    
    scalar tmp = 1 - (w.internalField()[celli]/250); 
    scalar delta = 2.61e-5 * tmp/(R_v*T.internalField()[celli]*16*(0.503*tmp*tmp + 0.497)); // Water vapour diffusion coefficient [s]

    K_pt.ref()[celli] = ( (delta*p_vsat*relhum)/(rho_l*R_v*pow(T.internalField()[celli],2)) ) * (rho_l*L_v - pc.internalField()[celli]);
}

//*********************************************************** //
