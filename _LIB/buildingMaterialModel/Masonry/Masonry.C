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

#include "Masonry.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace buildingMaterialModels
{
    defineTypeNameAndDebug(Masonry, 0);

    addToRunTimeSelectionTable
    (
        buildingMaterialModel,
		Masonry,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::buildingMaterialModels::Masonry::Masonry
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
void Foam::buildingMaterialModels::Masonry::update_w_C_cell(const volScalarField& pc, volScalarField& w, volScalarField& Crel, label& celli)
{
	List<scalar> reta_mortar; reta_mortar.setSize(2); reta_mortar[0] = -2.99e-8; reta_mortar[1] = -2.93e-7;
	List<scalar> retn_mortar; retn_mortar.setSize(2); retn_mortar[0] = 2.14; retn_mortar[1] = 3.24;
	List<scalar> retm_mortar; retm_mortar.setSize(2); retm_mortar[0] = 0.533; retm_mortar[1] = 0.691;
	List<scalar> retw_mortar; retw_mortar.setSize(2); retw_mortar[0] = 0.657; retw_mortar[1] = 0.343;
	
	List<scalar> reta_brick; reta_brick.setSize(2); reta_brick[0] = -2.1e-5; reta_brick[1] = -1.8e-5;
	List<scalar> retn_brick; retn_brick.setSize(2); retn_brick[0] = 4.1674; retn_brick[1] = 1.1698;
	List<scalar> retm_brick; retm_brick.setSize(2); retm_brick[0] = 0.76; retm_brick[1] = 0.915;
	List<scalar> retw_brick; retw_brick.setSize(2); retw_brick[0] = 0.29899; retw_brick[1] = 0.70101;
		
	scalar w_tmp_mortar = 0; scalar tmp_mortar = 0; scalar C_tmp_mortar = 0; scalar tmp2_mortar = 0;
	scalar w_tmp_brick = 0; scalar tmp_brick = 0; scalar C_tmp_brick = 0; scalar tmp2_brick = 0;
	
	for (int i = 0; i <= 1; i++)
	{
		tmp_mortar = pow((reta_mortar[i] * pc.internalField()[celli]), retn_mortar[i]);
		w_tmp_mortar = w_tmp_mortar + retw_mortar[i] / (pow((1 + tmp_mortar), retm_mortar[i]));
		tmp_brick = pow((reta_brick[i] * pc.internalField()[celli]), retn_brick[i]);
		w_tmp_brick = w_tmp_brick + retw_brick[i] / (pow((1 + tmp_brick), retm_brick[i]));		
		
		tmp2_mortar = pow((1 + tmp_mortar), retm_mortar[i]);
		C_tmp_mortar = C_tmp_mortar - retw_mortar[i] / tmp2_mortar * retm_mortar[i] * retn_mortar[i] * tmp_mortar / ((1 + tmp_mortar)*pc.internalField()[celli]);
		tmp2_brick = pow((1 + tmp_brick), retm_brick[i]);
		C_tmp_brick = C_tmp_brick - retw_brick[i] / tmp2_brick * retm_brick[i] * retn_brick[i] * tmp_brick / ((1 + tmp_brick)*pc.internalField()[celli]);		
	}
	
	scalar f_brick = 0.79;
	scalar f_mortar = 0.20;
	scalar f_crack = 0.01;
	
	w.ref()[celli] = ( (w_tmp_mortar * 97 * f_mortar) + (w_tmp_brick * 220 * f_brick) ) * (1-f_crack);
	Crel.ref()[celli] = mag( ( (C_tmp_mortar * 97 * f_mortar) + (C_tmp_brick * 220 * f_brick) ) * (1-f_crack) );
}

//- Correct the buildingMaterial liquid permeability (cell)
void Foam::buildingMaterialModels::Masonry::update_Krel_cell(const volScalarField& pc, const volScalarField& w, volScalarField& Krel, label& celli)
{
	scalar logpc = log10(-pc.internalField()[celli]);
	scalar logKl = 0;
	int i;
	List<scalar> logpc_M = { 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 
	                     3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 
	                     4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 
	                     5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 
	                     6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 
	                     7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 
	                     8.0, 8.1, 8.2, 8.3, 8.4 };
	
	List<scalar> logKl_M = { -8.6507724, -8.6507724, -8.6507724, -8.6507724, -8.6507724, -8.6507724, -8.6507724, -8.6507724, -8.6507724, -8.6507724, 
	                     -8.6507724, -8.6507724, -8.6507724, -8.6507724, -8.6507724, -8.6507724, -8.6507724, -8.6507724, -8.6507724, -8.6507724, 
	                     -8.6507724, -8.6507724, -8.6507724, -8.6507724, -8.6507724, -8.6507724, -8.6507724, -8.8798882, -9.6137666, -10.6223043, 
	                     -11.3272055, -11.7872541, -12.0896977, -12.3138365, -12.5080185, -12.6906902, -12.8663704, -13.0349722, -13.1952711, -13.3461295, 
	                     -13.4869154, -13.6175814, -13.7386963, -13.87408711, -14.00852236, -14.19588285, -14.40779296, -14.63630777, -14.85683619, -15.06708781, 
	                     -15.25560639, -15.42764354, -15.54980586, -15.66970026, -15.7942771, -15.92665579, -16.06639617, -16.21117188, -16.35908419, -16.50977421, 
						 -16.66444212, -16.82430118, -16.98779925, -17.14681363, -17.28937462 };

	if (logpc < logpc_M.first())
	{
		i = 0;
		logKl = logKl_M[i] + (((logKl_M[i + 1] - logKl_M[i]) / (logpc_M[i + 1] - logpc_M[i]))*(logpc - logpc_M[i]));
	}
	else if (logpc >= logpc_M.last())
	{
		i = logpc_M.size()-2;
		logKl = logKl_M[i] + (((logKl_M[i + 1] - logKl_M[i]) / (logpc_M[i + 1] - logpc_M[i]))*(logpc - logpc_M[i]));
	}
	else
	{
		for (i = 0; i < logpc_M.size()-1; ++i)
		{
			if ((logpc_M[i] <= logpc) && (logpc < logpc_M[i + 1]))
			{
				logKl = logKl_M[i] + (((logKl_M[i + 1] - logKl_M[i]) / (logpc_M[i + 1] - logpc_M[i]))*(logpc - logpc_M[i]));
				break;
			}
		}
	}
	Krel.ref()[celli] = pow(10, logKl);
}

//- Correct the buildingMaterial vapor permeability (cell)
void Foam::buildingMaterialModels::Masonry::update_Kvap_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_v, volScalarField& K_pt, label& celli)
{
	scalar rho_l = 1.0e3;
	scalar R_v = 8.31451 * 1000 / (18.01534);  
    scalar L_v = 2.5e6;
    	
	scalar a = 2.61e-5;
	scalar b = 307;
	scalar c = 0.503;
	scalar d = 0.497;
	
	scalar delta_a = 1.91e-10;

	scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3 / T.internalField()[celli] - 5.976*Foam::log(T.internalField()[celli])); // saturation vapour pressure [Pa]
	scalar relhum = Foam::exp(pc.internalField()[celli] / (rho_l*R_v*T.internalField()[celli])); // relative humidity [-]

	scalar f_brick = 0.79;
	scalar f_mortar = 0.20;
	scalar f_crack = 0.01;
	
	scalar wsat = ( (97 * f_mortar) + (220 * f_brick) ) * (1-f_crack);

    scalar tmp = 1 - (w.internalField()[celli]/wsat); 
	scalar delta_masonry = a * tmp / (R_v*T.internalField()[celli] * b * (c*tmp*tmp + d)); // Water vapour diffusion coefficient "for Masonry" [s]
    scalar delta_masonry_crack = delta_a*f_crack + delta_masonry*(1-f_crack);
	
    K_v.ref()[celli] = (delta_masonry_crack*p_vsat*relhum)/(rho_l*R_v*T.internalField()[celli]); 
	K_pt.ref()[celli] = ((delta_masonry_crack*p_vsat*relhum) / (rho_l*R_v*pow(T.internalField()[celli], 2))) * (rho_l*L_v - pc.internalField()[celli]);
}

//*********************************************************** //
