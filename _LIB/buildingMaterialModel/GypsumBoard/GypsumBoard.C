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

#include "GypsumBoard.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace buildingMaterialModels
{
    defineTypeNameAndDebug(GypsumBoard, 0);

    addToRunTimeSelectionTable
    (
        buildingMaterialModel,
		GypsumBoard,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::buildingMaterialModels::GypsumBoard::GypsumBoard
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
void Foam::buildingMaterialModels::GypsumBoard::update_w_C_cell(const volScalarField& pc, volScalarField& w, volScalarField& Crel, label& celli)
{
    List<scalar> reta; reta.setSize(1); reta[0]=-3.4e-6;
    List<scalar> retn; retn.setSize(1); retn[0]=1.62;
    List<scalar> retm; retm.setSize(1); retm[0]=0.62;
    List<scalar> retw; retw.setSize(1); retw[0]=1;
    scalar w_tmp = 0; scalar tmp = 0; scalar C_tmp = 0; scalar tmp2 = 0;    
	for (int i = 0; i <= 0; i++)
	{
		tmp = pow((reta[i] * pc.internalField()[celli]), retn[i]);
		w_tmp = w_tmp + retw[i] / (pow((1 + tmp), retm[i]));
		tmp2 = pow((1 + tmp), retm[i]);
		C_tmp = C_tmp - retw[i] / tmp2 * retm[i] * retn[i] * tmp / ((1 + tmp)*pc.internalField()[celli]);
	}
	w.ref()[celli] = w_tmp * 810;
	Crel.ref()[celli] = mag(C_tmp * 810);
}

//- Correct the buildingMaterial liquid permeability (cell)
void Foam::buildingMaterialModels::GypsumBoard::update_Krel_cell(const volScalarField& pc, const volScalarField& w, volScalarField& Krel, label& celli)
{
	scalar logpc = log10(-pc.internalField()[celli]);
	scalar logKl = 0;
	int i;
	List<scalar> logpc_M = { 3, 5.52026, 5.53615, 5.55196, 5.5677, 5.58341, 5.5991, 5.6148, 5.63052,
						5.64628, 5.6621, 5.67801, 5.69403, 5.71017, 5.72647, 5.74293, 5.7596, 
						5.77648, 5.80239, 5.82874, 5.84679, 5.86522, 5.88407, 5.90337, 5.92318,
						5.94354, 5.97552, 6.02041, 6.06885, 6.12176, 6.18035, 6.24642, 6.32263, 6.41335, 6.52644, 6.67848,
						6.91768, 7.15431, 7.48021, 7.6839, 7.8399, 7.97245, 8.09366, 8.21225 };	
	
	List<scalar> logKl_M = { -8.85398, -8.85398, -8.85398, -8.85398, -8.85398, -8.85398, -8.85398,
						-8.85398, -8.85398, -8.85398, -8.85398, -8.85398, -8.85398, -8.85398,
						-8.85398, -8.85398, -8.85398, -8.85398, -8.85398, -8.85398, -8.85398, -9.58278,
						-9.8215, -9.96441, -10.0695, -10.1569, -10.286, -10.4207, -10.561, -10.7225, -10.9058,
						-11.1138, -11.3534, -11.6288, -11.9632, -12.4104, -13.1461, -15.3938, -15.5518, -15.7333, -15.9512, -16.2285, -16.6545, -20 };

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
void Foam::buildingMaterialModels::GypsumBoard::update_Kvap_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_v, volScalarField& K_pt, label& celli)

{
	scalar rho_l = 1.0e3;
	scalar R_v = 8.31451 * 1000 / (18.01534);
    scalar L_v = 2.5e6;
	scalar a = 2.61e-5;
	scalar b = 5.6;
	scalar c = 0.503;
	scalar d = 0.497;

	scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3 / T.internalField()[celli] - 5.976*Foam::log(T.internalField()[celli])); // saturation vapour pressure [Pa]
	scalar relhum = Foam::exp(pc.internalField()[celli] / (rho_l*R_v*T.internalField()[celli])); // relative humidity [-]

	scalar tmp = 1 - (w.internalField()[celli] / 810);
	scalar delta = a * tmp / (R_v*T.internalField()[celli] * b * (c*tmp*tmp + d)); // Water vapour diffusion coefficient "for GypsumBoard" [s]

	K_v.ref()[celli] = (delta*p_vsat*relhum) / (rho_l*R_v*T.internalField()[celli]);
    K_pt.ref()[celli] = ( (delta*p_vsat*relhum)/(rho_l*R_v*pow(T.internalField()[celli],2)) ) * (rho_l*L_v - pc.internalField()[celli]);

}

//*********************************************************** //
