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

Class
    Foam::buildingMaterialModels::GlassFiberInsulation

Description
	buildingMaterial constitutive relationships as given in
	Heat-Moisture Transport Simulations according to Cameliet
	and Derome's research.

    References:

	Carmeliet, J., & Derome, D. (2011). Temperature Driven Inward Vapor Diffusion Under Constant
	and Cyclic Loading in Small-Scale Wall Assemblies: Part 2 Heat-Moisture Transport Simulations.
	Building and Environment. Vol. 47, pp. 161-169

SourceFiles
    GlassFiberInsulation.C

\*---------------------------------------------------------------------------*/

#ifndef GlassFiberInsulation_H
#define GlassFiberInsulation_H

#include "buildingMaterialModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace buildingMaterialModels
{

/*---------------------------------------------------------------------------*\
                           Class GlassFiberInsulation Declaration
\*---------------------------------------------------------------------------*/

class GlassFiberInsulation
:
    public buildingMaterialModel
{
    // Private data

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("GlassFiberInsulation");


    // Constructors

        //- Construct from components
		GlassFiberInsulation
        (
            const word& name, 
            const dictionary& buildingMaterialDict,
            const word& cellZoneModel
        );


    // Destructor

        ~GlassFiberInsulation()
        {}


    // Member Functions

        //- Correct the buildingMaterial moisture content (cell)
        void update_w_C_cell(const volScalarField& pc, volScalarField& w, volScalarField& Crel, label& celli);

        //- Correct the buildingMaterial liquid permeability (cell)
        void update_Krel_cell(const volScalarField& pc, const volScalarField& w, volScalarField& Krel, label& celli);

        //- Correct the buildingMaterial vapor permeability (cell)
        void update_Kvap_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_v, volScalarField& K_pt, label& celli);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace buildingMaterialModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
