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
    Foam::buildingMaterialModels::Soil

Description
    buildingMaterial constitutive relationship for soil, consisting of 
    a balanced proportion of sand, silt and clay based on 
    Janssen et al. 2002 and Schaap et al. 1998.
   
    References:
    Janssen H, Carmeliet J, Hens H. 2002. The influence of soil moisture in the un-
    saturated zone on the heat loss from buildings via the ground. 
    Journal of Building Physics, 25, 275–298.
    
    Schaap MG, Leij FJ, Van Genuchten MT. 1998. Neural network analysis for
    hierarchical prediction of soil hydraulic properties. Soil Science Society 
    of America Journal Soil Science Society of America, 62, 847–855.

SourceFiles
    Soil.C

\*---------------------------------------------------------------------------*/

#ifndef Soil_H
#define Soil_H

#include "buildingMaterialModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace buildingMaterialModels
{

/*---------------------------------------------------------------------------*\
                           Class Soil Declaration
\*---------------------------------------------------------------------------*/

class Soil
:
    public buildingMaterialModel
{
    // Private data

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("Soil");


    // Constructors

        //- Construct from components
        Soil
        (
            const word& name, 
            const dictionary& buildingMaterialDict,
            const word& cellZoneModel
        );


    // Destructor

        ~Soil()
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
