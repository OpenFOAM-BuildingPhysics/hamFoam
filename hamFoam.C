/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

Application
    hamFoam

Description
    Solves for coupled heat and moisture transport in porous building materials
    Written in 2014, ETH Zurich/Empa

    Contributions:
    Aytac Kubilay, aytac.kubilay@empa.ch, akubilay@ethz.ch
    Andrea Ferrari, andferra@ethz.ch

    Latest Update: 01.08.2019

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "buildingMaterialModel.H"
#include "mixedFvPatchFields.H"
#include "interpolationTable.H"  

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    #include "createMesh.H"
    #include "createFields.H"

    #include "readSolverControls.H"

    #include "updatebuildingMaterialsHeader.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    bool timeStepDecrease = false;

    Info<< "\nStarting time loop\n" << endl;
    TimeState tSCurrent(runTime);

    while (runTime.run())
    {
        Info << nl << "Time = " << runTime.timeName() << ", deltaT = " << runTime.deltaT().value() << endl;       

        #include "updatebuildingMaterials.H"
        
        //store values from previous timestep
        volScalarField ws_old = ws; //ws_old.boundaryField().updateCoeffs();
        volScalarField pc_old = pc; //pc_old.boundaryField().updateCoeffs();
        volScalarField Ts_old = Ts; //Ts_old.boundaryField().updateCoeffs();

        //store values from previous Picard iteration
        volScalarField ws_n = ws; //ws_n.boundaryField().updateCoeffs();
        volScalarField pc_n = pc; //pc_n.boundaryField().updateCoeffs();
        volScalarField Ts_n = Ts; //Ts_n.boundaryField().updateCoeffs();

        for (int nIterPicard=1; nIterPicard<=nIterPicardMax; nIterPicard++)
        {   
            //Moisture transfer////////////
            #include "pcEqn.H"    

            //Firstly, test if all pc values are valid
            if (gMax(pc) >= 0 || gMax(pc.boundaryField()) >= 0)
            {
                Info << "This is going to crash (pc)! Decreasing timestep and reverting fields..." << endl;
                Info << "Error: gMax(pc): " << gMax(pc) << ", gMax(pc.boundaryField()): " << gMax(pc.boundaryField()) << endl;
                timeStepDecrease = true;
                #include "setDeltaT.H"
                #include "revertValues.H"
                break;
            }
            pc.correctBoundaryConditions();
            ///////////////////////////////

            //Heat transfer////////////////
            #include "TsEqn.H" 
            //Firstly, test if all Ts values are valid
            if (gMin(Ts) <= 0 || gMin(Ts.boundaryField()) <= 0)
            {           
                Info << "This is going to crash (Ts)! Decreasing timestep and reverting fields..." << endl;
                timeStepDecrease = true;
                #include "setDeltaT.H"
                #include "revertValues.H"  
                break;                      
            }
            Ts.correctBoundaryConditions();
            ///////////////////////////////

            //Convergence test/////////////      
            #include "updatebuildingMaterials.H" //update values for convergence test
            
            #include "checkConvergence.H" 
            ///////////////////////////////
        }       

        runTime.write();    
        
        if (runTime.writeTime())
        {
        Info<< "Time = " << runTime.value() << " s"
            << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
        }
    
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

