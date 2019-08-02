/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
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

Global
    setDeltaT

Description
    Reset the timestep based on Picard iteration.
    Reduction of time-step is immediate, but increase is damped to avoid
    unstable oscillations.

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{
    if(timeStepDecrease == false && runTime.deltaT().value() < maxDeltaT) //time step is increased (slowly)
    {
      Info << "Increasing time step slowly..." 
           <<  min
               (
                 runTime.deltaT().value()*increase_factor,
                 maxDeltaT
               )
           << endl;

      runTime.setDeltaT
      (
          min
          (
              runTime.deltaT().value()*increase_factor,
              maxDeltaT
          )
      );
    }
    else if(timeStepDecrease == true) //time step is decreased
    {
      if (runTime.deltaT().value() - minDeltaT < VSMALL)
      {
        pc.write();
        Ts.write();
        ws.write();
        Crel.write();
    
        FatalError
          << "Cannot decrease time step further. Dumping fields and exiting..." 
          << abort(FatalError);
      }
      Info << "Decreasing time step..." 
           << max
              (
                  runTime.deltaT().value()*decrease_factor,
                  minDeltaT
              )
           << endl;
         

      runTime.setDeltaT
      (
          max
          (
              runTime.deltaT().value()*decrease_factor,
              minDeltaT
          )
      );
    }
    
    Info<< "deltaT = " <<  runTime.deltaT().value() << endl;
}
// ************************************************************************* //
