scalar maxChangews = gMax(mag(ws.primitiveField()-ws_n.primitiveField())); 
scalar maxChangeTs = gMax(mag(Ts.primitiveField()-Ts_n.primitiveField()));         

if(maxChangews < PicardTolerancews && maxChangeTs < PicardToleranceTs && nIterPicard>=2) //force at least 1 picard iteration
{
    pc_n = pc; //pc_n.boundaryFieldRef().updateCoeffs();
    ws_n = ws; //ws_n.boundaryFieldRef().updateCoeffs();
    Ts_n = Ts; //Ts_n.boundaryFieldRef().updateCoeffs();
    Info << "Total Picard iterations: " << nIterPicard << endl;
    if (debugFluxes)
    {
        #include "debugFluxes.H"
    }
    timeStepDecrease = false;
    runTime++;
    #include "setDeltaT.H" 
    break;
}
else if (nIterPicard == nIterPicardMax) //Picard iteration reached maximum
{
    Info
       << "Picard iteration didn't converge !"
       << endl
       << "maxChangews: " << maxChangews << endl
       << "maxChangeT: " << maxChangeTs << endl
       << "Total Picard iterations: " << nIterPicard << endl;
    timeStepDecrease = true;
    #include "setDeltaT.H"                    
    #include "revertValues.H"
    break;
}   
else //not converged nor reached the maximum iteration yet, continue
{
    pc_n = pc; //pc_n.boundaryFieldRef().updateCoeffs();
    ws_n = ws; //ws_n.boundaryFieldRef().updateCoeffs();
    Ts_n = Ts; //Ts_n.boundaryFieldRef().updateCoeffs();
    Info << "maxChangews: " << maxChangews << ", maxChangeTs: " << maxChangeTs << endl;
}

