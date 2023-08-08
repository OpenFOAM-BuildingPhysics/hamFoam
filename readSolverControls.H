word pcEqnForm =
    runTime.controlDict().lookupOrDefault<word>("pcEqnForm", "pc-based"); 

Switch adjustTimeStep =
    runTime.controlDict().lookupOrDefault<Switch>("adjustTimeStep",true);

scalar minCrel =
    runTime.controlDict().lookupOrDefault<scalar>("minCrel", VSMALL); 
    
bool debugFluxes =
    runTime.controlDict().lookupOrDefault<bool>("debugFluxes", false); 

//Picard iteration controls
const dictionary& picardDict = runTime.controlDict().subDict("Picard");

const int nIterPicardMax =
    picardDict.lookupOrDefault<int>("nIterPicardMax", 10);

scalar minDeltaT =
    picardDict.lookupOrDefault<scalar>("minDeltaT", 0.1);

scalar maxDeltaT =
    picardDict.lookupOrDefault<scalar>("maxDeltaT", 600.0);

scalar PicardTolerancews =
    picardDict.lookupOrDefault<scalar>("PicardTolerancews", 0.01);

scalar PicardToleranceTs =
    picardDict.lookupOrDefault<scalar>("PicardToleranceTs", 0.01); 

scalar increase_factor =
    picardDict.lookupOrDefault<scalar>("increase_factor", 1.2);

scalar decrease_factor =
    picardDict.lookupOrDefault<scalar>("decrease_factor", 0.5); 
