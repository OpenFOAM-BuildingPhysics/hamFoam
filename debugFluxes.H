wordReList patchSet =
    runTime.controlDict().lookupOrDefault("debugPatches", wordReList());
    
forAll(patchSet, i)
{   
    const label& patchi = mesh.boundaryMesh().findPatchID(patchSet[i]);
    const fvPatch& patch = mesh.boundary()[patchi];    

    const mixedFvPatchScalarField&
        pcFvPatch = refCast
            <const mixedFvPatchScalarField>
            (
                patch.lookupPatchField<volScalarField, scalar>("pc")
            ); 
            
    scalarField TsPatch(patch.size(), 0.0);
        TsPatch = patch.lookupPatchField<volScalarField, scalar>("Ts");
    scalarField pcPatch(patch.size(), 0.0);
        pcPatch = patch.lookupPatchField<volScalarField, scalar>("pc");     
    scalarField KrelPatch(patch.size(), 0.0);
        KrelPatch = patch.lookupPatchField<volScalarField, scalar>("Krel");
    scalarField K_vPatch(patch.size(), 0.0);
        K_vPatch = patch.lookupPatchField<volScalarField, scalar>("K_v");                 

    dictionary betaDict;
    betaDict.add("file", "$FOAM_CASE/0/" + patch.name() +  "/beta");
    dictionary pv_oDict;
    pv_oDict.add("file", "$FOAM_CASE/0/" + patch.name() +  "/pv_o");
    dictionary glDict;
    glDict.add("file", "$FOAM_CASE/0/" + patch.name() +  "/gl");        

    Function1s::Table<scalar> beta
    (
        "beta",
        betaDict
    ); 
    Function1s::Table<scalar> pv_o
    (
        "pv_o",
        pv_oDict
    ); 
    Function1s::Table<scalar> gl
    (
        "gl",
        glDict
    ); 
    
    scalar currentTime = runTime.value();
    scalar currentTimeStep = runTime.deltaT().value();
    
    scalarField pvsat_s = exp(6.58094e1-7.06627e3/TsPatch-5.976*log(TsPatch));
    scalar rhol=1.0e3; scalar Rv=8.31451*1000/(18.01534);
    scalarField pv_s = pvsat_s*exp((pcPatch)/(rhol*Rv*TsPatch));
    scalarField evap = beta.value(currentTime)*(pv_s-pv_o.value(currentTime));

    const scalarField& V = mesh.V();

    scalar A = 0; scalar B = 0; scalar C = 0;

    if(gMax(pcFvPatch.valueFraction())>0)
    {
        scalar gl_ = gl.value(currentTime);
        scalarField runoff = (gl_ - (KrelPatch+K_vPatch)*pcFvPatch.snGrad()) * pcFvPatch.valueFraction();
        A = ( gSum(runoff*patch.magSf())/gSum(patch.magSf()) )*currentTimeStep;
    }
    else
    {
        C = ( gSum(evap*patch.magSf())/gSum(patch.magSf()) )*currentTimeStep;
    }
    B = gSum(ws*V)/gSum(patch.magSf());

    Info << patch.name() << " | time, runoff, abs, evap: " << runTime.time().value() << " " << A << " " << B << " " << C << endl;
}
