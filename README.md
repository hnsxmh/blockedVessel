# blockedVessel
simulation of blocked vessel

Tips: workable in OpenFoam-v1906 and have to make some modifications in constant folder (rename file "turbulenceProperties" with momentumTransport) in OpenFoam-dev.

How it works:

1).Compile the tractionDisplacementPressureProfile with "wmake libso" in Terminal.

2). Run the case blockedVessel with script Allrun.

3). Use postprocessing utility with command "postprocess -func sampleDict" to sample the pressure profile on the patches blockFront, blockTop, blockBack.

4). Use the sample result in case "block" to intialize the boundary conditions on block (already done, its only for notation).

5). Run the case "block" with script Allrun.
