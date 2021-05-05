Occasionally changes to the test input files are required for compatibility with the latest HAZ code. Those changes are described here.

Test 1.10
The code now readjusts hzstep (hystep) to fit evenly into the seismogenic thickness for more stable results with respect to the sampling of depths for areal source zones. This necessitated a change to the seismic source input file to force the creation of a single layer of point sources at depth = 5 km. The seismogenic thickness was changed from 0.00001 km to 0.05 km.