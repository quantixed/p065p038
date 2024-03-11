# 3D Cluster analysis - LBR + MAPPER

A series of scripts is required to do this analysis.

`00_3dSeedFinderAll.ijm`

Finds all seeds for 3D object finding
`01_NucleusFinder.ijm`

LBR around the nucleus needs to be excluded from the analysis`02_3dClusterFinderAll.ijm`Find 3D objects. Using the results from 00 and 01`03_whatDoWeHave.ijm`Visualise the results of the 3D segmentation. Any cells in the corners of the image need to be removed from the analysis. Use the result of this script with 04 to remove objects outside the cell of interest.`04_quickClearII.ijm`Semi-automated removal of extraneous objects.`05_3dClusterFinderAllAlt.ijm`

Rerun of 02 following cleanup using 04.`06_simpleRetrieveClusterStats.ijm`Retrieve the stats on 3D objects detected in 05.`07_simpleRetrieveClusterStatsOverlap.ijm`

Examine colocalization.