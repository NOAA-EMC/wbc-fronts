# wbc-fronts

# Introduction
This repository is forked from NOAA EMC Ocean Verification branch frontal analysis written by Todd Spindler.
The purpose of this repository is to show comparison of observational analyses of three Western Boundary Currents:
the Gulf Stream, Loop CUrrent and Kurushio Current obtained from Naval Oceanographic Office for the Gulf Streamand the Naval Eastern Ocean Center for the Loop CUrrent and Kuroshio Current with the appropriate frontal location computed from the Global RTOFS model, overlaid on maps of the Global RTOFS Sea Surface Temperature, Sea Surface Height, and Surface Currents. Later on there will be new development to adapt this repository to regional MOM6 for HAFS basin HAT10 domain.

# Technical Details
The Gulf Stream and Loop Current location is defined as the intersection of the 12 degree celsius isotherm and 400 meter isobath.
The Kuroshio Current locaiton is defined as the intersection of the 14 degree celsius isotherm and the 200 meter isobath.

The Hausdorff distance between two sets of points (in this case the Gulf Stream North Wall Navy analyses and the North Wall from the Global model such as RTOFS) is a measure of how far apart two polygon shapes or set of points are and is computed as the greatest of all the distances from the points in one set to the closest point in the other set. The modified Hausdorff distance is the minimum of the Hausdorff distances computed by swapping the point sets. 

#How to Add Mom6 Capabilities
Copy the original RTOFS read_model and delete any lines relating to opening seperate .nc files for each variable(ssh, u, v ...). Open up single MOM6 .nc file and read in each variable, discarding the singleton dimensions and subsecting the appropriate layers. Next properly mask the data.