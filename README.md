# wbc-fronts


# Technical Details
The Gulf Stream and Loop Current location is defined as the intersection of the 12 degree celsius isotherm and 400 meter isobath.
The Kuroshio Current locaiton is defined as the intersection of the 14 degree celsius isotherm and the 200 meter isobath.

The Hausdorff distance between two sets of points (in this case the Gulf Stream North Wall Navy analyses and the North Wall from the Global model such as RTOFS) is a measure of how far apart two polygon shapes or set of points are and is computed as the greatest of all the distances from the points in one set to the closest point in the other set. The modified Hausdorff distance is the minimum of the Hausdorff distances computed by swapping the point sets. 

#How to Add Mom6 Capabilities on Orion
Step 1: Load your own Python or use the following:
             module use -a /home/cmartin/opt/modulefiles
     module load anaconda/anaconda3-2020.04.02

Step 2: 

Change line 42 imageDir='/work/noaa/marine/erob/soca-diag-eric/' to your desired image output directory 

Change line 44 dbFile='/home/erobin/fronts/fix/global_fronts.db'  Create a copy of this .db file and insert correct file path (or use this one)

Line 104- same change as dbfile line above

Line 241 change to file path to MOM6 data

Line 463+464 

Line 477 change to desired output file path

Line 492- change dbfile as above

Lines 558-560
Set desired output path for Hasudorff statistics

Line 607- set to Mom6 data file path

Step 3
python global_fronts 20200724