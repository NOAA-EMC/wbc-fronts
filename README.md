# wbc-fronts


# Technical Details
The Gulf Stream and Loop Current location is defined as the intersection of the 12 degree celsius isotherm and 400 meter isobath.
The Kuroshio Current locaiton is defined as the intersection of the 14 degree celsius isotherm and the 200 meter isobath.

The Hausdorff distance between two sets of points (in this case the Gulf Stream North Wall Navy analyses and the North Wall from the Global model such as RTOFS) is a measure of how far apart two polygon shapes or set of points are and is computed as the greatest of all the distances from the points in one set to the closest point in the other set. The modified Hausdorff distance is the minimum of the Hausdorff distances computed by swapping the point sets. 

#How to Run On Orion
Step 0: aquire data
MOM6
/work/noaa/marine/erob/soca-diag-eric/2020/2020072412/ctrl/ocn.bkg.2020072412.nc
Navy
/work/noaa/marine/erob/soca-diag-eric/2020/20200724/wtxtbul/gs206xx.mrf
/work/noaa/marine/erob/soca-diag-eric/2020/20200724/wtxtbul/gs207nw.sub
/work/noaa/marine/erob/soca-diag-eric/2020/20200724/wtxtbul/ms206xx.mrf
/work/noaa/marine/erob/soca-diag-eric/2020/20200724/wtxtbul/np206xx.mrf
/work/noaa/marine/erob/soca-diag-eric/2020/20200724/wtxtbul/wi206xx.mrf

RTOFS
/work/noaa/marine/erob/soca-diag-eric/2020_used/20200724/wtxtbul/rtofs_glo_2ds_n024_daily_diag.nc 
/work/noaa/marine/erob/soca-diag-eric/2020_used/20200724/wtxtbul/ rtofs_glo_3dz_n024_daily_3zsio.nc
/work/noaa/marine/erob/soca-diag-eric/2020_used/20200724/wtxtbul/ rtofs_glo_3dz_n024_daily_3zuio.nc
/work/noaa/marine/erob/soca-diag-eric/2020_used/20200724/wtxtbul/rtofs_glo_2ds_n024_daily_prog.nc
/work/noaa/marine/erob/soca-diag-eric/2020_used/20200724/wtxtbul/rtofs_glo_3dz_n024_daily_3ztio.nc
/work/noaa/marine/erob/soca-diag-eric/2020_used/20200724/wtxtbul/rtofs_glo_3dz_n024_daily_3zvio.nc


Step 1: Load your own Python or use the following:
  
     module use -a /home/cmartin/opt/modulefiles
     module load anaconda/anaconda3-2020.04.02

Step 2 Make the following changes to global_fronts_hat10.py: 


Read_navy.py
Line 8 change to location of DCOm files

Changes for global_fronts.py
Change line 42 imageDir='/work/noaa/marine/erob/soca-diag-eric/' to your desired image output directory 

Change line 43 dbFile='/home/erobin/fronts/fix/global_fronts.db'  Create a copy of this .db file and insert correct file path (or use this one)

Line 385 change to desired output file path

Line 466 change to desired output file path

Line 510  modelDir='/work/noaa/marine/erob/soca-diag-eric/2020_used/'+theDate.strftime('%Y%m%d')+'/wtxtbul'

Change to location of RTOFS data



Step 3
python global_fronts_hat10.py 20200724