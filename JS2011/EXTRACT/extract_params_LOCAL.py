from numpy import array

# standard path to utility function directory (used to source generic R functions)
utilFolder = '/home/pwg/map_utils/map_utils/'

# main input hdf5  file of simulated realisations of f
realisations_path = '/home/pwg/mbg-world/extraction/realizations/realizations_QRYPFPR010708_Africa_Run_9.10.2008_iterations_FILESTARTREL_FILEENDREL.hdf5'
#filename = '/home/pwg/Desktop/test_sim_KE.hdf5'

# location for export of raw extractions (as they come off each distributed instance)
exportPathDistributed_country = '/home/pwg/mbg-world/extraction/DistributedOutput_country/'
exportPathDistributed_perpixel = '/home/pwg/mbg-world/extraction/DistributedOutput_perpixel/'

# location for export of combined extractions (after distributed files joined by extract_combineDistribExtractions.py)
exportPathCombined_country = '/home/pwg/mbg-world/extraction/CombinedOutput_country/'
exportPathCombined_perpixel = '/home/pwg/mbg-world/extraction/CombinedOutput_perpixel/'

# input 1km salb raster of unique spatial IDs
salblim1km_path="/home/pwg/mbg-world/datafiles/auxiliary_data/salblim1km-e_AF.hdf5"
#salblim1km_path="/home/pwg/mbg-world/datafiles/auxiliary_data/salblim1km-e_ken.hdf5"

# input 1km and 5km raster of population per cell
gr001km_path="/home/pwg/mbg-world/datafiles/auxiliary_data/gr001km_AF.hdf5"
gr005km_path="/home/pwg/mbg-world/datafiles/auxiliary_data/gr005km_AF.hdf5"
#gr001km_path="/home/pwg/mbg-world/datafiles/auxiliary_data/gr001km_ken.hdf5"

# location of 5km hdf5  limits mask
lim5kmbnry_path="/home/pwg/mbg-world/datafiles/auxiliary_data/lim5kmbnry-e_y-x+_AF.hdf5"
#lim5kmbnry_path="/home/pwg/mbg-world/datafiles/auxiliary_data/lim5kmbnry-e_y-x+_ken.hdf5"

# files containing list of unique salb IDs in input raster and pixels per ID : generated as ouptut from FUNexamineSalb
uniqueSalb_path='/home/pwg/mbg-world/extraction/misc/uniqueSalb.txt'
pixelN_path='/home/pwg/mbg-world/extraction/misc/pixelN.txt'

# if string given, this is where the standard output and standard error will be sent, if 0, then none sent
STDOUTPUT = '/home/pwg/mbg-world/extraction/DistributedOutput/scriptoutput/'   

# class definition dictionaries
breaks_MBGW={"BREAKS":[0.,0.05,0.40,1.1],"BREAKNAMES":["lte05","gt05lte40","gt40lte100"],"NAME":"MBGW"}
breaks_MVI={"BREAKS":[0.,0.05,0.30,0.40,0.60,1.1],"BREAKNAMES":["lte05","gt05lte30","gt30lte40","gt40lte60","gt60lte100"],"NAME":"MVI"}
breaks_NOORKENYA={"BREAKS":[0.,0.01,0.05,0.4,1.1],"BREAKNAMES":["lte1","gt1lte0","gt5lte40","gt40lte100"],"NAME":"NOORKENYA"}

#breaks_LYSENKO={"BREAKS":[0.,0.1,0.5,0.75,1.1],"BREAKNAMES":["lte10","gt10lte50","gt50lte75","gt75lte100"],"NAME":"LYSENKO"}
#breaks_one={"BREAKS":[0.,1.1],"BREAKNAMES":["all"],"NAME":"ONE"}
#breaksDict={"MBGW":breaks_MBGW,"MVI":breaks_MVI,"LYSENKO":breaks_LYSENKO,"ONE":breaks_one}
breaksDict={"MBGW":breaks_MBGW}
#breaksDict={}

# ratio of high to low resolution
HiResLowResRatio=5

# parameter determining number of rows of the coarser grid (5km) to process at a time. Default is 1. (currently unsupported)
rowsInslice5km=5 

# what summary stats do we want to generate for meanPR, PAR etc in each country? If none, then set to None
summaryStats={'mean':'mean','SD':'SD','quantiles':array([0,0.025,0.25,0.5,0.75,0.975,1])}
#summaryStats=None



