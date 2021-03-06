from numpy import array

# define region code, used to bring in grid paramaters from appropriate lines of grid_lims.py
region = "EA"

# what year are we starting at (the january) and how many months are we predicting?
nmonths = 132
start_year = 2000
 
# standard path to utility function directory (used to source generic R functions)
utilFolder = '/mnt/pwg/map_utils/map_utils/'
 
# set path to file containg keys for amazon S3
keyPath = '/mnt/pwg/s3code.txt' 
 
# location and name of global 5km hdf5 stable limits mask
lim5kmbnry_path="/root/mbg-world/datafiles/auxiliary_data/st_mask5km-e_y-x+.hdf5"

# location and name of global 5km hdf5 urban indicator surface
urb5km_path="/root/mbg-world/datafiles/auxiliary_data/urb5km-e_y-x+.hdf5"

# location and name of global 5km hdf5 periurban indicator surface
periurb5km_path="/root/mbg-world/datafiles/auxiliary_data/periurb5km-e_y-x+.hdf5"

# location and name of trace file 
trace_path="/mnt/auxiliary_data/QRYPFPR010708_Africa_Run_9.10.2008.hdf5"

# location of folder to house realization to be generated
realizations_path = '/mnt/qrypfpr010708_africa_run_9.10.2008_try2/'

# how many points will we take from the block to condition the unconditioned field at the outside data locations?
NinThinnedBlock =10000

# set some other memory and kriging parameters
burn = 2000
memmax = 1.e8
thinning = 4
relp=1e-6

# have we merged the urban and periurban categories into a sinlge 'urban' category
merged_urb = False

# are we going to test each month has sensible range of max and min values ?
TESTRANGE = False

