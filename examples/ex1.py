import numpy as np

import os, sys
mfpytoolspth = '../'
if mfpytoolspth not in sys.path:
    sys.path.append(mfpytoolspth)
import mfpytools.binaryfile as bf
reload(bf)

#create the cell by cell flow budget object
hds = bf.HeadFile('ex1.hds')

#this will print a list of the records in the heads file
hds.list_records()

#you can get heads as a 3d array like
array3d = hds.get_data(kstp=1, kper=1)
print array3d.shape

#or as a 2d array by including the layer
array2d = hds.get_data(kstp=1, kper=1, ilay=1)
print array2d.shape

#print the times that data are available
print hds.times

#get a time series for a cell
ts = hds.get_ts(10, 1, 10)
print ts

#or get times series for multiple cells (layer, row, column)
ts = hds.get_ts([(1, 1, 1), (10, 1, 1)])
print ts

#these exact same utilities work for mt3dms concentrations also
#create the cell by cell flow budget object
ucn = bf.UcnFile('ex1.ucn')

#this will print a list of the records in the heads file
ucn.list_records()

#you can get heads as a 3d array like
array3d = ucn.get_data(kstp=1, kper=1)
print array3d.shape

#or as a 2d array by including the layer
array2d = ucn.get_data(kstp=1, kper=1, ilay=1)
print array2d.shape

#print the times that data are available
print ucn.times

#get a time series for a cell
ts = ucn.get_ts(10, 1, 10)
print ts

#or get times series for multiple cells (layer, row, column)
ts = ucn.get_ts([(1, 1, 1), (10, 1, 1)])
print ts
