import numpy as np

import os, sys
mfpytoolspth = '../'
if mfpytoolspth not in sys.path:
    sys.path.append(mfpytoolspth)
import mfpytools.binaryfile as bf
reload(bf)

#create the cell by cell flow budget object
cbb = bf.CellBudgetFile('ex2.cbc')

#this will print a list of the records in the cbc file
cbb.list_records()

#this shows the list of data types in the file
print 'types of data: ', cbb.textlist

#this is the list of unique times
print 'times: ', cbb.times

#this will get the data for each record in the file.
for idx in xrange(cbb.nrecords):
    data = cbb.get_data(idx=idx)

#you can also get the data for a specific budget item
#For compact list budgets, this is returned as a structured 
#numpy array.  See this page for more details:
# http://docs.scipy.org/doc/numpy/user/basics.rec.html
data = cbb.get_data(kstp=1, kper=1, text='SWR LEAKAGE')
print data['node']  #node numbers
print data['q']  #volumetric flux for each node

#for three dimensional arrays, this is just a 3d numpy
#array.
data = cbb.get_data(kstp=1, kper=1, text='FLOW RIGHT FACE')
print data

