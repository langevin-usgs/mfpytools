import numpy as np
from collections import OrderedDict


def binaryread_struct(file, vartype, shape=(1), charlen=16):
    '''Read text, a scalar value, or an array of values from a binary file.
       file is an open file object
       vartype is the return variable type: str, numpy.int32, numpy.float32, 
           or numpy.float64
       shape is the shape of the returned array (shape(1) returns a single value)
           for example, shape = (nlay, nrow, ncol)
       charlen is the length of the text string.  Note that string arrays cannot
           be returned, only multi-character strings.  Shape has no affect on strings.
    '''
    import struct
    import numpy as np
    
    #store the mapping from type to struct format (fmt)
    typefmtd = {np.int32:'i', np.float32:'f', np.float64:'d'}
        
    #read a string variable of length charlen
    if vartype is str:
        result = file.read(charlen*1)
        
    #read other variable types
    else:
        fmt = typefmtd[vartype]
        #find the number of bytes for one value
        numbytes = vartype(1).nbytes
        #find the number of values
        nval = np.core.fromnumeric.prod(shape)
        fmt = str(nval) + fmt
        s = file.read(numbytes * nval)
        result = struct.unpack(fmt, s)
        if nval == 1:
            result = vartype(result[0])
        else:
            result = np.array(result, dtype=vartype)
            result = np.reshape(result, shape)
    return result
    
def binaryread(file, vartype, shape=(1), charlen=16):
    '''uses numpy to read from binary file.  This
       was found to be faster than the struct
       approach and is used as the default.
    '''
    
    #read a string variable of length charlen
    if vartype is str:
        result = file.read(charlen*1)     
    else:
        #find the number of values
        nval = np.core.fromnumeric.prod(shape)
        result = np.fromfile(file,vartype,nval)
        if nval == 1:
            result = result[0]
        else:
            result = np.reshape(result, shape)
    return result


class HeadFile(object):
    '''
    The HeadFile class provides simple ways to retrieve 2d and 3d 
    head arrays from a MODFLOW binary head file and time series
    arrays for one or more cells.
    
    A HeadFile object is created as
    hdobj = HeadFile(filename, precision='single')
    
    The HeadFile class is built on an ordered dictionary consisting of 
    keys, which are tuples of the modflow header information
    (kstp, kper, pertim, totim, text, nrow, ncol, ilay)
    and long integers, which are pointers to first bytes of data for
    the corresponding data array.
    '''
    def __init__(self, filename, precision='single', verbose=False):
        self.filename = filename
        self.precision = precision
        self.verbose = verbose
        self.file = open(self.filename, 'rb')
        self.nrow = 0
        self.ncol = 0
        self.nlay = 0
        self.times = []
        self.kstpkper = []

        if precision is 'single':
            self.realtype = np.float32
        elif precision is 'double':
            self.realtype = np.float64
        else:
            raise Exception('Unknown precision specified: ' + precision)

        self.header_dtype = np.dtype([('kstp','i4'),('kper','i4'),('pertim','f4'),\
                                     ('totim','f4'),('text','a16'),\
                                     ('ncol','i4'),('nrow','i4'),('ilay','i4')])


        #read through the file and build the pointer index
        self._build_index()
        
        #allocate the head array
        self.head = np.empty( (self.nlay, self.nrow, self.ncol), 
                         dtype=self.realtype)

        return

    def _build_index(self):
        '''
        Build the ordered dictionary, which maps the MODFLOW header information
        to the position in the binary file.
        '''
        kstp, kper, pertim, totim, text, nrow, ncol, ilay = self.get_header()
        self.nrow = nrow
        self.ncol = ncol
        self.file.seek(0, 2)
        self.totalbytes = self.file.tell()
        self.file.seek(0, 0)        
        self.databytes = ncol * nrow * self.realtype(1).nbytes
        self.recorddict = OrderedDict()
        ipos = 0
        while ipos < self.totalbytes:
            kstp, kper, pertim, totim, text, nrow, ncol, ilay = self.get_header()
            self.nlay=max(self.nlay, ilay)
            if totim not in self.times:
                self.times.append(totim)
            self.kstpkper.append( (kstp, kper) )
            key = (kstp, kper, pertim, totim, text, nrow, ncol, ilay)
            ipos = self.file.tell()
            self.recorddict[key] = ipos
            self.file.seek(self.databytes, 1)
            ipos = self.file.tell()
        return

    def get_header(self):
        '''
        Read the MODFLOW header
        '''        
        header = binaryread(self.file,self.header_dtype,(1,))
        return header



    def list_records(self):
        '''
        Print a list of all of the records in the file
        hdobj.list_records()
        '''
        for key in self.recorddict.keys():
            print key
        return

    def _fill_head_array(self, kstp=0, kper=0, totim=-1, text='HEAD'):
        '''
        Fill the three dimensional head array, self.head, for the
        specified kstp and kper value or totim value.
        '''
        recordlist = []
        for key in self.recorddict.keys():
            if text.upper() not in key[4]: continue
            if kstp > 0 and kper > 0:
                if key[0] == kstp and key[1] == kper:
                    recordlist.append(key)
            elif totim >= 0.:
                if totim == key[3]:
                    recordlist.append(key)
            else:
                raise Exception('Data not found...')

        #initialize head with nan and then fill it
        self.head[:, :, :] = np.nan
        for key in recordlist:
            ipos = self.recorddict[key]
            self.file.seek(ipos, 0)
            ilay = key[7]
            self.head[ilay - 1, :, :] = binaryread(self.file, self.realtype, 
                shape=(self.nrow, self.ncol))
        return

    def get_data(self, kstp=0, kper=0, idx=0, totim=-1, ilay=0, text='HEAD'):
        '''
        Return a three dimensional head array for the specified kstp, kper
        pair or totim value, or return a two dimensional head array
        if the ilay argument is specified.
        '''
        self._fill_head_array(kstp, kper, totim, text)
        if ilay == 0:
            return self.head
        else:
            return self.head[ilay-1, :, :]
        return

    def get_ts(self, k=0, i=0, j=0, text='HEAD'):
        '''
        Create and return a time series array of size [ntimes, nstations].
        
        The get_ts method can be called as ts = hdobj.get_ts(1, 2, 3),
        which will get the time series for layer 1, row 2, and column 3.
        The time value will be in the first column of the array.
        
        Alternatively, the get_ts method can be called with a list like
        ts = hdobj.get_ts( [(1, 2, 3), (2, 3, 4)] ), which will return
        the time series for two different cells.
        
        '''
        if isinstance(k, list):
            kijlist = k
            nstation = len(kijlist)
        else:
            kijlist = [ (k, i, j) ]
            nstation = 1
        result = np.empty( (len(self.times), nstation + 1), dtype=self.realtype)
        result[:, :] = np.nan
        result[:, 0] = np.array(self.times)

        istat = 1
        for k, i, j in kijlist:
            recordlist = []
            ioffset = ((i - 1) * self.ncol + j - 1) * self.realtype(1).nbytes
            for key in self.recorddict.keys():
                if text not in key[4]: continue
                ilay = key[7]
                if k == ilay:
                    recordlist.append(key)
            for key in recordlist:
                ipos = self.recorddict[key]
                self.file.seek(ipos + np.long(ioffset), 0)
                itim = np.where(result[:, 0] == key[3])[0]
                result[itim, istat] = binaryread(self.file, np.float32)
            istat += 1
        return result
        
