import unittest
import numpy as np
import struct

import os, sys
mfpytoolspth = '../'
if mfpytoolspth not in sys.path:
    sys.path.append(mfpytoolspth)
import mfpytools.binaryfile as bf

class TestBinaryRead(unittest.TestCase):
    
    def setUp(self):

        self.ncol_write = np.int32(100)
        self.nrow_write = np.int32(200)
        self.arr_write = np.random.rand(self.nrow_write, 
            self.ncol_write).astype(np.float32)
        aname = 'array.bin'
        self.aname = aname
        f = open(aname,'wb')
        self.kstp_write = np.int32(99)
        self.kper_write = np.int32(100)
        self.pertim_write = np.float32(1.999)
        self.totim_write = np.float32(75.6)
        self.ilay_write = np.int32(1)
        self.txt_write = 'test string     '
        f.write(struct.pack('i', self.kstp_write))
        f.write(struct.pack('i', self.kper_write))
        f.write(struct.pack('f', self.pertim_write))
        f.write(struct.pack('f', self.totim_write))
        f.write(struct.pack('i', self.ncol_write))
        f.write(struct.pack('i', self.nrow_write))
        f.write(struct.pack('i', self.ilay_write))
        f.write(self.txt_write)
        self.arr_write.tofile(f)
        f.close()
        
        f = open(aname,'rb')
        self.kstp = bf.binaryread(f, np.int32)
        self.kper = bf.binaryread(f, np.int32)
        self.pertim = bf.binaryread(f, np.float32)
        self.totim = bf.binaryread(f, np.float32)
        self.ncol = bf.binaryread(f, np.int32)
        self.nrow = bf.binaryread(f, np.int32)
        self.ilay = bf.binaryread(f, np.int32)
        self.txt = bf.binaryread(f, str)
        self.arr = bf.binaryread(f, np.float32, shape=(self.nrow, self.ncol))
        f.close()
        return
        
    def test_read_int32(self):
        self.assertEqual(self.kstp, self.kstp_write)
        self.assertEqual(self.kper, self.kper_write)
        self.assertEqual(self.ncol, self.ncol_write)
        self.assertEqual(self.nrow, self.nrow_write)
        self.assertEqual(self.ilay, self.ilay_write)
        return

    def test_read_float32(self):
        self.assertEqual(self.pertim, self.pertim_write)
        self.assertEqual(self.totim, self.totim_write)
        return

    def test_read_text(self):
        self.assertEqual(self.txt, self.txt_write)

    def test_read_float32_array(self):
        for i in range(self.nrow):
            for j in range(self.ncol):
                self.assertEqual(self.arr[i, j], self.arr_write[i, j])
        return
        
    def tearDown(self):
        os.remove(self.aname)
        return
        
if __name__ == '__main__':
    unittest.main()

