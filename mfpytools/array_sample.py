import numpy as np
from scipy.ndimage.interpolation import map_coordinates
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from math import cos, sin

def modflow_global_coords(nrow, ncol, delr, delc, xoffset, yoffset, theta):
    '''
    Return xg and yg, the global x and y coordinates for a modflow grid.
    '''
    Ly = np.add.reduce(delc)
    Lx = np.add.reduce(delr)
    xlocal = np.add.accumulate(delr) - 0.5 * delr
    ylocal = Ly - ( np.add.accumulate(delc) - 0.5 * delc )
    X, Y = np.meshgrid(xlocal, ylocal)
    X, Y = np.meshgrid(xlocal, ylocal)
    xg, yg = rotate(X, Y, theta, xorigin=0., yorigin=Ly)
    xg += xoffset
    yg += yoffset
    return xg, yg

def rotate(x, y, theta, xorigin=0., yorigin=0.):

    '''
    Given x and y array-like values calculate the rotation about and arbitrary
    origin and then return the rotated coordinates'''
    xrot = xorigin + cos(theta) * (x - xorigin) - sin(theta) * (y - yorigin)
    yrot = yorigin + sin(theta) * (x - xorigin) + cos(theta) * (y - yorigin)
    return xrot, yrot

def array_sample(x,  y, arr, arr_extent, order=0, cval=-999., 
                      mode='nearest'):
    '''
    Function to interpolate a value from an array using x and y coordinates.
    The array can be mapped in space using an array extent.  For direct
    sampling without any interpolation, order should be set to zero.  If order
    is set to 1 then bilinear interpolation will be used.  This function
    should be relatively fast because it uses the ndimage.interpolation
    routine called map_coordinates.
    
    Note that arr[0, 0] value is assigned at xmin, ymax, so if using this with
    a finite difference grid, xmin should be set to x0 + 1/2 * dx
    '''
    nrow, ncol = arr.shape
    (xmin, xmax, ymin, ymax) = arr_extent
    dx = float((xmax - xmin) / (ncol - 1) )
    dy = float((ymax - ymin) / (nrow - 1) )
    icoord = (ymax - np.array([y])) / (ymax - ymin) * (nrow - 1)
    jcoord = (np.array([x]) - xmin) / (xmax - xmin) * (ncol - 1)
    coords = [icoord, jcoord]
    z = map_coordinates(arr, coords, order=order, cval=cval, mode=mode)
    return z
    
def modflow_array_sample(xg, yg, arr, xoffset=0., yoffset=0., delr=1., 
                             delc=1., rotation=0., order=0, cval=-999., 
                             mode='nearest'):
    '''
    This function is used to interpolate a value from an array given global
    coordinates, xg and yg.  xg and yg can be arrays.  The function works
    by rotating and then offsetting the global coordinates into local
    coordinates and then interpolating values from the specified array.
    '''
    nrow, ncol = arr.shape

    #calculate grid dimensions
    Ly = np.add.reduce(delc)
    Lx = np.add.reduce(delr)

    #rotate and offset the points
    xlocal, ylocal = rotate(xg, yg, -rotation, xorigin=xoffset, yorigin=yoffset)
    xlocal -= xoffset
    ylocal -= yoffset

    #get interpolated values
    arr_extent = (0 + 0.5*delr[0] , Lx - 0.5 * delr[-1], 
                  0 + 0.5*delc[-1], Ly - 0.5 * delc[0])
    z = array_sample(xlocal, ylocal, arr, arr_extent, order=order, 
                      cval=cval, mode=mode)
    return z
    
if __name__ == '__main__':

    #create an array and map it in space
    a = np.arange(9.).reshape(3,3)
    xmin = 0.5
    xmax = 2.5
    ymin = 0.5
    ymax = 2.5
    extent = (xmin, xmax, ymin, ymax)
    
    #create x and y points that will be used to sample the array
    x = np.arange(0, 3.25, 0.25) #[0., ]
    y = np.arange(0, 3.25, 0.25)
    X, Y = np.meshgrid(x, y)
    
    #sample the array and print results
    z = array_sample(X, Y, a, extent, order=1, mode='nearest')

    #plot the array and the interpolated values
    try:
        plt.close('all')
    except:
        pass
    plt.figure(1)
    plt.subplot(1, 1, 1, aspect='equal')
    #plt.imshow(a, interpolation='nearest', extent=extent)
    plt.pcolor(np.arange(0, 4, 1), np.arange(3, -1, -1), a)
    plt.plot(X, Y, 'bo')
    for i in xrange(X.shape[0]):
        for j in xrange(X.shape[1]):
            plt.text(X[i, j], Y[i, j], str(z[0, i, j]))
    
    #ticks
    majorLocator = MultipleLocator(1)
    plt.gca().xaxis.set_major_locator(majorLocator)
    plt.gca().yaxis.set_major_locator(majorLocator)
    plt.grid(b=True, which='major', color='b', linestyle='-')
    plt.xlim(0, 3)
    plt.ylim(0, 3)
    plt.show()


    #now test the modflow_array_interpolate function
    plt.figure(2)
    plt.subplot(1, 1, 1, aspect='equal')
    xoffset = 5.
    yoffset = 10.
    theta = 15. * 3.1415 / 180.
    nrow = 3
    ncol = 3
    delr = np.array([1, 2, 1])
    delc = np.array([1, 2, 3])

    #create an array that we want to sample and plot it with pcolor
    arr = np.arange(nrow*ncol).reshape( (nrow, ncol) )
    Xelocal = np.concatenate(([0.], np.add.accumulate(delr) ))
    Ly = np.add.reduce(delc)
    Yelocal = np.concatenate(([Ly], Ly - np.add.accumulate(delc) ) )
    Xe, Ye = np.meshgrid(Xelocal, Yelocal)
    Xe, Ye = rotate(Xe, Ye, theta, xorigin=0., yorigin=Ly)
    Xe += xoffset
    Ye += yoffset
    plt.pcolor(Xe, Ye, arr)
    plt.plot(xoffset, yoffset, 'rx')
    
    #find the global coordinates for the cell centers and plot them
    xlocal = np.add.accumulate(delr) - 0.5 * delr
    ylocal = Ly - ( np.add.accumulate(delc) - 0.5 * delc )
    X, Y = np.meshgrid(xlocal, ylocal)
    #plt.plot(X, Y, 'bo')
    X, Y = np.meshgrid(xlocal, ylocal)
    xg, yg = rotate(X, Y, theta, xorigin=0., yorigin=Ly)
    xg += xoffset
    yg += yoffset
    plt.plot(xg, yg, 'go')

    #sample the array using the global coordinates    
    z = modflow_array_sample(xg, yg, arr, xoffset=xoffset, yoffset=yoffset, 
                             delr=delr, delc=delc, rotation=theta, order=0, 
                             cval=-999., mode='nearest')

    #plot the results
    for i in xrange(xg.shape[0]):
        for j in xrange(yg.shape[1]):
            plt.text(xg[i, j], yg[i, j], str(z[0, i, j]))
    
    plt.show()
    