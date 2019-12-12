import os, re, datetime as dt 

from PIL import Image                                                       
import numpy as np
from netCDF4 import Dataset

from satlib_py.utils import solar

# Constants for xy to latlon calculations 
req = 6378137 #m
rpol = 6356752.31414 #m
H = 42164160 #m
lambda0 = -1.308996939 #rad
e = 0.0818191910435

def ignore(*args, **kwargs):
    pass

def xy2ab(x,y, log=ignore):
    """
    alpha and betta were defined on MSG. Keeping them for consistency of the 
    post-processing of the library
    """
    return -x,y

def ab2xy(alpha, beta, log=ignore):
    """
    alpha and betta were defined on MSG. Keeping them for consistency of the 
    post-processing of the library
    """
    return -alpha, beta

def xy2lonlat(xy, log=ignore):
    """
    Transforms sat view angles pairs (standarized coordinates x,y) to lats and lons
    Arguments
    ---------
    xy : np.array() - shape(n,2)
        paris of view angles from the satellitei in millirads
    Returns
    -------
    lonlat : np.array() - shape(n,2)
        pairs of latitudes and longitues in degrees
    """
    x = xy.T[0]
    y = xy.T[1]
    a = np.sin(x)**2 + np.cos(x)**2 * (np.cos(y)**2 + req**2 /rpol**2 * np.sin(y)**2 )
    b = -2 * H * np.cos(x) * np.cos(y)
    c = H **2 - req**2
    rs = (-b - np.sqrt(b**2 - 4 * a * c)) / (2 * a)
    sx = rs * np.cos(x) * np.cos(y)
    sy = -rs * np.sin(x)
    sz = rs * np.cos(x) * np.sin(y)
    lat, lon = ( np.arctan( req**2 / rpol**2 * sz / (np.sqrt((H - sx)**2 + sy**2)) ),
           lambda0 - np.arctan( sy / (H - sx) ) )
    log('lats and lons have been generated from x and y')
    lonlat = np.array([lon*180/np.pi, lat*180/np.pi]).T
    return lonlat

def lonlat2xy(lonlat, log=ignore):
    """
    Transforms lats and lons to sat view angles (standarized coordinates x,y) 
    Arguments
    ---------
    lonlat : np.array() - shape(n,2)
        pairs of lat and longitues in degrees
    Returns
    -------
    xy : np.array() - shape(n,2)
        pairs of view angles from the satellite in millirads
    """
    lon = lonlat.T[0]
    lat = lonlat.T[1]
    lat *= np.pi / 180
    lon *= np.pi / 180
    phic = np.arctan(rpol**2 / req**2 * np.tan(lat)) # geocentric latitued
    rc = rpol / np.sqrt(1 - e**2 * np.cos(phic)**2) # geocentric distance
    sx = H - rc * np.cos(phic) * np.cos(lon - lambda0)
    sy = -rc * np.cos(phic) * np.sin(lon - lambda0)
    sz = rc * np.sin(phic)
    x, y = (np.arcsin(-sy / np.sqrt(sx**2 + sy**2 + sz**2)),
            np.arctan(sz / sx))
    log('x and y (sat view angles) have been generated from lats and lons')
    not_in_goesr_view = H * (H - sx) < sy**2 + req**2 / rpol**2 * sz**2
    x[not_in_goesr_view] = np.nan
    y[not_in_goesr_view] = np.nan
    xy = np.array([x, y]).T
    return xy
    
def conus_extraction_fulldisck (x, y, var, log=ignore):
    """
    Extract the section conus using the viewing anlges (standarased coordinates) x and y 
    Arguments
    ---------
    x : np.array() 
        x axis "lon viewing anlges" for the whole disk in millirads
    y : np.array() 
        y axis "lat viewing anlges" for the whole disk in millirads
    var : np.array(2d) - shape=(y.shape, x.shape)
        variable to be extracted for the conus section
    Returns
    -------
    xconus : np.array()
        x axis "lon viewing anlges" of the conus section
    yconus : np.array()
        y axis "lat viewing anlges" of the conus section
    varconus : np.array() - shape=(yconus.shape, xconus.shape)
        extracted section for the conus section of the variable
    """
    #-- valeus for conus section from PUG-L1b-vol3.pdf (GOES16 data description)
    idx = np.where((x > -0.101360) & (x < 0.038640))[0]
    idy = np.where((y < 0.128240) & (y > 0.044240))[0]
    xconus = x[idx]
    yconus = y[idy]
    varconus = var[idy[0] : idy[-1] +1, idx[0] : idx[-1]+1]
    log('The conus section has been extracted from the FullDisk')
    return xconus, yconus, varconus

def gen_latlon_grid(x, y):
    '''
    get lat and lon grids (2d) for an image havving viewing anlge axis x and y (1d)

    Arguments
    ---------
    x, y : numpy.array
        axis of the viewing angles x and y (from the satellite info)
    Returns
    -------
    lat_grid, lon_grid : numpy.array - shape=(len(y), len(x))
        latitude and longitude grids
    '''
    mx, my = np.meshgrid(x, y, indexing = 'xy')  
    lin,col = mx.shape
    mxy = np.array([mx.reshape((lin * col,)), my.reshape((lin * col,))]).T
    lonlat_grid = xy2lonlat(mxy)
    lon_grid = lonlat_grid.T[0].reshape((lin, col))
    lat_grid = lonlat_grid.T[1].reshape((lin, col))
    return lat_grid, lon_grid

def xyAxis2lonlatMat(x_axis, y_axis, log=ignore):
    """
    Creates the lons lats Matrix from the x and y axis
    Arguments
    ---------
    x_axis : np.array() - shape(a,)
        Axis of the alpha sat view anlgles in millirads
    y_axis : np.array() - shape(b,)
        Axis of the beta sat view anlgles in millirads
    Returns
    -------
    np.array() - shape(b,a,2)
        Matrix containing the lon,lat pairs in degrees for each alpha,beta pairs
    """
    mx, my = np.meshgrid(x_axis, y_axis, indexing = 'xy')  
    lins, cols = mx.shape
    mxy = np.array([
        mx.reshape(lins * cols, ), 
        my.reshape(lins * cols, )
    ]).T
    lonlatMat = xy2lonlat(mxy).reshape(lins,cols,2)
    log('lonlat Matrix has been calculated from x and y axis')
    return lonlatMat

def abAxis2lonlatMat(alpha_axis, beta_axis, log=ignore):
    """
    wrapper of xyAxis2lonlatMat using alpha_axis and beta_axis
    """
    x_axis, y_axis = ab2xy(alpha_axis, beta_axis)
    return xyAxis2lonlatMat(x_axis, y_axis, log=log)
    
def rad2tiff_deflate(rad, dst, log=ignore):
    """
    Write tif file from rad variable using deflate compression.
    To read the tif file use :
        arr = np.array(Image.open(src))

    Arguments
    ---------
    rad : numpy.array 
        rad variable to save
    dst : str
        destination path
    log : func
        log function to use
    """
    #--Normalie image array to uint8
    if type(rad) is np.ma.core.MaskedArray: 
        rad = rad.filled(1023) # space values to 1023
    rad[rad < 0] = 0 # night values to 0
    if rad.max() < 1024:
        rad_norm = rad / 4 # to normalize int10 (0-1023) to int8(0-255)
    else:
        raise ValueError ('Not possible to Normalized! '
            'The asumption that the maximun value of the image is < 1024 is not valid. '
            'Rad.max() = {val}'.format(val = rad.max()))
    rad_norm = rad_norm.astype(np.uint8)
    Image.fromarray(rad_norm, mode='L').save(dst,compression='tiff_deflate')
    log('The Rad variable has been written into {file}'.format(file=dst))
    return rad_norm

def var2tiff(var, dst):
    Image.fromarray(var).save(dst)

class G16nc(object):
    """
    Class to handel netCDF files from goes16
    """
    
    def __init__(self, filepath, log=ignore):
        self.filepath = filepath
        self.filenameinfo = self.get_filenameinfo(log)
    
    def load_dataset(self):
        self.dataset = Dataset(self.filepath)
    
    def get_filenameinfo(self, log=ignore):
        _regex_ncdf_file = 'OR_ABI-L1b-RadF-M(?P<abi_num>[\d])C(?P<channel>[\d]{2})_G16_s(?P<dt_start>[\d]{11})[\d]{3}_e(?P<dt_end>[\d]{11})[\d]{3}_c(?P<dt_create>[\d]{11})[\d]{3}(?P<last_num>-[\d]{6}_[\d])?\.nc'
        basename = os.path.basename(self.filepath)
        research = re.search(_regex_ncdf_file, basename)
        log('The info in the file is: \n {dic}'.format(dic=research.groupdict()))
        return research.groupdict()
    
    def read_var(self, var):
        return self.dataset.variables[var][:]
    
    def gen_tiff_name(self):
        datetime = dt.datetime.strptime(self.filenameinfo['dt_start'], '%Y%j%H%M') 
        return '{date:%Y%m%d%H%M}G16C{channel}.tif'.format(date=datetime,**self.filenameinfo)
    
    def close(self):
        self.dataset.close()
        self.dataset = None

    def get_BRF(self, rad):
        kappa = self.read_var('kappa0')
        BRF =  kappa * rad 
        # esun = self.read_var('esun')
        # es_dist = self.read_var('earth_sun_distance_anomaly_in_AU')
        # BRF = (np.pi * es_dist**2)/esun * rad
        return BRF

    def get_BRF_using_defaults(self, rad, esun=1631.3351):
        date = dt.datetime.strptime(self.get_filenameinfo()['dt_start'],'%Y%j%H%M')
        es_dist = solar.earth_sun_dist(solar.solartime(date.timestamp())['eccentricity_factor'])
        BRF = (np.pi * es_dist**2)/esun * rad
        return BRF

