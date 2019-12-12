import re 

import datetime as dt, numpy as np

def ignore(*args, **kwargs):
    pass

class LastNdays(object): 
    def __init__(self, n): 
        self.n=n 
     
    def get_dates(self, date): 
        start = date - dt.timedelta(days = self.n) 
        for day in range(0,self.n): 
            yield (start + dt.timedelta(days = day))
     
    __call__ = get_dates          

def yield_daterange(s):
    """ 
    yileds datetime objects for the datetime range specified on 'yyyymmddHHMM:yyyymmddHHMM:N[d/h/m/s - default m]'
        example : 201906011200:201906051200:10  --> default = minutes
        example : 201906011200:201906051200:10m 
    """
    resolutions={
        'd':'days',
        'h':'hours',
        'm':'minutes',
        's':'seconds',
    }
    strings = s.split(':')
    start = dt.datetime.strptime(strings[0], '%Y%m%d%H%M')

    if (np.array(list(map(len, strings)))[0:2]!=12).any():
        raise ValueError('input must be of the form "yyyymmddHHMM:yyyymmddHHMM:N[d/h/m/s]"')

    if len(strings) == 1:
        yield start
    elif len(strings) == 3:
        end = dt.datetime.strptime(strings[1], '%Y%m%d%H%M')
        if end <= start:
            raise RuntimeError('end date {end} must be grater than start date {start}'.format(
                                start=strings[0], end=strings[1]))
        try :
            res_dic = re.search('^(?P<num>[\d]*)(?P<res>[dhms])?$', strings[2]).groupdict() 
            if not(res_dic['res']):
                res_dic['res'] = 'm' ## default in minutes - to test other codes
        except AttributeError as e :
            raise RuntimeError('resolution "{res}" not valid'.format(res=strings[2]))
        timedelta = dt.timedelta(**{resolutions[res_dic['res']]:int(res_dic['num'])})
        for num in range(int((end-start)/timedelta) + 1):
            yield(start + num * timedelta)
    else :
        raise ValueError('input must be of the form "yyyymmddHHMM:yyyymmddHHMM:N[d/h/m/s]"')

def lonlatMat2bin_annette(lat_grid, lon_grid, dst, log=ignore):
    """
    In order to use the C code (Annette Hammer)
    wirtes [lat1, lon1, lat2, lon2, ..., latN, lonN] in binary file in dtype int16
        -point (1,1) = northeast grid point
        -point (N,N) = southwest grid point
        -nan = 99
    lat/lon values are rounded to 2 decimals and multiplied by 100 to fit in int16 
        -e.g. "-75.324" is the int16 "-7532" 
        -e.g  nan is the int16 "9900"
    endianess = Least Significan Byte First
    -
    To read the file
        np.fromfile(src, dtype='int16')
    """
    shape_grid = lat_grid.shape

    lat_grid[np.isnan(lat_grid)] = 99 
    lon_grid[np.isnan(lon_grid)] = 99

    lats = (lat_grid.reshape(shape_grid[0] * shape_grid[1], ).round(2)*100).astype(np.int16)
    lons = (lon_grid.reshape(shape_grid[0] * shape_grid[1], ).round(2)*100).astype(np.int16)
    latlon_arr = np.array([lats, lons]).T.reshape(len(lats)*2, )
    latlon_arr.tofile(dst) 
    log('The binary latlon file {file} has been written'.format(file=dst))
