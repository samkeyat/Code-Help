
import os, datetime as dt, glob, re, traceback, sys

import numpy as np, h5py 
import matplotlib.pyplot as plt

from satlib_py.proc import datatype, heliosat, error
from satlib_py.utils import solar, etc
from satlib_py.satellites import goes16

class Func_BRF(object):
    def __init__(self, datadir):
        self.datadir = datadir

    def get_BRF(self, datetime, channel, section):
        ESUN=1631.3351
        try:
            File = h5py.File(
                    os.path.join(self.datadir, '{datetime:%Y%m%d%H%M}-{section}.h5'.format(
                                 datetime=datetime, section=section
                    ))
            , 'r')
        except OSError as e:
            return None 

#        import pdb; pdb.set_trace()
        rad = File['counts'][:]
        x = File['counts'].attrs['x_0_step_len']
        y = File['counts'].attrs['y_0_step_len']
        alpha = [-x[0], -x[1], x[2]]
        beta = y
        es_dist = solar.earth_sun_dist(solar.solartime(datetime.timestamp())['eccentricity_factor'])
        kappa0 = (np.pi * es_dist**2) / ESUN 
        BRF = kappa0 * rad * 100
        return datatype.BRFImg(datetime, channel, section)(
                        values = BRF,
                        alpha_info = [alpha[0], alpha[0] +alpha[1] * alpha[2] - alpha[1], alpha[2]],
                        beta_info = [beta[0], beta[0] +beta[1] * beta[2] - beta[1], beta[2]],
                        sat_name = 'GOES16',
                        )

    __call__ = get_BRF

class Func_albedo(object):
    def __init__(self, percentile, h5u_BRF, days_hist=30):
        self.h5u_BRF = h5u_BRF
        self.percentile = percentile
        self.days_hist = days_hist

    def get_albedo(self, *args, **kwargs):
        alb, attrs = heliosat.Albedo_percentile(self.percentile, self.h5u_BRF, self.days_hist, log=print)(*args, **kwargs) 
        return datatype.Albedo(*args, **kwargs)(
                values = alb,
                alpha_info = attrs['alpha_info'],
                beta_info = attrs['beta_info'],
                sat_name = attrs['sat_name']
                ) 

    __call__ = get_albedo

class Func_cloud_index(object):
    def __init__(self, h5u_BRF, h5u_alb, rhomax_func):
        self.h5u_BRF = h5u_BRF
        self.h5u_alb = h5u_alb
        self.rhomax_func = rhomax_func

    def get_cloud_index(self, *args, **kwargs):
        cloud_index, attrs = heliosat.Cloud_index(self.h5u_BRF, self.h5u_alb, self.rhomax_func, log=print)(*args, **kwargs) 
        return datatype.Cloud_index(*args, **kwargs)(
                values = cloud_index,
                alpha_info = attrs['alpha_info'],
                beta_info = attrs['beta_info'],
                sat_name = attrs['sat_name'],
                rhomax = attrs['rhomax']
                ) 

    __call__ = get_cloud_index
#
class Func_ghi_ground(object):
    def __init__(self, h5u_cld_idx):
        self.h5u_cld_idx = h5u_cld_idx

    def get_ghi_ground(self, *args, **kwargs):
        ghi_ground, attrs = heliosat.Ghi_ground(self.h5u_cld_idx, goes16.abAxis2lonlatMat, log=print)(*args, **kwargs)
        return datatype.Ghi(*args, **kwargs)(
                values = ghi_ground,
                alpha_info = attrs['alpha_info'],
                beta_info = attrs['beta_info'],
                sat_name = attrs['sat_name']
                )

    __call__ = get_ghi_ground

if __name__ == '__main__':
    ## Inputs ##
    ## daterange =
    ##    'yyyymmddHHMM' --> only one timestamp
    ##      ex= '2018051700'
    ##    'yyyymmddHHmm:yyyyddmmHHMM:N[d,h,m,s]' --> range of time stamps
    ##      ex= '201805010000:201805050000:15m'
#    DATERANGE='201707011315'
    DATERANGE='201707011315:201708011315:15m'
    ## do not change below##
    ## Inputs ##
    DATADIR='/user/auzo8195/work/samkeyat/goes_data'
    channel = 2
    section = 'ecuador'
    filefmt =  '/user/auzo8195/work/samkeyat/ghi/{section}/{type}/{date:%Y%m%d%H%M}.h5'
    days_hist = 30
    percentile_alb = 5
    percentile_rho = 95
    #############

    h5u_BRF = datatype.H5Update(
        datatype = datatype.BRFImg,
        filefmt = filefmt,
        func = Func_BRF(DATADIR),
        log=print
    )
    
    h5u_alb =  datatype.H5Update(
        datatype = datatype.Albedo,
        filefmt = filefmt,
        func = Func_albedo(percentile_alb, h5u_BRF, days_hist),
        log=print
    )
    
    h5u_cld_idx =  datatype.H5Update(
        datatype = datatype.Cloud_index,
        filefmt = filefmt,
        func = Func_cloud_index(
            h5u_BRF, 
            h5u_alb, 
            heliosat.Rhomax_percentile(h5u_BRF, percentile_rho, days_hist, only_values_over=50, log=print)
        ),
        log=print
    )

    h5u_ghi = datatype.H5Update(
        datatype = datatype.Ghi,
        filefmt = filefmt,
#        h5u_cld_idx = 0.8<h5u_cld_idx<=1.05, 
        func = Func_ghi_ground(h5u_cld_idx),
        log = print
    )

 
    for datetime in etc.yield_daterange(DATERANGE): 
        print('for: ', datetime)
        try:
            #BRF = h5u_BRF(datetime, channel, section)
            #albedo = h5u_alb(datetime, channel, section)
            #cld_idx = h5u_cld_idx(datetime, channel, section)
            ghi = h5u_ghi(datetime, channel, section)
            pass
        except Exception as e:
            traceback.print_exc()
            



