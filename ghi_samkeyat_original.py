
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
    DATADIR='/sat/san/tempfs/jorge/samkeyat/test_ghi/data/'
    channel = 2
    section = 'ecuador'
    filefmt = '/sat/san/tempfs/jorge/samkeyat/test_ghi/{section}/{type}/{date:%Y%m%d%H%M}.h5'
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
        func = Func_ghi_ground(h5u_cld_idx),
        log = print
    )

    for datetime in etc.yield_daterange('201805301700:201805311700:1d'): 
#    for datetime in etc.yield_daterange('201708011700:201708311700:1d'): 
#    for datetime in etc.yield_daterange('201708041700'): 
        print('for: ', datetime)
        try:
#            BRF = h5u_BRF(datetime, channel, section)
#            albedo = h5u_alb(datetime, channel, section)
#            cld_idx = h5u_cld_idx(datetime, channel, section)
            ghi = h5u_ghi(datetime, channel, section)
        except Exception as e:
            traceback.print_exc()
            





#    for file in glob.glob(os.path.join(DATADIR, '201708??1700-ecuador.h5')):
#    for file in glob.glob(os.path.join(DATADIR, '*-ecuador.h5')):
#        filename = '{date:%Y%m%d%H%M}-ecuador.h5'.format(date=date)
#        datetime = dt.datetime.strptime(re.search('([\d]{12})-ecuador.h5', filename).groups()[0],'%Y%m%d%H%M')
#        print(np.max(ghi.values))

#    data = ghi
#    plt.imshow(data.values.astype(np.float))
#    plt.show()
#
#    print(np.max(data.values))




        
    
#    
#    alb = h5u_alb(date, channel, section)
#    rhomax = heliosat.Rhomax_percentile(h5u_BRF, 95)(date, channel, section)
#    cld_idx = h5u_cld_idx(date, channel, section)
#    ghi = h5u_ghi(date, channel, section)
    
#    fig, ax = plt.subplots(2,1)
#    ax[0].imshow(cld_idx.values.astype(np.float))
#    ax[1].imshow(ghi.values.astype(np.float))
#    plt.savefig('/sat/san/tempfs/jorge/py_figs/ghi.png')

#    for date in np.arange('2019-04-01', '2019-05-27', dtype='datetime64[D]') + np.timedelta64(60*11,'m'):
#        BRF = h5u_BRF(date.item(), channel, section)

#    for date in np.arange('2019-05-01', '2019-05-19', dtype='datetime64[D]') + np.timedelta64(60*11,'m'):
#        print(date)
##        try:
##        BRF = h5u_BRF(date.item(), channel, section)
##        alb = h5u_alb(date.item(), channel, section)
##        cld_idx = h5u_cld_idx(date.item(), channel, section)
#        print(h5u_cld_idx(date.item(), channel, section).attrs['rhomax'])
##        ghi = h5u_ghi(date.item(), channel, section)
##        except Exception as e:
##            print(repr(e))
#    #--------------------------------------------
#    #--------------------------------------------
#    ## comparison with Annettes Images ##
#
#    roll_lin_col = (4,3) 
#    max_ci = 3
#    min_ci = -1
#    
#    ci=cld_idx.values.astype(np.float)
#    ci=np.roll(ci,roll_lin_col,(0,1))
#    ##TODO limits uppper/bottom?
##    ci[ci>max_ci] = max_ci
##    ci[ci<min_ci] = min_ci
#    header, data = xpif.read_xpif(
#        '/sat/san/products/msg/{section}/alb_ci/{date:%Y}/{date:%m}/{date:%Y%m%d%H%M}_{section}_vishrv.CI.XPIF'.format(date=date, section=section))
#    ci2 = 1.4/1023 * data - 0.2
#
#    ci = ci[roll_lin_col[0]:, roll_lin_col[1]:]
#    ci2 = ci2[roll_lin_col[0]:, roll_lin_col[1]:]
#    fig, ax = plt.subplots(3,1)
#    im = ax[0].imshow(ci.astype(np.float))
#    clim = im.properties()['clim']
#    ax[1].imshow(ci2, clim=clim)
#    ax[2].imshow(ci - ci2, clim=clim)
#    at = AnchoredText("mean bias = {mb:.5f}\nRMSE = {rmse:.5f}".format(
#        mb= (error.img_bias(ci, ci2)), rmse= (error.img_RMSE(ci, ci2))), prop=dict(size=8), frameon=True, loc=2)
#    ax[2].add_artist(at)                                                                             
#    fig.colorbar(im, ax=ax.ravel().tolist(), shrink=0.5)
#    plt.show()
#







#    fig, ax = plt.subplots(2,1)
#    im = ax[0].imshow(ci.astype(np.float))
#    ax0 = ax[0]; ax0.set_ylim([0,30]); ax0.set_xlim([660,768])
#    clim = im.properties()['clim']
#    ax[1].imshow(ci2, clim=clim)
#    ax1 = ax[1]; ax1.set_ylim([0,30]); ax1.set_xlim([660,768])
#    fig.colorbar(im, ax=ax.ravel().tolist(), shrink=0.5)
#    plt.show()
    
