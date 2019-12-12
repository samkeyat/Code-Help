
import numpy as np

from satlib_py.utils.etc import LastNdays
from satlib_py.utils import solar
#from satlib_py.satellites import msg, goes16
from satlib_py.proc import datatype

def ignore(*args, **kwargs):
    pass

class Albedo_percentile(object):
    
    def __init__(self, percentile, h5u_BRF, days_hist=30, log=ignore):
        self.days_hist = days_hist 
        self.percentile = percentile
        self.h5u_BRF = h5u_BRF
        self.log = log
    
    def gen_albedo(self, datetime, channel, section):
        imgs = []
        for past_date in LastNdays(self.days_hist)(datetime):
            try:
                BRF = self.h5u_BRF(past_date, channel, section)
                imgs.append(BRF.values)
                attrs = BRF.attrs
            except AttributeError:
                self.log('No BRF image found for the section {section} for channel {channel} on {date:%Y%m%d%H%M}'.format(
                    date=past_date, section=section, channel=channel))
        return np.percentile(imgs, self.percentile, axis=0), attrs

    __call__ = gen_albedo


class Rhomax_percentile(object):
    
    def __init__(self, h5u_BRF, percentile=95 , days_hist=30, only_values_over=50, log=ignore):
        ## TODO default value for only_values_over (using dichotomie on canary images?)
        ## TODO rhomas tables
        ## TODO sensitivity on percentile (compared to Xpifs)
        ## TODO find min number of points needed for the percentile to work -> if not default
        ## TODO default rawmax
        self.days_hist = days_hist 
        self.percentile = percentile
        self.h5u_BRF = h5u_BRF
        self.only_values_over = only_values_over
        self.log = log
    
    def gen_rhomax(self, datetime, channel, section):
        imgs=[]
        for past_date in LastNdays(self.days_hist)(datetime):
            try:
                imgs.append(self.h5u_BRF(past_date, channel, section).values)
            except AttributeError:
                self.log('No BRF image found for the section {section} for channel {channel} on {date:%Y%m%d%H%M}'.format(
                    date=past_date, section=section, channel=channel))
#        imgs2 = np.array([self.h5u_BRF(past_date, channel, section).values for past_date in LastNdays(self.days_hist)(datetime)])
        imgs = np.array(imgs)
        return np.percentile(imgs[imgs>self.only_values_over], self.percentile)

    __call__ = gen_rhomax

class Cloud_index(object):

    def __init__(self, h5u_BRF, h5u_alb, rhomax_func, log=ignore):
        self.h5u_BRF = h5u_BRF
        self.h5u_alb = h5u_alb
        self.rhomax_func = rhomax_func
        self.log = log

    def gen_cloud_index(self, *args, **kwargs):
        datetime = args[0]
        channel = args[1]
        section = args[2]
        BRF = self.h5u_BRF(*args, **kwargs)
        if not BRF:
            raise RuntimeError('no cloud index img generated for the section {section} for channel {channel}'
                               ' on {datetime:%Y%m%d%H%M} : No BRF img found'.format(section=section, 
                               channel=channel, datetime=datetime))
        alb = self.h5u_alb(*args, **kwargs)
        if not alb:
            raise RuntimeError('no cloud index img generated for the section {section} for channel {channel}'
                               ' on {datetime:%Y%m%d%H%M} : No albedo img found'.format(section=section, 
                               channel=channel, datetime=datetime))
        rhomax = self.rhomax_func(*args, **kwargs)
        if not rhomax:
            raise RuntimeError('no cloud index img generated for the section {section} for channel {channel}'
                               ' on {datetime:%Y%m%d%H%M} : No rhomax value found'.format(section=section, 
                               channel=channel, datetime=datetime))
#        cloud_index = (self.h5u_BRF(*args, **kwargs).values - self.h5u_alb(*args, **kwargs).values)/(
#                       rhomax - self.h5u_alb(*args, **kwargs).values)
        cloud_index = (BRF.values - alb.values) / (rhomax - alb.values)
        cloud_index[cloud_index<0.05]=0.05
        cloud_index[cloud_index>1]=1
        attrs = BRF.attrs
        attrs['rhomax'] = rhomax
        return cloud_index, attrs 

    __call__ = gen_cloud_index


class Ghi_ground(object):

    def __init__(self, h5u_cld_idx, abAxis2lonlatMat_func, log=ignore):
        self.h5u_cld_idx = h5u_cld_idx
        self.ab2lonlatMat = abAxis2lonlatMat_func 
        self.log = log

    def gen_ghi(self, datetime, channel, section):
        cld_idx = self.h5u_cld_idx(datetime, channel, section)
        lonlatMat =self.ab2lonlatMat(
                alpha_axis = datatype.axis_info_2_values(cld_idx.attrs['alpha_info']),
                beta_axis = datatype.axis_info_2_values(cld_idx.attrs['beta_info']),
        )
        rolled_lonlatMat = np.rollaxis(lonlatMat,2,0)
        geo_map = solar.Map(lats = rolled_lonlatMat[1], lons = rolled_lonlatMat[0])
        ghi_cs = geo_map.get_clearsky_dumortier(datetime.timestamp())
        ghi = (1-cld_idx.values) * ghi_cs['ghi']
        return ghi, cld_idx.attrs
    
#    def gen_ghi_goes(self, datetime, channel, section):
#        cld_idx = self.h5u_cld_idx(datetime, channel, section)
#
#        lonlatMat = goes16.abAxis2lonlatMat(
#                alpha_axis = datatype.axis_info_2_values(cld_idx.attrs['alpha_info']),
#                beta_axis = datatype.axis_info_2_values(cld_idx.attrs['beta_info']),
#        )
#        rolled_lonlatMat = np.rollaxis(lonlatMat,2,0)
#        geo_map = solar.Map(lats = rolled_lonlatMat[1], lons = rolled_lonlatMat[0])
#        ghi_cs = geo_map.get_clearsky_dumortier(datetime.timestamp())
#        ghi = (1-cld_idx.values) * ghi_cs['ghi']
#        return ghi, cld_idx.attrs

    __call__ = gen_ghi

