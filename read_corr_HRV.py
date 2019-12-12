#import functools
import os, pickle as pkl
import datetime as dt

import numpy as np

from satlib_py.op import msg_op_daterange as msg_op

#class HDict(dict):
#    def __hash__(self):
#        return hash(frozenset(self.items()))
#def hash_dict(func):
#    """Transform mutable dictionnary
#    Into immutable
#    Useful to be compatible with cache
#    """
#    class HDict(dict):
#        def __hash__(self):
#            return hash(frozenset(self.items()))
#
#    @functools.wraps(func)
#    def wrapped(*args, **kwargs):
#        args = tuple([HDict(arg) if isinstance(arg, dict) else arg for arg in args])
#        kwargs = {k: HDict(v) if isinstance(v, dict) else v for k, v in kwargs.items()}
#        return func(*args, **kwargs)
#    return wrapped

#@hash_dict
#@functools.lru_cache(maxsize=4)


def RMSE(Ddata, Dref):
    return np.sqrt(np.sum((Ddata - Dref) ** 2) / (len(Ddata)-1))

def rMB(Ddata,Dref):
    return (np.sum(Ddata - Dref) / len(Ddata)) / np.mean(Dref) * 100

def get_chan(input):
    chan = getattr(msg_op.extract_and_tif(**input), 'channel')
    if 'VIS' in input['channel']['name']:
        chan.gen_hr_from_lr_section(input['section'])
        input['section'] = input['section'] + '_HR'
    return chan.extractions[input['section']]['img'].reshape(1263*1593,) 

def get_data(datetime,cache_file_fmt):
    try:
        with open (cache_file_fmt, 'br') as fobj:
            dic = pkl.load(fobj)
            v6=dic['v6']
            v8=dic['v8']
            hrv=dic['hrv']
        print('reading from data cached in {file}'.format(file=cache_file_fmt))
    except:
        vis006_info=dict(
                datetime = datetime,
                decomp_bin = '/sat/san/hrit/PublicDecompWT/2.06/xRITDecompress/xRITDecompress',
                datapath = '/sat/san/dvbhrit/dvb/VIS006/{date:%Y}/{date:%m}/{date:%d}/',
                alt_EPIpath = '/sat/san/dvbhrit/dvb/EPI/{date:%Y}/{date:%m}/{date:%d}/',
                alt_PROpath = '/sat/san/dvbhrit/dvb/PRO/{date:%Y}/{date:%m}/{date:%d}/',
                cfg_sections = '/home/user/jorge/repos/satlib_py/src/satlib_py/op/cfg_msg_vis_sections.json',
                channel = msg_op.parse_channel('VIS006'),
                calibration = 'counts',
                output_tif_dir = None,
                norm255 = True,
                verbose = True,
                section = 'creturkis'
                )
        vis008_info=dict(
                datetime = datetime,
                decomp_bin = '/sat/san/hrit/PublicDecompWT/2.06/xRITDecompress/xRITDecompress',
                datapath = '/sat/san/dvbhrit/dvb/VIS008/{date:%Y}/{date:%m}/{date:%d}/',
                alt_EPIpath = '/sat/san/dvbhrit/dvb/EPI/{date:%Y}/{date:%m}/{date:%d}/',
                alt_PROpath = '/sat/san/dvbhrit/dvb/PRO/{date:%Y}/{date:%m}/{date:%d}/',
                cfg_sections = '/home/user/jorge/repos/satlib_py/src/satlib_py/op/cfg_msg_vis_sections.json',
                channel = msg_op.parse_channel('VIS008'),
                calibration = 'counts',
                output_tif_dir = None,
                norm255 = True,
                verbose = True,
                section = 'creturkis'
                )
        hrv_info=dict(
                datetime = datetime,
                decomp_bin = '/sat/san/hrit/PublicDecompWT/2.06/xRITDecompress/xRITDecompress',
                datapath = '/sat/san/dvbhrit/dvb/HRV/{date:%Y}/{date:%m}/{date:%d}/',
                alt_EPIpath = '/sat/san/dvbhrit/dvb/EPI/{date:%Y}/{date:%m}/{date:%d}/',
                alt_PROpath = '/sat/san/dvbhrit/dvb/PRO/{date:%Y}/{date:%m}/{date:%d}/',
                cfg_sections = '/home/user/jorge/repos/satlib_py/src/satlib_py/op/cfg_msg_hrv_sections.json',
                channel = msg_op.parse_channel('HRV'),
                calibration = 'counts',
                output_tif_dir = None,
                norm255 = True,
                verbose = True,
                section = 'creturkis'
                )
        v6 = get_chan(vis006_info)
        v8 = get_chan(vis008_info)
        hrv = get_chan(hrv_info)
        dic = {'v6':v6, 'v8':v8, 'hrv':hrv}
        with open (cache_file_fmt, 'bw') as fobj:
            pkl.dump(dic, fobj)

    return v6, v8, hrv    

def gen_hrv_syn(v6,v8,hrv):
    
    minc, maxc = int(np.min([v8,hrv,v6])) , int(np.max([v8,hrv,v6]))
    
    hist_v6, bins_v6 = np.histogram(v6, bins=range(maxc+1)) 
    hist_v8, bins_v8 = np.histogram(v8, bins=range(maxc+1)) 
    hist_hrv, bins_hrv = np.histogram(hrv, bins=range(maxc+1)) 
    bins = bins_v6[1:]

    cs_hist_v6 = np.cumsum(hist_v6) 
    cs_hist_v8 = np.cumsum(hist_v8) 
    cs_hist_hrv = np.cumsum(hist_hrv) 
    
    norm_val = maxc 
    max_hist = np.max([cs_hist_v6,cs_hist_v8,cs_hist_hrv]) 
    norm_cs_hist_v6=np.round(cs_hist_v6/max_hist*norm_val).astype(np.int)
    norm_cs_hist_v8=np.round(cs_hist_v8/max_hist*norm_val).astype(np.int)
    norm_cs_hist_hrv=np.round(cs_hist_hrv/max_hist*norm_val).astype(np.int)
    
    FDv6 = {} 
    FDv8 = {} 
    FDhrv = {}
    for num,bin in enumerate(bins): 
        FDv6[bin]=norm_cs_hist_v6[num]
        FDv8[bin]=norm_cs_hist_v8[num]
        FDhrv[bin]=norm_cs_hist_hrv[num]

    ## -- max(bins) == max(norm_cs_hist) . both axes are the same
    FDhrv_inv ={0:0}
    j=1
    for bin in bins:
        while(FDhrv[j] < bin):
            j+=1
        FDhrv_inv[bin]=j

    LUT={}
    for bin in bins:
        LUT[bin]=FDhrv_inv[FDv6[bin]]

    hrv_syn = np.array(list(map(lambda x : LUT[x], v6.astype(np.int))))

    return hrv_syn

##--- main  ---
#datetime = dt.datetime(2018,10,24,9,0)
datetime = dt.datetime(2018,10,1,15,0)
cache_file_fmt = '/tmp/{date:%Y%m%d%H%M}.pkl'.format(date=datetime)
v6, v8, hrv = get_data(datetime,cache_file_fmt)
#hrv_syn = gen_hrv_syn(v6,v8,hrv)

	


