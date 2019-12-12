import os

import h5py, numpy as np

import matplotlib.pyplot as plt

def ignore(*args, **kwargs):
    pass

def axis_values_2_info(ax):
    if len(ax): 
        return np.array([ax[0], ax[-1], len(ax)])
    else:
        return ax

def axis_info_2_values(info):
    if len(info):
        return np.linspace(info[0], info[1], info[2], endpoint=True)
    else:
        return info


class H5Update(object):
    def __init__(self, datatype, filefmt, func, log=ignore):
        self.datatype = datatype
        self.filefmt = filefmt
        self.func = func
        self.log = log
    
    def update_data(self, *args, **kwargs):
        try:
            return self.h5_load(*args, **kwargs)
        except OSError:
            pass
        ##TODO are there other cases where load is not done other than OSError
        return self.h5_store(*args, **kwargs)

    __call__ = update_data
    
    def h5_load(self, *args, **kwargs):
        instance = self.datatype(*args, **kwargs) ##initializing datatype class
        filepath = self.filefmt.format(**instance.__h5context__())
#        filepath = instance.get_filepath()
        data = h5py.File(filepath, 'r', swmr=True)
        if data.attrs.get('__datatype__') != self.datatype.__name__:
            type_h5 =data.attrs.get('__datatype__')
            data.close()
            raise ValueError('Invalid __datatype__ "{h5}" in cache file "{file}", expecting "{dtype}"'.format(
                h5=type_h5, file=filepath, dtype=self.datatype.__name__)
            )
        instance.__h5load__(data)
        self.log('data has been loaded from {filename}'.format(filename=filepath))
        return instance
    
    def h5_store(self, *args, **kwargs):
        instance = self.func(*args, **kwargs) ## Calling func object
        if not instance:
            return None
        filepath = self.filefmt.format(**instance.__h5context__())
#        filepath = instance.get_filepath()
        if instance.__class__.__name__ != self.datatype.__name__:
            raise ValueError('__datatype__ missmatch: __datatype__ from argument datatype ({data})'
                ' should coincided with __datatype__ from argument func ({func})'.format(
                    data=self.datatype.__name__, func=instance.__class__.__name__
                )
            )
        dirname = os.path.dirname(filepath)
        if dirname:
            os.makedirs(dirname, exist_ok=True)
        # FIXME Handle cases where file opened twice.
        h5_ds = h5py.File(filepath, 'w', swmr=True)
        h5_ds.attrs['__datatype__'] = self.datatype.__name__
        instance.__h5store__(h5_ds)
        self.log('{filename} has been created'.format(filename=filepath))
        return instance


class SatImg(object):
    """ General class for satellite related images. Most general case """

    def __init__(self, datetime, channel='', section='', type='SatImg', dtype=np.float):
        self._datetime = datetime
        self._channel = channel
        self._section = section
        self._type = type
        self._dtype = dtype

    def __h5context__(self):
        info =  dict(
            date=self._datetime, 
            channel=self._channel, 
            section=self._section,
            type = self._type
        )
        keys2pop = [key for key in info if info[key] == '']
        for key in keys2pop:
            info.pop(key)
        return info

    def __h5store__(self, h5_ds):
        values_ds = h5_ds.create_dataset('values', data=self.values, compression='gzip')
        ##TODO: optimize dtype to int ?
        values_ds.attrs.update(self.attrs)
        return self

    def __h5load__(self, h5_ds):
        self.attrs = dict(h5_ds['values'].attrs)
        if '{date:%Y%m%d%H%M}'.format(date=self._datetime) != self.attrs['datetime_str']:
            h5_ds.close()
            raise ValueError('Mismatching timestamps')
        self.values = h5_ds['values'][:].astype(self._dtype)

    def set_data(self, values, **kwargs):
        self.values = values.astype(self._dtype)
        self.attrs = dict(
                datetime_str = '{date:%Y%m%d%H%M}'.format(date=self._datetime),
                channel = self._channel,
                section = self._section,
                type = self._type,
                **kwargs
                )
        return self

    def show_img(self):
        plt.imshow(self.values.astype(np.float)) ## imshow does not accept np.float16 type, so forced np.float64
        plt.show()

    __call__ = set_data

class BRFImg(SatImg):
    """ child class to force consistency on BRF images """
    
    def __init__(self, *args, **kwargs):
        kwargs['type'] = 'BRF'
        super().__init__(*args, **kwargs)
    
    def set_data(self, values, 
            sat_name='',
            cal_slope='', 
            cal_offset='', 
            alpha_info='', 
            beta_info=''):
        self.attrs = {
                'type':self._type,
                'datetime_str':'{date:%Y%m%d%H%M}'.format(date=self._datetime),
                'sat_name': sat_name,
                'channel':self._channel,
                'section':self._section,
                'cal_slope':cal_slope,
                'cal_offset':cal_offset,
                'alpha_info':alpha_info,
                'beta_info':beta_info
                }
        self.values = values.astype(self._dtype)
        return self

    __call__ = set_data

class Albedo(SatImg):
    """ child class to force consistency on Albedo images """
    
    def __init__(self, *args, **kwargs):
        kwargs['type'] = 'albedo'
        super().__init__(*args, **kwargs)
    
    def set_data(self, values, 
            sat_name='',
            alpha_info='', 
            beta_info=''):
        self.attrs = {
                'type':self._type,
                'datetime_str':'{date:%Y%m%d%H%M}'.format(date=self._datetime),
                'sat_name':sat_name,
                'channel':self._channel,
                'section':self._section,
                'alpha_info':alpha_info,
                'beta_info':beta_info
                }
        self.values = values.astype(self._dtype)
        return self

    __call__ = set_data

class Cloud_index(SatImg):
    """ child class to force consistency on cloud index images """

    def __init__(self, *args, **kwargs):
        kwargs['type'] = 'cloud_index'
        super().__init__(*args, **kwargs)
    
    def set_data(self, values, 
            sat_name='',
            alpha_info='', 
            beta_info='',
            rhomax=''):
        self.attrs = {
                'type':self._type,
                'datetime_str':'{date:%Y%m%d%H%M}'.format(date=self._datetime),
                'sat_name':sat_name,
                'channel':self._channel,
                'section':self._section,
                'alpha_info':alpha_info,
                'beta_info':beta_info,
                'rhomax':rhomax
                }
        self.values = values.astype(self._dtype)
        return self

    __call__ = set_data

class Ghi(SatImg):
    """ child class to force consistency on GHI images """
    
    def __init__(self, *args, **kwargs):
        kwargs['type'] = 'ghi'
        super().__init__(*args, **kwargs)
    
    def set_data(self, values, 
            sat_name='',
            alpha_info='', 
            beta_info=''):
        self.attrs = {
                'type':self._type,
                'datetime_str':'{date:%Y%m%d%H%M}'.format(date=self._datetime),
                'sat_name':sat_name,
                'channel':self._channel,
                'section':self._section,
                'alpha_info':alpha_info,
                'beta_info':beta_info
                }
        self.values = values.astype(self._dtype)
        return self

    __call__ = set_data
