import os
from functools import lru_cache

import numpy as np, h5py
from scipy.interpolate import RegularGridInterpolator

SOLAR = 1367
JULIAN_2000 = 2451544.5
DAY_SECONDS = 60 * 60 * 24
YEAR_JULIANS = 365.2425
UTC_TSTAMP_2000 = 946684800.0


def utc_tstamps2julians(utc_tstamps):
        return JULIAN_2000 + (utc_tstamps - UTC_TSTAMP_2000) / 86400 

def solartime(utc_tstamps):
    ##TODO eot has small differences with robust models (spa or Uni oregon sun calc)
    """
    Computes equation of time 
    .. note::
        Formulas from Spencer (1972) and can be also found in "Solar energy
        fundamentals and modeling techniques" from Zekai Sen
    """
    try:
        len(utc_tstamps)
    except TypeError:
        utc_tstamps = np.array([utc_tstamps])
    # Bind names for convenience.
    radians, pi, sin, cos, = (np.radians, np.pi, np.sin, np.cos)
    # Compute julian dates relative to 2000-01-01 00:00.
    julians = JULIAN_2000 + (utc_tstamps - UTC_TSTAMP_2000) / 86400 
    julians_2000 = np.asarray(julians, dtype=np.float) - JULIAN_2000
    # Compute fractional year (gamma) in radians
    gamma = 2 * pi * (julians_2000 % YEAR_JULIANS) / YEAR_JULIANS
    cos_gamma = cos(gamma), cos(gamma * 2), cos(gamma * 3)
    sin_gamma = sin(gamma), sin(gamma * 2), sin(gamma * 3)
    day_time = (julians_2000 % 1) * 24
    # Eccentricity factor: correction factor of the earth's orbit.
    E_null = (1.00011 + 0.034221 * cos_gamma[0] + 0.001280 * sin_gamma[0] +
            0.000719 * cos_gamma[1] + 0.000077 * sin_gamma[1])
    # declination.
    declination = (0.006918 - 0.399912 * cos_gamma[0] +
            0.070257 * sin_gamma[0] -
            0.006758 * cos_gamma[1] + 0.000907 * sin_gamma[1] -
            0.002697 * cos_gamma[2] + 0.001480 * sin_gamma[2])
    # Equation of time (difference between standard time and solar time).
    eot = (0.000075 + 0.001868 * cos_gamma[0] - 0.032077 * sin_gamma[0] -
            0.014615 * cos_gamma[1]  - 0.040849 * sin_gamma[1]) * 229.18
    return {
            'eccentricity_factor':E_null,
            'declination':declination,
            'eot':eot,
    }


def local_solar_params(utc_tstamps, lats, lons):
    """
    computes solar position parameters (apparent zenith [radians], apparent azimuth[radinas], 
    ghi_ext[w.m-2]) relative to a location on earth for the given points in time.

    arguments
    =========
    utc_tstamps : numpy.array shape=(n,) or float
    lats : numpy.array shape=(n,m) or float
    lons : numpy.array shape=(n,m) or float

    .. note::
        * case Timeseries : when utc_tstamps shape is (n,), shape of lats and lons must be (1,) or float
        * case Map : when lats and lons shape is (n,), shape of utc_tstamps must be (1,) or float
        * formulas from spencer (1972) and can be also found in "solar energy
        fundamentals and modeling techniques" from zekai sen
    """
    try:
        len(utc_tstamps)
    except TypeError:
        utc_tstamps = np.array([utc_tstamps])
    # Bind names for convenience.
    radians, pi, sin, cos, arcsin, arccos = (
        np.radians, np.pi, np.sin, np.cos, np.arcsin, np.arccos)
    # Compute julian dates relative to 2000-01-01 00:00.
    julians = JULIAN_2000 + (utc_tstamps - UTC_TSTAMP_2000) / 86400 
    julians_2000 = np.asarray(julians, dtype=np.float) - JULIAN_2000
    lats, lats_deg = radians(lats), lats
    lons, lons_deg = radians(lons), lons
    day_time = (julians_2000 % 1) * 24
    # solartime variables
    stime =  solartime(utc_tstamps) 
    eccentricity_factor, declination, eot = (
            stime['eccentricity_factor'], stime['declination'], stime['eot']
    )
    # True local time
    tlt = (day_time + lons_deg / 15 + eot / 60) % 24 - 12
    # Solar hour angle
    ha = radians(tlt * 15)
    # Calculate sun elevation.
    sin_sun_elevation = (
        sin(declination) * sin(lats) + cos(declination) * cos(lats) * cos(ha)
    )
    # Compute the sun's elevation and zenith angle.
    el = arcsin(sin_sun_elevation)
    zenith = pi / 2 - el
    # Compute the sun's azimuth angle.
    y = -(sin(lats) * sin(el) - sin(declination)) / (cos(lats) * cos(el))
    azimuth = arccos(y)
    # Convert azimuth angle from 0-pi to 0-2pi.
    tlt_filter = 0 <= tlt
    azimuth[tlt_filter] = 2 * pi - azimuth[tlt_filter]
    # Calculate the extraterrestrial radiation.
    ghi_ext = 1367.7 * sin_sun_elevation * eccentricity_factor
    return {
        'zenith': zenith,
        'azimuth_rad': azimuth,
        'ghi_ext': ghi_ext,
    }

def AM_rozenberg(zenith):
    """
    airmas calculation rozenberg

    arguemnts
    =========
    zenith : np.array shape=(n,m)
        
    """
    try:
        len(zenith)
        zenith_am = zenith.copy()
        zenith_am[zenith_am > np.pi/2] = np.pi/2
    except TypeError:
        if zenith > np.pi/2:
            zenith_am = np.pi/2
        else:
            zenith_am = zenith
    return 1 / (np.cos(zenith_am) + 0.025 * np.exp(-11*np.cos(zenith_am)))

def earth_sun_dist(eccentricity_factor):
    return np.sqrt(1/eccentricity_factor)

@lru_cache(maxsize=1)
def remund_map_interpolator(filename):
    try:
        with h5py.File(filename, 'r', swmr=True) as db:
            data = db['map'][:]
            ## Append january data to allow smooth interpolation over the year
            ## so that it is a continues value
            data = np.append(data, data[:1], axis=0)
    except OSError as e:
        raise OSError(
            'Unable to open file \'{filename}\''.format(filename=filename))
    julian_of_year = np.linspace(0, YEAR_JULIANS, data.shape[0], endpoint=True)
    lats = np.linspace(-90, 90, data.shape[1], endpoint=True)
    lons = np.linspace(-180, 180, data.shape[2], endpoint=True)
    return RegularGridInterpolator(
        (julian_of_year, lats, lons), data, bounds_error=False,
        fill_value=np.nan)

default_turbidity_h5 = os.path.join(__file__[:__file__.find('utils/solar')], 'data/remund.h5')
def turbidity_remund(utc_tstamps, lats, lons, turbidity_filename=None):
    """
    turbidity map from meteonorm
    
    arguments
    =========
    utc_tstamps : numpy.array shape=(n,) or float
    lats : numpy.array shape=(n,) or float
    lons : numpy.array shape=(n,) or float
    turbidity_filename : string
        path to the file containing turbidity info

    .. note::
        * case Timeseries : when utc_tstamps shape is (n,), shape of lats and lons must be (1,) or float
        * case Map : when lats and lons shape is (n,), shape of utc_tstamps must be (1,) or float
    """
    julians = JULIAN_2000 + (utc_tstamps - UTC_TSTAMP_2000) / 86400 
    if turbidity_filename is None:
        turbidity_filename = default_turbidity_h5
    return remund_map_interpolator(turbidity_filename)(
        ((julians - JULIAN_2000) % YEAR_JULIANS, lats, lons))


@lru_cache(maxsize=1)
def height_map_interpolator(filename):
    with h5py.File(filename, 'r', swmr=True) as db:
        data = db['map'][:]
        attr = dict(db['map'].attrs)
    lats = np.linspace(attr['y1'], attr['y2'], data.shape[0], endpoint=True)
    lons = np.linspace(attr['x1'], attr['x2'], data.shape[1], endpoint=True)
    return RegularGridInterpolator((lats, lons), data, bounds_error=False,
            fill_value=np.nan)
    
default_heightmap_h5 = os.path.join(__file__[:__file__.find('utils/solar')], 'data/heightmap.h5')
def height_map(lats, lons, height_filename=None):
    """
    get heigth for any lat,lon point in the world 
    
    arguments
    =========
    lats : numpy.array shape=(n,) or float
    lons : numpy.array shape=(n,) or float
    turbidity_filename : string
        path to the file containing turbidity info
    """
    if height_filename is None:
        height_filename = default_heightmap_h5
    return height_map_interpolator(height_filename)((lats, lons))

####################################################
### All conditions been equal,
### for a minute resolution time series for one day:
### clearsky_dumortier = 2.7 ms
### pvlib inichien =  58 ms
### pvlib simplified solis = 48 ms
#####################################################
def clearsky_dumortier(zenith, eccentricity_factor, turbidity, height):
    """
    Clearsky irradiation following Dumortier Model

    arguments
    =========
        zenith : np.array shape=(n,1) or float
            same shape as the utc_timestamps studied
        eccentricity_factor: np.array shape=(n,1) or float
            shape as the utc_timestamps studied
        turbidity:  np.array shape=(n,m) or float
            for a Map same shape of lats, lons and for a timeseries same shape of utc_timestamps
        height : numpy.array shape=(n,m) or floats
            same shape of lats, lons
    """
    with np.errstate(invalid='ignore'):
        coszen = np.cos(zenith)
        elevation_degree = 90 - 180 / np.pi * zenith
        m = (
            (1 - height / 10000) /
            (coszen + 0.50572 * (6.07995 + elevation_degree) ** -1.6364)
        )
        selection = m <= 20
        delta = np.empty_like(m)
        m_s = m[selection]
        delta[selection] = (
            1 / (
                6.6296 + 1.7513 * m_s - 0.1202 * m_s ** 2 + 0.0065 * m_s ** 3 -
                0.00013 * m_s ** 4
            )
        )
        delta[~selection] = 1 / (10.4 + 0.718 * m[~selection])
        f1 = coszen * np.exp(-0.8662 * turbidity * m * delta)
        f2 = (
            0.0065 + coszen * (-0.045 + 0.0646 * turbidity) -
            coszen ** 2 * (-0.014 + 0.0327 * turbidity)
        )
        I_cl = SOLAR * (f1 + f2) * eccentricity_factor
        k_dir = f1 / (f1 + f2)
        selection = coszen <= 0
        I_cl[selection] = 0
        k_dir[selection] = 0
    return {'ghi': I_cl, 'dni_fraction': k_dir}



################################################
## Location : one location / many utc_tstamps ##
################################################
class Location(object):
    def __init__(self, lat, lon, height=None, turbidity_func=turbidity_remund, height_func=height_map, 
            am_func=AM_rozenberg):
        self.lat = lat
        self.lon = lon
        self.__height =height  
        self.turbidity_func = turbidity_func
        self.height_func = height_func
        self.am_func = am_func
    
    def get_turbidities(self, utc_tstamps):
        return  self.turbidity_func(utc_tstamps=utc_tstamps, lats=self.lat, lons=self.lon)

    def get_height(self):
        if type(self.__height) != type(None):
            return self.__height
        return self.height_func(lats=self.lat, lons=self.lon)

    def get_solar_params(self, utc_tstamps):
        return local_solar_params(utc_tstamps=utc_tstamps, lats=self.lat, lons=self.lon)

    def get_airmass(self, utc_tstamps):
        return self.am_func(self.get_solar_params(utc_tstamps)['zenith'])

    def get_clearsky_dumortier(self, utc_tstamps):
        return clearsky_dumortier(
                zenith = self.get_solar_params(utc_tstamps)['zenith'],
                eccentricity_factor = solartime(utc_tstamps)['eccentricity_factor'],
                turbidity = self.get_turbidities(utc_tstamps),
                height = self.get_height()
        )

###########################################
## MAP : many locations / one utc_tstamp ##
###########################################
class Map(object):
    def __init__(self, lats, lons, heights=None, turbidity_func=turbidity_remund, height_func=height_map, 
            am_func=AM_rozenberg):
        self.lats = lats
        self.lons = lons
        self.__heights =heights  
        self.turbidity_func = turbidity_func
        self.height_func = height_func
        self.am_func = am_func
    
    def get_turbidities(self, utc_tstamp):
        return  self.turbidity_func(utc_tstamps=utc_tstamp, lats=self.lats, lons=self.lons)

    def get_heights(self):
        if type(self.__heights) != type(None):
            return self.__heights
        return self.height_func(lats=self.lats, lons=self.lons)

    def get_solar_params(self, utc_tstamp):
        return local_solar_params(utc_tstamps=utc_tstamp, lats=self.lats, lons=self.lons)

    def get_airmass(self, utc_tstamp):
        return self.am_func(self.get_solar_params(utc_tstamp)['zenith'])

    def get_clearsky_dumortier(self, utc_tstamp):
        return clearsky_dumortier(
                zenith = self.get_solar_params(utc_tstamp)['zenith'],
                eccentricity_factor = solartime(utc_tstamp)['eccentricity_factor'],
                turbidity = self.get_turbidities(utc_tstamp),
                height = self.get_heights()
        )

