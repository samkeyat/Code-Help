# solar-python — Solar data calculation and prediction library for Python
# Copyright © 2014, 2015  Mattias Andrée (maandree@member.fsf.org)
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.



SOLAR_APPARENT_RADIUS = 32.0 / 60.0
'''
:float  Approximate apparent size of the Sun in degrees
'''


SOLAR_ELEVATION_PRESUNSET_POSTSUNRISE = 32.0 / 60.0
'''
:float  The Sun's elevation the beginning of sunset and
        end of sunrise, measured in degrees
'''

SOLAR_ELEVATION_SUNSET_SUNRISE = -32.0 / 60.0
'''
:float  The Sun's elevation the (end of) at sunset and
        (beginning of) sunrise, measured in degrees
'''

SOLAR_ELEVATION_CIVIL_DUSK_DAWN = -6.0
'''
:float  The Sun's elevation at civil dusk and civil
        dawn, measured in degrees
'''

SOLAR_ELEVATION_NAUTICAL_DUSK_DAWN = -12.0
'''
:float  The Sun's elevation at nautical dusk and
        nautical dawn, measured in degrees
'''

SOLAR_ELEVATION_AMATEUR_ASTRONOMICAL_DUSK_DAWN = -15.0
'''
:float  The Sun's elevation at amateur astronomical dusk
        and amateur astronomical dawn, measured in degrees
'''

SOLAR_ELEVATION_ASTRONOMICAL_DUSK_DAWN = -18.0
'''
:float  The Sun's elevation at astronomical dusk
        and astronomical dawn, measured in degrees
'''

SOLAR_ELEVATION_RANGE_TWILIGHT = (-18.0, 0.0)
'''
:(float, float)  The Sun's lowest and highest elevation during
                 all periods of twilight, measured in degrees
'''

SOLAR_ELEVATION_RANGE_CIVIL_TWILIGHT = (-6.0, SOLAR_ELEVATION_SUNSET_SUNRISE)
'''
:(float, float)  The Sun's lowest and highest elevation
                 during civil twilight, measured in degrees
'''

SOLAR_ELEVATION_RANGE_NAUTICAL_TWILIGHT = (-12.0, SOLAR_ELEVATION_SUNSET_SUNRISE)
'''
:(float, float)  The Sun's lowest and highest elevation
                 during nautical twilight, measured in degrees
'''

SOLAR_ELEVATION_RANGE_ASTRONOMICAL_TWILIGHT = (-18.0, SOLAR_ELEVATION_SUNSET_SUNRISE)
'''
:(float, float)  The Sun's lowest and highest elevation during
                 astronomical twilight, measured in degrees
'''

SOLAR_ELEVATION_RANGE_AMATEUR_ASTRONOMICAL_TWILIGHT = (-18.0, -15.0)
'''
:(float, float)  The Sun's lowest and highest elevation during
                 amateur astronomical twilight, measured in degrees
'''

SOLAR_ELEVATION_RANGE_GOLDEN_HOUR = (-4.0, 6.0)
'''
:(float, float)  The Sun's lowest and highest elevation during
                 the golden "hour" (also known as magic hour),
                 measured in degrees. These elevations are
                 approximate.
'''

SOLAR_ELEVATION_RANGE_BLUE_HOUR = (-6.0, -4.0)
'''
:(float, float)  The Sun's lowest and highest elevation during
                 the blue "hour", measured in degrees. These
                 elevations are approximate.
'''



# The following functions are used to calculate the result for `sun`
# (most of them) but could be used for anything else. There name is
# should tell you enough, `t` (and `noon`) is in Julian Centuries
# except for in the convertion methods.


def julian_day_to_epoch(t):
    '''
    Converts a Julian Day timestamp to a POSIX time timestamp
    
    @param   t:float  The time in Julian Days
    @return  :float   The time in POSIX time
    '''
    return (t - 2440587.5) * 86400.0


def epoch_to_julian_day(t):
    '''
    Converts a POSIX time timestamp to a Julian Day timestamp
    
    @param   t:float  The time in POSIX time
    @return  :float   The time in Julian Days
    '''
    return t / 86400.0 + 2440587.5


def julian_day_to_julian_centuries(t):
    '''
    Converts a Julian Day timestamp to a Julian Centuries timestamp
    
    @param   t:float  The time in Julian Days
    @return  :float   The time in Julian Centuries
    '''
    return (t - 2451545.0) / 36525.0


def julian_centuries_to_julian_day(t):
    '''
    Converts a Julian Centuries timestamp to a Julian Day timestamp
    
    @param   t:float  The time in Julian Centuries
    @return  :float   The time in Julian Days
    '''
    return t * 36525.0 + 2451545.0


def epoch_to_julian_centuries(t):
    '''
    Converts a POSIX time timestamp to a Julian Centuries timestamp
    
    @param   t:float  The time in POSIX time
    @return  :float   The time in Julian Centuries
    '''
    return julian_day_to_julian_centuries(epoch_to_julian_day(t))


def julian_centuries_to_epoch(t):
    '''
    Converts a Julian Centuries timestamp to a POSIX time timestamp
    
    @param   t:float  The time in Julian Centuries
    @return  :float   The time in POSIX time
    '''
    return julian_day_to_epoch(julian_centuries_to_julian_day(t))


def epoch():
    '''
    Get current POSIX time
    
    @return  :float  The current POSIX time
    '''
    import time
    return time.time()


def julian_day():
    '''
    Get current Julian Day time
    
    @return  :float  The current Julian Day time
    '''
    return epoch_to_julian_day(epoch())


def julian_centuries():
    '''
    Get current Julian Centuries time (100 Julian days since J2000)
    
    @return  :float  The current Julian Centuries time
    '''
    return epoch_to_julian_centuries(epoch())


def radians(deg):
    '''
    Convert an angle from degrees to radians
    
    @param   deg:float  The angle in degrees
    @return  :float     The angle in radians
    '''
    import math
    return deg * math.pi / 180


def degrees(rad):
    '''
    Convert an angle from radians to degrees
    
    @param   rad:float  The angle in radians
    @return  :float     The angle in degrees
    '''
    import math
    return rad * 180 / math.pi


def sun_geometric_mean_longitude(t):
    '''
    Calculates the Sun's geometric mean longitude
    
    @param   t:float  The time in Julian Centuries
    @return  :float   The Sun's geometric mean longitude in radians
    '''
    return radians((0.0003032 * t ** 2 + 36000.76983 * t + 280.46646) % 360)
    # CANNIBALISERS:
    #     The result of this function should always be positive, this
    #     means that after division modulo 360 but before `radians`,
    #     you will need to add 360 if the value is negative. This can
    #     only happen if `t` is negative, which can only happen for date
    #     times before 2000-(01)Jan-01 12:00:00 UTC par division modulo
    #     implementations with the signess of at least the left operand.
    #     More precively, it happens between circa 1970-(01)Jan-11
    #     16:09:02 UTC and circa 374702470660351740 seconds before
    #     January 1, 1970 00:00 UTC, which is so far back in time
    #     it cannot be reliable pinned down to the right year, but it
    #     is without a shadow of a doubt looooong before the Earth
    #     was formed, is right up there with the age of the Milky Way
    #     and the universe itself.


def sun_geometric_mean_anomaly(t):
    '''
    Calculates the Sun's geometric mean anomaly
    
    @param   t:float  The time in Julian Centuries
    @return  :float   The Sun's geometric mean anomaly in radians
    '''
    return radians(-0.0001537 * t ** 2 + 35999.05029 * t + 357.52911)


def earth_orbit_eccentricity(t):
    '''
    Calculates the Earth's orbit eccentricity
    
    @param   t:float  The time in Julian Centuries
    @return  :float   The Earth's orbit eccentricity
    '''
    return -0.0000001267 * t ** 2 - 0.000042037 * t + 0.016708634


def sun_equation_of_centre(t):
    '''
    Calculates the Sun's equation of the centre, the difference between
    the true anomaly and the mean anomaly
    
    @param   t:float  The time in Julian Centuries
    @return  :float   The Sun's equation of the centre, in radians
    '''
    import math
    a = sun_geometric_mean_anomaly(t)
    rc = math.sin(1 * a) * (-0.000014 * t ** 2 - 0.004817 * t + 1.914602)
    rc += math.sin(2 * a) * (-0.000101 * t + 0.019993)
    rc += math.sin(3 * a) * 0.000289
    return radians(rc)


def sun_real_longitude(t):
    '''
    Calculates the Sun's real longitudinal position
    
    @param   t:float  The time in Julian Centuries
    @return  :float   The longitude, in radians
    '''
    rc = sun_geometric_mean_longitude(t)
    return rc + sun_equation_of_centre(t)


def sun_apparent_longitude(t):
    '''
    Calculates the Sun's apparent longitudinal position
    
    @param   t:float  The time in Julian Centuries
    @return  :float   The longitude, in radians
    '''
    import math
    rc = degrees(sun_real_longitude(t)) - 0.00569
    rc -= 0.00478 * math.sin(radians(-1934.136 * t + 125.04))
    return radians(rc)


def mean_ecliptic_obliquity(t):
    '''
    Calculates the mean ecliptic obliquity of the Sun's
    apparent motion without variation correction
    
    @param   t:float  The time in Julian Centuries
    @return  :float   The uncorrected mean obliquity, in radians
    '''
    rc = 0.001813 * t ** 3 - 0.00059 * t ** 2 - 46.815 * t + 21.448
    rc = 26 + rc / 60
    rc = 23 + rc / 60
    return radians(rc)


def corrected_mean_ecliptic_obliquity(t):
    '''
    Calculates the mean ecliptic obliquity of the Sun's
    apparent motion with variation correction
    
    @param   t:float  The time in Julian Centuries
    @return  :float   The mean obliquity, in radians
    '''
    import math
    rc = -1934.136 * t + 125.04
    rc = 0.00256 * math.cos(radians(rc))
    rc += degrees(mean_ecliptic_obliquity(t))
    return radians(rc)


def solar_declination(t):
    '''
    Calculates the Sun's declination
    
    @param   t:float  The time in Julian Centuries
    @return  :float   The Sun's declination, in radians
    '''
    import math
    rc = math.sin(corrected_mean_ecliptic_obliquity(t))
    rc *= math.sin(sun_apparent_longitude(t))
    return math.asin(rc)


def equation_of_time(t):
    '''
    Calculates the equation of time, the discrepancy
    between apparent and mean solar time
    
    @param   t:float  The time in Julian Centuries
    @return  :float   The equation of time, in degrees
    '''
    import math
    l = sun_geometric_mean_longitude(t)
    e = earth_orbit_eccentricity(t)
    m = sun_geometric_mean_anomaly(t)
    y = corrected_mean_ecliptic_obliquity(t)
    y = math.tan(y / 2) ** 2
    rc = y * math.sin(2 * l)
    rc += (4 * y * math.cos(2 * l) - 2) * e * math.sin(m)
    rc -= 0.5 * y ** 2 * math.sin(4 * l)
    rc -= 1.25 * e ** 2 * math.sin(2 * m)
    return 4 * degrees(rc)


def hour_angle_from_elevation(latitude, declination, elevation):
    '''
    Calculates the solar hour angle from the Sun's elevation
    
    @param   longitude:float    The longitude in degrees eastwards
                                from Greenwich, negative for westwards
    @param   declination:float  The declination, in radians
    @param   elevation:float    The Sun's elevation, in radians
    @return  :float             The solar hour angle, in radians
    '''
    import math
    if elevation == 0:
        return 0
    rc = math.cos(abs(elevation))
    rc -= math.sin(radians(latitude)) * math.sin(declination)
    rc /= math.cos(radians(latitude)) * math.cos(declination)
    rc = math.acos(rc)
    return -rc if (rc < 0) == (elevation < 0) else rc


def elevation_from_hour_angle(latitude, declination, hour_angle):
    '''
    Calculates the Sun's elevation from the solar hour angle
    
    @param   longitude:float    The longitude in degrees eastwards
                                from Greenwich, negative for westwards
    @param   declination:float  The declination, in radians
    @param   hour_angle:float   The solar hour angle, in radians
    @return  :float             The Sun's elevation, in radians
    '''
    import math
    rc = math.cos(radians(latitude))
    rc *= math.cos(hour_angle) * math.cos(declination)
    rc += math.sin(radians(latitude)) * math.sin(declination)
    return math.asin(rc)


def time_of_solar_noon(t, longitude):
    '''
    Calculates the time of the closest solar noon
    
    @param   t:float          A time close to the sought time,
                              in Julian Centuries
    @param   longitude:float  The longitude in degrees eastwards from
                              Greenwich, negative for westwards
    @return  :float           The time, in Julian Centuries,
                              of the closest solar noon
    '''
    t, rc = julian_centuries_to_julian_day(t), longitude
    for (k, m) in ((-360, 0), (1440, -0.5)):
        rc = julian_day_to_julian_centuries(t + m + rc / k)
        rc = 720 - 4 * longitude - equation_of_time(rc)
    return rc


def time_of_solar_elevation(t, noon, latitude, longitude, elevation):
    '''
    Calculates the time the Sun has a specified apparent
    elevation at a geographical position
    
    @param   t:float          A time close to the sought time,
                              in Julian Centuries
    @param   noon:float       The time of the closest solar noon
    @param   latitude:float   The latitude in degrees northwards from
                              the equator, negative for southwards
    @param   longitude:float  The longitude in degrees eastwards from
                              Greenwich, negative for westwards
    @param   elevation:float  The solar elevation, in radians
    @return  :float           The time, in Julian Centuries,
                              of the specified elevation
    '''
    rc = noon
    rc, et = solar_declination(rc), equation_of_time(rc)
    rc = hour_angle_from_elevation(latitude, rc, elevation)
    rc = 720 - 4 * (longitude + degrees(rc)) - et
    
    rc = julian_day_to_julian_centuries(julian_centuries_to_julian_day(t) + rc / 1440)
    rc, et = solar_declination(rc), equation_of_time(rc)
    rc = hour_angle_from_elevation(latitude, rc, elevation)
    rc = 720 - 4 * (longitude + degrees(rc)) - et
    return rc


def solar_elevation_from_time(t, latitude, longitude):
    '''
    Calculates the Sun's elevation as apparent
    from a geographical position
    
    @param   t:float          The time in Julian Centuries
    @param   latitude:float   The latitude in degrees northwards from
                              the equator, negative for southwards
    @param   longitude:float  The longitude in degrees eastwards from
                              Greenwich, negative for westwards
    @return  :float           The Sun's apparent elevation at the specified
                              time as seen from the specified position,
                              measured in radians
    '''
    rc = julian_centuries_to_julian_day(t)
    rc = (rc - float(int(rc + 0.5)) - 0.5) * 1440
    rc = 720 - rc - equation_of_time(t)
    rc = radians(rc / 4 - longitude)
    return elevation_from_hour_angle(latitude, solar_declination(t), rc)


def solar_elevation(latitude, longitude, t = None):
    '''
    Calculates the Sun's elevation as apparent
    from a geographical position
    
    @param   latitude:float   The latitude in degrees northwards from
                              the equator, negative for southwards
    @param   longitude:float  The longitude in degrees eastwards from
                              Greenwich, negative for westwards
    @param   t:float?         The time in Julian Centuries, `None`
                              for the current time
    @return  :float           The Sun's apparent elevation at the specified
                              time as seen from the specified position,
                              measured in degrees
    '''
    rc = julian_centuries() if t is None else t
    rc = solar_elevation_from_time(rc, latitude, longitude)
    return degrees(rc)



def have_sunrise_and_sunset(latitude, t = None):
    '''
    Determine whether solar declination currently is
    so that there can be sunrises and sunsets. If not,
    you either have 24-hour daytime or 24-hour nighttime.
    
    @param   latitude:float  The latitude in degrees northwards from
                             the equator, negative for southwards
    @param   t:float?        The time in Julian Centuries, `None`
                             for the current time
    @return                  Whether there can be sunrises and
                             sunsets where you are located
    '''
    t = julian_centuries() if t is None else t
    d = degrees(solar_declination(t))
    ## Convert everything to the Northern hemisphere
    latitude = abs(latitude)
    if d >= 0:
        ## Northern summer
        return -90 + d < latitude < 90 - d
    else:
        ## Northern winter
        return -90 - d < latitude < 90 + d


def is_summer(latitude, t = None):
    '''
    Determine whether it is summer
    
    @param   latitude:float  The latitude in degrees northwards from
                             the equator, negative for southwards
    @param   t:float?        The time in Julian Centuries, `None`
                             for the current time
    @return                  Whether it is summer on the hemisphere
                             you are located on
    '''
    t = julian_centuries() if t is None else t
    d = solar_declination(t)
    return (d > 0) == (latitude > 0)


def is_winter(latitude, t = None):
    '''
    Determine whether it is winter
    
    @param   latitude:float  The latitude in degrees northwards from
                             the equator, negative for southwards
    @param   t:float?        The time in Julian Centuries, `None`
                             for the current time
    @return                  Whether it is winter on the hemisphere
                             you are located on
    '''
    t = julian_centuries() if t is None else t
    d = solar_declination(t)
    return not ((d > 0) == (latitude > 0))



def solar_prediction(delta, requested, fun, epsilon = 0.000001, span = 0.01, t = None):
    '''
    Predict the time point of the next or previous
    time an arbitrary condition is meet
    
    @param   delta:float          Iteration step size, negative for past
                                  event, positive for future event
    @param   requested:float      The value returned by `fun` for which to
                                  calculate the time point of occurrence
    @param   fun:(t:float)→float  Function that calculate the data of interest
    @param   epsilon:float        The tolerance for the result
    @param   span:float           The number of Julian Centuries (0,01 for
                                  one year) to restrict the search to
    @param   t:float?             The time in Julian Centuries, `None` for
                                  the current time
    @return  :float?              The calculated time point, `None` if none
                                  were found within the specified time span
    '''
    fun_ = lambda t : fun(t) - requested
    t = julian_centuries() if t is None else t
    t1 = t2 = t
    v1 = v0 = fun_(t)
    
    # Predicate time point to within a small time span
    while True:
        if abs(t2 - t) > span:
            return None
        t2 += delta
        v2 = fun_(t2)
        if (v1 <= 0 <= v2) or ((0 >= v1 >= v2) and (0 <= v0)):
            break
        if (v1 >= 0 >= v2) or ((0 <= v1 <= v2) and (0 >= v0)):
            break
        t1 = t2
        v2 = v1
    
    # Binary search the small time span for the exact time point
    for _itr in range(1000):
        tm = (t1 + t2) / 2
        v1 = fun_(t1)
        v2 = fun_(t2)
        vm = fun_(tm)
        if abs(v1 - v2) < epsilon:
            return tm
        if v1 < v2:
            if 0 < vm:
                t2 = tm
            else:
                t1 = tm
        elif v1 > v2:
            if 0 > vm:
                t2 = tm
            else:
                t1 = tm
    return tm



def future_past_equinox(delta, t = None):
    '''
    Predict the time point of the next or previous equinox
    
    @param   delta:float  Iteration step size, negative for
                          past event, positive for future event
    @param   t:float?     The time in Julian Centuries, `None`
                          for the current time
    @return  :float       The calculated time point
    '''
    return solar_prediction(delta, 0, solar_declination, t = t)


def future_equinox(t = None):
    '''
    Predict the time point of the next equinox
    
    @param   t:float?  The time in Julian Centuries, `None`
                       for the current time
    @return  :float    The calculated time point
    '''
    return future_past_equinox(0.01 / 2000, t)
    

def past_equinox(t = None):
    '''
    Predict the time point of the previous equinox
    
    @param   t:float?  The time in Julian Centuries, `None`
                       for the current time
    @return  :float    The calculated time point
    '''
    return future_past_equinox(0.01 / -2000, t)



def future_past_solstice(delta, t = None):
    '''
    Predict the time point of the next or previous solstice
    
    @param   delta:float  Iteration step size, negative for
                          past event, positive for future event
    @param   t:float?     The time in Julian Centuries, `None`
                          for the current time
    @return  :float       The calculated time point
    '''
    e = 0.00001
    fun = solar_declination
    dfun = lambda t : (fun(t + e) - fun(t - e)) / (2 * e)
    return solar_prediction(delta, 0, dfun, t = t)


def future_solstice(t = None):
    '''
    Predict the time point of the next solstice
    
    @param   t:float?  The time in Julian Centuries,
                       `None` for the current time
    @return  :float    The calculated time point
    '''
    return future_past_solstice(0.01 / 2000, t)
    

def past_solstice(t = None):
    '''
    Predict the time point of the previous solstice
    
    @param   t:float?  The time in Julian Centuries,
                       `None` for the current time
    @return  :float    The calculated time point
    '''
    return future_past_solstice(0.01 / -2000, t)



def future_past_elevation(delta, latitude, longitude, elevation, t = None):
    '''
    Predict the time point of the next or previous time
    the Sun reaches or reached a specific elevation
    
    @param   delta:float      Iteration step size, negative for past
                              event, positive for future event
    @param   latitude:float   The latitude in degrees northwards from
                              the equator, negative for southwards
    @param   longitude:float  The longitude in degrees eastwards from
                              Greenwich, negative for westwards
    @param   elevation:float  The elevation of interest
    @param   t:float?         The time in Julian Centuries, `None`
                              for the current time
    @return  :float?          The calculated time point, `None` if
                              none were found within a year
    '''
    fun = lambda t : solar_elevation(latitude, longitude, t)
    return solar_prediction(delta, elevation, fun, t = t)


def future_elevation(latitude, longitude, elevation, t = None):
    '''
    Predict the time point of the next time the Sun
    reaches a specific elevation
    
    @param   latitude:float   The latitude in degrees northwards from
                              the equator, negative for southwards
    @param   longitude:float  The longitude in degrees eastwards from
                              Greenwich, negative for westwards
    @param   elevation:float  The elevation of interest
    @param   t:float?         The time in Julian Centuries, `None`
                              for the current time
    @return  :float?          The calculated time point, `None` if
                              none were found within a year
    '''
    return future_past_elevation(0.01 / 2000, latitude, longitude, elevation, t)
    

def past_elevation(latitude, longitude, elevation, t = None):
    '''
    Predict the time point of the previous time the Sun
    reached a specific elevation
    
    @param   latitude:float   The latitude in degrees northwards from
                              the equator, negative for southwards
    @param   longitude:float  The longitude in degrees eastwards from
                              Greenwich, negative for westwards
    @param   elevation:float  The elevation of interest
    @param   t:float?         The time in Julian Centuries, `None`
                              for the current time
    @return  :float?          The calculated time point, `None` if
                              none were found within a year
    '''
    return future_past_elevation(0.01 / -2000, latitude, longitude, elevation, t)



def future_past_elevation_derivative(delta, latitude, longitude, derivative, t = None):
    '''
    Predict the time point of the next or previous time the
    Sun reaches or reached a specific elevation derivative
    
    @param   delta:float       Iteration step size, negative for past
                               event, positive for future event
    @param   latitude:float    The latitude in degrees northwards from
                               the equator, negative for southwards
    @param   longitude:float   The longitude in degrees eastwards from
                               Greenwich, negative for westwards
    @param   derivative:float  The elevation derivative value of interest
    @param   t:float?          The time in Julian Centuries, `None`
                               for the current time
    @return  :float?           The calculated time point, `None` if
                               none were found within a year
    '''
    e = 0.00001
    fun = lambda t : solar_elevation(latitude, longitude, t)
    dfun = lambda t : (fun(t + e) - fun(t - e)) / (2 * e)
    return solar_prediction(delta, derivative, dfun, t = t)


def future_elevation_derivative(latitude, longitude, derivative, t = None):
    '''
    Predict the time point of the next time the
    Sun reaches a specific elevation derivative
    
    @param   latitude:float    The latitude in degrees northwards from
                               the equator, negative for southwards
    @param   longitude:float   The longitude in degrees eastwards from
                               Greenwich, negative for westwards
    @param   derivative:float  The elevation derivative value of interest
    @param   t:float?          The time in Julian Centuries, `None`
                               for the current time
    @return  :float?           The calculated time point, `None` if
                               none were found within a year
    '''
    return future_past_elevation_derivative(0.01 / 2000, latitude, longitude, derivative, t)
    

def past_elevation_derivative(latitude, longitude, derivative, t = None):
    '''
    Predict the time point of the previous time
    the Sun reached a specific elevation derivative
    
    @param   latitude:float    The latitude in degrees northwards from
                               the equator, negative for southwards
    @param   longitude:float   The longitude in degrees eastwards from
                               Greenwich, negative for westwards
    @param   derivative:float  The elevation derivative value of interest
    @param   t:float?          The time in Julian Centuries, `None`
                               for the current time
    @return  :float?           The calculated time point, `None`
                               if none were found within a year
    '''
    return future_past_elevation_derivative(0.01 / -2000, latitude, longitude, derivative, t)



# TODO: This algorithm is imprecise, gives an incorrent sunrise and I do not fully know its behaviour
def sunrise_equation(latitude, longitude, t = None):
    import math
    # Calculate Julian Cycle
    j_cent = julian_centuries() if t is None else t
    j_date = julian_centuries_to_julian_day(j_cent)
    j_cycle = int(j_date - 2451545.0009 - longitude / 360 + 0.5)
    
    # Calculate approximate solar noon and solar man anomaly
    approx_solar_noon = 451545.0009 + longitude / 360 + j_cycle
    solar_mean_anomaly = int(357.5291 + 0.98560028 * (j_cycle - 2451545)) % 360
    
    # Calculate solar equation of centre
    equation_of_centre  = 1.9148 * math.sin(1 * solar_mean_anomaly)
    equation_of_centre += 0.0200 * math.sin(2 * solar_mean_anomaly)
    equation_of_centre += 0.0003 * math.sin(3 * solar_mean_anomaly)
    
    # Calculate solar ecliptic longitude
    ecliptic_longitude = (solar_mean_anomaly + 102.9372 + equation_of_centre + 180) % 360
    
    # Calculate solar transit
    solar_transit  = approx_solar_noon + 0.0053 * math.sin(solar_mean_anomaly)
    solar_transit -= 0.0069 * math.sin(2 * ecliptic_longitude)
    
    # Calculate solar declination
    declination = math.asin(math.sin(ecliptic_longitude) * math.sin(radians(23.45)))
    
    # Calculate solar hour angle
    hour_angle  = math.sin(radians(-0.83))
    hour_angle -= math.sin(latitude) * math.sin(declination)
    hour_angle /= math.cos(latitude) * math.cos(declination)
    hour_angle = degrees(math.acos(hour_angle))
    
    # Calculate time of sunset and sunrise
    sunset  = 2451545.0009 + (hour_angle + longitude) / 360
    sunset += j_cycle + solar_transit - approx_solar_noon
    sunrise = 2 * solar_transit - sunset
    
    # Convert to Julian Centuries
    return (julian_day_to_julian_centuries(sunset),
            julian_day_to_julian_centuries(sunrise))

