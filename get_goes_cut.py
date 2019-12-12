
import datetime as dt, glob, os, argparse

import numpy as np, h5py

from satlib_py.satellites import goes16

def ignore(*args, **kwargs):
    pass

def parse_date(s):
    return dt.datetime.strptime(s,'%Y%m%d%H%M') 

def parse_domain(s):
    arr = list(map(float, s.split(',')))
    return np.array([
            [arr[0], arr[1]], 
            [arr[2], arr[3]]
    ])


class Section_from_lonlat(object):
    def __init__(self, datapath, lonlat_domain, sec_name, outpath, log=ignore):
        self.datapath = datapath
        self.lonlat_domain = lonlat_domain
        self.sec_name = sec_name
        self.dst = outpath
        self.log = log

    def get_cut(self, datetime):
        self.datetime = datetime
        path2file = glob.glob(self.datapath.format(date=self.datetime))[0]
        g16 = goes16.G16nc(path2file, log=self.log)
        g16.load_dataset()
        rad = g16.read_var('Rad')
        x = g16.read_var('x')
        y = g16.read_var('y')
        xy_domain = goes16.lonlat2xy(self.lonlat_domain)
        idx = np.where((x > xy_domain[0][0]) & (x < xy_domain[1][0]))[0]
        idy = np.where((y < xy_domain[0][1]) & (y > xy_domain[1][1]))[0]
        cut = rad[idy[0]:idy[-1]+1, idx[0]:idx[-1]+1]
        self.store_h5(cut, x[idx], y[idy])
#        goes16.rad2tiff_deflate(cut, 
#            os.path.join(self.dst, '{date:%Y%m%d%H%M}-{name}.tif'.format(date=self.datetime, name=self.sec_name)))
        return cut, x[idx], y[idy]

    def store_h5(self, cut, x, y):
        xinfo = [x[0], x[1]-x[0], len(x)]
        yinfo = [y[0], y[1]-y[0], len(y)]
        F = h5py.File(
                os.path.join(self.dst, '{date:%Y%m%d%H%M}-{name}.h5'.format(date=self.datetime, name=self.sec_name)),
                'w')
        ds = F.create_dataset('counts',data=cut, dtype = np.float16)
        ds.attrs['datetime']='{date:%Y%m%d%H%M}'.format(date=self.datetime)
        ds.attrs['x_0_step_len'] = xinfo
        ds.attrs['y_0_step_len'] = yinfo
        ds.attrs['sat_name'] = 'goes16'
        ds.attrs['calibration'] = 'radiance'
        F.close()
    
    __call__ = get_cut


def call():
    parser = argparse.ArgumentParser(prog='PROG',
            description='cut lonlat domain in goes16 image')
    parser.add_argument(
            '-d', '--date', type=parse_date, required=True, 
            help='datetime in the form of yyyymmddHHMM')
    parser.add_argument(
            '-p', '--path', type=str,
            default = '/sat/san/raw/goesreast/conus/netcdf/{date:%Y}/{date:%m}/{date:%d}/OR_ABI-L1b-RadF-M3C02_G16_s{date:%Y%j%H%M}*nc',
            help='path template to netCDF file. Accepts "{date:<format>}"')
    parser.add_argument(
            '-m', '--domain', action='store',  type=parse_domain, required=True, 
            help='domain to cut in the form: "lon1,lat1,lon2,lat2". Point1 is the top-left corner and point2 is the rigth-bottom corner')
    parser.add_argument(
            '-n', '--sec_name', type=str, 
            help='name to be given to the cutted section')
    parser.add_argument(
            '-o', '--out_dir',  type=str, required=True, 
            help='make dryrun (default is False)')
    args = parser.parse_args()
    return Section_from_lonlat(args.path, args.domain, args.sec_name, args.out_dir)(args.date)

if __name__ == '__main__':
    call()
#    datetime = dt.datetime(2019,1,22,19,0)
#    path_fmt = '/sat/san/raw/goesreast/conus/netcdf/{date:%Y}/{date:%m}/{date:%d}/OR_ABI-L1b-RadF-M3C02_G16_s{date:%Y%j%H%M}*nc'
#    lonlat_domain_ecuador = np.array([[-82.705078, 2.195996],[-72.861328, -7.101619]])
#    cut, x, y = Section_from_lonlat(path_fmt, lonlat_domain_ecuador, 'ecuador', '/tmp/')(datetime)



