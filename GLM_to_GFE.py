import sys
import os
import numpy as np
from astropy.coordinates import EarthLocation, AltAz, ICRS
from astropy.table import Table
from astropy.time import Time
import astropy.units as u
import pandas as pd


# very rough calculation based on pixel size (8km)
GLM_angular_uncertainty = np.rad2deg(np.sin(8/36e3))/2

# Earth radii
R_eq = 6378e3
R_pol = 6357e3
wgs84_flatenning = 1/298.257223563

lightning_height_eq = 16e3
lightning_height_pol = 6e3

lightning_eq_R = R_eq + lightning_height_eq

lightning_ellipsoid_flatenning = (lightning_eq_R - R_pol - lightning_height_pol) / lightning_eq_R

def ellipsoid_radius(lat, r=R_eq, f=wgs84_flatenning):
    # from https://au.mathworks.com/help/aeroblks/radiusatgeocentriclatitude.html
    return (r**2 / (1+(1/(1-f)**2-1)*np.sin(np.deg2rad(lat))**2))**0.5

def lightning_height_above_wgs84(lat):
    '''
    determine top of cloud height (aka. lightning ellipsoid) for a particular latitude
    '''
    return ellipsoid_radius(lat, r=lightning_eq_R, f=lightning_ellipsoid_flatenning) - ellipsoid_radius(lat)


def parse_header(fname):
    '''
    parse useful field in CSV header of GLM file
    '''
    header = {}
    #Satellite (lat/lon/height):,0.0 / -75.199997 / 35786.02km
    with open(fname) as f:
        head = [next(f).replace("\n", "") for x in range(9)]
    sat_loc_line = head[4].replace(' ','')
    items = re.split(',|/', sat_loc_line)
    # ['#Satellite(lat', 'lon', 'height):', '0.0', '-75.199997', '35786.02km']
    header['satellite_latitude'] = float(items[3])
    header['satellite_longitude'] = float(items[4])
    header['satellite_altitude'] = float(items[5].replace('km', ''))
    
    
    '#Platform:,G16',
    sat_platform_line = head[3].replace(' ','')
    
    header['satellite_platform'] = re.split(',', sat_platform_line)[1]
    
    return header

def ECEF2ENU(lon, lat): # slightly inaccurate -> should use enu_matrix(obs_LLH, t_jd)?
    """
    # convert to local ENU coords
    # http://www.navipedia.net/index.php/Transformations_between_ECEF_and_ENU_coordinates
    # Title     Transformations between ECEF and ENU coordinates
    # Author(s)   J. Sanz Subirana, J.M. Juan Zornoza and M. Hernandez-Pajares, Technical University of Catalonia, Spain.
    # Year of Publication     2011 

    # use long to make greenwich mean turn to meridian: A clockwise rotation over the z-axis by and angle to align the east-axis with the x-axis
    # use lat to rotate z to zenith
    """
    
    res = np.array([[-np.sin(lon)         , np.cos(lon)            , 0], 
        [-np.cos(lon)*np.sin(lat) , -np.sin(lon) * np.sin(lat), np.cos(lat)],
        [np.cos(lon) * np.cos(lat), np.sin(lon) * np.cos(lat) , np.sin(lat)]])

    return res


def alt_az(sat_lat, sat_lon, sat_alt,
            obj_lat, obj_lon):

    # following https://gssc.esa.int/navipedia//index.php/Transformations_between_ECEF_and_ENU_coordinates
    
    pos_sat = EarthLocation(lat=sat_lat*u.deg, lon=sat_lon*u.deg, height=sat_alt * u.m)
    pos_sat_XYZ = np.array([pos_sat.x.value, pos_sat.y.value, pos_sat.z.value])
    
    # determine assumed height of lightning
    obj_alt = lightning_height_above_wgs84(obj_lat)
    
    pos_obj = EarthLocation(lat=obj_lat*u.deg, lon=obj_lon*u.deg, height=obj_alt * u.m)
    pos_obj_XYZ = np.array([pos_obj.x.value, pos_obj.y.value, pos_obj.z.value])

    #  line of sight unit vector 
    p = (pos_obj_XYZ - pos_sat_XYZ) / np.linalg.norm(pos_obj_XYZ - pos_sat_XYZ)

    enu = ECEF2ENU(np.deg2rad(sat_lon), np.deg2rad(sat_lat))

    elevation = np.rad2deg(np.arcsin((np.dot(p, enu[2]))))
    azimuth = np.rad2deg(np.arctan(np.dot(p, enu[0])/np.dot(p, enu[1]))) % 360

    return azimuth, elevation



def main(ifile):
    
    header = parse_header(ifile)
    t = Table.from_pandas(pd.read_csv(ifile, skiprows=10))
    
    t.meta['obs_latitude'] = header['satellite_latitude']
    t.meta['obs_longitude'] =  header['satellite_longitude']
    t.meta['obs_elevation'] = header['satellite_altitude']*1e3
    t.meta['location'] = 'Geostationary'
    t.meta['telescope'] = header['satellite_platform']
    t.meta['camera_id'] = header['satellite_platform']
    t.meta['instrument'] = 'GLM'
    t.meta['cx'] = 1372
    t.meta['cy'] = 1300
    t.meta['obs_az'] = 0.
    t.meta['obs_ev'] = -90.
    t.meta['obs_rot'] = 0.
    t.meta['fov_horiz'] = 8.
    t.meta['fov_vert'] = 8.
    t.meta['photometric_band'] = '777.4 nm'

    
    
    # Datetime is UNIX millisecond
    datetime_array = Time(t['time (ms)']/1e3, format='unix')
    t['datetime'] = datetime_array.isot
    
    # convert navigated lat/longs to alt-az
    t['azimuth'] = np.nan
    t['altitude'] = np.nan
    for i in range(len(t)):
        azimuth, elevation = alt_az(t.meta['obs_latitude'], t.meta['obs_longitude'], t.meta['obs_elevation'],
            t['latitude'][i], t['longitude'][i])
        t['altitude'][i] = elevation
        t['azimuth'][i] = azimuth
        
    t['err_minus_altitude'] = GLM_angular_uncertainty
    t['err_plus_altitude'] = GLM_angular_uncertainty
    t['err_minus_azimuth'] = GLM_angular_uncertainty
    t['err_plus_azimuth'] = GLM_angular_uncertainty

        
    station = EarthLocation(lat=t.meta['obs_latitude']*u.deg,
                        lon=t.meta['obs_longitude']*u.deg,
                        height=t.meta['obs_elevation']*u.meter)
    

    # Convert horizontal to equatorial coordinates
    horizontal_data = AltAz(t['azimuth']*u.deg,t['altitude']*u.deg,obstime=datetime_array,location=station)
    equatorial_data = horizontal_data.transform_to(ICRS)
    
    ra = equatorial_data.ra.value
    dec = equatorial_data.dec.value
    
    t['ra'] = ra
    t['dec'] = dec
    
    
    # rename the energy column
    t['energy'] = t['energy (joules)'] *u.joule
    
    
    t.meta['file_standard'] = 'GFE_1.2'
    t.remove_columns(['time (ms)','longitude','latitude','energy (joules)'])
    
    o_fname = '_'.join([(t['datetime'][0][:19]).replace(':','_'),header['satellite_platform'], t.meta['location'].replace('_','')]) + '.ecsv'
    o_fname = os.path.join(os.path.dirname(ifile), o_fname)
    
    
    # write the file to disk
    t.write(o_fname, format='ascii.ecsv', delimiter=',', overwrite=False)
    print(f'table has been written to {o_fname}')


def test_data():
    # 6414f5447199482c98f554cb/6414f5447199482c98f554cb_OR_GLM-L2-LCFA_G16_s20230611020400_e20230611021000_c20230611021020.nc.csv
    sat_lat = 0.0
    sat_lon = -75.199997
    sat_alt = 35786.02e3

    obj_lat = 35.30725860595703
    obj_lon = -123.09264373779297

    print(alt_az(sat_lat, sat_lon, sat_alt,
                obj_lat, obj_lon))





if __name__ == "__main__":
    main(sys.argv[1])

