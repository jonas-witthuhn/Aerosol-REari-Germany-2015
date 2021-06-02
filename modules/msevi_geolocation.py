
from pyproj import Proj


lres = 3000.40316582
hres = 1000.1343886066667

proj_param = { 
    'proj':  'geos',
    'h':     35785831.0, 
    'a':      6378169.0, 
    'b':      6356583.8, 
    'lon_0':        0.0
}

def latlon2xy(lat, lon, res=lres):
    ''' 
    Convert geographic coordinates lat/lon to
    satellite pixel coordinates
    '''
    proj = Proj(**proj_param)
    (x,y) = proj(lon,lat)
    x = 1856.0+(x/res) #1856
    y = 1856.0-(y/res)
    return (x,y)

def xy2latlon(x,y,res=lres):
    ''' 
    Convert satellite pixel coordinates
    to geographic coordinates
    '''    
    proj = Proj(**proj_param)
    x = (x-1856.0)*res
    y = (1856.0-y)*res
    (lon,lat) = proj(x,y,inverse=True)
    return (lat,lon)
