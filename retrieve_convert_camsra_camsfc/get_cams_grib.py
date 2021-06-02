#!/usr/bin/python
import os
import numpy as np
import jstyleson as json


def read_json_file( fn ):
    with open(fn,'r') as f:
        return json.loads(f.read())
    return None

def get_camsfc_date(start,area,grid,targetfile,pfx='',
                    levs=['sfc','ml','sfc_chem','sfc_met','sfc_rad'],
                    debug=False):
    """
    Request and download one day of data from cams reanalysis.
    
    Parameters:
    -----------
    start : datetime object
        Day to request
    area : list
        Coordinates box of [N, W, S, E], of latitude (N,S) and longitude (W,E) in degrees N and degrees E
    grid : list
        grid resolution for [latitude,longitude] in degrees, e.g. [0.25,0.25]
    targetfile : string
        path for output file
    pfx : string
        tag in the output filename
    levs : list
        list of key names in 'camsra_variables_streams.json'
    debug: bool
        switch for debug mode (no download) and normal mode, default False (downloads will be issued)
        
    
    """
    if debug:
        print(f"Checking {start:%Y-%m-%d} to download cams forecast tables: {levs}")
    else:
        from ecmwfapi import ECMWFDataServer
    
    source = {
                "dataset" : "cams_nrealtime",
                "class"   : "mc",
                "expver"  : "1",
                "stream"  : "oper",
                "type"    : "fc",
             }
    product = 'cams-fc'
    param = "camsra_variables_streams.json"

    
    req = {'dataset': 'cams_nrealtime',
           'class': 'mc',
           'expver': '0001',
           'stream': 'oper',
           'type': 'fc',
           'date': f"{start:%Y-%m-%d}",
           'grid': '/'.join([str(a) for a in grid]),
           'area': '/'.join([str(a) for a in area]),
           'time': '00:00:00/12:00:00', # forcast from time ..
           'step' : '3/6/9/12', # time steps from forcast time 
           'levelist': '/'.join([str(a) for a in (np.arange(60)+1).astype('U2')])
          }
    
    # separate ml and sfc file
    req_param = read_json_file(param)
    for rq in req_param:
        if not rq['levtype'] in levs:
            continue
        target = targetfile.format(start=start,
                                   product=product,
                                   pfx=pfx,
                                   levtype=rq['levtype'])
        if os.path.exists(target):
            if debug:
                print("  >> Request already exists or already done. Skip!")
            continue
        else:
            # make a dummy file as place holde to avoid double requests
            os.system("touch {}".format(target))
            
            
        r = req.copy()
        if rq['levtype'] != 'ml':
            r.pop('levelist')
        r.update({'variable':'/'.join([str(a) for a in rq['param_ID']]),
                  'levtype':rq['levtype'].split('_')[0],
                  'target':target}) 

        if debug:
            print("  >> send request:")
            print(json.dumps(r,indent=4))
        else:
            client = ECMWFDataServer()
            client.retrieve(r)
    return 0

def get_camsra_date(start,area,grid,targetfile,pfx='',
                    levs=['sfc','ml','sfc_chem','sfc_met','sfc_rad'],
                    debug=False):
    """
    Request and download one day of data from cams reanalysis.
    
    Parameters:
    -----------
    start : datetime object
        Day to request
    area : list
        Coordinates box of [N, W, S, E], of latitude (N,S) and longitude (W,E) in degrees N and degrees E
    grid : list
        grid resolution for [latitude,longitude] in degrees, e.g. [0.25,0.25]
    targetfile : string
        path for output file
    pfx : string
        tag in the output filename
    levs : list
        list of key names in 'camsra_variables_streams.json'
    debug: bool
        switch for debug mode (no download) and normal mode, default False (downloads will be issued)
        
    
    """
    
    if debug:
        print(f"Checking {start:%Y-%m-%d} to download cams reanalysis tables: {levs}")
    else:
        import cdsapi
    
    
    product = 'cams-ra'
    param = "camsra_variables_streams.json"
    source = "cams-global-reanalysis-eac4"

    req = dict( date = f"{start:%Y-%m-%d}",
                format = 'grib',
                model_level = list((np.arange(60)+1).astype('U2')),# model level 1 to 60
                time = ['00:00','03:00','06:00','09:00',
                        '12:00','15:00','18:00','21:00'],
                area = area,
                grid = grid)

    
    # separate ml and sfc file
    req_param = read_json_file(param)
    for rq in req_param:
        print(rq)
        if not rq['levtype'] in levs:
            continue
        target = targetfile.format(start=start,
                                   product=product,
                                   pfx=pfx,
                                   levtype=rq['levtype'])
        if os.path.exists(target):
            if debug:
                print("  >> Request already exists or already done. Skip!")
            continue
        else:
            # make a dummy file as place holde to avoid double requests
            os.system("touch {}".format(target))
        r = req.copy()
        if rq['levtype'] != 'ml':
            r.pop('model_level')
        r.update({'variable':rq['param_name']})
        if debug:
            print("  >> send request:")
            print(json.dumps(r,indent=4))
            print(f"  >> target: {target}")
        else:
            client = cdsapi.Client()
            client.retrieve(source,r,target)
    return 0
