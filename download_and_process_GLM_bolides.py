#!/usr/bin/env python

import os
import requests
import json
import wget
import datetime

import pandas as pd


from astropy.table import Table


DATA_DIR = 'glm_bolides_events'

URL = "https://neo-bolide.ndc.nasa.gov/service/event/public"
BASE_CSV_URL = 'https://neo-bolide.ndc.nasa.gov/csv/'



def process_attachements(attachments, basedir, overwrite=False):
    for att in attachments:
        identifier = att['_id']
        dirname = os.path.join(basedir, str(identifier))
        os.makedirs(dirname, exist_ok=True)
        netCdf_filename = att['netCdfFilename']
        csv_filename = netCdf_filename + '.csv'
        dest_csv = os.path.join(dirname, csv_filename)
        if os.path.isfile(dest_csv) and not overwrite:
            print('WARNING: {} exists, not downloading again'.format(dest_csv))
            return
        csv_url = BASE_CSV_URL + csv_filename
        try:
            wget.download(csv_url, out=dest_csv)
        except Exception as e:
            print(e)
            print(f'Could not download {csv_url}')

def main():
    os.makedirs(DATA_DIR, exist_ok=True)
    
    r = requests.get(URL)
    js = r.json()
    json_write = json.dumps(js)

    date = datetime.datetime.utcnow().isoformat()[:10]
    json_fname = os.path.join(DATA_DIR, 'nasa_bolides_GLM_asof_' + date + '.json')
    with open(json_fname, "w") as f:
        f.write(json_write)

    df  = pd.DataFrame.from_dict(js['data'])



    tab = Table.from_pandas(df)
    
    tab.write(os.path.join(DATA_DIR, 'nasa_bolides_GLM_asof_' + date + '.xml'), format='votable', overwrite=True)

    for e in tab:
        identifier = e['_id']
        if identifier!=''
        dirname = os.path.join(DATA_DIR, str(identifier))
        os.makedirs(dirname, exist_ok=True)
        process_attachements(e['attachments'], dirname)


        
        



if __name__ == "__main__":
    main()
