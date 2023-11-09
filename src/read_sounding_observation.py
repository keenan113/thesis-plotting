import pandas
import os
from datetime import datetime,timedelta
import requests
from bs4 import BeautifulSoup
import argparse


def main():
    args = _parse_args()
    valid_time = datetime.strptime(args.valid_time,"%Y%m%d%H")
    output_filename = write_sounding_text_file(args.outdir,valid_time,args.station_id)
    print(read_fixed_width_sounding(output_filename))

def read_fixed_width_sounding(filename):
    column_names=['pressure_hPa','height_m','temperature_C','dewpoint_C',
    'realtive_humidity','mixing_ratio_gkg','wind_direction',
    'wind_speed_kts','potential_temperature',
    'equivalent_potential_temperature','virtual_potential_temperature']
    data = pandas.read_fwf(filename,names=column_names,skiprows=5)
    return data

def write_sounding_text_file(outdir,valid_time,station_id):
    response_text = query_wyoming_state_sounding_archive(valid_time,station_id)
    data_string = get_data_from_wyoming_state_html(response_text)
    output_filename = os.path.join(outdir,f'{station_id}_{valid_time:%Y%m%d_%H00Z}')
    with open(output_filename,'w') as f:
        f.write(data_string)
    return output_filename

def get_data_from_wyoming_state_html(html):
    bs = BeautifulSoup(html)
    data = bs.find('pre')
    return data.string
    
def query_wyoming_state_sounding_archive(valid_time,station_index):
    toplevelurl = 'https://weather.uwyo.edu/cgi-bin/sounding?region=naconf&TYPE=TEXT%3ALIST' 
    full_url = f'{toplevelurl}&YEAR={valid_time:%Y}&MONTH={valid_time:%m}&FROM={valid_time:%d%H}&TO={valid_time:%d%H}&STNM={station_index}'
    print(full_url)
    response = requests.get(full_url,verify=False)
    return response.text

def _parse_args():
    parser = argparse.ArgumentParser(
        description = 'Write and then read into a dataframe a sounding text file obtained from the Wyoming state archive'
    )
    parser.add_argument('--station-id',help='a 5 digit integer indicating the station to request data for')
    parser.add_argument('--valid-time',help='A valid-time indicating when a profile should be requested from. Formatted like YYYYMMDDHH')
    parser.add_argument('--outdir',help='full path to output directory')
    return parser.parse_args()

if __name__ == "__main__":
    main()