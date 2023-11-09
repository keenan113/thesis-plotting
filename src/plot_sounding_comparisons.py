
from ensemble_datasets import Dataset
from iris_grib_mappings import _field_keyval_mapping
from glob import glob
import xarray
import numpy
import seaborn as sns
import pandas
from read_sounding_observation import read_fixed_width_sounding
import os
import argparse
from matplotlib import pyplot

sounding_location_meta = {
    "OKX":{
        "full_name":"Upton, NY",
        "station_id": 72501,
        "lat":40.87,
        "lon":-72.86
    },
    "CHH":{
        "full_name":"Upton, NY",
        "station_id": 74494,
        "lat":41.66,
        "lon":-69.96,
    },
}

def main():
    args = parse_args()
    single_member_filenames = glob(os.path.join(args.grib_directory,"YesradRun/wrfout_201302*0001.grib2"))
    single_member_filenames.sort()

    observed_sounding_00Z = read_fixed_width_sounding(
        os.path.join(
            '/Users/keenanfryer/projects/thesis-plotting/sounding_observations',
            f"{sounding_location_meta[args.station_name]['station_id']}_20130209_0000Z")
    )
    if args.station_name == "OKX":
        observed_sounding_12Z = read_fixed_width_sounding(
            os.path.join(
                '/Users/keenanfryer/projects/thesis-plotting/sounding_observations',
                f"{sounding_location_meta[args.station_name]['station_id']}_20130209_1200Z")
        )


    for fn_time in single_member_filenames:
        #get filename globs for building the datasets
        yesrad_fn_glob = fn_time.replace("0001.grib2","0002.grib2")
        print(yesrad_fn_glob)
        norad_fn_glob = yesrad_fn_glob.replace("YesradRun","NoradRun")
        
        #load the two datasets
        yesrad_ds = Dataset(yesrad_fn_glob,ensemble_data=True)
        norad_ds = Dataset(norad_fn_glob,ensemble_data=True)

        #get that single point profile
        yesrad_profile = yesrad_ds.get_isobaric_profile_dataset(
            numpy.array(sounding_location_meta[args.station_name]['lat']),
            numpy.array(sounding_location_meta[args.station_name]['lon'])
        )
        norad_profile = norad_ds.get_isobaric_profile_dataset(
            numpy.array(sounding_location_meta[args.station_name]['lat']),
            numpy.array(sounding_location_meta[args.station_name]['lon'])
        )

        #convert datasets to dataframes for plotting
        yes_plot_profile = yesrad_profile.to_dataframe()
        no_plot_profile = norad_profile.to_dataframe()

        #add a celsius temperature column to each dataframe
        yes_plot_profile['temperature_C'] = yes_plot_profile['temperature'] - 273.15
        no_plot_profile['temperature_C'] = no_plot_profile['temperature'] - 273.15

        #get valid-time for plotting 
        valid_time = pandas.to_datetime(yes_plot_profile.time.unique()[0])
        print(f'Data loaded and now plotting ensemble data for {valid_time:%Y%m%d_%H%M}')

        ax = sns.lineplot(data=subset_profile_by_pressure(yes_plot_profile,500),x='temperature_C',y='pressure',orient='y',color='red',**{'label':'yes-radar'})
        ax = sns.lineplot(data=subset_profile_by_pressure(no_plot_profile,500),x='temperature_C',y='pressure',orient='y',color='blue',**{'label':'no-radar'})
        ax = sns.lineplot(data=subset_profile_by_pressure(observed_sounding_00Z,500,pressure_keyname='pressure_hPa'),x='temperature_C',y='pressure_hPa',orient='y',**{'color':'black','linestyle':'solid'})
        if args.station_name == "OKX":
            ax = sns.lineplot(data=subset_profile_by_pressure(observed_sounding_12Z,500,pressure_keyname='pressure_hPa'),x='temperature_C',y='pressure_hPa',orient='y',**{'color':'black','linestyle':'dashed'})
        ax.invert_yaxis()
        ax.figure.savefig(os.path.join(args.outdir,f"Sounding_comparison_{args.station_name}_{valid_time:%Y%m%dT%H%M}.png"))
        pyplot.close('all')

def subset_profile_by_pressure(indf,pressure_ceiling,pressure_keyname='pressure'):
    return indf.query(f'{pressure_keyname} >= @pressure_ceiling').copy()


def parse_args():
    parser = argparse.ArgumentParser(
        description='Plot sounding comparisons of all members - compared to observed sounding profiles'
    )
    parser.add_argument('--station-name',help='sounding location name. Either CHH or OKX')
    parser.add_argument('--grib-directory',help='Absolute path to input grib files. Should have subdirs like YesradRun, NoradRun, etc.')
    parser.add_argument('--outdir',help='absolute path to output directory')
    return parser.parse_args()

if __name__ == "__main__":
    main()