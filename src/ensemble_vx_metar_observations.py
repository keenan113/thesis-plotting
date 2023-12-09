import ensemble_datasets
import xarray
import util
import argparse
import os

def main():
    args = _parse_args()
    
    ensemble_ds = ensemble_datasets.Dataset(
        args.ensemble_file_glob,
        ensemble_data=True,
        dataset_name=args.dataset_name
        )
    ensemble_data = rename_fields(ensemble_ds.data.mean(dim='member'))
    for time in ensemble_ds.valid_time:
        forecast_ds = ensemble_data[args.field_to_verify].sel(time=time)
        metar_ds = load_metar_ds(args.metar_file_path,time)
        interp_forecast = interpolate_forecast_data_to_metar(forecast_ds,metar_ds,ensemble_ds.cartopy.projection)
        errors = calculate_differences(interp_forecast,metar_ds,args.field_to_verify)
        rmse = (errors**2).mean(dim="index")**0.5
        print(f"For the valid-time: {time:%Y%m%d %H:%M}; the {args.dataset_name}:{args.field_to_verify} RMSE is: {rmse}")

def calculate_differences(forecast,observations,field_name):
    assert f"temperatureDD" in observations.variables
    qc_inds = observations.where(observations[f"temperatureDD"] == b"V",drop=True)
    print(forecast)
    print(observations)
    differences = forecast.sel(index=qc_inds.index.values) - observations.sel(index=qc_inds.index.values)
    return differences

def interpolate_forecast_data_to_metar(fcst_ds,metar_ds,crs):
    projected_metar_coords = util.get_x_and_y_proj_coords_from_latlon(
        metar_ds.latitude.to_numpy(),
        metar_ds.longitude.to_numpy(),
        crs
    )

    interp_ds = fcst_ds.sel(
        projection_x_coordinate=projected_metar_coords['x'].to_xarray(),
        projection_y_coordinate=projected_metar_coords['y'].to_xarray(),
        method='nearest'
    )
    return interp_ds

def load_metar_ds(metar_dir,valid_time):
    metar_fn = os.path.join(
        metar_dir,
        f"{valid_time:%Y%m%d_%H%M}"
    )
    metar_ds = xarray.open_dataset(metar_fn)
    metar_ds = metar_ds.assign_coords({"recNum":metar_ds.recNum})
    metar_ds = metar_ds.rename({"recNum":"index"})
    metar_ds = rename_fields(metar_ds)
    return metar_ds

def rename_fields(inds):
    rename_dict = {
        'air_temperature_surface':'temperature-2m',
        'temperature':'temperature-2m', 
    }
    rename_dict = {k:v for k,v in rename_dict.items() if k in inds.variables}
    return inds.rename(rename_dict)

def _parse_args():
    parser = argparse.ArgumentParser(
        description = 'Calculate the RMSE of a WRF-based ensemble forecast'
    )
    parser.add_argument(
        "--ensemble-file-glob",
        help = "the file glob that grabs all the necessary WRF grib files",
        required = True
    )
    parser.add_argument(
        "--dataset-name",
        help = "human-readable name to assign to dataset",
        required = True
    )
    parser.add_argument(
        "--metar-file-path",
        help = "absolute path to the METAR netCDF data",
        default = '/F0/keenan/ThesisProject/Feb2013Case/MADISobs/point/metar/netcdf'
    )
    parser.add_argument(
        "--field-to-verify",
        help = 'The field name that we want to verify',
        required=True,
        choices = ['temperature-2m']
    )
    return parser.parse_args()


if __name__  == "__main__":
    main()
