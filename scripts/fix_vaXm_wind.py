import xarray as xr
import numpy as np
import sys
import numpy as np
import os
import shutil

def correct_va_component(ua_file_path, va_wrong_file_path, geo_em_path, xoffset, yoffset):
    """
    Corrects the va component at the specified height (m). 
    It uses ua from one file and va_wrong from another, and corrects va 
    using the rotation info from the geo_em file (COSALPHA and SINALPHA), 
    masking out the relaxation zone as the geo_em file give the complete domain.

    Parameters:
        ua_file_path (str): Path to the postprocessed NetCDF file with correct ua.
        va_wrong_file_path (str): Path to the file with incorrect va.
        geo_em_path (str): Path to the geo_em file.
        xoffset (int): Number of the grids in x direction in the relaxation zone
        yoffset (int): Number of the grids in y direction the relaxation zone
        
    Returns:
        va_corrected (xarray.DataArray): Corrected va variable.
        
    To run:
        python fix_vaXm_wind.py ua_file va_wrong_file geo_em xoffset {yoffset}
        NOTE: yoffset not obligatory if xoffset=yoffset
    """
        
    # Load ua from correct file
    ds_ua = xr.open_dataset(ua_file_path)
    dirname, basename = os.path.split(ua_file_path)
    ua_var = basename.split("_", 1)
    ua = ds_ua[f'{ua_var[0]}']
    lat = ds_ua['lat']
    lon = ds_ua['lon']

    # Load va_wrong from the broken file
    ds_va = xr.open_dataset(va_wrong_file_path)
    dirname, basename = os.path.split(va_wrong_file_path)
    va_var = basename.split("_", 1)
    va_wrong = ds_va[f'{va_var[0]}']

    # Load geo_em fields
    ds_geo = xr.open_dataset(geo_em_path)
    if xoffset > 0 and yoffset > 0:
        sinalpha = ds_geo['SINALPHA'].squeeze()[yoffset:-(yoffset+2), xoffset:-(xoffset+2)]
        cosalpha = ds_geo['COSALPHA'].squeeze()[yoffset:-(yoffset+2), xoffset:-(xoffset+2)]
        xlat = ds_geo['XLAT_M'].squeeze()[yoffset:-(yoffset+2), xoffset:-(xoffset+2)]
        xlon = ds_geo['XLONG_M'].squeeze()[yoffset:-(yoffset+2), xoffset:-(xoffset+2)]
    elif xoffset > 0 and yoffset == 0:
        sinalpha = ds_geo['SINALPHA'].squeeze()[yoffset:-(yoffset+2), :]
        cosalpha = ds_geo['COSALPHA'].squeeze()[yoffset:-(yoffset+2), :]
        xlat = ds_geo['XLAT_M'].squeeze()[yoffset:-(yoffset+2), :]
        xlon = ds_geo['XLONG_M'].squeeze()[yoffset:-(yoffset+2), :]
    elif xoffset == 0 and yoffset > 0:
        sinalpha = ds_geo['SINALPHA'][:, xoffset:xoffset-2]
        cosalpha = ds_geo['COSALPHA'][:, xoffset:xoffset-2]
        xlat = ds_geo['XLAT_M'][:, xoffset:xoffset-2]
        xlon = ds_geo['XLONG_M'][:, xoffset:xoffset-2]
    else:
        sinalpha = ds_geo['SINALPHA'].squeeze()
        cosalpha = ds_geo['COSALPHA'].squeeze()
        xlat = ds_geo['XLAT_M'].squeeze()
        xlon = ds_geo['XLONG_M'].squeeze()
        
    # Check if raw and postprocessed domain mach
    if (lat.shape != xlat.shape) or (lon.shape != xlon.shape) \
        or not np.isclose(lat.values[0, 0], xlat.values[0, 0], atol=1e-4) \
        or not np.isclose(lon.values[0, 0], xlon.values[0, 0], atol=1e-4):
        print(lat.values[0, 0])
        print(xlat.values)
        
        print("Error: Shape or value mismatch between raw and postprocessed domain.")
        print("Check if the numebr of grids in the relaxation zone is correct.")
        print(f"{lat.shape},{xlat.shape},{lat.values[0, 0]},{xlat.values[0, 0]}, {lon.values[0, 0]},{xlon.values[0, 0]}")
        sys.exit(1)  # Exit the script with an error code
        
    else:
        print("lat, lon, xlat, xlon shapes match and values at (0, 0) are equal. Proceeding with the process.")

    # Fill 2D variables with 3rd dimension, to fit the shapes of ua and va_wrong
    sinalpha = np.broadcast_to(sinalpha, va_wrong.shape)
    cosalpha = np.broadcast_to(sinalpha, va_wrong.shape)

    # Compute corrected components
    denom = sinalpha + cosalpha
    uava2 = va_wrong / denom
    uava1 = (ua + va_wrong * sinalpha / denom) / cosalpha
    va_corrected = uava1 * sinalpha + uava2 * cosalpha
    
    return va_corrected

if __name__ == "__main__":
    # Arguments
    if len(sys.argv) < 5:
        print("Usage: python fix_vaXm_wind.py [ua_file] [va_wrong_file] [geo_em_file] [xoffset]")
        print("If xoffset != yoffset: python fix_vaXm_wind [ua_file] [va_wrong_file] [geo_em_file] [xoffset] [yoffset]")
        sys.exit(1)
    elif len(sys.argv) == 5:
        xoffset = int(sys.argv[4])
        yoffset = xoffset
    else:
        xoffset = int(sys.argv[4])
        yoffset = int(sys.argv[5])       

    ua_file = sys.argv[1]
    va_wrong_file = sys.argv[2]
    geo_em_file = sys.argv[3]

    # Fix va component
    va_fixed = correct_va_component(ua_file, va_wrong_file, geo_em_file, xoffset, yoffset)
    
    # Backup the wrong file. Comment out if not necessary
    dirname, fname = os.path.split(va_wrong_file)
    va_var = fname.split("_", 1)    
    if len(va_var) == 2:
        new_fname = f"{va_var[0]}wrong_{va_var[1]}"
    else:
        new_fname = f"{basename}.backup"    
    backup_path = os.path.join(dirname, new_fname)
    shutil.copy(va_wrong_file, backup_path)
    print(f"📝 Backed up original file to: {backup_path}")
    
    # Save into the file named as orginal
    ds_va = xr.open_dataset(va_wrong_file)
    ds_va.attrs['comment1'] = f'Derotation of {va_var[0]} corrected.'  # Add comment that the file if fixed. If not necessary commnet it out
    ds_va[f'{va_var[0]}'].data[:] = va_fixed.data                      # Replacing the values
    os.remove(va_wrong_file)                                           # Remove the original file to avoid permission problems
    ds_va.to_netcdf(va_wrong_file)                                     # Save to the original netcdf file
    print(f"✅ Corrected va written to: {va_wrong_file}")
    ds_va.close()
