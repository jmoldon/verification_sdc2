import os
import sys
import shutil
import re
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
import yaml
from astropy.io import fits
from astropy.wcs import WCS
from shutil import which

data_parameters = './parameters/data.yml'
param_development_small = './parameters/sofia_dev_small.par'

results_path = 'results'

dev_small_cat = './results/development_small/developement_small_cat.txt'
final_cat = './results/development_small/final_dev_smal.csv'

with open(data_parameters, "r") as f:
    data_yml = yaml.load(f, Loader=yaml.FullLoader)

data_path = data_yml['data_path']
fitsfile = os.path.join(data_path, 'development_small/sky_dev.fits')



# Functions
def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""
    return which(name) is not None


def download_data(data_parameters, force=False):
    if not os.path.isdir(data_path):
        os.mkdir(data_path)
        
    for dataset in data_yml['download_locations'].keys():
        print(dataset)
        dataset_dir = os.path.join(data_path, dataset)
        if force == True:
            shutil.rmtree(dataset_dir)

        if not os.path.isdir(dataset_dir):
            pathname = data_yml['download_locations'][dataset]['path']
            os.mkdir(dataset_dir)
            for filename in data_yml['download_locations'][dataset]['files']:
                command = f'wget --no-check-certificate "{pathname}download?path=%2F&files={filename}" -O {dataset_dir}/{filename}'
                print(command)
                os.system(command)

def run_sofia(parameters, outputdir):
    """Only executed if directory does not exist"""
    if not os.path.isdir(results_path):
        os.mkdir(results_path)
        
    if not os.path.isfile(os.path.join(results_path, outputdir,'developement_small_cat.txt')):
        try:
            os.mkdir(os.path.join(results_path, outputdir))
        except FileExistsError:
            pass
        command = f"sofia {parameters}"
        print(command)
        if is_tool('sofia'):
            os.system(command)
        else:
            print('sofia not available. Please install Sofia-2')
            sys.exit(1)
    return

def read_sofia_header(filename):
    with open(filename, 'r') as f:
        head_line = f.readlines()[10]
    head = re.split('\s+', head_line.strip('\n'))[1:] # 1: to remove #
    return head
        
def sofia2cat(catalog):
    head = read_sofia_header(catalog)
    raw_cat = pd.read_csv(catalog, delim_whitespace=True, header=None, names=head, comment='#')
    raw_cat.sort_values(by='f_sum', ascending=False, inplace=True)
    raw_cat_filtered = raw_cat[raw_cat['kin_pa']>0]
    print('Sofia raw catalog filtered:')
    print(raw_cat_filtered[['x', 'y', 'ell_maj', 'ell_min', 'f_sum', 'freq', 'kin_pa', 'w20']])
    return raw_cat_filtered

def pix2coord(wcs, x, y):
    coord = wcs.pixel_to_world(x, y, 1)
    #print('coord')
    #print(coord)
    return coord[0].ra.deg, coord[0].dec.deg

def compute_inclination(bmaj, bmin):
    # returns an angle in degrees
    return np.arctan2(bmin, bmaj)*180./np.pi

def convert_units(raw_cat, fitsfile):
    f = fits.open(fitsfile)
    wcs=WCS(f[0].header)
    f.close()
    # Convert x,y in pixels to R.A.,Dec. in deg
    ra_deg, dec_deg = pix2coord(wcs, raw_cat['x'], raw_cat['y'])
    # Get pixel size
    pix2arcsec = wcs.wcs.get_cdelt()[1]*3600. # This assumes same pixel size in both direction
    return ra_deg, dec_deg, pix2arcsec

def process_catalog(raw_cat, fitsfile):
    # Unit conversion
    ra_deg, dec_deg, pix2arcsec = convert_units(raw_cat, fitsfile)
    hi_size = raw_cat['ell_maj']*pix2arcsec
    # Estimate inclination based on fitted ellipsoid, assuming the galaxy is intrinsically circular
    inclination = compute_inclination(raw_cat['ell_maj'], raw_cat['ell_min'])

    # Construct the output catalog
    processed_cat = pd.DataFrame()
    processed_cat['id'] = raw_cat['id']
    processed_cat['ra'] = ra_deg
    processed_cat['dec'] = dec_deg
    processed_cat['hi_size'] = hi_size
    processed_cat['line_flux_integral'] = raw_cat['f_sum']  # we need to clarify if this is the right magnitude and the right units
    processed_cat['central_freq'] =  raw_cat['freq'] # we need to clarify if what sofia gives is the central freq
    processed_cat['pa'] = raw_cat['kin_pa']  # we need to clarify if Sofia kinematic angle agrees with their P.A. 
    processed_cat['i'] = inclination
    processed_cat['w20'] = raw_cat['w20'] # we need to clarify if the units and the definition is the same 
    processed_cat.reset_index(drop=True, inplace=True)
    return processed_cat


def main():
    download_data(data_parameters, force=False)
    run_sofia(parameters=param_development_small,
              outputdir='development_small')
    raw_cat = sofia2cat(catalog=dev_small_cat)
    processed_cat = process_catalog(raw_cat, fitsfile)
    print(processed_cat)
    print(f'This catalog is being saved in: {final_cat}')
    processed_cat.to_csv(final_cat, sep=' ', index=False, float_format="%.4f")

if __name__ == "__main__":                
    main()

