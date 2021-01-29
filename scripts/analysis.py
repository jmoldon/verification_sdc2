import os
import shutil
import re
from datetime import datetime
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
import yaml
from astropy.io import fits
from astropy.wcs import WCS


#set_type ='evaluation'
#set_type ='development_small'
set_type= 'debug'
run_fat = False
type = {'evaluation': {'dir': 'evaluation', 'file': 'sky_eval'},
        'development_small': {'dir': 'development_small', 'file': 'sky_dev'},
        'debug': {'dir': 'debug', 'file': 'test'},
}


main_dir = '/home/peter/SDC2'

data_parameters = './parameters/data.yml'
param_development_small = './parameters/sofia_dev_small.par'
fitsfile = f'''{main_dir}/{type[set_type]['dir']}/{type[set_type]['file']}.fits'''

data_path = f'''{main_dir}/{type[set_type]['dir']}'''
results_path = f'''{main_dir}/results'''

dev_small_cat = f'''{results_path}/{set_type}_small_cat.txt'''
final_cat = f'''{results_path}/final_catalogue_{set_type}.csv'''

def download_data(data_parameters, type = 'debug', force=False):
    with open(data_parameters, "r") as f:
        data_yml = yaml.load(f, Loader=yaml.FullLoader)

    if force == True:
            shutil.rmtree(data_path)

    if not os.path.isdir(data_path):
        os.mkdir(data_path)
    else:
        if os.path.isfile(fitsfile):
            print(f'There is no need to download {fitsfile} as it already exists')
            # This could be expanded to check for the readme and continuum
            return

    for filename in data_yml[type]['files']:
        pathname = data_yml[type]['path']
        command = f'wget --no-check-certificate "{pathname}download?path=%2F&files={filename}" -O {data_path}/{filename}'
        print(command)
        os.system(command)


def run_sofia(parameters, outputdir):
    """Only executed if the output catalog  does not exist"""

    #It makes sense to not run this when the results exist but maybe a check on an existing catalog is better
    if not os.path.isfile(os.path.join(results_path, outputdir,f'{outputdir}_cat.txt')):
        if not os.path.isdir(os.path.join(results_path, outputdir)):
            os.mkdir(os.path.join(results_path, outputdir))
        # I guess the 2 is because of my dual installation of SoFiA versions we should implement a version check
        command = f"sofia2 {parameters}"
        print(command)
        os.system(command)
        command = f'mv {parameters} {os.path.join(results_path, outputdir)}/sofia_input_parameters.par'
        print(command)
        os.system(command)
    else:
        print(f'''We have already found the catalogue {os.path.join(results_path, outputdir,f'{outputdir}_cat.txt')}, continuing to process.''' )

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
    if 'freq' in raw_cat_filtered:
        print(raw_cat_filtered[['x', 'y', 'ell_maj', 'ell_min', 'f_sum', 'freq', 'kin_pa', 'w20']])
    elif 'v_app' in raw_cat_filtered:
        print(raw_cat_filtered[['x', 'y', 'ell_maj', 'ell_min', 'f_sum', 'v_app', 'kin_pa', 'w20']])
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
    pix2freq = f[0].header['CDELT3']
    return ra_deg, dec_deg, pix2arcsec,pix2freq

def frequency_to_vel(freq, invert=False):
    f0 = 1420405751.786           #Hz
    c = 299792.458 #km/s
    if not invert:
        return c*((f0**2-freq**2)/(f0**2+freq**2))
    else:
        return f0*np.sqrt((1-freq/c)/(1+freq/c))

# Convert the frequency axis of a cube
def convert_frequency_axis(filename, outname, velocity_req = 'radio'):
    print(filename)
    cube = fits.open(filename)
    hdr = cube[0].header
    # Check we have a proper third axis
    if hdr['CTYPE3'].lower() != 'freq' or hdr['NAXIS'] < 3:
        print('We can not convert this axis as it is not a frequency axis')
        return
    f0 = 1420405751.786           #Hz
    c = 299792458 #m/s
    # get central values
    crpix = float(hdr['CRPIX3'])
    crval = float(hdr['CRVAL3'])
    naxis_len = float(hdr['NAXIS3'])
    # make sure the central pixel is rather central else large errors are introduce in both vrad and rel
    if naxis_len/2.-5 < crpix <  naxis_len/2.+5:
            hdr_wcs = WCS(hdr)
            centralx,centraly, new_freq = hdr_wcs.pix2world([hdr['CRPIX1'],hdr['CRPIX2'],naxis_len/2.],1)
            hdr['CRPIX3'] = new_pix
            crval = new_freq
    #Now convert
    if velocity_req == 'radio':
          # convert from frequency to radio velocity
            cdelt_vel = -c*float(hdr['CDELT3'])/f0
            crval_vel = c*(1-crval/f0)
            # https://fits.gsfc.nasa.gov/standard40/fits_standard40aa-le.pdf
            hdr['CTYPE3'] = 'VRAD'
    elif velocity_req == 'relativistic':
        # This should always only ever be used for cubes with small velocity range
        crval_vel = frequency_to_vel(crval)
        freq_step = float(hdr['CDELT3'])
        central_two = frequency_to_vel(crval+freqstep)
        lower_one = frequency_to_vel(crval-(naxis_len/2.)*freqstep)
        lower_two = frequency_to_vel(crval-(naxis_len/2.+1)*freqstep)
        upper_one = frequency_to_vel(crval+(naxis_len/2.-1.)*freqstep)
        upper_two = frequency_to_vel(crval+(naxis_len/2.)*freqstep)
        cdelt_vel = np.mean([central_two-crval_vel,lower_two-lower_one,upper_two-upper_one])
        if cdelt_vel*naxis_len > 1e6:
            print('This cube is too big for a relativistic conversion')
            return
        hdr['CTYPE3'] = 'VELO'
    else:
        print('We dont do those things here.')
        return
    hdr['CDELT3'] = cdelt_vel
    hdr['CRVAL3'] = crval_vel

    if 'CUNIT3' in hdr:
        # delete cunit3 because we adopt the default units = m/s
        del hdr['CUNIT3']
    fits.writeto(outname,cube[0].data,hdr,overwrite = True)





def process_catalog(raw_cat, fitsfile):
    # Unit conversion
    ra_deg, dec_deg, pix2arcsec,pix2vel = convert_units(raw_cat, fitsfile)
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
    if 'freq' in raw_cat:
        processed_cat['central_freq'] =  raw_cat['freq']
        processed_cat['central_velocity'] = frequency_to_vel(raw_cat['freq'])
        processed_cat['w20'] = frequency_to_vel(raw_cat['freq']-raw_cat['w20']/2.*pix2vel)-frequency_to_vel(raw_cat['freq']+raw_cat['w20']/2.*pix2vel) # we need to clarify if the units and the definition is the same
    else:
        processed_cat['central_velocity'] =  raw_cat['v_app']
        processed_cat['central_freq'] = frequency_to_vel(raw_cat['v_app'],invert=True)
        processed_cat['w20'] = raw_cat['w20']*pix2vel
         # we need to clarify if what sofia gives is the central freq
    processed_cat['pa'] = raw_cat['kin_pa']  # we need to clarify if Sofia kinematic angle agrees with their P.A.
    processed_cat['i'] = inclination
    processed_cat.reset_index(drop=True, inplace=True)
    return processed_cat

def prepare_parameters(parameters=param_development_small, type ='debug'):
    parameters_in = read_sofia_parameters(param_development_small)
    parameters_in['input.data'] = f'{fitsfile}'
    parameters_in['output.directory'] = f'{results_path}/{type}'
    parameters_in['output.filename'] = f'{type}'
    if not os.path.isdir(results_path):
        os.mkdir(results_path)
    write_sofia_parameters(parameters_in, f'{results_path}/sofia_settings.par')


def write_sofia_parameters(template,name, debug = False):
    with open(name,'w') as file:
        for key in template:
            if key[0] == 'E' or key [0] == 'H':
                file.write(template[key])
            else:
                file.write(f"{key} = {template[key]}\n")

def read_sofia_parameters(filename,debug = False):
    with open(filename,'r') as f:
        template = f.readlines()
    result = {}
    counter = 0
    counter2 = 0
    # Separate the keyword names
    for line in template:
        key = str(line.split('=')[0].strip())
        if key == '':
            result[f'EMPTY{counter}'] = line
            counter += 1
        elif key[0] == '#':
            result[f'HASH{counter2}'] = line
            counter2 += 1
        else:
            result[key] = str(line.split('=')[1].strip())
    return result

def organize_sofia(catalog,convert= True, type='debug'):
    fat_catalog = {'id': ['number'], 'dist': ['Distance'], 'dir': ['Directoryname'], 'cube': ['Cubename']}
    #sofia_output = ['spec.txt','chan.fits','mom0.fits','mom1.fits','mom2.fits','mask.fits','cube.fits']o
    sofia_output = ['cube.fits']
    for source in catalog['id']:
        if not os.path.isdir(f'{results_path}/fat/sofia_{source}'):
            os.mkdir(f'{results_path}/fat/sofia_{source}')
        #Move all sofia out put to a proper directory
        for file in sofia_output:
            if convert:
                convert_frequency_axis(f'{results_path}/{type}/{type}_cubelets/{type}_{source}_{file}',\
                f'{results_path}/fat/sofia_{source}/{type}_{source}_{file}')
            else:
                command= f'cp -f {results_path}/{type}/{type}_cubelets/{type}_{source}_{file} {results_path}/fat/sofia_{source}/{type}_{source}_{file}'
                os.system(command)
        fat_catalog['id'].append(source)
        fat_catalog['dist'].append('-1')
        fat_catalog['dir'].append(f'sofia_{source}')
        fat_catalog['cube'].append(f'{type}_{source}_cube')
    with open(f'{results_path}/fat/fit_catalogue.txt','w') as f:
        for i in range(len(fat_catalog['id'])):
            f.write(f'''{fat_catalog['id'][i]}|{fat_catalog['dist'][i]}|{fat_catalog['dir'][i]}|{fat_catalog['cube'][i]}\n''')

def fat_configuration(filename,type='debug'):
    with open(filename,'r') as f:
        template = f.readlines()

    with open(f'{results_path}/fat/FAT_INPUT.config','w') as f:
        f.write(f'#This is the configuration file for fit {type} at {datetime.now()} \n')

    with open(f'{results_path}/fat/FAT_INPUT.config','a') as f:
        for line in template:
            setting = line.split('=')[0].strip()
            if setting == 'catalogue':
                line = f'catalogue = {results_path}/fat/fit_catalogue.txt \n'
            elif setting == 'maindir':
                line = f'maindir = {results_path}/fat/ \n'
            elif setting == 'outputcatalogue':
                line = f'outputcatalogue={results_path}/fat/fat_results.txt \n'
            elif setting == 'outputlog':
                line = f'outputlog = log.txt \n'
            f.write(line)


def main():
    download_data(data_parameters, type= set_type, force=False)
    prepare_parameters(parameters=param_development_small, type = set_type)
    run_sofia(parameters=f'{results_path}/sofia_settings.par',
              outputdir= set_type)

    raw_cat = sofia2cat(catalog=os.path.join(results_path, set_type,f'{set_type}_cat.txt'))
    processed_cat = process_catalog(raw_cat, fitsfile)
    print(processed_cat)
    print(f'This catalog is being saved in: {final_cat}')
    processed_cat.to_csv(final_cat, sep=' ', index=False, float_format="%.4f")

    if run_fat:
        if not os.path.isdir(f'{results_path}/fat'):
            os.mkdir(f'{results_path}/fat')
        convert = False
        if 'freq' in raw_cat:
            convert = True
        organize_sofia(processed_cat,convert= convert, type=set_type)
        fat_configuration('./parameters/FAT_INPUT.config',type=set_type)
        command = f'pyFAT -t -c {results_path}/fat/FAT_INPUT.config'
        print(command)
        os.system(command)


if __name__ == "__main__":
    main()
