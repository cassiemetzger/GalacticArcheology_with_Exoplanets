import os 
print("Installing dependencies...")
os.system("pip install -r requirements.txt")
print("Done!")
import time
import numpy as np
import matplotlib.pyplot as plt
import astropy 
from astropy.io import fits, ascii
from astropy.table import Table
import astroquery 
import glob
from astroquery.mast import Observations, Mast, MastMissions, MastMissionsClass, Catalogs
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.xmatch import XMatch
from astropy.timeseries import BoxLeastSquares
from astroquery.sdss import SDSS
import pandas as pd
import multiprocessing
import time
import ipywidgets as widgets
__all__ = [
    "np", 
    "plt", 
    "fits", 
    "ascii", 
    "Table",
    "Observations",
    "u", 
    "SkyCoord", 
    "XMatch", 
    "BoxLeastSquares", 
    "SDSS", 
    "pd", 
    "widgets"

]
# configure matplotlib settings
params = {'axes.labelsize': 12, 'xtick.labelsize': 12, 'ytick.labelsize': 12, 
          'text.usetex': False, 'lines.linewidth': 1,
          'axes.titlesize': 18, 'font.family': 'serif', 'font.size': 12}
plt.rcParams.update(params)

current_files = glob.glob("*")
if("allStar-dr17-synspec_rev1.fits" not in current_files): 
    print("Downloading APOGEE allStar...")
    Observations.download_file('mast:SDSS/apogee/allStar-dr17-synspec_rev1.fits')
else: 
    print("APOGEE allStar found!")
     
def check_file_download(name): 
    files = glob.glob('tess_data/*')
    files = [f[10:] for f in files]
    if(name in files): 
        return True 
    else: 
        return False 

# function for downloading tess data 
def tess_data_download(num_obj = None, order_random = False, target_name = None, search_radius = 0.1): 
    transiting_exo = ascii.read('transiting_exoplanets.csv')
    os.system("mkdir tess_data")
    if(target_name is not None): 
        target_index = np.where(transiting_exo['pl_name'] == target_name.strip())[0]
        if(len(target_index) < 0): 
            print("Could not resolve target name.")
        else: 
            print("Target name resolved.")
            row = transiting_exo[target_index]
            ra = row['ra'][0]
            dec = row['dec'][0]
            obj_name = row['pl_name'][0]
            obj_name=obj_name.replace(" ", "_")
            obj_downloaded=check_file_download(obj_name)
            if(obj_downloaded): 
                print("File has already been downloaded")
            else: 
                observations = Observations.query_region(f"{ra} {dec}", radius=f"{search_radius} deg")
                if(len(observations)<0):
                    print("No APOGEE match found")
                else:  
                    print("Observations found.")
                    obs_wanted = ((observations['dataproduct_type'] == 'timeseries') &
                    (observations['obs_collection'] == 'TESS'))
                    if(len(observations[obs_wanted])>0): 
                        print("Retrieving data products...")
                        data_products = Observations.get_product_list(observations[obs_wanted][:5])
                        products_wanted = Observations.filter_products(data_products, 
                                                        productSubGroupDescription=["LC"])
                        os.system(f"mkdir tess_data/{obj_name}")
                        print(f"Downloading data for {obj_name}...")
                        download = Observations.download_products(products_wanted)
                        for file in download['Local Path']: 
                            os.system(f"cp {file} ./tess_data/{obj_name}")
                        os.system(f"rm -r mastDownload")
                    else: 
                        print("No APOGEE match found")
    else: 
        if(num_obj is not None): 
            obj_limit = num_obj
        else: 
            obj_limit = len(transiting_exo)
        if(order_random): 
            for i in range(obj_limit): 
                obs_found = False 
                while obs_found == False: 
                    idx = np.random.randint(0, len(transiting_exo))
                    row = transiting_exo[idx]
                    ra = row['ra']
                    dec = row['dec']
                    obj_name = row['pl_name']
                    obj_name=obj_name.replace(" ", "_")
                    obj_downloaded=check_file_download(obj_name)
                    if(obj_downloaded): 
                        print("File has already been downloaded")
                    else: 
                        observations = Observations.query_region(f"{ra} {dec}", radius=f"{search_radius} deg")
                        if(len(observations)>0): 
                            obs_wanted = ((observations['dataproduct_type'] == 'timeseries') &
                            (observations['obs_collection'] == 'TESS'))
                            if(len(observations[obs_wanted])>0):    
                                obs_found= True 
                                data_products = Observations.get_product_list(observations[obs_wanted])
                                products_wanted = Observations.filter_products(data_products, 
                                                                productSubGroupDescription=["LC"])
                                os.system(f"mkdir tess_data/{obj_name}")
                                print(f"Downloading data for {obj_name}...")
                                download = Observations.download_products(products_wanted)
                                for file in download['Local Path']: 
                                    os.system(f"cp {file} ./tess_data/{obj_name}")
                                os.system(f"rm -r mastDownload")
        else: 
            num_objects = 0 
            i = 0
            while (num_objects < obj_limit) and (i < len(transiting_exo)):
                row = transiting_exo[i]
                ra = row['ra']
                dec = row['dec']
                obj_name = row['pl_name']
                obj_name=obj_name.replace(" ", "_")
                obj_downloaded=check_file_download(obj_name)
                if(obj_downloaded): 
                    print("File has already been downloaded")
                else: 
                    observations = Observations.query_region(f"{ra} {dec}", radius=f"{search_radius} deg")
                    if(len(observations)> 0): 
                        obs_wanted = ((observations['dataproduct_type'] == 'timeseries') &
                                        (observations['obs_collection'] == 'TESS'))
                        if(len(observations[obs_wanted])>0): 
                            obs_found = True 
                            num_objects+=1
                        
                            data_products = Observations.get_product_list(observations[obs_wanted])
                            products_wanted = Observations.filter_products(data_products, 
                                                        productSubGroupDescription=["LC"])
                            os.system(f"mkdir tess_data/{obj_name}")
                            print(f"Downloading data for {obj_name}...")
                            download = Observations.download_products(products_wanted)
                            for file in download['Local Path']: 
                                os.system(f"cp {file} ./tess_data/{obj_name}")
                            os.system(f"rm -r mastDownload")
                        i+=1
                    else: 
                        print("No source found.") 
                        i+=1 
                        print(f"Now trying {row['pl_name']}")               
            else: 
                print("No source found.")
                i+=1 
                print(f"Now trying {row['pl_name'][i]}")
                
## functions for modeling an exoplanet 
def orbit(R, T, time):
    """
    Defines a circular orbit of radius R and period T over time 
    """
    # time is expected in days 
    # radius in parsec 
    T = T * 86400
    time = time * 86400 
    R = R* 3.086e+18
    w = 2*np.pi/T 
    return R*np.cos(w*time), R*np.sin(w*time)
def target_star(x1, y1, vt, time): 
    """
    Defines the motion of some star relative to the observing planet 
    Looks at the trpansverse motion of the star only
    """
    # vt is expected to be in km/s 
    # time in days 
    # x1 and y1 are in coordinate units 
    time = time*86400
    vt = vt* 100000
    x_init = x1 * 3E18 
    y_init = y1 * 3E18 
    xlist = []
    ylist= []
    for val in time: 
        d = np.sqrt(x_init**2 + y_init**2)
        vtx= vt * (x_init / d)
        vty = vt * (y_init / d)
        x = x_init + vtx*val 
        y = y_init + vty*val 
        x_init = x 
        y_init = y 
        xlist.append(x)
        ylist.append(y)
    return xlist, ylist
def distance(x1, y1, x2, y2): 
    """
    Computes the distance between the planet and target star
    """
    return np.sqrt((x1-x2)**2 + (y1-y2)**2)
def flux(d, Lsun): 
    """
    Computes the observed flux as a function of distance
    """
    return Lsun/(4*np.pi *d**2)
