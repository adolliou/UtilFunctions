from glob import glob
import numpy as np
from astropy.io import fits
import astropy.time
import os
import astropy.units as u
from tqdm import tqdm
import copy
import warnings


def create_dict_file(suffix: str,  path_instrument: str=None, name_list_txt: str = None,  window: int=None, sort_dict: bool=True,
                     verbose: int=1):
    """
    Create a dict_file dictionnary allowing to deal with the paths of the files, and easily access some of their properties :
    - "paths": list of paths to the FITS files. (list)
    - "date-avg", parameter DATE-AVG, time of the middle of the exposure. (list)
    - "dsun-obs", Distance in meter of the spacecraft to the sun center. (list)

    Args:
        suffix (str): _description_Suffix to of files to look for. typical suffix is suffix="*fsi174*.fits"
        path_instrument (str, optional): _description_. Defaults to None. Path to the folder where the instruments FITS files are, or 
        where is the name_list_txt.txt file is.
        name_list_txt (str, optional): _description_. Defaults to None. Path to a text files where the absolute paths to the FITS files are written.
        The latter is obtained through Selektor prior.  
        window (int, optional): _description_. Defaults to None. Window of the hdu list where the header is extracted 
        sort_dict (bool, optional): _description_. Defaults to True. If True, then the paths are sorted with time in the output list.
        verbose (int, optional): level of printing you want

    Returns:
        _type_: _description_
    """
    data_dict = {}

    if (name_list_txt is not None) and (path_instrument is not None):
        paths = []
        with open(os.path.join(path_instrument, name_list_txt), "r") as f:
            paths = f.read().splitlines()

    elif path_instrument is not None:
        paths = glob(os.path.join(path_instrument, suffix))
    else: 
        raise NotImplementedError("either path_instrument and path_list_txt, or path_instrument alone parameter is necessary")
    if verbose>0:
        print("create dictionary file with files in the folder %s " % path_instrument)
    data_dict['path_instrument'] = path_instrument
    data_dict['path'] = paths
    data_dict['date-avg'] = []
    data_dict['dsun-obs'] = []

    for kk, path in enumerate(tqdm(data_dict['path'], desc="Adding files to dict")):
        f = fits.open(path)

        if window is None:
            if len(f) > 1:
                idx = 1
            else:
                idx = 0
        else:
            idx = window
        if "DATE-AVG" in f[idx].header:
            data_dict['date-avg'].append(astropy.time.Time(f[idx].header['DATE-AVG']))
        elif ("DATE_AVG" in f[idx].header):
            data_dict['date-avg'].append(astropy.time.Time(f[idx].header['DATE_AVG']))

        if "DSUN_OBS" in f[idx].header:
            data_dict['dsun-obs'].append(f[idx].header['DSUN_OBS'])
        elif "DSUN-OBS" in f[idx].header:
            data_dict['dsun-obs'].append(f[idx].header['DSUN-OBS'])

        # data_dict['telescop'].append(f[idx].header['TELESCOP'])
        if ("DATE-AVG" not in f[idx].header) & ("DATE_AVG" not in f[idx].header):
            warnings.warn("DATE-AVG not found in header, manually compute it.")
            if f[idx].header["TELESCOP"] == "SDO/HMI":
                warnings.warn("use date-obs for HMI")
                cad = f[idx].header["TRECSTEP"]
                unit = f[idx].header["TRECUNIT"] 
                if unit == 'secs':
                    data_dict['date-avg'].append(astropy.time.Time(f[idx].header['DATE-OBS']) + 0.5*cad*u.s)
                else:
                    raise NotImplementedError
            else:
                raise NotImplementedError("The code below does not work for SPICE files. Better to raise an error.")

    data_dict['path'] = np.array(data_dict['path'])
    data_dict['date-avg'] = np.array(data_dict['date-avg'])
    data_dict['dsun-obs'] = np.array(data_dict['dsun-obs'])
    if len(data_dict['path']) == 0:
        raise ValueError(f"could construct dict data at {path_instrument=} with {suffix=}")
    if verbose>0:
        print("%i FITS files added in the dict_file." % len(data_dict["path"]))

    if sort_dict:
        return _sort_dict_file(data_dict)
    else:
        return data_dict


def _sort_dict_file(dict_file: dict):
    """Sort the dict file by FITS file date. 

    Args:
        dict_file (dict): input dict file

    Returns:
        dict : output sorted dict file
    """    
    ref_time = dict_file["date-avg"][0]
    d = copy.deepcopy(dict_file)
    time = [(n - ref_time).to(u.s).value for n in dict_file["date-avg"]]
    sort = np.argsort(time)

    d["path"] = dict_file["path"][sort]
    d["date-avg"] = dict_file["date-avg"][sort]
    d["dsun-obs"] = dict_file["dsun-obs"][sort]
    return d


def select_time_interval(dict_file: dict, date_start=None, date_stop=None, ):

    """ In a dict file, select only the FITS files over a specific time interval.

    """
    selection1 = np.ones(len(dict_file["date-avg"]), dtype="bool")
    selection2 = np.ones(len(dict_file["date-avg"]), dtype="bool")
    d = copy.deepcopy(dict_file)

    if date_start is not None:
        selection1 = (np.array(dict_file["date-avg"]) >= date_start)
    if date_stop is not None:
        selection2 = (np.array(dict_file["date-avg"]) <= date_stop)
    selection = selection1 & selection2

    d["path"] = dict_file["path"][selection]
    d["date-avg"] = dict_file["date-avg"][selection]
    d["dsun-obs"] = dict_file["dsun-obs"][selection]

    return d


def remove_paths_with_str(dict_file: dict, str_to_remove: str, verbose=1):
    """ 
    In a dict file, remove fits files with a specific str in their path string. 
    """    
    d = copy.deepcopy(dict_file)

    selection_to_rm = np.zeros(len(dict_file["path"]), dtype="bool")
    for ii, path in enumerate(dict_file["path"]):
        selection_to_rm[ii] = str_to_remove in os.path.basename(path)
    selection_to_keep = ~selection_to_rm

    d["path"] = dict_file["path"][selection_to_keep]
    d["date-avg"] = dict_file["date-avg"][selection_to_keep]
    d["dsun-obs"] = dict_file["dsun-obs"][selection_to_keep]
    if verbose:
        print(f"removed {selection_to_rm.sum()} files in dict")

    return d


def select_path_with_str(dict_file: dict, str_to_keep: str, additional_key: list = None) -> dict:
    """Select FITS files with a specific str in their filepath

    Args:
        dict_file (dict): input dict file
        str_to_keep (str): str to keep in their file path
        additional_key (list, optional): List of header keys to shortlist as well in dict. Defaults to None.

    Returns:
        _type_: a corrected dict
    """    
    d = copy.deepcopy(dict_file)

    selection_to_keep = np.zeros(len(dict_file["path"]), dtype="bool")
    for ii, path in enumerate(dict_file["path"]):
        selection_to_keep[ii] = str_to_keep in os.path.basename(path)

    d["path"] = dict_file["path"][selection_to_keep]
    d["date-avg"] = dict_file["date-avg"][selection_to_keep]
    d["dsun-obs"] = dict_file["dsun-obs"][selection_to_keep]
    if additional_key is not None:
        for key in additional_key:
            d[key] = dict_file[key][selection_to_keep]

    return d

def write_txt(dict_file: dict, path_to_txt: str):
    """Save the paths of all the FITS files in a dict into a text file

    Args:
        dict_file (dict): dict file to save 
        path_to_txt (str): Name of the .text file to create
    """    
    with open(path_to_txt, "w") as f:
        for ii, path in enumerate(dict_file["path"]):
            if ii != (len(dict_file["path"]) - 1):
                f.write(path + "\n")
            else:
                f.write(path)
 
def select_values_in_header_keywords(dict_file: str, window: int | str, keyword: str, values: list) -> dict:
    """Sort paths by selecting only those with keyword values in a given values list. For instance to select a list of SOOP names

    Args:
        dict_file (str): input dict file
        window (int | str): HDULIST window to use to get the header 
        keyword (str): keyword where to make the selection
        values (list): list of values that are selected.

    Returns:
        dict: output sorted dict file.
    """    
    d = copy.deepcopy(dict_file)
    selection_to_keep = np.zeros(len(dict_file["path"]), dtype="bool")

    for ii, path in enumerate(dict_file["path"]):
        with fits.open(path) as hdul:
            hdu = hdul[window]
            header = hdu.header
            val_header = header[keyword]
            val_header = val_header.replace(" ", "")
            if val_header in values:
                selection_to_keep[ii] = True

    d["path"] = dict_file["path"][selection_to_keep]
    d["date-avg"] = dict_file["date-avg"][selection_to_keep]
    d["dsun-obs"] = dict_file["dsun-obs"][selection_to_keep]

    return d

            


