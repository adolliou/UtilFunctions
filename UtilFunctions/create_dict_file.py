from glob import glob
import numpy as np
from astropy.io import fits
import astropy.time
import os
import astropy.units as u
from tqdm import tqdm
import copy
import warnings


def create_dict_file(path_instrument: str, suffix: str, window=None, sort_dict=True):
    data_dict = {}

    paths = glob(os.path.join(path_instrument, suffix))
    print("create dictionary file with files in the folder %s " % path_instrument)
    data_dict['path_instrument'] = path_instrument
    data_dict['path'] = paths
    data_dict['date-avg'] = []
    # data_dict['date-beg'] = []
    data_dict['dsun-obs'] = []
    # data_dict['telescop'] = []

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
        # if "DATE-BEG" in f[idx].header:
        #     data_dict['date-beg'].append(astropy.time.Time(f[idx].header['DATE-BEG']))
        # elif:
        #     data_dict['date-beg'].append(astropy.time.Time(f[idx].header['DATE-OBS']))


        data_dict['dsun-obs'].append(f[idx].header['DSUN_OBS'])
        # data_dict['telescop'].append(f[idx].header['TELESCOP'])
        if "DATE-AVG" not in f[idx].header:
            warnings.warn("DATE-AVG not found in header, manually compute it.")
            if "EXPTIME" in f[idx].header:
                data_dict['date-avg'].append(astropy.time.Time(f[idx].header['DATE-OBS']) +
                                         0.5*u.Quantity(f[idx].header["EXPTIME"], "s"))
            elif "CADENCE" in f[idx].header:
                print("create HMI date-avg")
                data_dict['date-avg'].append(astropy.time.Time(f[idx].header['DATE-OBS']) +
                                             0.5 * u.Quantity(f[idx].header["CADENCE"], "s"))
            else:
                data_dict['date-avg'].append(astropy.time.Time(f[idx].header['DATE-OBS']))
    # print(data_dict["date-avg"])
    data_dict['path'] = np.array(data_dict['path'])
    data_dict['date-avg'] = np.array(data_dict['date-avg'])
    # data_dict['date-beg'] = np.array(data_dict['date-beg'])
    data_dict['dsun-obs'] = np.array(data_dict['dsun-obs'])
    # data_dict['telescop'] = np.array(data_dict['telescop'])

    print("%i FITS files added in the dict_file." % len(data_dict["path"]))

    if sort_dict:
        return _sort_dict_file(data_dict)
    else:
        return data_dict


def _sort_dict_file(dict_file: dict):
    ref_time = dict_file["date-avg"][0]
    d = copy.deepcopy(dict_file)
    time = [(n - ref_time).to(u.s).value for n in dict_file["date-avg"]]
    sort = np.argsort(time)

    d["path"] = dict_file["path"][sort]
    d["date-avg"] = dict_file["date-avg"][sort]
    # d["date-beg"] = dict_file["date-beg"][sort]
    d["dsun-obs"] = dict_file["dsun-obs"][sort]
    # d["telescop"] = dict_file["telescop"][sort]
    return d


def select_time_interval(dict_file: dict, date_start=None, date_stop=None, ):
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
    # d["date-beg"] = dict_file["date-beg"][selection]
    d["dsun-obs"] = dict_file["dsun-obs"][selection]
    # d["telescop"] = dict_file["telescop"][selection]

    return d


def remove_paths_with_str(dict_file: dict, str_to_remove: str):
    d = copy.deepcopy(dict_file)

    selection_to_rm = np.zeros(len(dict_file["path"]), dtype="bool")
    for ii, path in enumerate(dict_file["path"]):
        selection_to_rm[ii] = str_to_remove in os.path.basename(path)
    selection_to_keep = ~selection_to_rm

    d["path"] = dict_file["path"][selection_to_keep]
    d["date-avg"] = dict_file["date-avg"][selection_to_keep]
    # d["date-beg"] = dict_file["date-beg"][selection_to_keep]
    d["dsun-obs"] = dict_file["dsun-obs"][selection_to_keep]
    # d["telescop"] = dict_file["telescop"][selection_to_keep]
    print(f"removed {selection_to_rm.sum()} files in dict")

    return d
