from glob import glob
import numpy as np
from astropy.io import fits
import astropy.time
import os
import astropy.units as u


def create_dict_file(path_instrument: str, suffix: str, window=None, sort_dict=True):
    data_dict = {}

    paths = glob(os.path.join(path_instrument, suffix))
    data_dict['path_instrument'] = path_instrument
    data_dict['path'] = paths
    data_dict['date-avg'] = []
    data_dict['dsun-obs'] = []
    data_dict['telescop'] = []

    for kk, path in enumerate(data_dict['path']):
        f = fits.open(path)
        if window is None:
            if len(f) > 1:
                idx = 1
            else:
                idx = 0
        else:
            idx = window
        data_dict['date-avg'].append(astropy.time.Time(f[idx].header['DATE-AVG']))
        data_dict['dsun-obs'].append(f[idx].header['DSUN_OBS'])
        data_dict['telescop'].append(f[idx].header['TELESCOP'])

    if sort_dict:
        return _sort_dict_file(data_dict)
    else:
        return data_dict


def _sort_dict_file(dict_file: dict):
    ref_time = dict_file["date-avg"][0]
    time = [(n - ref_time).to(u.s).value for n in dict_file["date-avg"]]
    sort = np.argsort(time)
    dict_file["path"] = dict_file["path"][sort]
    dict_file["date-avg"] = dict_file["date-avg"][sort]
    dict_file["dsun-obs"] = dict_file["dsun-obs"][sort]

    dict_file["telescop"] = dict_file["telescop"][sort]
    return dict_file


def select_time_interval(dict_file: dict, date_start=None, date_stop=None):
    selection1 = np.ones(len(dict_file["date-avg"]), dtype="bool")
    selection2 = np.ones(len(dict_file["date-avg"]), dtype="bool")

    if date_start is not None:
        selection1 = (np.array(dict_file["date-avg"]) >= date_start)
    if date_stop is not None:
        selection2 = (np.array(dict_file["date-avg"]) <= date_stop)
    selection = selection1 & selection2

    dict_file["path"] = dict_file["path"][selection]
    dict_file["date-avg"] = dict_file["date-avg"][selection]
    dict_file["dsun-obs"] = dict_file["dsun-obs"][selection]
    dict_file["telescop"] = dict_file["telescop"][selection]

    return dict_file
