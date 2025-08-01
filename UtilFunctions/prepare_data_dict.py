from .create_dict_file import create_dict_file, select_time_interval
from UtilFunctions.folder_manager import ResultFolderManager, InputFolderManager, OutputFolderManager


def prepare_data_dict(files, data_folder,sequence_folder_name, results_folder=None):

    files["data_folder"] = data_folder
    files["sequence_folder_name"] = sequence_folder_name
    params_results = {
        "results_folder": results_folder,
        "sequence_folder_name": sequence_folder_name,
        "name_function": "write_tables"
    }
    folderman = InputFolderManager(files)
    if results_folder is not None:
        folderman["res"] = ResultFolderManager(params_results)["res"]
    if "out_folder" in files:
        folderman["out"] = OutputFolderManager(files)["out"]


    dict_files = create_dict_file(path_instrument=folderman["in"]["path"], suffix=folderman["in"]["suffix"],
                                     window=folderman["in"]["window"], name_list_txt = folderman["in"]["name_list_txt"])
    dict_files = select_time_interval(dict_file=dict_files,
                                         date_start=files["date_start"], date_stop=files["date_stop"])