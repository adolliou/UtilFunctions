import numpy as np
import os
from pathlib import Path
from datetime import datetime


class FolderManager(dict):
    def __init__(self, dict_input, list_needed_keys: list):
        super().__init__()
        self._assert_correct_dict(dict_input, list_needed_keys)

    @staticmethod
    def _assert_correct_dict(dict_input, list_needed_keys: list):
        assertion = np.array([(n in dict_input) for n in list_needed_keys], dtype="bool")
        if assertion.sum() != len(list_needed_keys):
            print(f"{dict_input=}")
            print(f"{list_needed_keys=}")
            raise ValueError("Not correct dict keys in the input")


class InputFolderManager(FolderManager):
    def __init__(self, dict_input):
        list_needed_keys = ["data_folder", "sequence_folder_name",
                            "in_folder", "in_subfolder", "in_level"]
        super().__init__(dict_input, list_needed_keys)
        self._initialize_input_folder(dict_input)

    def _initialize_input_folder(self, dict_input):
        path_spice_old_level = os.path.join(dict_input["data_folder"], dict_input["sequence_folder_name"],
                                            dict_input["in_folder"], dict_input["in_subfolder"])
        self["in"] = {"path": path_spice_old_level,
                      "level": dict_input["in_level"]}
        if "in_suffix" in dict_input:
            self["in"]["suffix"] = dict_input["in_suffix"]
        if "in_window" in dict_input:
            self["in"]["window"] = dict_input["in_window"]
        if "in_read_list" in dict_input:
            self["in"]["read_list"] = dict_input["in_read_list"]
        else:
            self["in"]["read_list"] = False


class OutputFolderManager(FolderManager):
    def __init__(self, dict_input):
        list_needed_keys = ["data_folder", "sequence_folder_name",
                            "out_folder", "out_subfolder", "out_level"]
        super().__init__(dict_input, list_needed_keys)
        self._initialize_output_folder(dict_input)

    def _initialize_output_folder(self, dict_input):
        path_spice_new_level = os.path.join(dict_input["data_folder"], dict_input["sequence_folder_name"],
                                            dict_input["out_folder"], dict_input["out_subfolder"], )

        Path(os.path.join(dict_input["data_folder"], dict_input["sequence_folder_name"],
                          dict_input["out_folder"])).mkdir(parents=False, exist_ok=True)

        Path(path_spice_new_level).mkdir(parents=False, exist_ok=True)
        self["out"] = {"path": path_spice_new_level,
                       "level": dict_input["out_level"]}
        if "out_window" in dict_input:
            self["out"]["window"] = dict_input["out_window"]



class ResultFolderManager(FolderManager):
    def __init__(self, dict_input):
        list_needed_keys = ["results_folder", "sequence_folder_name", "name_function", ]
        super().__init__(dict_input, list_needed_keys)
        self._initialise_result_folder(dict_input)

    def _initialise_result_folder(self, dict_input):
        results_folder = os.path.join(dict_input["results_folder"])
        Path(results_folder).mkdir(parents=False, exist_ok=True)

        results_folder = os.path.join(results_folder, dict_input["sequence_folder_name"])
        Path(results_folder).mkdir(parents=False, exist_ok=True)

        results_folder = os.path.join(results_folder, dict_input['name_function'])
        Path(results_folder).mkdir(parents=False, exist_ok=True)
        if "personalize_name" in dict_input:
            path_save = os.path.join(results_folder, dict_input["personalize_name"])
            Path(path_save).mkdir(parents=False, exist_ok=True)
        else:
            now = datetime.now()
            dt_string = now.strftime("D_%d_%m_%Y_T_%H_%M_%S")
            path_save = os.path.join(results_folder, dt_string)
            Path(path_save).mkdir(parents=False, exist_ok=True)

        self["res"] = {"path": path_save}
