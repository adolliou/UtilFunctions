import numpy as np
import os
from pathlib import Path
from datetime import datetime


class FolderManager(dict):
    def __init__(self, dict_input, list_needed_keys: list):
        self._assert_correct_dict(dict_input, list_needed_keys)

    @staticmethod
    def _assert_correct_dict(dict_input, list_needed_keys: list):
        assertion = np.array([(n in dict_input) for n in list_needed_keys], dtype="bool")
        if assertion.sul != len(assertion):
            raise ValueError("Not correct dict keys in the input")


class OutputFolderManager(FolderManager):
    def __init__(self, dict_input):
        list_needed_keys = ["data_folder", "sequence_folder_name",
                            "old_folder", "old_subfolder", "old_level",
                            "new_folder", "new_subfolder", "new_level"]
        super().__init__(dict_input, list_needed_keys)
        self._initialize_output_folder(dict_input)
        self._initialize_input_folder(dict_input)

    def _initialize_input_folder(self, dict_input):
        path_spice_old_level = os.path.join(dict_input["data_folder"], dict_input["sequence_folder_name"],
                                            dict_input["old_folder"], dict_input["old_subfolder"])
        self["old"] = {"path": path_spice_old_level,
                       "level": dict_input["old_level"]}

    def _initialize_output_folder(self, dict_input):
        path_spice_new_level = os.path.join(dict_input["data_folder"], dict_input["sequence_folder_name"],
                                            dict_input["new_folder"], dict_input["new_subfolder"], )

        Path(os.path.join(dict_input["data_folder"], dict_input["sequence_folder_name"],
                          dict_input["new_folder"])).mkdir(parents=False, exist_ok=True)

        Path(path_spice_new_level).mkdir(parents=False, exist_ok=True)
        self["new"] = {"path": path_spice_new_level,
                       "level": dict_input["new_level"]}


class ResultFolderManager(FolderManager):
    def __init__(self, dict_input):
        list_needed_keys = ["results_folder", "sequence_folder_name", "name_function", ]
        super().__init__(dict_input, list_needed_keys)

    def _initialise_result_folder(self, dict_input):
        results_folder = os.path.join(dict_input["results_folder"])
        Path(results_folder).mkdir(parents=False, exist_ok=True)

        results_folder = os.path.join(results_folder, dict_input["sequence_folder_name"])

        results_folder = os.path.join(results_folder, dict_input['alignement_pixels'])
        Path(results_folder).mkdir(parents=False, exist_ok=True)

        now = datetime.now()
        dt_string = now.strftime("D_%d_%m_%Y_T_%H_%M_%S")
        path_save = os.path.join(results_folder, dt_string)
        Path(path_save).mkdir(parents=False, exist_ok=True)

        self["res"] = {"path": path_save}
