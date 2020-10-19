"""
This module defines tasks that acts as a glue between other vasp Firetasks to allow communication
between different Firetasks and Fireworks. This module also contains tasks that affect the control
flow of the workflow, e.g. tasks to check stability or the gap is within a certain range.
"""

import os
import glob
import shutil
from monty.io import zopen
from monty.os.path import zpath
from monty.shutil import gzip_dir

from pymatgen.io.cp2k.inputs import Coord, Cell
from pymatgen.io.cp2k.outputs import Cp2kOutput

from fireworks import explicit_serialize, FiretaskBase

from atomate.utils.utils import get_logger
from atomate.common.firetasks.glue_tasks import (
    get_calc_loc,
    CopyFiles,
)

logger = get_logger(__name__)

__author__ = "Nicholas Winner"
__email__ = "nwinner@berkeley.edu"


@explicit_serialize
class UpdateStructureFromPrevCalc(FiretaskBase):
    """
    Using the location of a previous calculation. The CP2K output parser will
    get the final structure from a previous calculation and update this FW's
    cp2k_input_set using it.
    """

    required_params = ["prev_calc_loc"]
    optional_params = ["cp2k_output_file", "perturb"]

    def run_task(self, fw_spec):
        calc_loc = (
            get_calc_loc(self.get("prev_calc_loc"), fw_spec["calc_locs"])
            if self.get("prev_calc_loc")
            else {}
        )

        if self.get("cp2k_output_file"):
            out = Cp2kOutput(
                zpath(
                    os.path.join(calc_loc["path"], self.get("cp2k_output_file"))
                )
            )
        else:
            out = Cp2kOutput(glob.glob(calc_loc["path"] + "/cp2k.out*")[0])

        ci = out.input
        out.parse_structures()
        struc = out.final_structure
        if self.get('perturb', False):
            struc.perturb(0.01)

        ci["FORCE_EVAL"]["SUBSYS"]["COORD"] = Coord(struc)
        ci["FORCE_EVAL"]["SUBSYS"]["CELL"] = Cell(struc.lattice)
        fw_spec["cp2k_input_set"] = ci.as_dict()


@explicit_serialize
class CopyCp2kOutputs(CopyFiles):
    """
    Copy CP2K outputs from from one location to another, unzipping them if necessary.
    Unlike VASP, which might have CONTCAR copied to POSCAR in order to continue a calculation,
    Cp2k doesn't exactly use this file system. What will generally be used for is to copy
    a .wfn file in order to restart a calculation or as an initial guess for a hybrid calculation.


    Note that you must specify either "calc_loc" or "calc_dir" to indicate
    the directory containing the previous run.

    Required params:
        (none) - but you must specify either "calc_loc" OR "calc_dir"

    Optional params:
        calc_loc (str OR bool): if True will set most recent calc_loc. If str
            search for the most recent calc_loc with the matching name
        calc_dir (str): path to dir that contains output files.
        filesystem (str): remote filesystem. e.g. username@host
    """

    optional_params = ["files_to_copy", "calc_loc", "calc_dir", "filesystem"]

    def run_task(self, fw_spec):

        calc_loc = (
            get_calc_loc(self["calc_loc"], fw_spec["calc_locs"])
            if self.get("calc_loc")
            else {}
        )

        files_to_copy = self.get("files_to_copy", [])
        # setup the copy
        self.setup_copy(
            self.get("calc_dir", None),
            filesystem=self.get("filesystem", None),
            files_to_copy=files_to_copy,
            from_path_dict=calc_loc,
        )
        # do the copying
        self.copy_files()

    def copy_files(self):
        all_files = self.fileclient.listdir(self.from_dir)
        # start file copy
        for f in self.files_to_copy:
            prev_path_full = os.path.join(self.from_dir, f)
            dest_fname = f
            dest_path = os.path.join(self.to_dir, dest_fname)

            # detect .gz extension if needed - note that monty zpath() did not seem useful here
            gz_ext = ""
            if not f in all_files:
                for possible_ext in [".gz", ".GZ"]:
                    if (f + possible_ext) in all_files:
                        gz_ext = possible_ext

            if not (f + gz_ext) in all_files:
                raise ValueError("Cannot find file: {}".format(f))

            # copy the file (minus the relaxation extension)
            try:
                self.fileclient.copy(prev_path_full + gz_ext, dest_path + gz_ext)
            except shutil.SameFileError:
                pass

            # unzip the .gz if needed
            if gz_ext in [".gz", ".GZ"]:
                # unzip dest file
                try:
                    f = zopen(dest_path + gz_ext, "rt")
                    file_content = f.read()
                except (UnicodeDecodeError, AttributeError):
                    f = zopen(dest_path + gz_ext, "rb")
                    file_content = f.read()
                if isinstance(file_content, (bytes, bytearray)):
                    with open(dest_path, "wb") as f_out:
                        f_out.write(file_content)
                else:
                    with open(dest_path, "w") as f_out:
                        f_out.write(file_content)

                f.close()
                os.remove(dest_path + gz_ext)


@explicit_serialize
class GzipDir(FiretaskBase):
    """
    Task to gzip the current directory.
    """

    required_params = []
    optional_params = []

    def run_task(self, fw_spec=None):
        cwd = os.getcwd()
        gzip_dir(cwd)


@explicit_serialize
class DeleteFiles(FiretaskBase):
    """
    Delete files
    Uses glob to search for files so any pattern it can accept can be used
    Required params:
        files: list of files to remove
    """

    required_params = ["files"]

    def run_task(self, fw_spec=None):
        cwd = os.getcwd()

        for file in self.get("files", []):
            for f in glob.glob(os.path.join(cwd, file)):
                if os.path.isdir(f):
                    shutil.rmtree(f)
                else:
                    os.remove(f)


@explicit_serialize
class DeleteFilesPrevFolder(DeleteFiles):
    """
    Can delete files, also from a previous folder in the wf if one of the optional parameters are given
    Required params:
        files: list of files to remove
    Optional params:
        calc_dir: directory to delete the files from
        calc_loc (str OR bool): if True will set most recent calc_loc. If str
             search for the most recent calc_loc with the matching name
    """

    required_params = ["files"]
    optional_params = ["calc_dir", "calc_loc"]

    def run_task(self, fw_spec=None):

        calc_dir = self.get("calc_dir", None)
        calc_loc = (
            get_calc_loc(self["calc_loc"], fw_spec["calc_locs"])
            if self.get("calc_loc")
            else {}
        )

        base_folder = os.getcwd()
        if calc_loc:
            base_folder = calc_loc["path"]
        elif calc_dir is not None:
            base_folder = calc_dir
        for file in self.get("files", []):
            for f in glob.glob(os.path.join(base_folder, file)):
                if os.path.isdir(f):
                    shutil.rmtree(f)
                else:
                    os.remove(f)
