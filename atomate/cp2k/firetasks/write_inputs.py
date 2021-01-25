"""
This module defines tasks for writing vasp input sets for various types of vasp calculations
"""

import os
from importlib import import_module
from monty.os.path import zpath
from monty.serialization import dumpfn

from pymatgen.io.cp2k.sets import Cp2kInputSet
from pymatgen.io.cp2k.outputs import Cp2kOutput
from pymatgen.alchemy.materials import TransformedStructure
from pymatgen.alchemy.transmuters import StandardTransmuter

from fireworks import FiretaskBase, explicit_serialize
from atomate.utils.utils import load_class
from atomate.common.firetasks.glue_tasks import get_calc_loc

__author__ = "Nicholas Winner"
__email__ = "nwinner@berkeley.edu"


@explicit_serialize
class WriteCp2kFromIOSet(FiretaskBase):
    """
    Create CP2K input files using the pymatgen Cp2kInputSet object or dict representation.

    Required params:
        structure (Structure): structure
        cp2k_input_set (Cp2kInputSet or dict): Either a Cp2kInputSet object or a dict representation of one

    Optional params:
        cp2k_input_params (dict): When using a string name for CP2K input set, use this as a dict
            to specify kwargs for instantiating the input set parameters. For example, if you wan
    """

    required_params = ["structure", "cp2k_input_set"]
    optional_params = ["cp2k_input_params"]

    def run_task(self, fw_spec):

        if isinstance(self["cp2k_input_set"], dict):
            self['cp2k_input_set'].update(self.get('cp2k_input_params', {}))
            cis = load_class(
                self["cp2k_input_set"]["@module"],
                self["cp2k_input_set"]["@class"]
            ).from_dict(self["cp2k_input_set"])
        elif isinstance(self['cp2k_input_set'], str):
            cis_cls = load_class(
                "pymatgen.io.cp2k.sets", self["cp2k_input_set"]
            )
            cis = cis_cls(
                self["structure"], **self.get("cp2k_input_params", {})
            )
        else:
            cis = self["cp2k_input_set"]
        cis.verbosity(False)
        cis.write_file(input_filename="cp2k.inp", output_dir=".")


@explicit_serialize
class WriteCp2kWithStrucUpdate(FiretaskBase):
    """
    Using the location of a previous calculation. The CP2K output parser will
    get the final structure from a previous calculation and update this FW's
    cp2k_input_set using it.
    """

    required_params = ["cp2k_input_set", "prev_calc_loc"]
    optional_params = ["cp2k_output_file"]

    def run_task(self, fw_spec):
        calc_loc = (
            get_calc_loc(self.get("prev_calc_loc"), fw_spec["calc_locs"])
            if self.get("prev_calc_loc")
            else {}
        )
        out = Cp2kOutput(
            zpath(
                os.path.join(
                    calc_loc["path"],
                    self.get("cp2k_output_file", 'cp2k.out'))
            )
        )

        out.parse_structures()
        struc = out.final_structure
        cis = self['cp2k_input_set']
        cis.create_subsys(struc)
        cis.write_file('cp2k.inp')


@explicit_serialize
class WriteTransmutedStructureIOSet(FiretaskBase):
    """
    Apply the provided transformations to the input structure and write the
    input set for that structure. Note that if a transformation yields many
    structures from one, only the last structure in the list is used.

    Required params:
        transformations (list): list of names of transformation classes as
            defined in the modules in pymatgen.transformations
        cp2k_input_set (Cp2kInputSet): CP2K input set.

    Optional params:
        structure (Structure): input structure if this is the first FW
        transformation_params (list): list of dicts where each dict specifies
            the input parameters to instantiate the transformation class in the
            transformations list.
        override_default_params (dict): additional user input settings.
        prev_calc_dir: path to previous calculation if using structure from
            another calculation.
    """

    required_params = ["transformations", "cp2k_input_set"]
    optional_params = [
        "structure"
        "prev_calc_loc",
        "transformation_params",
        "override_default_params",
        'cp2k_output_file'
    ]

    def run_task(self, fw_spec):

        transformations = []
        transformation_params = self.get(
            "transformation_params",
            [{} for _ in range(len(self["transformations"]))],
        )
        for t in self["transformations"]:
            found = False
            t_cls = None
            for m in [
                "advanced_transformations",
                "defect_transformations",
                "site_transformations",
                "standard_transformations",
            ]:
                mod = import_module("pymatgen.transformations.{}".format(m))

                try:
                    t_cls = getattr(mod, t)
                    found = True
                    continue
                except AttributeError:
                    pass

            if not found:
                raise ValueError("Could not find transformation: {}".format(t))

            t_obj = t_cls(**transformation_params.pop(0))
            transformations.append(t_obj)

        # If transmuter comes mid-wf, use the output of the previous FW
        if self.get('prev_calc_loc'):
            calc_loc = (
                get_calc_loc(self.get("prev_calc_loc"), fw_spec["calc_locs"])
            )
            out = Cp2kOutput(
                zpath(
                    os.path.join(
                        calc_loc["path"],
                        self.get("cp2k_output_file") if self.get("cp2k_output_file") else 'cp2k.out')
                )
            )
            out.parse_structures()
            structure = out.final_structure
        else:
            structure = self['structure']

        ts = TransformedStructure(structure)
        transmuter = StandardTransmuter([ts], transformations)
        final_structure = transmuter.transformed_structures[
            -1
        ].final_structure.copy()

        # TODO Temporary way to deal with perturb, which removes charges
        for ts in transformations:
            if ts in [
                "PerturbStructureTransformation",
                "DeformStructureTransformation",
                "ConventionalCellTransformation",
                "PrimitiveCellTransformation",
                "SupercellTransformation"

            ]:
                transformation_is_charge_preserving = True
            else:
                transformation_is_charge_preserving = False

            if self.get('structure'):
                if transformation_is_charge_preserving:
                    final_structure.set_charge(self['structure'].charge)

        cis_orig = self["cp2k_input_set"]
        cis_dict = cis_orig.as_dict()
        cis_dict["structure"] = final_structure.as_dict()
        cis_dict.update(self.get("override_default_params", {}) or {})
        cis = cis_orig.__class__.from_dict(cis_dict)
        cis.verbosity(False)
        cis.write_file("cp2k.inp")

        dumpfn(transmuter.transformed_structures[-1], "transformations.json")
