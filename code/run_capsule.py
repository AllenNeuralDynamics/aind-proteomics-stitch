"""
Run the capsule
This script is used to run the capsule for stitching images using BigStitcher.
"""

import os
from pathlib import Path

from aind_proteomics_stitch import bigstitcher
from aind_proteomics_stitch.utils import utils


def run():
    """Function that runs image stitching with BigStitcher"""
    data_folder = Path(os.path.abspath("../data/HCR_785631_2025-04-19_02-00-00"))
    results_folder = Path("../results")  # os.path.relpath(
    # scratch_folder = Path(os.path.abspath("../scratch"))

    # It is assumed that these files
    # will be in the data folder

    required_input_elements = [
        f"{data_folder}/SPIM/derivatives/processing_manifest.json",
        f"{data_folder}/SPIM/derivatives/all_channel_tile_metadata.json",
        f"{data_folder}/data_description.json",
        f"{data_folder}/acquisition.json",
    ]

    missing_files = utils.validate_capsule_inputs(required_input_elements)

    if len(missing_files):
        raise ValueError(
            f"We miss the following files in the capsule input: {missing_files}"
        )

    pipeline_config, proteomics_dataset_name, acquisition_dict = utils.get_data_config(
        data_folder=data_folder,
        processing_manifest_path="SPIM/derivatives/processing_manifest.json",
        data_description_path="data_description.json",
        acquisition_path="acquisition.json",
    )

    voxel_resolution = utils.get_resolution(acquisition_dict)
    stitching_channel = pipeline_config["pipeline_processing"]["stitching"]["channel"]

    output_json_file = results_folder.joinpath(
        f"{proteomics_dataset_name}_tile_metadata.json"
    )

    # Computing image transformations with bigtstitcher
    path_to_tile_metadata = required_input_elements[1]
    relative_data_folder = os.path.relpath(
        "/data/HCR_785631_2025-04-19_02-radially-corrected"
    )

    bigstitcher.proteomics_main(
        data_folder=relative_data_folder,
        channel_wavelength=stitching_channel,
        path_to_tile_metadata=path_to_tile_metadata,
        voxel_resolution=voxel_resolution,
        output_json_file=output_json_file,
        results_folder=results_folder,
        proteomics_dataset_name=proteomics_dataset_name,
        res_for_transforms=(0.19, 0.19, 0.85),
    )


if __name__ == "__main__":
    run()
