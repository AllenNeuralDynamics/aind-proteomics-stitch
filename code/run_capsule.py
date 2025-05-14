"""
Run the capsule
This script is used to run the capsule for stitching images using BigStitcher.
"""

from pathlib import Path

from aind_proteomics_stitch import bigstitcher
from aind_proteomics_stitch.utils import utils
import subprocess

def run():
    """Function that runs image stitching with BigStitcher"""
    data_folder = Path("../data")
    results_folder = Path("../results")  # os.path.relpath(
    # scratch_folder = Path(os.path.abspath("../scratch"))
    
    data_folder = Path(data_folder)
    #data_in_data_folder = [str(l) for l in list(data_folder.glob("*"))]
    #print("Data in data folder: ", data_in_data_folder)
    #result = subprocess.run(['ls', '-al'], capture_output=True, text=True, cwd=data_folder)
    #print("ls -al command: ", result.stdout)
    
    # It is assumed that these files
    # will be in the data folder

    required_input_elements = [
        f"{data_folder}/processing_manifest.json",
        f"{data_folder}/all_channel_tile_metadata.json",
        f"{data_folder}/data_description.json",
        f"{data_folder}/acquisition.json",
        f"{data_folder}/processed_data_description.json",
        f"{data_folder}/radial_correction_parameters.json",
    ]

    missing_files = utils.validate_capsule_inputs(required_input_elements)

    if len(missing_files):
        raise ValueError(
            f"We miss the following files in the capsule input: {missing_files}"
        )

    pipeline_config, proteomics_dataset_name, acquisition_dict = utils.get_data_config(
        data_folder=data_folder,
        processing_manifest_path="processing_manifest.json",
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
    # print("Contents data folder: ", list(data_folder.glob("*")))

    processed_data_description = utils.read_json_as_dict(required_input_elements[4])
    radial_parameters = utils.read_json_as_dict(required_input_elements[5])

    processed_asset_name = processed_data_description.get("name", None)
    bucket_name = radial_parameters.get('bucket_name', None)

    if processed_asset_name is None or bucket_name is None:
        raise ValueError("Stitching requires S3 paths in Code Ocean at the moment.")

    path_to_data = f"s3://{bucket_name}/{processed_asset_name}/image_radial_correction"

    bigstitcher.main(
        path_to_data=path_to_data,
        channel_wavelength=stitching_channel,
        path_to_tile_metadata=path_to_tile_metadata,
        voxel_resolution=voxel_resolution,
        output_json_file=output_json_file,
        results_folder=results_folder,
        proteomics_dataset_name=proteomics_dataset_name,
        res_for_transforms=(0.76, 0.76, 3.4),
        scale_for_transforms=2,
        # If this is provided, res for
        # transforms is ignored
    )


if __name__ == "__main__":
    run()
