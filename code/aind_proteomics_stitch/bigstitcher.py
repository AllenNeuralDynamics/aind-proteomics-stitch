"""
Computes stitching transformations using
bigstitcher for proteomics data structure
"""

import json
from pathlib import Path
from time import time
from typing import List, Optional, Tuple

from aind_data_schema.core.processing import DataProcess, ProcessName

from . import (__maintainers__, __pipeline_version__, __version__,
               bigstitcher_utilities)
from .utils import utils


def validate_capsule_inputs(input_elements: List[str]) -> List[str]:
    """
    Validates input elemts for a capsule in
    Code Ocean.

    Parameters
    -----------
    input_elements: List[str]
        Input elements for the capsule. This
        could be sets of files or folders.

    Returns
    -----------
    List[str]
        List of missing files
    """

    missing_inputs = []
    for required_input_element in input_elements:
        required_input_element = Path(required_input_element)

        if not required_input_element.exists():
            missing_inputs.append(str(required_input_element))

    return missing_inputs


def get_data_config(
    data_folder: str,
    processing_manifest_path: str = "processing_manifest.json",
    data_description_path: str = "data_description.json",
    acquisition_path: str = "acquisition.json",
) -> Tuple:
    """
    Returns the first proteomics dataset found
    in the data folder

    Parameters
    -----------
    data_folder: str
        Path to the folder that contains the data

    processing_manifest_path: str
        Path for the processing manifest

    data_description_path: str
        Path for the data description

    Returns
    -----------
    Tuple[Dict, str]
        Dict: Empty dictionary if the path does not exist,
        dictionary with the data otherwise.

        Str: Empty string if the processing manifest
        was not found
    """

    # Returning first proteomics dataset found
    # Doing this because of Code Ocean, ideally we would have
    # a single dataset in the pipeline

    processing_manifest_path = Path(f"{data_folder}/{processing_manifest_path}")
    data_description_path = Path(f"{data_folder}/{data_description_path}")

    if not processing_manifest_path.exists():
        raise ValueError(
            f"Please, check processing manifest path: {processing_manifest_path}"
        )

    if not data_description_path.exists():
        raise ValueError(
            f"Please, check data description path: {data_description_path}"
        )

    derivatives_dict = utils.read_json_as_dict(str(processing_manifest_path))
    data_description_dict = utils.read_json_as_dict(str(data_description_path))
    acquisition_dict = utils.read_json_as_dict(f"{data_folder}/{acquisition_path}")

    proteomics_dataset = data_description_dict["name"]

    return derivatives_dict, proteomics_dataset, acquisition_dict


def get_stitching_dict(
    specimen_id: str, dataset_xml_path: str, downsample: Optional[int] = 2
) -> dict:
    """
    A function that writes a stitching dictioonary that will be used for
    creating a json file that gives parmaters to bigstitcher sittching run

    Parameters
    ----------
    specimen_id: str
        Specimen ID
    dataset_xml_path: str
        Path where the xml is located
    downsample: Optional[int] = 2
        Image multiscale used for stitching

    Returns
    -------
    dict
        Dictionary with the stitching parameters
        used for bigstitcher
    """
    # assert pathlib.Path(dataset_xml_path).exists()

    stitching_dict = {
        "session_id": str(specimen_id),
        "memgb": 100,
        "parallel": utils.get_code_ocean_cpu_limit(),
        "dataset_xml": str(dataset_xml_path),
        "do_phase_correlation": True,
        "do_detection": False,
        "do_registrations": False,
        "phase_correlation_params": {
            "downsample": downsample,
            "min_correlation": 0.6,
            "max_shift_in_x": 100,
            "max_shift_in_y": 100,
            "max_shift_in_z": 100,
        },
    }
    return stitching_dict


def get_estimated_downsample(
    voxel_resolution: List[float], phase_corr_res: Tuple[float] = (8.0, 8.0, 4.0)
) -> int:
    """
    Get the estimated multiscale based on the provided
    voxel resolution. This is used for image stitching.

    e.g., if the original resolution is (1.8. 1.8, 2.0)
    in XYZ order, and you provide (3.6, 3.6, 4.0) as
    image resolution, then the picked resolution will be
    1.

    Parameters
    ----------
    voxel_resolution: List[float]
        Image original resolution. This would be the resolution
        in the multiscale "0".
    phase_corr_res: Tuple[float]
        Approximated resolution that will be used for bigstitcher
        in the computation of the transforms. Default: (8.0, 8.0, 4.0)
    """

    downsample_versions = []
    for idx in range(len(voxel_resolution)):
        downsample_versions.append(phase_corr_res[idx] // voxel_resolution[idx])

    downsample_res = int(min(downsample_versions) - 1)
    return downsample_res


def main(
    data_folder,
    channel_wavelength,
    path_to_tile_metadata,
    voxel_resolution,
    output_json_file,
    results_folder,
    proteomics_dataset_name,
    res_for_transforms=(0.19, 0.19, 0.85),
):
    """
    Computes image stitching with BigStitcher using Phase Correlation

    Parameters
    ----------
    channel_wavelength: str
        Channel wavelength
    voxel_resolution: Tuple[float]
        Voxel resolution in order XYZ
    output_json_file: str
        Path where the json file will be written
    results_folder: Path
        Results folder
    proteomics_dataset_name: str
        Proteomics dataset name
    """
    start_time = time()
    metadata_folder = results_folder.joinpath("metadata")
    utils.create_folder(str(metadata_folder))

    tile_metadata = utils.read_json_as_dict(path_to_tile_metadata)
    channel_metadata = []

    for t in tile_metadata:
        t["file"] = f"{Path(t['file']).stem}.ome.zarr"
        if int(channel_wavelength) == int(t["channel_wavelength"]):
            channel_metadata.append(t)

    utils.save_dict_as_json(filename=output_json_file, dictionary=channel_metadata)

    tree = bigstitcher_utility.parse_json(
        output_json_file, str(data_folder), microns=True
    )
    output_big_stitcher_xml = f"/{results_folder}/{proteomics_dataset_name}_stitching_channel_{channel_wavelength}.xml"

    bigstitcher_utility.write_xml(tree, output_big_stitcher_xml)

    estimated_downsample = get_estimated_downsample(
        voxel_resolution=voxel_resolution, phase_corr_res=res_for_transforms
    )

    # print(f"Voxel resolution: {voxel_resolution} - Estimating transforms in res: {res_for_transforms} - Scale: {estimated_downsample}")

    proteomics_stitching_params = get_stitching_dict(
        specimen_id=proteomics_dataset_name,
        dataset_xml_path=output_big_stitcher_xml,
        downsample=estimated_downsample,
    )
    end_time = time()

    output_big_stitcher_json = f"{results_folder}/{proteomics_dataset_name}_stitch_channel_{channel_wavelength}_params.json"

    data_processes = []
    data_processes.append(
        DataProcess(
            name=ProcessName.IMAGE_TILE_ALIGNMENT,
            software_version="1.2.11",
            start_date_time=start_time,
            end_date_time=end_time,
            input_location=str(proteomics_dataset_name),
            output_location=str(output_big_stitcher_json),
            outputs={"output_file": str(output_big_stitcher_json)},
            code_url="",
            code_version=__version__,
            parameters=proteomics_stitching_params,
            notes="Creation of stitching parameters",
        )
    )

    utils.generate_processing(
        data_processes=data_processes,
        dest_processing=metadata_folder,
        processor_full_name=__maintainers__[0],
        pipeline_version=__pipeline_version__,
    )

    with open(output_big_stitcher_json, "w") as f:
        json.dump(proteomics_stitching_params, f, indent=4)

    # Printing to get output on batch script
    print(output_big_stitcher_json)


if __name__ == "__main__":
    main()
