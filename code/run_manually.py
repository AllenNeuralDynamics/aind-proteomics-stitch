
from aind_proteomics_stitch import bigstitcher
import json
from pathlib import Path
def run_offpipeline(): 
    scale_for_transforms = 2
    results_folder = '/results'
    channel_wavelength = 405
    proteomics_dataset_name = f'HCR_807074_2025-08-27_11-45-00_processed_2025-10-02_06-30-17'
    # xml_path = f'/root/capsule/data/{proteomics_dataset_name}_bigstitcher_pc.xml'

    
    output_big_stitcher_xml = f"{results_folder}/{proteomics_dataset_name}_stitching_channel_{channel_wavelength}.xml"

    # Creating XML
    # bigstitcher.create_nominal_positions.create_xml_from_acquisition(
    #     acquisition_json_path=acquisition_path,
    #     output_xml_path=output_big_stitcher_xml,
    #     zarr_base_path=path_to_data,
    #     stitching_channel=channel_wavelength,
    # )

    # zarr_path_xml = tree.find("SequenceDescription").find("ImageLoader").find("zarr")
    # if not zarr_path_xml.text.startswith("s3://"):
    #     zarr_path_xml.text = os.path.abspath(zarr_path_xml.text)

    # if scale_for_transforms is None:
    #     scale_for_transforms = get_estimated_downsample(
    #         voxel_resolution=voxel_resolution, phase_corr_res=res_for_transforms
    #     )

    # scale_for_transforms = int(scale_for_transforms)

    # proteomics_stitching_params = bigstitcher.get_stitching_dict(
    #     specimen_id=proteomics_dataset_name,
    #     dataset_xml_path=output_big_stitcher_xml,
    #     downsample=scale_for_transforms,
    # )
    # output_big_stitcher_json = f"{results_folder}/{proteomics_dataset_name}_stitch_channel_{channel_wavelength}_params.json"
    # with open(output_big_stitcher_json, "w") as f:
    #     json.dump(proteomics_stitching_params, f, indent=4)

    # # Printing to get output on batch script
    # print(output_big_stitcher_json)

def run_bigstitcher(): 
        processed_asset_name = f'HCR_807074_2025-08-26_15-45-00_processed_2025-10-02_22-49-01'
        stitching_channel=405
        path_to_data = f"s3://aind-open-data/{processed_asset_name}/image_radial_correction"
        voxel_resolution = (1, 0.3880046677791278, 0.3880046677791278)
        results_folder = Path('/results')
        acquisition_path = f"/data/{processed_asset_name}/acquisition.json"

        bigstitcher.main(
        path_to_data=path_to_data,
        channel_wavelength=stitching_channel,
        acquisition_path=acquisition_path,
        voxel_resolution=voxel_resolution,
        results_folder=results_folder,
        proteomics_dataset_name=processed_asset_name,
        res_for_transforms=(0.76, 0.76, 3.4),
        scale_for_transforms=2,
        # If this is provided, res for
        # transforms is ignored
    )

def combine_all_xmls():
    from aind_proteomics_stitch.utils.xml_utils import transfer_stitching_to_multichannel, split_multichannel_xml
    proteomics_datasets = [
        "HCR_807074_2025-07-30_13-45-00_processed_2025-10-02_06-30-15", 
        "HCR_807074_2025-08-26_15-45-00_processed_2025-10-02_22-49-01",
        "HCR_807074_2025-08-27_11-45-00_processed_2025-10-02_06-30-17"
    ]
    for dataset in proteomics_datasets:
        output_big_stitcher_xml = f'../data/{dataset}_bigstitcher_pc_405.xml'
        CAMERA_ALIGNED_XML_PATH =  f"../data/{dataset}/image_tile_alignment/stitching_cam_alignment_forward_transform_spot_channels.xml"
        COMBINED_XML_PATH = f"../results/{dataset}_combined_stitching_pc_405_cam_alignment_all_channels.xml"
        try: 
            transfer_stitching_to_multichannel(single_channel_xml = output_big_stitcher_xml, 
            multichannel_xml = CAMERA_ALIGNED_XML_PATH, 
            output_xml = COMBINED_XML_PATH)
        except: 
            print(f'Error combining xmls')

if __name__ == "__main__":
    run_bigstitcher()
    # combine_all_xmls()