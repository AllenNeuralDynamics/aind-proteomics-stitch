from . import __maintainers__, __pipeline_version__, __version__
from .utils import create_nominal_positions, utils
from .utils.xml_utils import transfer_stitching_to_multichannel

def combine_xmls():
    output_big_stitcher_xml = '/results/bigstitcher.xml'
    CAMERA_ALIGNED_XML_PATH =  "/data/stitching_cam_alignment_spot_channels.xml"
    COMBINED_XML_PATH = "/results/combined_stitching_cam_alignment_all_channels.xml"
    transfer_stitching_to_multichannel(single_channel_xml = output_big_stitcher_xml, 
    multichannel_xml = CAMERA_ALIGNED_XML_PATH, 
    output_xml = COMBINED_XML_PATH)

if __name__ == "__main__": 
    combine_xmls()