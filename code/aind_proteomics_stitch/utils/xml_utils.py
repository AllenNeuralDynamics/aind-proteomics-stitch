import xml.etree.ElementTree as ET


# The approach we will be taking is collecting the affine transformations from the single -channel stitching results,
# confirming the tile identity, and applying it to the appropriate tiles in the multichannel camera-corrected xml. 
# the multichannel camera aligned xml will then be the main xml used for qc. We can also include a utility to break it up into 
# individual channels if that would be helpful. 

import boto3
import re
import json
import xmltodict
from collections import OrderedDict
from typing import Optional, Dict, List, Tuple, Union
import xml.etree.ElementTree as ET
from xml.dom import minidom
import copy
from pathlib import Path


class BigStitcherXMLManager:
    """
    Manager class for handling BigStitcher XML transformations between 
    single-channel stitched and multichannel camera-aligned datasets.
    """
    
    def __init__(self):
        self.s3_client = boto3.client('s3')
    
    def load_xml(self, xml_path: str) -> OrderedDict:
        """
        Load XML from local file or S3 path.
        
        Parameters
        ----------
        xml_path : str
            Path to XML file. Can be local path or S3 path (s3://bucket/path/file.xml).
            
        Returns
        -------
        OrderedDict
            Parsed XML content
        """
        if xml_path.startswith('s3://'):
            # Parse S3 path
            s3_path = xml_path.replace("s3://", "")
            bucket_name, key = s3_path.split("/", 1)
            
            # Read from S3
            try:
                response = self.s3_client.get_object(Bucket=bucket_name, Key=key)
                xml_content = response['Body'].read().decode('utf-8')
            except Exception as e:
                raise Exception(f"Could not read XML file from S3: {e}")
        else:
            # Read from local file
            with open(xml_path, "r") as file:
                xml_content = file.read()
        
        # Parse XML content
        data = xmltodict.parse(xml_content)
        return data
    
    def save_xml(self, data: OrderedDict, output_path: str, pretty: bool = True) -> None:
        """
        Save XML data to file or S3.
        
        Parameters
        ----------
        data : OrderedDict
            XML data to save
        output_path : str
            Output path (local or S3)
        pretty : bool
            Whether to format XML with indentation
        """
        # Convert to XML string
        xml_str = xmltodict.unparse(data, pretty=pretty, indent='  ')
        
        if output_path.startswith('s3://'):
            # Save to S3
            s3_path = output_path.replace("s3://", "")
            bucket_name, key = s3_path.split("/", 1)
            
            try:
                self.s3_client.put_object(
                    Bucket=bucket_name,
                    Key=key,
                    Body=xml_str.encode('utf-8'),
                    ContentType='application/xml'
                )
                print(f"Saved XML to S3: {output_path}")
            except Exception as e:
                raise Exception(f"Could not save XML to S3: {e}")
        else:
            # Save to local file
            Path(output_path).parent.mkdir(exist_ok=True, parents=True)

            with open(output_path, 'w') as f:
                f.write(xml_str)
            print(f"Saved XML to local file: {output_path}")
    
    def get_tile_id_from_name(self, data: dict, tilename: str) -> int:
        """
        Get tile ID from tile name in XML data.
        
        Parameters
        ----------
        data : dict
            Parsed XML data
        tilename : str
            Name of the tile to find
            
        Returns
        -------
        int
            Tile ID number
        """
        viewsetups = data["SpimData"]["SequenceDescription"]["ViewSetups"]["ViewSetup"]
        
        # Ensure viewsetups is a list
        if not isinstance(viewsetups, list):
            viewsetups = [viewsetups]
        
        matching_viewsetup = [v for v in viewsetups if tilename in v.get('name', '')]
        
        if not matching_viewsetup:
            raise ValueError(f"No ViewSetup found for tile name: {tilename}")
        
        # Get tile_number from this viewsetup
        matching_tile_number = matching_viewsetup[0]['attributes']['tile']
        
        return int(matching_tile_number)
    
    def extract_tile_name_from_viewsetup(self, viewsetup: dict) -> str:
        """
        Extract tile name from ViewSetup entry.
        
        Parameters
        ----------
        viewsetup : dict
            ViewSetup dictionary from XML
            
        Returns
        -------
        str
            Extracted tile name
        """
        name = viewsetup.get('name', '')
        # Remove channel information if present
        # Common patterns: Tile_X_0001_Y_0002_Z_0003_ch_488
        pattern = r"(Tile_X_\d+_Y_\d+_Z_\d+)"
        match = re.search(pattern, name)
        if match:
            return match.group(1)
        return name
    
    def get_channel_from_viewsetup(self, viewsetup: dict) -> Optional[int]:
        """
        Extract channel wavelength from ViewSetup.
        
        Parameters
        ----------
        viewsetup : dict
            ViewSetup dictionary from XML
            
        Returns
        -------
        int or None
            Channel wavelength if found
        """
        name = viewsetup.get('name', '')
        pattern = r"(ch|CH)_(\d+)"
        match = re.search(pattern, name)
        if match:
            return int(match.group(2))
        return None
    
    def extract_transforms_for_tile(self, data: dict, tile_id: int) -> List[dict]:
        """
        Extract all transforms for a specific tile ID.
        
        Parameters
        ----------
        data : dict
            Parsed XML data
        tile_id : int
            Tile ID to extract transforms for
            
        Returns
        -------
        list
            List of transform dictionaries
        """
        view_registrations = data["SpimData"]["ViewRegistrations"]["ViewRegistration"]
        
        if not isinstance(view_registrations, list):
            view_registrations = [view_registrations]
        
        for view_reg in view_registrations:
            if int(view_reg.get("@setup", -1)) == tile_id:
                transforms = view_reg.get("ViewTransform", [])
                if not isinstance(transforms, list):
                    transforms = [transforms]
                return transforms
        
        return []
    
    def transfer_stitching_transforms(
        self,
        single_channel_xml_path: str,
        multichannel_xml_path: str,
        output_xml_path: str,
        transform_index: int = 1,
        channels_to_apply: Optional[List[int]] = None
    ) -> None:
        """
        Transfer stitching transforms from single-channel XML to multichannel XML.
        
        Parameters
        ----------
        single_channel_xml_path : str
            Path to single-channel stitched XML
        multichannel_xml_path : str
            Path to multichannel camera-aligned XML
        output_xml_path : str
            Path for output XML with combined transforms
        transform_index : int
            Which transform to extract from single-channel (default: 1 for second transform)
        channels_to_apply : list, optional
            List of channel wavelengths to apply transforms to. If None, applies to all.
        """
        print("Loading single-channel stitched XML...")
        single_channel_data = self.load_xml(single_channel_xml_path)
        
        print("Loading multichannel camera-aligned XML...")
        multichannel_data = self.load_xml(multichannel_xml_path)
        
        # Create a copy for modification
        output_data = copy.deepcopy(multichannel_data)
        
        # Build mapping of tile positions to transforms from single-channel
        print("Extracting transforms from single-channel XML...")
        single_channel_transforms = self._build_tile_transform_map(
            single_channel_data, 
            transform_index
        )
        
        print(f"Found transforms for {len(single_channel_transforms)} tile positions")
        
        # Apply transforms to multichannel data
        print("Applying transforms to multichannel XML...")
        tiles_updated = self._apply_transforms_to_multichannel(
            output_data,
            single_channel_transforms,
            channels_to_apply
        )
        
        print(f"Updated {tiles_updated} tiles in multichannel XML")
        
        # Save the result
        self.save_xml(output_data, output_xml_path)
        print(f"Successfully saved combined XML to: {output_xml_path}")
    
    def _build_tile_transform_map(
        self, 
        data: dict, 
        transform_index: int
    ) -> Dict[str, dict]:
        """
        Build a mapping of tile positions to their transforms.
        
        Parameters
        ----------
        data : dict
            Parsed XML data
        transform_index : int
            Which transform to extract (0-based index)
            
        Returns
        -------
        dict
            Mapping of tile position strings to transform dictionaries
        """
        transform_map = {}
        
        viewsetups = data["SpimData"]["SequenceDescription"]["ViewSetups"]["ViewSetup"]
        if not isinstance(viewsetups, list):
            viewsetups = [viewsetups]
        
        for viewsetup in viewsetups:
            tile_id = self.get_tile_id_from_name(data, viewsetup['name'])
            tile_name = self.extract_tile_name_from_viewsetup(viewsetup)
            
            # Get transforms for this tile
            transforms = self.extract_transforms_for_tile(data, tile_id)
            
            if len(transforms) > transform_index:
                # Store the specified transform
                if transforms[transform_index]['Name'] == "Stitching Transform":
                    transform_map[tile_name] = transforms[transform_index]
                elif transforms[0]['Name'] == "Stitching Transform":
                    transform_map[tile_name] = transforms[0]
                print(f"  Extracted transform for {tile_name}")
            else:
                print(f"  Warning: No transform at index {transform_index} for {tile_name}")
        
        return transform_map
    
    def _apply_transforms_to_multichannel(
        self,
        data: dict,
        transform_map: Dict[str, dict],
        channels_to_apply: Optional[List[int]] = None
    ) -> int:
        """
        Apply transforms from map to multichannel XML data.
        
        Parameters
        ----------
        data : dict
            Multichannel XML data (modified in place)
        transform_map : dict
            Mapping of tile positions to transforms
        channels_to_apply : list, optional
            Channels to apply transforms to
            
        Returns
        -------
        int
            Number of tiles updated
        """
        tiles_updated = 0
        
        viewsetups = data["SpimData"]["SequenceDescription"]["ViewSetups"]["ViewSetup"]
        if not isinstance(viewsetups, list):
            viewsetups = [viewsetups]
        
        view_registrations = data["SpimData"]["ViewRegistrations"]["ViewRegistration"]
        if not isinstance(view_registrations, list):
            view_registrations = [view_registrations]
            data["SpimData"]["ViewRegistrations"]["ViewRegistration"] = view_registrations
        
        for viewsetup in viewsetups:
            tile_id = self.get_tile_id_from_name(data, viewsetup['name'])
            tile_name = self.extract_tile_name_from_viewsetup(viewsetup)
            channel = self.get_channel_from_viewsetup(viewsetup)
            
            # Check if we should apply to this channel
            if channels_to_apply and channel not in channels_to_apply:
                continue
            
            # Check if we have a transform for this tile position
            if tile_name in transform_map:
                # Find the corresponding ViewRegistration
                for view_reg in view_registrations:
                    if int(view_reg.get("@setup", -1)) == tile_id:
                        # Get existing transforms
                        existing_transforms = view_reg.get("ViewTransform", [])
                        if not isinstance(existing_transforms, list):
                            existing_transforms = [existing_transforms]
                        
                        # Append the stitching transform
                        new_transform = copy.deepcopy(transform_map[tile_name])
                        new_transform["Name"] = f"Stitching Transform from Single Channel"
                        
                        # Add to transform list
                        existing_transforms.append(new_transform)
                        view_reg["ViewTransform"] = existing_transforms
                        
                        tiles_updated += 1
                        print(f"  Applied transform to {viewsetup.get('name', 'unknown')} (channel {channel})")
                        break
        
        return tiles_updated
    
    def split_multichannel_xml(
        self,
        multichannel_xml_path: str,
        output_dir: str,
        output_prefix: str = "channel"
    ) -> Dict[int, str]:
        """
        Split a multichannel XML into individual channel XMLs.
        
        Parameters
        ----------
        multichannel_xml_path : str
            Path to multichannel XML
        output_dir : str
            Directory for output XMLs
        output_prefix : str
            Prefix for output files (default: "channel")
            
        Returns
        -------
        dict
            Mapping of channel wavelengths to output file paths
        """
        print("Loading multichannel XML...")
        data = self.load_xml(multichannel_xml_path)
        
        # Group ViewSetups by channel
        channel_groups = self._group_viewsetups_by_channel(data)
        
        print(f"Found {len(channel_groups)} channels to split")
        
        output_files = {}
        
        for channel, viewsetup_ids in channel_groups.items():
            print(f"\nProcessing channel {channel} with {len(viewsetup_ids)} tiles...")
            
            # Create a copy for this channel
            channel_data = copy.deepcopy(data)
            
            # Filter ViewSetups and ViewRegistrations
            # self._filter_xml_to_channel(channel_data, viewsetup_ids)
            reindex_tiles = False            
            # Build ID mapping if reindexing
            id_mapping = {}
            if reindex_tiles:
                for new_id, old_id in enumerate(sorted(viewsetup_ids)):
                    id_mapping[old_id] = new_id
                print(f"  Reindexing {len(id_mapping)} tiles starting from 0")
            else:
                # Identity mapping
                id_mapping = {vid: vid for vid in viewsetup_ids}
            self._filter_xml_to_channel_complete(channel_data, viewsetup_ids, id_mapping = id_mapping, channel = int(channel))
            
            # Save the channel-specific XML
            output_path = f"{output_dir}/{output_prefix}_{channel}.xml"
            self.save_xml(channel_data, output_path)
            output_files[channel] = output_path
            
            print(f"  Saved channel {channel} to: {output_path}")
        
        return output_files
    
    def _group_viewsetups_by_channel(self, data: dict) -> Dict[int, List[int]]:
        """
        Group ViewSetup IDs by channel.
        
        Parameters
        ----------
        data : dict
            Parsed XML data
            
        Returns
        -------
        dict
            Mapping of channel wavelengths to lists of ViewSetup IDs
        """
        channel_groups = {}
        
        viewsetups = data["SpimData"]["SequenceDescription"]["ViewSetups"]["ViewSetup"]
        if not isinstance(viewsetups, list):
            viewsetups = [viewsetups]
        
        for viewsetup in viewsetups:
            tile_id = int(viewsetup.get("id", -1))
            channel = self.get_channel_from_viewsetup(viewsetup)
            
            if channel:
                if channel not in channel_groups:
                    channel_groups[channel] = []
                channel_groups[channel].append(tile_id)
        
        return channel_groups
    
    def _filter_xml_to_channel(self, data: dict, viewsetup_ids: List[int]) -> None:
        """
        Filter XML data to only include specified ViewSetup IDs.
        
        Parameters
        ----------
        data : dict
            XML data (modified in place)
        viewsetup_ids : list
            List of ViewSetup IDs to keep
        """
        # Filter ViewSetups
        viewsetups = data["SpimData"]["SequenceDescription"]["ViewSetups"]["ViewSetup"]
        if not isinstance(viewsetups, list):
            viewsetups = [viewsetups]
        
        filtered_viewsetups = [
            vs for vs in viewsetups 
            if int(vs.get("id", -1)) in viewsetup_ids
        ]
        
        data["SpimData"]["SequenceDescription"]["ViewSetups"]["ViewSetup"] = filtered_viewsetups
        
        # Filter ViewRegistrations
        view_registrations = data["SpimData"]["ViewRegistrations"]["ViewRegistration"]
        if not isinstance(view_registrations, list):
            view_registrations = [view_registrations]
        
        filtered_registrations = [
            vr for vr in view_registrations 
            if int(vr.get("@setup", -1)) in viewsetup_ids
        ]
        
        data["SpimData"]["ViewRegistrations"]["ViewRegistration"] = filtered_registrations
    
    def _filter_xml_to_channel_complete(
        self, 
        data: dict, 
        viewsetup_ids: List[int],
        id_mapping: Dict[int, int],
        channel: int
    ) -> None:
        """
        Comprehensively filter XML data to only include data for a single channel.
        This includes ViewSetups, ViewRegistrations, ImageLoader paths, and all references.
        
        Parameters
        ----------
        data : dict
            XML data (modified in place)
        viewsetup_ids : list
            List of ViewSetup IDs to keep
        id_mapping : dict
            Mapping from old IDs to new IDs (can be identity mapping)
        channel : int
            Channel wavelength being extracted
        """
        # 1. Filter and optionally reindex ViewSetups
        viewsetups = data["SpimData"]["SequenceDescription"]["ViewSetups"]["ViewSetup"]
        if not isinstance(viewsetups, list):
            viewsetups = [viewsetups]
        
        filtered_viewsetups = []
        for vs in viewsetups:
            old_id = int(vs.get("id", -1))
            if old_id in viewsetup_ids:
                # Update ID if reindexing
                vs["id"] = str(id_mapping[old_id])
                
                # Update tile attribute if reindexing
                if "attributes" in vs and "tile" in vs["attributes"]:
                    vs["attributes"]["tile"] = str(id_mapping[old_id])
                
                filtered_viewsetups.append(vs)
        
        data["SpimData"]["SequenceDescription"]["ViewSetups"]["ViewSetup"] = filtered_viewsetups
        
        # 2. Filter and update ViewRegistrations
        view_registrations = data["SpimData"]["ViewRegistrations"]["ViewRegistration"]
        if not isinstance(view_registrations, list):
            view_registrations = [view_registrations]
        
        filtered_registrations = []
        for vr in view_registrations:
            old_setup = int(vr.get("@setup", -1))
            if old_setup in viewsetup_ids:
                # Update setup reference
                vr["@setup"] = str(id_mapping[old_setup])
                
                filtered_registrations.append(vr)
        
        data["SpimData"]["ViewRegistrations"]["ViewRegistration"] = filtered_registrations
        
        # 3. Update ImageLoader section to only include relevant zarr paths
        if "SequenceDescription" in data["SpimData"] and "ImageLoader" in data["SpimData"]["SequenceDescription"]:
            image_loader = data["SpimData"]["SequenceDescription"]["ImageLoader"]
            
            # Handle Zarr format
            if "format" in image_loader and "zarr" in image_loader["format"].lower():
                if "zarr" in image_loader:
                    zarr_data = image_loader["zarr"]
                    
                    # Filter dataset entries
                    if "dataset" in zarr_data:
                        datasets = zarr_data["dataset"]
                        if not isinstance(datasets, list):
                            datasets = [datasets]
                        
                        filtered_datasets = []
                        for ds in datasets:
                            # Check if this dataset is for our channel
                            if "path" in ds and f"ch_{channel}" in ds["path"]:
                                # Update setup reference if present
                                if "@setup" in ds:
                                    old_setup = int(ds["@setup"])
                                    if old_setup in viewsetup_ids:
                                        ds["@setup"] = str(id_mapping[old_setup])
                                        filtered_datasets.append(ds)
                                else:
                                    # If no setup attribute, include if path matches channel
                                    filtered_datasets.append(ds)
                        
                        zarr_data["dataset"] = filtered_datasets
        
        # 4. Filter ViewInterestPoints if present
        if "ViewInterestPoints" in data["SpimData"]:
            vip = data["SpimData"]["ViewInterestPoints"]
            if "ViewInterestPoint" in vip:
                view_interest_points = vip["ViewInterestPoint"]
                if not isinstance(view_interest_points, list):
                    view_interest_points = [view_interest_points]
                
                filtered_vips = []
                for vip_entry in view_interest_points:
                    old_setup = int(vip_entry.get("@setup", -1))
                    if old_setup in viewsetup_ids:
                        vip_entry["@setup"] = str(id_mapping[old_setup])
                        filtered_vips.append(vip_entry)
                
                if filtered_vips:
                    data["SpimData"]["ViewInterestPoints"]["ViewInterestPoint"] = filtered_vips
                else:
                    # Remove empty ViewInterestPoints section
                    del data["SpimData"]["ViewInterestPoints"]
        
        # 5. Update base path if it contains channel information
        if "BasePath" in data["SpimData"]["SequenceDescription"]:
            base_path = data["SpimData"]["SequenceDescription"]["BasePath"]
            # You might want to update this to reflect single-channel output
            # For now, keep it as is or add channel suffix
            # base_path_new = f"{base_path}_ch_{channel}"
            # data["SpimData"]["SequenceDescription"]["BasePath"] = base_path_new
        
        # 6. Update any Attributes that might reference multiple channels
        if "Attributes" in data["SpimData"]["SequenceDescription"]["ViewSetups"]:
            attributes = data["SpimData"]["SequenceDescription"]["ViewSetups"]["Attributes"]
            if "Channel" in attributes[1]:
                channels = attributes[1]["Channel"]
                if not isinstance(channels, list):
                    channels = [channels]
                
                # Find the channel entry that matches our wavelength
                matching_channel = None
                for ch in channels:
                    if "name" in ch and str(channel) in ch["name"]:
                        matching_channel = ch
                        break
                    elif "id" in ch:
                        # Check associated ViewSetups
                        ch_id = int(ch["id"])
                        if ch_id in id_mapping.values():
                            matching_channel = ch
                            break
                
                if matching_channel:
                    # Keep only this channel
                    matching_channel["id"] = str(channel)  # Single channel gets ID channel str
                    attributes[1]["Channel"] = matching_channel
        
        # 7. Clean up the total number of setups/timepoints if specified
        if "SequenceDescription" in data["SpimData"]:
            seq_desc = data["SpimData"]["SequenceDescription"]
            
            # Update Timepoints if it lists specific setups
            if "Timepoints" in seq_desc:
                timepoints = seq_desc["Timepoints"]
                if "@type" in timepoints and timepoints["@type"] == "range":
                    # Update range to match filtered setups
                    if "first" in timepoints and "last" in timepoints:
                        new_first = min(id_mapping.values())
                        new_last = max(id_mapping.values())
                        timepoints["first"] = str(new_first)
                        timepoints["last"] = str(new_last)
        
        print(f"    Filtered data to {len(filtered_viewsetups)} tiles for channel {channel}")
        print(f"    Updated {len(filtered_registrations)} view registrations")
        
        # Log what was filtered
        if "ImageLoader" in data["SpimData"]["SequenceDescription"]:
            if "zarr" in data["SpimData"]["SequenceDescription"]["ImageLoader"]:
                zarr_data = data["SpimData"]["SequenceDescription"]["ImageLoader"]["zarr"]
                if "dataset" in zarr_data:
                    datasets = zarr_data["dataset"]
                    if not isinstance(datasets, list):
                        datasets = [datasets]
                    print(f"    Retained {len(datasets)} zarr dataset entries")

    def validate_transform_transfer(
        self,
        original_xml_path: str,
        updated_xml_path: str
    ) -> Dict[str, any]:
        """
        Validate that transforms were correctly transferred.
        
        Parameters
        ----------
        original_xml_path : str
            Path to original multichannel XML
        updated_xml_path : str
            Path to updated XML with transferred transforms
            
        Returns
        -------
        dict
            Validation results
        """
        original_data = self.load_xml(original_xml_path)
        updated_data = self.load_xml(updated_xml_path)
        
        original_viewregs = original_data["SpimData"]["ViewRegistrations"]["ViewRegistration"]
        updated_viewregs = updated_data["SpimData"]["ViewRegistrations"]["ViewRegistration"]
        
        if not isinstance(original_viewregs, list):
            original_viewregs = [original_viewregs]
        if not isinstance(updated_viewregs, list):
            updated_viewregs = [updated_viewregs]
        
        results = {
            "total_tiles": len(updated_viewregs),
            "tiles_with_new_transforms": 0,
            "transform_differences": []
        }
        
        for updated_vr in updated_viewregs:
            tile_id = int(updated_vr.get("@setup", -1))
            
            # Find corresponding original
            original_vr = None
            for orig_vr in original_viewregs:
                if int(orig_vr.get("@setup", -1)) == tile_id:
                    original_vr = orig_vr
                    break
            
            if original_vr:
                original_transforms = original_vr.get("ViewTransform", [])
                updated_transforms = updated_vr.get("ViewTransform", [])
                
                if not isinstance(original_transforms, list):
                    original_transforms = [original_transforms]
                if not isinstance(updated_transforms, list):
                    updated_transforms = [updated_transforms]
                
                if len(updated_transforms) > len(original_transforms):
                    results["tiles_with_new_transforms"] += 1
                    results["transform_differences"].append({
                        "tile_id": tile_id,
                        "original_count": len(original_transforms),
                        "updated_count": len(updated_transforms)
                    })
        
        return results


# Convenience functions for common operations
def transfer_stitching_to_multichannel(
    single_channel_xml: str,
    multichannel_xml: str,
    output_xml: str,
    channels: Optional[List[int]] = None
):
    """
    Convenience function to transfer stitching transforms.
    
    Parameters
    ----------
    single_channel_xml : str
        Path to single-channel stitched XML
    multichannel_xml : str
        Path to multichannel camera-aligned XML
    output_xml : str
        Output path for combined XML
    channels : list, optional
        Channels to apply transforms to
    """
    manager = BigStitcherXMLManager()
    manager.transfer_stitching_transforms(
        single_channel_xml,
        multichannel_xml,
        output_xml,
        transform_index=1,  # Second transform (stitching)
        channels_to_apply=channels
    )


def split_multichannel_xml(xml_path: str, output_dir: str):
    """
    Convenience function to split multichannel XML.
    
    Parameters
    ----------
    xml_path : str
        Path to multichannel XML
    output_dir : str
        Directory for output files
        
    Returns
    -------
    dict
        Mapping of channels to output files
    """
    manager = BigStitcherXMLManager()
    return manager.split_multichannel_xml(xml_path, output_dir)

    


# Example usage
if __name__ == "__main__":
    # Example 1: Transfer stitching transforms
    print("=" * 60)
    print("Example 1: Transfer stitching transforms")
    print("=" * 60)
    
    transfer_stitching_to_multichannel(
        single_channel_xml="s3://aind-open-data/HCR_000000-s43_2025-07-24_13-00-00_processed_2025-08-28_22-50-35/image_tile_alignment/bigstitcher.xml",
        multichannel_xml="s3://aind-open-data/HCR_000000-s43_2025-07-24_13-00-00_processed_2025-08-28_22-50-35/image_tile_alignment/stitching_cam_alignment_spot_channels.xml",
        output_xml="/scratch/multichannel_with_stitching.xml",
        # channels=[488, 561, 647]  # Optional: only apply to specific channels
    )
    
    # Example 2: Split multichannel XML
    print("\n" + "=" * 60)
    print("Example 2: Split multichannel XML")
    print("=" * 60)
    
    output_files = split_multichannel_xml(
        xml_path="/scratch/multichannel_with_stitching.xml",
        output_dir="/scratch/single_channel_xmls"
    )
    print(f"Created {len(output_files)} channel-specific XMLs")
    
    # Example 3: Validate transform transfer
    # print("\n" + "=" * 60)
    # print("Example 3: Validate transform transfer")
    # print("=" * 60)
    
    # manager = BigStitcherXMLManager()
    # validation = manager.validate_transform_transfer(
    #     original_xml_path="multichannel_aligned.xml",
    #     updated_xml_path="multichannel_with_stitching.xml"
    # )
    # print(f"Validation results:")
    # print(f"  Total tiles: {validation['total_tiles']}")
    # print(f"  Tiles with new transforms: {validation['tiles_with_new_transforms']}")