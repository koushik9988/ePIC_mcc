import h5py
import os
import sys
from os.path import join as pjoin

# Get the path to the directory containing the HDF5 file from command-line arguments
path = sys.argv[1]
file_name = "result.h5"

# Open the HDF5 file
file_path = pjoin(path, file_name)
with h5py.File(file_path, 'r') as file:
    # Access and display the '/metadata' group if it exists
    if '/metadata' in file:
        metadata_group = file['/metadata']
        print("Content of /metadata group:\n")
        
        # Display datasets and attributes inside the group
        for key, value in metadata_group.items():
            print(f"Dataset: {key}, Value: {value[:]}")
        
        for attr_name, attr_value in metadata_group.attrs.items():
            print(f"{attr_name} : {attr_value}")
    
    # Access and display the '/metadata_species' group if it exists
    if 'metadata_species' in file:
        species_group = file['metadata_species']
        print("Content of /metadata_species group:\n")
        
        # Iterate over each species in the group
        for species_name, species_data in species_group.items():
            print(f"Species: {species_name}")
            # Display datasets for the species
            for key, value in species_data.items():
                print(f"{key}: {value[:]}")
            
            # Display attributes for the species
            for attr_name, attr_value in species_data.attrs.items():
                print(f"{attr_name} : {attr_value}")
            print()
