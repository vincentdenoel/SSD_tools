import os
import numpy as np
import h5py
import re

"""
Galeria del Vento di Politecnico di Milano (GPVM)
This module contains the functions to load the data from the GVPM HDF5 files.
The data storage corresponds to the format used at the wind tunnel of Politecnico di Milano.

"""
def load_gvpm_hdf5(f, p):
    # Load the HDF5 file into a structured data format
    HDFdata = HDF2Struct(os.path.join(p, f))
    
    # Initialize the data dictionary with the required fields
    data = {}
    data['Configuration'] = HDFdata['']['Configuration']
    data['Angle'] = HDFdata['']['Exposure_Angle']
    data['fsamp'] = HDFdata['']['Sample_Frequency']
    data['q'] = np.array(HDFdata['flowData']['q'], dtype=float)
    data['U'] = np.array(HDFdata['flowData']['U'], dtype=float)
    data['tapsCoords'] = np.array(HDFdata['metadata']['tapsCoords'], dtype=float)
    data['raw'] = {
        'Names': HDFdata['metadata']['tapsNames'],
        'pData': np.array(HDFdata['pressureData']['pressureData'], dtype=float)
    }
    data['dataTypes'] = ['raw']
    data['dataNames'] = ['Raw data']
    data['time'] = np.array(HDFdata['flowData']['time'], dtype=float)
    data['Tile'] = HDFdata['metadata']['tapsNames'][0][0]
    data['qref'] = np.mean(HDFdata['flowData']['q'])
    data['Cp'] = data['raw']['pData'] / float(np.mean(HDFdata['flowData']['q']))

    # Check if 'reference_pressure' field exists and assign it, otherwise set to 0
    if 'reference_pressure' in HDFdata['pressureData']:
        data['pRef'] = HDFdata['pressureData']['reference_pressure']
    else:
        data['pRef'] = 0
    
    return data

def HDF2Struct(f, verbose=False):
    data = {}

    def loadcontent(pathStr):
        # Gets info of current group (or root)
        with h5py.File(f, 'r') as hdf:
            group = hdf[pathStr]

            # Loading variables (Datasets)
            for var_name in group:
                if isinstance(group[var_name], h5py.Dataset):
                    fields = pathStr.strip('/').split('/')
                    fieldsName = validate_field_name(fields)
                    varName = validate_field_name([var_name])
                    setfield(data, fieldsName, varName, group[var_name][:])

            # Loading attributes
            for att_name in group.attrs:
                fields = pathStr.strip('/').split('/')
                fieldsName = validate_field_name(fields)
                attName = validate_field_name([att_name])
                setfield(data, fieldsName, attName, group.attrs[att_name])

            # Loading groups (recursively calls loadcontent for each subgroup)
            for group_name in group:
                if isinstance(group[group_name], h5py.Group):
                    loadcontent(pathStr + '/' + group_name)

    def validate_field_name(names):
        valid_names = []
        for name in names:
            if not re.match(r'^[a-zA-Z_]\w*$', name):
                if verbose:
                    print(f'Warning: "{name}" is not a valid field name.')
                name = re.sub(r'\W|^(?=\d)', '_', name)
                if verbose:
                    print(f'Changed to "{name}"')
            valid_names.append(name)
        return valid_names

    def setfield(data, fieldsName, varName, value):
        current_dict = data
        for key in fieldsName:
            if key not in current_dict:
                current_dict[key] = {}
            current_dict = current_dict[key]
        current_dict[varName[0]] = value

    # Start loading content from the root
    loadcontent('/')

    return data
