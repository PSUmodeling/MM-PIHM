#!/usr/bin/env python3

import os
import struct
import numpy as np
from datetime import datetime

def read_mesh(simulation):

    # Full file name
    fname = 'input/' + simulation + '/' + simulation + '.mesh'

    # Read mesh file into an array of strings with '#' and leading white spaces
    # removed
    meshstr = []
    with open(fname) as file:
        [meshstr.append(line.strip()) for line in file
         if line.strip() and line.strip()[0] != '#']

    # Read number of elements
    num_elem = int(meshstr[0].split()[1])

    # Read nodes information
    tri = []
    for nline in range(num_elem):
        strs = meshstr[2 + nline].split()[1:4]
        tri.append([int(i) - 1 for i in strs])

    # Read number of nodes
    num_nodes = int(meshstr[2 + num_elem].split()[1])

    # Read X, Y, ZMIN, and ZMAX
    x = []
    y = []
    zmin = []
    zmax = []
    for nline in range(num_nodes):
        strs = meshstr[4 + num_elem + nline].split()[1:]
        x.append(float(strs[0]))
        y.append(float(strs[1]))
        zmin.append(float(strs[2]))
        zmax.append(float(strs[3]))

    return num_elem, num_nodes, tri, x, y, zmin, zmax


def read_river(simulation):

    # Full file name
    fname = 'input/' + simulation + '/' + simulation + '.riv'

    # Read river file into an array of strings with '#' and leading white spaces
    # removed
    riverstr = []
    with open(fname) as file:
        [riverstr.append(line.strip()) for line in file
         if line.strip() and line.strip()[0] != '#']

    # Read number of river segments
    num_rivers = int(riverstr[0].split()[1])

    # Read nodes information
    from_node = []
    to_node = []
    for kriver in range(num_rivers):
        strs = riverstr[2 + kriver].split()[1:3]
        from_node.append(int(strs[0]) - 1)
        to_node.append(int(strs[1]) - 1)

    return num_rivers, from_node, to_node


def read_output(simulation, outputdir, ext):

    # Read number of river segments and elements from input files
    num_rivers, _, _ = read_river(simulation)
    num_elem, _, _, _, _, _, _ = read_mesh(simulation)

    # Determine output dimension, variable name and unit from extension
    if ext == 'surf':
        dim = num_elem
        varname = 'Surface water'
        unit = 'm'
    elif ext == 'unsat':
        dim = num_elem
        varname = 'Unsatureted zone storage'
        unit = 'm'
    elif ext == 'gw':
        dim = num_elem
        varname = 'Groundwater storage'
        unit = 'm'
    elif ext == 'stage':
        dim = num_rivers
        varname = 'River stage'
        unit = 'm'
    elif ext == 'rivgw':
        dim = num_rivers
        varname = 'River groundwater storage'
        unit = 'm'
    elif ext == 'snow':
        dim = num_elem
        varname = 'Water equivalent snow depth'
        unit = 'm'
    elif ext == 'is':
        dim = num_elem
        varname = 'Interception storage'
        unit = 'm'
    elif ext == 'infil':
        dim = num_elem
        varname = 'Infiltration'
        unit = 'm s$^{-1}$'
    elif ext == 'recharge':
        dim = num_elem
        varname = 'Recharge'
        unit = 'm s$^{-1}$'
    elif ext == 'ec':
        dim = num_elem
        varname = 'Canopy evaporation'
        unit = 'm s$^{-1}$'
    elif ext == 'ett':
        dim = num_elem
        varname = 'Transpiration'
        unit = 'm s$^{-1}$'
    elif ext == 'edir':
        dim = num_elem
        varname = 'Soil evaporation'
        unit = 'm s$^{-1}$'
    elif ext[0:6] == 'rivflx':
        dim = num_rivers
        varname = 'River flux ' + ext[6:]
        unit = 'm$^3$ s$^{-1}$'
    elif ext[0:6] == 'subflx':
        dim = num_elem
        varname = 'Subsurface flux ' + ext[6]
        unit = 'm$^3$ s$^{-1}$'
    elif ext[0:7] == 'surfflx':
        dim = num_elem
        varname = 'Surface flux ' + ext[6]
        unit = 'm$^3$ s$^{-1}$'
    elif ext == 't1':
        dim = num_elem
        varname = 'Land surface temperature'
        unit = 'K'
    elif ext[0:3] == 'stc':
        dim = num_elem
        varname = 'Soil temperature (Layer ' + str(int(ext[3:]) + 1) + ')'
        unit = 'K'
    elif ext[0:3] == 'smc':
        dim = num_elem
        varname = 'Soil moisture content (Layer ' + str(int(ext[3:]) + 1) + ')'
        unit = 'm$^3$ m$^{-3}$'
    elif ext[0:3] == 'swc':
        dim = num_elem
        varname = 'Soil water content (Layer ' + str(int(ext[3:]) + 1) + ')'
        unit = 'm$^3$ m$^{-3}$'
    elif ext == 'snowh':
        dim = num_elem
        varname = 'Snow depth'
        unit = 'm'
    elif ext == 'iceh':
        dim = num_elem
        varname = 'Ice depth'
        unit = 'm'
    elif ext == 'albedo':
        dim = num_elem
        varname = 'Albedo'
        unit = '-'
    elif ext == 'le':
        dim = num_elem
        varname = 'Latent heat flux'
        unit = 'W m$^{-2}$'
    elif ext == 'sh':
        dim = num_elem
        varname = 'Sensible heat flux'
        unit = 'W m$^{-2}$'
    elif ext == 'g':
        dim = num_elem
        varname = 'Ground heat flux'
        unit = 'W m$^{-2}$'
    elif ext == 'etp':
        dim = num_elem
        varname = 'Potential evaporation'
        unit = 'W m$^{-2}$'
    elif ext == 'esnow':
        dim = num_elem
        varname = 'Snow sublimation'
        unit = 'W m$^{-2}$'
    elif ext == 'rootw':
        dim = num_elem
        varname = 'Root zone vailable soil moisture'
        unit = '-'
    elif ext == 'soilm':
        dim = num_elem
        varname = 'Total soil moisture'
        unit = 'm'
    elif ext == 'solar':
        dim = num_elem
        varname = 'Solar radiation'
        unit = 'W m$^{-2}$'
    elif ext == 'ch':
        dim = num_elem
        varname = 'Surface exchange coefficient'
        unit = 'm s$^{-1}$'
    elif ext == 'lai':
        dim = num_elem
        varname = 'LAI'
        unit = 'm$^2$ m$^{-2}$'
    elif ext == 'npp':
        dim = num_elem
        varname = 'NPP'
        unit = 'kgC m$^{-2}$ day$^{-1}$'
    elif ext == 'nep':
        dim = num_elem
        varname = 'NEP'
        unit = 'kgC m$^{-2}$ day$^{-1}$'
    elif ext == 'nee':
        dim = num_elem
        varname = 'NEE'
        unit = 'kgC m$^{-2}$ day$^{-1}$'
    elif ext == 'gpp':
        dim = num_elem
        varname = 'GPP'
        unit = 'kgC m$^{-2}$ day$^{-1}$'
    elif ext == 'mr':
        dim = num_elem
        varname = 'Maintenace respiration'
        unit = 'kgC m$^{-2}$ day$^{-1}$'
    elif ext == 'gr':
        dim = num_elem
        varname = 'Growth respiration'
        unit = 'kgC m$^{-2}$ day$^{-1}$'
    elif ext == 'hr':
        dim = num_elem
        varname = 'Heterotrophic respiration'
        unit = 'kgC m$^{-2}$ day$^{-1}$'
    elif ext == 'fire':
        dim = num_elem
        varname = 'Fire losses'
        unit = 'kgC m$^{-2}$ day$^{-1}$'
    elif ext == 'litfallc':
        dim = num_elem
        varname = 'Litter fall'
        unit = 'kgC m$^{-2}$ day$^{-1}$'
    elif ext == 'vegc':
        dim = num_elem
        varname = 'Vegetation carbon'
        unit = 'kgC m$^{-2}$'
    elif ext == 'agc':
        dim = num_elem
        varname = 'Aboveground carbon'
        unit = 'kgC m$^{-2}$'
    elif ext == 'litrc':
        dim = num_elem
        varname = 'Litter carbon'
        unit = 'kgC m$^{-2}$'
    elif ext == 'soilc':
        dim = num_elem
        varname = 'Soil carbon'
        unit = 'kgC m$^{-2}$'
    elif ext == 'totalc':
        dim = num_elem
        varname = 'Total carbon'
        unit = 'kgC m$^{-2}$'
    elif ext == 'sminn':
        dim = num_elem
        varname = 'Soil mineral nitrogen'
        unit = 'kgN m$^{-2}$'
    elif ext == 'deepunsat':
        dim = num_elem
        varname = 'Deep layer unsaturated storage'
        unit = 'm'
    elif ext == 'deepgw':
        dim = num_elem
        varname = 'Deep groundwater storage'
        unit = 'm'
    elif ext == 'deepinfil':
        dim = num_elem
        varname = 'Deep layer infiltration'
        unit = 'm s$^{-1}$'
    elif ext == 'deeprechg':
        dim = num_elem
        varname = 'Deep layer recharge'
        unit = 'm s$^{-1}$'
    elif ext[0:8] == 'deepflow':
        dim = num_elem
        varname = 'Deep layer lateral flow ' + ext[8:]
        unit = 'm$^3$ s$^{-1}$'
    elif ext[0:4]== 'conc':
        dim = num_elem
        varname = ext[8:] + ' concentration'
        unit = 'mol L$^{-1}$'
    elif ext[0:9] == 'deep_conc':
        dim = num_elem
        varname = 'Deep zone ' + ext[15:] + ' concentration'
        unit = 'mol L$^{-1}$'
    elif ext[0:10] == 'river_conc':
        dim = num_rivers
        varname = 'Stream ' + ext[11:] + ' concentration'
        unit = 'mol L$^{-1}$'
    elif ext[0:15] == 'river_discharge':
        dim = num_rivers
        varname = 'River ' + ext[16:] + ' flux'
        unit = 'kmol s$^{-1}$'

    # Full file name (binary file)
    fname = 'output/' + outputdir + '/' + simulation + '.' + ext + '.dat'

    # Check size of output file
    fsize = int(os.path.getsize(fname) / 8)

    with open(fname, 'rb') as binfile:
        # Read binary output file
        data_str = binfile.read()
        data_tuple = struct.unpack('%dd' %(fsize), data_str)

        # Rearrange read values to numpy array
        data_array = np.resize(data_tuple, (int(fsize / (dim + 1)), dim + 1))

        # Output values
        sim_val = data_array[:, 1:]

        # Convert simulation time
        sim_time = []
        [sim_time.append(datetime.utcfromtimestamp(data_array[i, 0]))
         for i in range(int(fsize / (dim + 1)))]

    return sim_time, sim_val, varname, unit
