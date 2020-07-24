#!/usr/bin/env python3

import os
import struct
import numpy as np
from datetime import datetime

def read_mesh(simulation):

    # Full file name
    fname = 'input/' + simulation + '/' + simulation + '.mesh'

    # Read mesh file into an array of strings with leading white spaces removed
    # Line starting with "#" are not read in
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

    return (num_elem, num_nodes, np.array(tri), np.array(x), np.array(y),
            np.array(zmin), np.array(zmax))


def read_river(simulation):

    # Full file name
    fname = 'input/' + simulation + '/' + simulation + '.riv'

    # Read river file into an array of strings with leading white spaces removed
    # Line starting with "#" are not read in
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

    return (num_rivers, np.array(from_node), np.array(to_node))


def read_output(simulation, outputdir, ext):

    # Read number of river segments and elements from input files
    num_rivers, _, _ = read_river(simulation)
    num_elem, _, _, _, _, _, _ = read_mesh(simulation)

    # Determine output dimension, variable name and unit from extension
    if ext[0:6] == 'river.':
        dim = num_elem
    else:
        dim = num_rivers

    if ext == 'surf':
        varname = 'Surface water'
        unit = 'm'
    elif ext == 'unsat':
        varname = 'Unsaturated zone storage'
        unit = 'm'
    elif ext == 'gw':
        varname = 'Groundwater storage'
        unit = 'm'
    elif ext == 'river.stage':
        varname = 'River stage'
        unit = 'm'
    elif ext == 'snow':
        varname = 'Water equivalent snow depth'
        unit = 'm'
    elif ext == 'is':
        varname = 'Interception storage'
        unit = 'm'
    elif ext == 'infil':
        varname = 'Infiltration'
        unit = 'm s$^{-1}$'
    elif ext == 'recharge':
        varname = 'Recharge'
        unit = 'm s$^{-1}$'
    elif ext == 'ec':
        varname = 'Canopy evaporation'
        unit = 'm s$^{-1}$'
    elif ext == 'ett':
        varname = 'Transpiration'
        unit = 'm s$^{-1}$'
    elif ext == 'edir':
        varname = 'Soil evaporation'
        unit = 'm s$^{-1}$'
    elif ext[0:9] == 'river.flx':
        varname = 'River flux ' + ext[9:]
        unit = 'm$^3$ s$^{-1}$'
    elif ext[0:6] == 'subflx':
        varname = 'Subsurface flux ' + ext[6]
        unit = 'm$^3$ s$^{-1}$'
    elif ext[0:7] == 'surfflx':
        varname = 'Surface flux ' + ext[7]
        unit = 'm$^3$ s$^{-1}$'
    elif ext == 't1':
        varname = 'Land surface temperature'
        unit = 'K'
    elif ext[0:3] == 'stc':
        varname = 'Soil temperature (Layer ' + str(int(ext[3:]) + 1) + ')'
        unit = 'K'
    elif ext[0:3] == 'smc':
        varname = 'Soil moisture content (Layer ' + str(int(ext[3:]) + 1) + ')'
        unit = 'm$^3$ m$^{-3}$'
    elif ext[0:3] == 'swc':
        varname = 'Soil water content (Layer ' + str(int(ext[3:]) + 1) + ')'
        unit = 'm$^3$ m$^{-3}$'
    elif ext == 'snowh':
        varname = 'Snow depth'
        unit = 'm'
    elif ext == 'iceh':
        varname = 'Ice depth'
        unit = 'm'
    elif ext == 'albedo':
        varname = 'Albedo'
        unit = '-'
    elif ext == 'le':
        varname = 'Latent heat flux'
        unit = 'W m$^{-2}$'
    elif ext == 'sh':
        varname = 'Sensible heat flux'
        unit = 'W m$^{-2}$'
    elif ext == 'g':
        varname = 'Ground heat flux'
        unit = 'W m$^{-2}$'
    elif ext == 'etp':
        varname = 'Potential evaporation'
        unit = 'W m$^{-2}$'
    elif ext == 'esnow':
        varname = 'Snow sublimation'
        unit = 'W m$^{-2}$'
    elif ext == 'rootw':
        varname = 'Root zone vailable soil moisture'
        unit = '-'
    elif ext == 'soilm':
        varname = 'Total soil moisture'
        unit = 'm'
    elif ext == 'solar':
        varname = 'Solar radiation'
        unit = 'W m$^{-2}$'
    elif ext == 'ch':
        varname = 'Surface exchange coefficient'
        unit = 'm s$^{-1}$'
    elif ext == 'lai':
        varname = 'LAI'
        unit = 'm$^2$ m$^{-2}$'
    elif ext == 'npp':
        varname = 'NPP'
        unit = 'kgC m$^{-2}$ day$^{-1}$'
    elif ext == 'nep':
        varname = 'NEP'
        unit = 'kgC m$^{-2}$ day$^{-1}$'
    elif ext == 'nee':
        varname = 'NEE'
        unit = 'kgC m$^{-2}$ day$^{-1}$'
    elif ext == 'gpp':
        varname = 'GPP'
        unit = 'kgC m$^{-2}$ day$^{-1}$'
    elif ext == 'mr':
        varname = 'Maintenace respiration'
        unit = 'kgC m$^{-2}$ day$^{-1}$'
    elif ext == 'gr':
        varname = 'Growth respiration'
        unit = 'kgC m$^{-2}$ day$^{-1}$'
    elif ext == 'hr':
        varname = 'Heterotrophic respiration'
        unit = 'kgC m$^{-2}$ day$^{-1}$'
    elif ext == 'fire':
        varname = 'Fire losses'
        unit = 'kgC m$^{-2}$ day$^{-1}$'
    elif ext == 'litfallc':
        varname = 'Litter fall'
        unit = 'kgC m$^{-2}$ day$^{-1}$'
    elif ext == 'vegc':
        varname = 'Vegetation carbon'
        unit = 'kgC m$^{-2}$'
    elif ext == 'agc':
        varname = 'Aboveground carbon'
        unit = 'kgC m$^{-2}$'
    elif ext == 'litrc':
        varname = 'Litter carbon'
        unit = 'kgC m$^{-2}$'
    elif ext == 'soilc':
        varname = 'Soil carbon'
        unit = 'kgC m$^{-2}$'
    elif ext == 'totalc':
        varname = 'Total carbon'
        unit = 'kgC m$^{-2}$'
    elif ext == 'sminn':
        varname = 'Soil mineral nitrogen'
        unit = 'kgN m$^{-2}$'
    elif ext[0:11] == 'grain_yield':
        varname = 'Grain yield ' + ext[12:]
        unit = 'Mg ha$^{-1}$'
    elif ext[0:12] == 'forage_yield':
        varname = 'Forage yield ' + ext[13:]
        unit = 'Mg ha$^{-1}$'
    elif ext[0:5] == 'shoot':
        varname = 'Shoot biomass ' + ext[6:]
        unit = 'Mg ha$^{-1}$'
    elif ext[0:4] == 'root':
        varname = 'Root biomass ' + ext[5:]
        unit = 'Mg ha$^{-1}$'
    elif ext[0:8] == 'radintcp':
        varname = 'Radiation interception ' + ext[9:]
        unit = '-'
    elif ext[0:7] == 'wstress':
        varname = 'Water stress ' + ext[8:]
        unit = '-'
    elif ext[0:7] == 'nstress':
        varname = 'N stress ' + ext[8:]
        unit = '-'
    elif ext[0:6] == 'transp':
        varname = 'Transpiration ' + ext[7:]
        unit = 'mm day$^{-1}$'
    elif ext[0:9] == 'pottransp':
        varname = 'Potential transpiration ' + ext[10:]
        unit = 'mm day$^{-1}$'
    elif ext == 'NO3':
        varname = 'Soil profile NO$_3$'
        unit = 'Mg ha$^{-1}$'
    elif ext == 'NH4':
        varname = 'Soil profile NH$_4$'
        unit = 'Mg ha$^{-1}$'
    elif ext == 'river_NO3':
        varname = 'River NO$_3$'
        unit = 'Mg ha$^{-1}$'
    elif ext == 'river_NH4':
        varname = 'River NH$_4$'
        unit = 'Mg ha$^{-1}$'
    elif ext == 'denitrif':
        varname = 'Denitrification'
        unit = 'Mg ha$^{-1}$'
    elif ext == 'NO3_leaching':
        varname = 'NO$_3$ leaching'
        unit = '0.1 kg s$^{-1}$'
    elif ext == 'NH4_leaching':
        varname = 'NH$_4$ leaching'
        unit = '0.1 kg s$^{-1}$'
    elif ext == 'deep.unsat':
        varname = 'Deep zone unsaturated storage'
        unit = 'm'
    elif ext == 'deep.gw':
        varname = 'Deep groundwater storage'
        unit = 'm'
    elif ext == 'deep.infil':
        varname = 'Deep zone infiltration'
        unit = 'm s$^{-1}$'
    elif ext == 'deep.rechg':
        varname = 'Deep zone recharge'
        unit = 'm s$^{-1}$'
    elif ext[0:9] == 'deep.flow':
        varname = 'Deep layer lateral flow ' + ext[9:]
        unit = 'm$^3$ s$^{-1}$'
    elif ext[0:4]== 'conc':
        varname = ext[5:] + ' concentration'
        unit = 'mol L$^{-1}$'
    elif ext[0:9] == 'deep.conc':
        varname = 'Deep zone ' + ext[10:] + ' concentration'
        unit = 'mol L$^{-1}$'
    elif ext[0:10] == 'river.conc':
        varname = 'Stream ' + ext[11:] + ' concentration'
        unit = 'mol L$^{-1}$'
    elif ext[0:11] == 'river.chflx':
        varname = 'River ' + ext[12:] + ' flux'
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

    return (sim_time, sim_val, varname, unit)
