#!/usr/bin/env python3

import argparse
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from matplotlib.collections import LineCollection
from pihm import read_mesh
from pihm import read_river
from pihm import read_output

def main():
    '''
    An example Python script that plots MM-PIHM output.
    Requires PIHM-utils Python package, which can be installed using:
        pip install PIHM-utils

    For variables of triangular elements, two plots will be produced:
        1. A spatial plot showing the temporal average of variables, and
        2. A time series plot showing the spatial average of variables.
    For variables of river segments, only the time series plot will be produced.
    Note that not all river segments will be plotted, but only the outlet
    segment(s).
    '''

    parser = argparse.ArgumentParser(description='Plot MM-PIHM output')
    parser.add_argument(
        '-s',
        '--simulation',
        type=str,
        #required=True,
        default='ShaleHills',
        help='Name of simulation')
    parser.add_argument(
        '-o',
        '--output_dir',
        type=str,
        #required=True,
        default='ShaleHillsTestRun',
        help='Name of output directory')
    parser.add_argument(
        '-v',
        '--var',
        type=str,
        #required=True,
        default='gw',
        help='Output variable')

    args = parser.parse_args()

    _, from_node, to_node, outlets = read_river('.', args.simulation)
    sim_time, var, vname, unit = read_output('.', args.simulation, args.output_dir, args.var)

    # Set font size
    matplotlib.rcParams.update({'font.size': 14})

    # Plot spatial distribution
    if (args.var[0:6] != 'river.'):
        # Read mesh file and model output. Mesh file is needed to plot spatial distribution
        _, _, tri, x, y, _, _ = read_mesh('.', args.simulation)

        # Create line collection of river segments for plotting
        lines = [[(x[node1], y[node1]), (x[node2], y[node2])] for (node1, node2) in zip(from_node, to_node)]
        river_segments = LineCollection(lines,
            linewidths=5,
            colors='blue',
            alpha=0.7,
            linestyle='solid',
        )

        fig = plt.figure()

        # Test the shape of domain to determine if portrait or landscape direction should be used
        domain_xy_ratio = (np.max(x) - np.min(x)) / (np.max(y) - np.min(y))
        ax = fig.add_axes([0, 0, 0.8, 1]) if domain_xy_ratio < 1.33 else fig.add_axes([0, 0.2, 1, 0.8])

        # Plot output in space
        tpc = ax.tripcolor(x, y, tri,
            facecolors=np.average(var, axis=0),     # Average in time
            edgecolors='k',
            lw=0.5,
            cmap='RdBu'
        )
        ax.add_collection(river_segments)
        ax.set_aspect('equal')

        # Add color bar with variable name and unit
        if domain_xy_ratio < 1.33:
            cbaxes = fig.add_axes([0.8, 0.2, 0.025, 0.6])
            cbar = fig.colorbar(tpc, cax=cbaxes)
        else:
            cbaxes = fig.add_axes([0.2, 0.15, 0.6, 0.04])
            cbar = fig.colorbar(tpc, orientation='horizontal', cax=cbaxes)
        cbar.set_label(label='%s (%s)' % (vname, unit), size=16)
        ax.axis('off')

    # Plot average in time
    fig, ax = plt.subplots()
    if (args.var[0:6] != 'river.'):
        ax.plot(sim_time, np.average(var, axis=1), '-o')
    else:
        ax.plot(sim_time, var[:, outlets], '-o')

    # Set x and y axis labels
    ax.set_xlabel('Date', fontsize=16)
    ax.set_ylabel('%s (%s)' % (vname, unit), fontsize=16)

    # Clean up the x axis dates
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m/%y"))
    fig.tight_layout()  # Otherwise y-label could be slightly clipped

    plt.show()


if __name__ == '__main__':
    main()
