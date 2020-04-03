#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from pihm_func import read_mesh
from pihm_func import read_output

def main():

    '''
    An example python script that uses functions in pihm_func to plot PIHM
    output results
    Assuming a ShaleHills simulation has been executed and the output directory
    is ShaleHillsTestRun inside the output directory.
    '''

    matplotlib.rcParams.update({'font.size': 14})

    # Plot spatial distribution of groundwater
    _, _, tri, x, y, _, _ = read_mesh('ShaleHills')
    sim_time, gw, vname, unit = read_output('ShaleHills', 'ShaleHillsTestRun',
                                            'gw')

    sim_time, unsat, vname, unit = read_output('ShaleHills', 'ShaleHillsTestRun',
                                            'unsat')
    #domain_xy_ratio = (np.max(x) - np.min(x)) / (np.max(y) - np.min(y))

    ## Average groundwater storage in time
    #var = np.average(gw, axis=0)

    #fig = plt.figure(figsize=(8, 6))
    #if domain_xy_ratio < 1.33:
    #    ax = fig.add_axes([0, 0, 0.8, 1])
    #else:
    #    ax = fig.add_axes([0, 0.2, 1, 0.8])

    ## Plot groundwater in space
    #tpc = ax.tripcolor(x, y, tri, facecolors=var, vmin = 0, vmax = 1.5,
    #                   edgecolors='k', lw = 0.5, cmap='RdBu')
    #ax.set_aspect('equal')

    ## Add color bar with variable name and unit
    #if domain_xy_ratio < 1.33:
    #    cbaxes = fig.add_axes([0.8, 0.2, 0.025, 0.6])
    #    cbar = fig.colorbar(tpc, shrink=0.8, extend='max', cax=cbaxes)
    #else:
    #    cbaxes = fig.add_axes([0.2, 0.15, 0.6, 0.04])
    #    cbar = fig.colorbar(tpc, orientation='horizontal', shrink=0.8,
    #                        extend='max', cax=cbaxes)
    #cbar.set_label(label=vname + '(' + unit + ')', size=18)
    #ax.axis('off')

    # Plot average groundwater storage in time
    fig, ax = plt.subplots(figsize=(8, 6))

    totalw = gw + unsat

    [ax.plot(sim_time, np.average(totalw, 1), '-o') for i in range(9)]

    # Set x and y axis labels
    ax.set_xlabel('Date', fontsize=18)
    ax.set_ylabel(vname + ' (' + unit + ')', fontsize=18)

    ## Clean up the x axis dates
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m/%y"))

    #yticks = ax.get_yticks()
    #ax.set_yticklabels(["$%.2f$" % y for y in yticks])

    #ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
    #color = 'tab:red'
    #ax2.set_ylabel(var1 + ' (' + unit1 + ')', color=color, fontsize=18)
    #ax2.plot(sim_time1, sim_val1[:, 19], color=color)
    #ax2.tick_params(axis='y', labelcolor=color)
    #fig.tight_layout()  # otherwise the right y-label is slightly clipped

    plt.show()


if __name__ == '__main__':
    main()
