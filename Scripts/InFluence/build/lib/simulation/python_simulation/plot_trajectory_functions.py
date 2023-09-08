import time
import random
from numba.typed import List
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def PlotTrajectories(vector_coordinates_x, vector_coordinates_y, vector_coordinates_z, points_per_trajectory):
    
    
    path_to_trajectory_plot = "./3D_plot.png"
    path_to_trajectory_plot_svg = "./3D_plot.svg"
    path_to_trajectory_plot_v2 = "./3D_plot_v2.png"
    path_to_trajectory_plot_v2_svg = "./3D_plot_v2.svg"
    path_to_trajectory_histogram_plot = './histogram_plot.png'
    path_to_trajectory_histogram_plot_svg = './histogram_plot.svg'

    
    start = time.time()
    print('\nCalculating electron trajectories...')
    random.seed()
    x_axis, y_axis, z_axis, p= vector_coordinates_x, vector_coordinates_y, vector_coordinates_z, points_per_trajectory
    print('\nPlotting electron trajectories...')
    del p[0]
    del x_axis[0]
    del y_axis[0]
    del z_axis[0]
    p = List(map(int, p))

    distance_from_centre = [] #for plotting lateral spread bar chart

    #suppress plotting to use savefig



    #Plot 3D spread
    fig = plt.figure(figsize=(12, 12), dpi= 100) #change resolution of plot, size of plot
    ax = plt.axes(projection = '3d') #initialise axes
    ax.view_init(elev = 185, azim = 135) #view down z axis #185, 135 fro diaganal view
    plt.xlabel('x direction [$\mu$m]', fontsize = 24, labelpad= 40) #um because multiplied by 10000
    # later
    plt.ylabel('y direction [$\mu$m]', fontsize = 24, labelpad= 40)
    ax.set_zlabel('z direction [$\mu$m]', fontsize = 24, labelpad = 40)
    #ax.w_zaxis.line.set_lw(0.)  #make z linewidth zero and remove labels -> effectively make z axis invisible
    #ax.set_zticks([])
    ax.xaxis.set_tick_params(which='major', size=22, width=2, direction='inout', pad = 15)
    ax.xaxis.set_tick_params(which='minor', size=17, width=2, direction='inout', pad = 15)
    ax.yaxis.set_tick_params(which='major', size=22, width=2, direction='inout', pad = 15)
    ax.yaxis.set_tick_params(which='minor', size=17, width=2, direction='inout', pad = 15)
    ax.zaxis.set_tick_params(which='major', size=22, width=2, direction='inout', pad = 15)
    ax.zaxis.set_tick_params(which='minor', size=17, width=2, direction='inout', pad = 15)
    #ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
    #ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('w')
    ax.yaxis.pane.set_edgecolor('w')
    ax.zaxis.pane.set_edgecolor('w')
    num = len(p)
    beginning = 0
    for i in range(0, num):
            end = beginning  + p[i]
            x_coords = 10000*np.array(x_axis[beginning:end]) #multiply to change to um
            y_coords = 10000*np.array(y_axis[beginning:end])
            z_coords = 10000*np.array(z_axis[beginning:end])


            #Get info. for distance from centre
            distance = (x_coords[-1]**2 +y_coords[-1]**2)**0.5
            distance_from_centre.append(distance)

            #Plot 3D spread
            ax.plot(x_coords, y_coords, z_coords, color=(random.uniform(0, 1), random.uniform(0, 1), random.uniform(
             0, 1)))

            beginning = end

    #Plot 3D spread
    plt.tight_layout()
    plt.autoscale()
    lim_upper = np.max([abs(max(plt.xlim())), abs(max(plt.ylim())), abs(min(plt.xlim())), abs(min(plt.ylim()))])
    plt.xlim(-lim_upper, lim_upper)
    plt.ylim(-lim_upper, lim_upper) #these limits ensure that x and y axes are the same size, ie. square box. Size is
                                     #biggest distance away from centre




    #filename = r'I:\Monte Carlo stuff\ImagesForPaper\sideview60keV.svg'
    plt.savefig(path_to_trajectory_plot, dpi = 300,transparent = True, bbox_inches = 'tight', pad_inches = 0)
    plt.savefig(path_to_trajectory_plot_svg, dpi=300, transparent=True, bbox_inches='tight', pad_inches=0)
    ax.view_init(elev=90, azim=0)
    ax.zaxis.set_ticklabels([])
    ax.set_zlabel('')
    plt.savefig(path_to_trajectory_plot_v2, dpi=300, transparent=True, bbox_inches='tight', pad_inches=0)
    plt.savefig(path_to_trajectory_plot_v2_svg, dpi=300, transparent=True, bbox_inches='tight', pad_inches=0)
    plt.close(fig)
    #print('time elapsed: ', time.clock() - start)

    # Plot histogram

    def bins_labels(bins, **kwargs):
        bin_w = ((max(bins) - min(bins)) / (len(bins) - 1))
        lbls = np.around(bins[1:] - bin_w / 2, decimals=2)
        plt.xticks(np.arange(min(bins) + bin_w / 2, max(bins), bin_w), lbls,
                   **kwargs)


    from matplotlib.ticker import PercentFormatter

    distance_from_centre = np.around(np.array(distance_from_centre), decimals=2)
    distance_from_centre = distance_from_centre[~np.isnan(distance_from_centre)]
    # Plot distance from centre
    fig_hist = plt.figure(figsize=(14, 7), dpi=100)
    plt.style.use('seaborn-whitegrid')
    n, bins, patches = plt.hist(distance_from_centre, bins=10, facecolor='#2ab0ff', edgecolor='black',
                                linewidth=1.5, weights=np.ones(len(distance_from_centre)) / len(distance_from_centre))
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    # plt.title('Histogram of the final distance travelled by each electron relative to the centre of the pixel',
    #        fontsize=16)
    plt.xlabel('Final distance of an electron from its point of incidence [$\mu$m]', fontsize=14, labelpad=25)
    plt.ylabel('Counts', fontsize=14, labelpad=25)
    for count, patch in enumerate(patches):
        #patch.set_facecolor(colors[np.random.randint(100) % nbins])
        x = patch.get_x() + patch.get_width() / 2
        y = patch.get_height()
        #print(y)
        plt.annotate('{}%'.format(np.around(y*100, decimals = 2)), (x, y+0.0075), ha='center')
    #plt.grid()
    #plt.box(on=None)
    bins_labels(bins, rotation=0, fontsize=10)
    plt.ylim()
    bottom, top = plt.ylim()
    plt.ylim(top = 1.2*top)

    plt.savefig(path_to_trajectory_histogram_plot, dpi = 300,transparent = True, bbox_inches = 'tight', pad_inches = 0 )
    plt.savefig(path_to_trajectory_histogram_plot_svg, dpi = 300,transparent = True, bbox_inches = 'tight', pad_inches = 0 )
    plt.close(fig_hist)
    print('\n\nRun time: ', time.time() - start)