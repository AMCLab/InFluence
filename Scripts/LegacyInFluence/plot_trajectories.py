from params import *

N = scc.physical_constants['Avogadro constant'][0]
material_thickness = pixel_dimensions[2]


@jit(nopython = True)
def MCSA_Loop(E_i, minimum_energy, N, Z, A, p, material_thickness, dose):
    #counters
    number_transmitted = 0
    number_backscattered = 0
    number_stopped = 0
    vector_coordinates_x = List([np.float64(0)])
    vector_coordinates_y = List([np.float64(0)])
    vector_coordinates_z = List([np.float64(0)])
    points_per_trajectory = List([np.float64(0)])
    for _ in range(1, dose+1):
        point_counter = 1
        # calculate initial conditions for all electrons
        alpha = evaluate_alpha(E_i, Z)
        cross = evaluate_cross_section(E_i, Z, alpha)
        path_length = evaluate_path_length(A, N, p, cross)
        RND_step = RND(a=0.000001, b=0.999999)
        step = evaluate_step(path_length, RND=RND_step)
        ip = initialise_postions(step = step, d = 27.5*10**-7 ) #units a are in cm
        cx = ip[0]
        cy = ip[1]
        cz = ip[2]
        z0 = ip[3]
        y0 = ip[4]
        x0 = ip[5]
        vector_coordinates_x.append(x0)
        vector_coordinates_y.append(y0)
        vector_coordinates_z.append(z0)
        condition = True
        E = E_i
        while condition:

            if z0 < 0.01:
                number_backscattered  += 1
                condition = False

            if E <= minimum_energy:
                number_stopped  += 1
                condition = False
            if z0 > material_thickness:
                number_transmitted  += 1
                condition = False

            RND_phi = RND(a=0, b=1)  # generate random number for the phi angle
            RND_step = RND(a = 0, b = 1)
            RND_pho = RND(a = 0, b = 1)
            alpha = evaluate_alpha(E, Z) #calc screening constant, function of previous energy
            cross = evaluate_cross_section(E, Z, alpha) #calc cross section
            path_length = evaluate_path_length(A, N, p, cross) #calc mean free path length
            step = evaluate_step(path_length = path_length, RND=RND_step) #calculate step of this iteration
            E = E + step*p*evaluate_energy_loss_rate(E, Z, A) #calc new energy

            phi = evaluate_phi(RND = RND_phi, alpha = alpha) #calc scattering angle
            psi = evaluate_pho(RND = RND_pho) #calc other scattering angle

            ca = evaluate_direction_cosine_a(phi, psi, cx, cy, cz) #calc direction cosines
            cb = evaluate_direction_cosine_b(phi, psi, cx, cy, cz)
            cc = evaluate_direction_cosine_c(phi, psi, cz)

            x0 = x0 + step*ca #find and reset to new positions
            y0 = y0 + step*cb
            z0 = z0 + step*cc
            cx = ca #reset direction cosines
            cy = cb
            cz = cc
            vector_coordinates_x.append(x0)
            vector_coordinates_y.append(y0)
            vector_coordinates_z.append(z0)
            point_counter += 1

        points_per_trajectory.append(point_counter)

    return vector_coordinates_x, vector_coordinates_y, vector_coordinates_z, points_per_trajectory


if __name__ == '__main__':
    start = time.time()
    print('\nCalculating electron trajectories...')
    random.seed()
    x_axis, y_axis, z_axis, p= MCSA_Loop(E_i, minimum_energy, N, Z, A, p, material_thickness, dose)
    print('\nPlotting electron trajectories...')
    del p[0]
    del x_axis[0]
    del y_axis[0]
    del z_axis[0]
    p = List(map(int, p))

    distance_from_centre = [] #for plotting lateral spread bar chart

    #suppress plotting to use savefig

    import matplotlib as mpl
    #mpl.use('Agg')
    import matplotlib.font_manager as fm
    import matplotlib.pyplot as plt
    # Rebuild the matplotlib font cache
    fm._rebuild()
    mpl.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.linewidth'] = 2

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
    fig.show()
    plt.pause(1e-30)



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
    plt.pause(1e-30)
    bottom, top = plt.ylim()
    plt.ylim(top = 1.2*top)

    plt.savefig(path_to_trajectory_histogram_plot, dpi = 300,transparent = True, bbox_inches = 'tight', pad_inches = 0 )
    plt.savefig(path_to_trajectory_histogram_plot_svg, dpi = 300,transparent = True, bbox_inches = 'tight', pad_inches = 0 )
    plt.close(fig_hist)
    print('\n\nRun time: ', time.time() - start)
