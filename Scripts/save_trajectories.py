from params import *
import simplejson
import time
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
    energy_list = List([np.float64(0)])
    points_per_trajectory = List([np.float64(0)])
    for _ in range(1, dose+1):
        point_counter = 1
        # calculate initial conditions for all electrons
        alpha = evaluate_alpha(E_i, Z)
        cross = evaluate_cross_section(E_i, Z, alpha)
        path_length = evaluate_path_length(A, N, p, cross)
        RND_step = RND(a=0.000001, b=0.999999) #
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
        energy_list.append(E)
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
            energy_list.append(E)
            point_counter += 1

        points_per_trajectory.append(point_counter)

    return vector_coordinates_x, vector_coordinates_y, vector_coordinates_z, energy_list, points_per_trajectory


greater_list = []
distance_from_centre = []
if __name__ == '__main__':
    start = time.time()
    print('\nCalculating electron trajectories...')
    random.seed()
    x_axis, y_axis, z_axis, energies_list, p= MCSA_Loop(E_i, minimum_energy, N, Z, A, p, material_thickness, dose)
    print('\nSplitting electron trajectories...')
    del p[0]
    del x_axis[0]
    del y_axis[0]
    del z_axis[0]
    del energies_list[0]
    p = List(map(int, p))
    num = len(p)
    beginning = 0

    # indent=2 is not needed but makes the file human-readable
    for i in range(0, num):
                end = beginning  + p[i]
                x_coords = np.float64(np.array(x_axis[beginning:end])) #multiply to change to um
                #x_coords = x_coords[~np.isnan(x_coords)]
                y_coords = np.float64(np.array(y_axis[beginning:end]))
                #y_coords = y_coords[~np.isnan(y_coords)]
                z_coords = np.float64(np.array(z_axis[beginning:end]))
                #z_coords = z_coords[~np.isnan(z_coords)]
                energies_list_traj = np.float64(np.array(energies_list[beginning:end]))
                list_1 = [list(x_coords), list(y_coords), list(z_coords), list(energies_list_traj)]
                #distance = (x_coords[-1] ** 2 + y_coords[-1] ** 2)** 0.5
                #print(x_coords[-1], y_coords[-1])
                #distance_from_centre.append(distance)
                greater_list.append(list_1)
                #print(len(list_1[0]))


                beginning = end
    print('\nSaving electron trajectories...')
    #print(np.array(distance_from_centre).mean())
    with open(sample_distribution_file, 'w', encoding='utf-8') as f:
                simplejson.dump(greater_list, f, indent=2, ensure_ascii=True)


