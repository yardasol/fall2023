import numpy as np
from numpy.linalg import norm
from scipy.special import roots_legendre
import argparse

def run():
    threads, i = parse_arguments()

    # Get cross sections, regions, boundaries
    N, cross_sections, sources, regions, global_bounds = read_input_file(i)

    # Construct cross section arrays:
    xs = np.linspace(global_bounds[0], global_bounds[1], N)
    Sigma_s_array, Sigma_t_array, Q_array = get_constants_arrays(cross_sections, sources, regions, xs)


    # Get positive and negative ordinates
    mus, weights = roots_legendre(2)

    positive_mus = mus[mus > 0]
    negative_mus = mus[mus < 0]

    # Iteration loop
    dx = (N - 1)**-1
    indices = np.arange(0, N, 1)

    # Initial guess
    phi_m = np.ones((2, N))
    phi_mp1 = iteration_step(phi_m, mus, indices, weights, Sigma_s_array, Sigma_t_array, dx)
    while ...

def iteration_step(phi_m, mus, indices, weights, Sigma_s_array, Sigma_t_array, dx):
    for j, mu in enumerate(mus):
        phi_mp1 = np.zeros((2, N))
        # Negative directions
        if mu < 0:
            idxs = np.flip(indices)
            phi_iphalf = 0 # from vaccum bc
            for i in idxs:
                qn = q_n(i, weights, phi_m, Sigma_s_array, Q_array)
                phi_i, phi_imhalf = negtaive_mu_iteration(i, Sigma_t_array, dx, mu, phi_iphalf)

                phi_mp1[j, i] = phi_i
                phi_iphalf = phi_imhalf
        if mu > 0:
            idxs = indices
            phi_imhalf = 0 # from vaccum bc
            for i in idxs:
                qn = q_n(i, weights, phi_m, Sigma_s_array, Q_array)
                phi_i, phi_iphalf = positive_mu_iteration(i, Sigma_t_array, dx, mu, phi_imhalf)

                phi_mp1[j, i] = phi_i
                phi_imhalf = phi_iphalf
    return phi_mp1

def get_constants_arrays(cross_sections, sources, regions, xs):
    Sigma_s_array = []
    Sigma_t_array = []
    Q_array = []
    for x in xs:
        region_id = 0
        for region, bounds in regions.items():
            if x in bounds:
                region_id = region
                break
        Sigma_s_array += cross_sections[region_id]['Sigma_s']
        Sigma_t_array += cross_sections[region_id]['Sigma_t']
        Q_array += sources[region_id]

    return Sigma_s_array, Sigma_t_array, Q_array

def q_n(i, weights, phi_m, Sigma_s_array, Q_array):
    return Q_array[i]/2 + np.dot(phi_m[:,i], weights) * Sigma_s_array[i] / 2

def negtaive_mu_iteration(i, Sigma_t_arr, dx, mu, phi_iphalf):
    alpha_i = 1 + (Sigma_t_array[i] * dx) / (2 * np.abs(mu))
    beta_i = dx / (2 * np.abs(mu))
    phi_i = (qn * beta_i / alpha_i) + phi_iphalf / alpha_i
    phi_imhalf = 2 * phi_i - phi_iphalf
    return phi_i, phi_imhalf

def positive_mu_iteration(i, Sigma_t_arr, dx, mu, phi_imhalf):
    alpha_i = 1 + (Sigma_t_array[i] * dx) / (2 * np.abs(mu))
    beta_i = dx / (2 * np.abs(mu))
    phi_i = (qn * beta_i / alpha_i) + phi_imhalf / alpha_i
    phi_iphalf = 2 * phi_i - phi_imhalf
    return phi_i, phi_iphalf



def read_input_file(main_inp_file):
    """Read JSON input file and get
    to get cross section and spatial bounds
    in each region.

    Parameters
    ---------
    f : str
        Path to JSON input file

    Returns
    -------
    N : int
        Number of points for spatial discretizations

    xs : dict of dict of str to float
        Ditionary mapping regions to
        cross sections

    regions : dict of int to 2-tuple
        Dictionary mapping regions to
        region bounds.

    global_bounds :  2-tuple
        Global geometry bounds

    """

    with open(main_inp_file) as f:
        input_parameters = json.load(f)
        # TODO: Validate input with JSON schema

    N = input_parameters['N']

    xs = {}
    regions = {}
    sources = {}
    xmins = []
    xmaxes = []
    # Get XS, regions for each region
    for region_id, values in input_parameters['regions'].items():
        xs[region_id] = {}
        Sigma_a = values['Sigma_a']
        Sigma_s = values['Sigma_s']
        Sigma_t = Sigma_a + Sigma_s

        Q = values['Q']

        xs[region_id]['Sigma_a'] = Sigma_a
        xs[region_id]['Sigma_s'] = Sigma_s
        xs[region_id]['Sigma_t'] = Sigma_t

        sources[region_id] = Q

        regions[region_id] = [values['xmin'], values['xmax']]

        # TODO: check that xmin < xmax
        xmins += [values['xmin']]
        xmaxes += [values['xmax']]

    # TODO: check that regions do not overlap and that
    # there is no void space
    global_bounds = [np.min(xmins), np.max(xmaxes)]

    return N, xs, sources, regions, global_bounds



def region
