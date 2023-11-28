import numpy as np
from numpy.linalg import norm
from scipy.special import roots_legendre
import json
import argparse

def parse_arguments():
    """Parses arguments from command line.

    Returns
    -------
    i : str
        Path and name of main SaltProc input file (json format).

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',      # main input file
                        type=str,
                        default=None,
                        help='Path and name of JSON input file')
    args = parser.parse_args()
    return args.i

def run(f, N = None, S = None):
    S_override = S
    N_override = N

    # Get cross sections, regions, boundaries
    N, S, cross_sections, sources, regions, global_bounds = read_input_file(f)
    if S_override is None:
        pass
    else:
        S = S_override

    if N_override is None:
        pass
    else:
        N = N_override

    # Construct cross section arrays:
    xs = np.linspace(global_bounds[0], global_bounds[1], N)
    Sigma_s_array, Sigma_t_array, Q_array = get_constants_arrays(cross_sections, sources, regions, xs)

    # Get ordinates and weights
    mus, weights = roots_legendre(S)

    # Initial guess
    phi_m = np.ones((S, N))
    phi_mp1 = iteration_step(S, N, phi_m, mus, weights, Sigma_s_array, Sigma_t_array, Q_array)
    while norm(phi_mp1 - phi_m, 1) > 0.001:
        phi_m = phi_mp1
        phi_mp1 = iteration_step(S, N, phi_m, mus, weights, Sigma_s_array, Sigma_t_array, Q_array)

    return xs, mus, weights, phi_mp1

def iteration_step(S, N, phi_m, mus, weights, Sigma_s_array, Sigma_t_array, Q_array):
    indices = np.arange(0, N, 1)
    dx = (N - 1)**-1
    phi_mp1 = np.zeros((S, N))
    for j, mu in enumerate(mus):
        # Negative directions
        if mu < 0:
            idxs = np.flip(indices)
            phi_iphalf = 0 # from vaccum bc
            for i in idxs:
                qn = q_n(i, weights, phi_m, Sigma_s_array, Q_array)
                phi_i, phi_imhalf = negtaive_mu_iteration(i, qn, Sigma_t_array, dx, mu, phi_iphalf)

                phi_mp1[j, i] = phi_i
                phi_iphalf = phi_imhalf
        if mu > 0:
            idxs = indices
            phi_imhalf = 0 # from vaccum bc
            for i in idxs:
                qn = q_n(i, weights, phi_m, Sigma_s_array, Q_array)
                phi_i, phi_iphalf = positive_mu_iteration(i, qn, Sigma_t_array, dx, mu, phi_imhalf)
                phi_mp1[j, i] = phi_i
                phi_imhalf = phi_iphalf
    return phi_mp1

#def iteration_step(S, N, phi_m, mus, weights, Sigma_s_array, Sigma_t_array, Q_array):
#    idxs = np.arange(0, N, 1)
#    dx = (N - 1)**-1
#    phi_mp1 = np.zeros((S, N))
#    for j, mu in enumerate(mus[:int(len(mus)/2)]):
#        # Negative directions
#        phi_neg_iphalf = 0 # from vaccum bc
#        phi_pos_imhalf = 0 # from vaccum bc
#        for i in idxs:
#            # pass over i from the left to the right
#            i_neg = N-1-i
#            # pass over j from the left to the right
#            j_pos = S-1-j
#            qn = q_n(i_neg, weights, phi_m, Sigma_s_array, Q_array)
#
#            # mu < 0
#            phi_i, phi_neg_imhalf = negtaive_mu_iteration(i_neg, qn, Sigma_t_array, dx, mu, phi_neg_iphalf)
#
#            phi_mp1[j, i_neg] = phi_i
#            phi_neg_iphalf = phi_neg_imhalf
#
#            # mu > 0
#            phi_i, phi_pos_iphalf = positive_mu_iteration(i, qn, Sigma_t_array, dx, -mu, phi_pos_imhalf)
#
#            phi_mp1[j_pos, i] = phi_i
#            phi_pos_imhalf = phi_pos_iphalf
#    return phi_mp1

def get_constants_arrays(cross_sections, sources, regions, xs):
    Sigma_s_array = []
    Sigma_t_array = []
    Q_array = []
    for x in xs:
        for region, bounds in regions.items():
            if x >= bounds[0] and x <= bounds[1]:
                region_id = region
                break
        Sigma_s_array += [cross_sections[region_id]['Sigma_s']]
        Sigma_t_array += [cross_sections[region_id]['Sigma_t']]
        Q_array += [sources[region_id]]

    return Sigma_s_array, Sigma_t_array, Q_array

def q_n(i, weights, phi_m, Sigma_s_array, Q_array):
    return Q_array[i]/2 + np.dot(phi_m[:,i], weights) * Sigma_s_array[i] / 2

def negtaive_mu_iteration(i, qn, Sigma_t_array, dx, mu, phi_iphalf):
    alpha_i = 1 + (Sigma_t_array[i] * dx) / (2 * np.abs(mu))
    beta_i = dx / (2 * np.abs(mu))
    phi_i = (qn * beta_i / alpha_i) + phi_iphalf / alpha_i
    phi_imhalf = 2 * phi_i - phi_iphalf
    return phi_i, phi_imhalf

def positive_mu_iteration(i, qn, Sigma_t_array, dx, mu, phi_imhalf):
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

    S : int
        Number of oridnates for angular discretization. Must be an even number

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
    S = input_parameters['S']
    try:
        assert S % 2 == 0
    except AssertionError:
        print('S must be an even number')

    xs = {}
    regions = {}
    sources = {}
    xmins = []
    xmaxes = []
    # Get XS, regions for each region
    for region_id, values in input_parameters['regions'].items():
        xs[region_id] = {}
        Sigma_a = values['Sigma_a']
        Sigma_t = values['Sigma_t']
        Sigma_s = Sigma_t - Sigma_a

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

    return N, S, xs, sources, regions, global_bounds

#i = parse_arguments()
#run(i)


