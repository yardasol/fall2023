import numpy as np
import argparse
import json
import cProfile

from particle import Particle

def run():
    threads, i = parse_arguments()

    # Get XS, regions, boundaries
    xs, xs_bins, regions, global_bounds, batches, inactive_batches, particles_per_batch = read_input_file(i)

    # Particle tracking loop
    source_bank = []
    old_source_bank = []
    # TODO: implement machinery to allow simulating neutrons on a process
    for b in range(0, batches):
        if b < inactive_batches:
            tally = False
        else:
            tally = True
            track_length_arr= []
            collision_arr= []
        # Source bank neutrons
        for r in old_source_bank:
            track_length, collisions, fission_neutrons, fission_r = \
                simulate_particle(regions, global_bounds, xs, xs_bins, r=r, tally=tally)
            if tally:
                track_length_arr += [track_length]
                collision_arr += [collisions]
            if fission_neutrons:
                source_bank += fission_neutrons * [fission_r]

        # Base neutrons
        for i in range(0, particles_per_batch):
     #       print(f'Particle {i}')
            track_length, collisions, fission_neutrons, fission_r = \
                simulate_particle(regions, global_bounds, xs, xs_bins, r=None, tally=tally)
            if tally:
                track_length_arr += [track_length]
                collision_arr += [collisions]
            if fission_neutrons:
                source_bank += fission_neutrons * [fission_r]

        N = particles_per_batch + len(old_source_bank)
        if b > 0:
            print(f'k = {N / N_old} for generation {b}')

        N_old = N
        old_source_bank = source_bank
        source_bank = []

def parse_arguments():
    """Parses arguments from command line.

    Returns
    -------
    s : int
        Number of threads to use for shared-memory parallelism.
    i : str
        Path and name of main SaltProc input file (json format).

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--threads',
                        type=int,
                        default=None,
                        help='Number of threads to use for shared-memory \
                        parallelism.')
    parser.add_argument('-i',      # main input file
                        type=str,
                        default=None,
                        help='Path and name of JSON input file')
    args = parser.parse_args()
    return args.threads, args.i


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
    xs : dict of dict of str to float
        Ditionary mapping regions to
        cross sections

    xs_bins : dict of dict of tuple to str
        Dictionary mapping reigons to reaction
        bins.

    regions : dict of int to 2-tuple
        Dictionary mapping regions to
        region bounds.

    global_bounds :  2-tuple
        Global geometry bounds

    batches : int
        Total number of batches to simulate

    inactive_batches : int
        Number of batches to simulate without tallies

    particles_per_batch : int
        Number of particles to simulate per batch
    """

    with open(main_inp_file) as f:
        input_parameters = json.load(f)
        # TODO: Validate input with JSON schema

    batches = input_parameters['batches']
    inactive_batches = input_parameters['inactive_batches']
    particles_per_batch = input_parameters['particles_per_batch']
    xs = {}
    xs_bins = {}
    regions = {}
    xmins = []
    xmaxes = []
    # Get XS, regions for each region
    for region_id, values in input_parameters['regions'].items():
        xs[region_id] = {}
        xs_bins[region_id] = {}
        Sigma_a = values['Sigma_a']
        Sigma_s = values['Sigma_s']
        nuSigma_f = values['nuSigma_f']
        Sigma_t = Sigma_a + Sigma_s

        xs[region_id]['Sigma_a'] = Sigma_a
        xs_bins[region_id]['Sigma_a'] = (0,Sigma_a/Sigma_t)

        xs[region_id]['Sigma_s'] = Sigma_s
        xs_bins[region_id]['Sigma_s'] = (Sigma_a/Sigma_t, 1)

        xs[region_id]['nuSigma_f'] = nuSigma_f
        xs[region_id]['Sigma_t'] = Sigma_t

        regions[region_id] = [values['xmin'], values['xmax']]

        # TODO: check that xmin < xmax
        xmins += [values['xmin']]
        xmaxes += [values['xmax']]

    # TODO: check that regions do not overlap and that
    # there is no void space
    global_bounds = [np.min(xmins), np.max(xmaxes)]

    return xs, xs_bins, regions, global_bounds, batches, inactive_batches, particles_per_batch

def simulate_particle(regions, global_bounds, xs, xs_bins, r=None, tally=True):
    """Simulate particle from birth til death

    Parameters
    ----------


    Returns
    -------
    """
    p = Particle(regions, global_bounds, r=r)
    track_length = 0
    fission_neutrons = 0
    fission_r = None
    collisions = 0
    while p.alive:
        d_b = p.distance_to_boundary(regions)
        d_c = p.distance_to_collision(xs)

        if d_b <= d_c:
           p.translate(d_b)
           if tally:
               track_length += d_b
           p.get_current_region(regions, global_bounds)
        else:
           p.translate(d_c)
           if tally:
               track_length += d_c
               collisions += 1
           rxn = p.sample_reaction(xs, xs_bins)
           if rxn == 'Sigma_s':
               p.sample_direction()
           elif rxn == 'Sigma_a':
               fission_r, fission_neutrons = p.sample_fission_neutrons(xs)

    return track_length, collisions, fission_neutrons, fission_r

#cProfile.run('run()')
run()
