import numpy as np
import uncertainties
from uncertainties import unumpy as unp
import argparse
import json
import pandas as pd
import cProfile

from particle import Particle

def run():
    threads, i = parse_arguments()

    # Get XS, regions, boundaries
    xs, xs_bins, regions, volumes, global_bounds, batches, inactive_batches, particles_per_batch = read_input_file(i)

    # Particle tracking loop

    # Global arrays
    keff = []
    batch_track_lengths = {'global': []}
    batch_collisions = {'global': []}
    for region_id in regions.keys():
        batch_track_lengths[region_id] = []
        batch_collisions[region_id] = []

    # Loop arrays
    source_bank = []
    old_source_bank = []
    # TODO: implement machinery to allow simulating neutrons on a process
    for b in range(0, batches):
        # Determine if inactive or active batches
        if b < inactive_batches:
            tally = False
        else:
            tally = True
            track_length_dict = {'global': []}
            collision_dict = {'global': []}
            for region_id in regions.keys():
                track_length_dict[region_id] = []
                collision_dict[region_id] = []

        # Source bank neutrons
        #for r in old_source_bank:
        #    track_length, collisions, N_fission_neutrons, r_fission = \
        #        simulate_particle(regions, global_bounds, xs, xs_bins, r=r, tally=tally)
        #    if tally:
        #        track_length_arr += [track_length]
        #        collision_arr += [collisions]
        #    if N_fission_neutrons:
        #        source_bank += N_fission_neutrons * [r_fission]

        # Base neutrons
        for i in range(0, particles_per_batch):
            # Spawn neutrons from locations in old_source_bank
            # until we run out or we reach the number of neutrons
            # per batch
            if len(old_source_bank):
                r_x = old_source_bank.pop(0)
            else:
                r_x = None
            track_length, collisions, N_fission_neutrons, r_fission_x = \
                simulate_particle(regions, global_bounds, xs, xs_bins, r_x=r_x, tally=tally)
            if tally:
                for region_id, t_length in track_length.items():
                    track_length_dict[region_id] += [t_length]
                for region_id, col in collisions.items():
                    collision_dict[region_id] += [col]
            if N_fission_neutrons:
                source_bank += N_fission_neutrons * [r_fission_x]

        # Get number of particle produced from fission
        N = particles_per_batch + len(old_source_bank)
        if b > 0:
            k = N / N_old
            print(f'k = {k} for generation {b}')
            keff += [k]

        if tally:
            for region_id, track_length_arr in track_length_dict.items():
                batch_track_lengths[region_id] += track_length_arr
            for region_id, collision_arr in collision_dict.items():
                batch_collisions[region_id] += collision_arr

        N_old = N
        old_source_bank = source_bank
        source_bank = []

    print("SIMULATION COMPLETE\n"
          "-------------------")
    print(f"k_eff: {np.average(keff)} +/- {np.var(keff)}")

    vol_avg_sigma, total_volume = get_vol_avg_sigma(xs, volumes)

    print(f"Collision estimate of global flux per unit : {np.average(batch_collisions['global']) / (vol_avg_sigma * total_volume)}")
    print(f"Tracklength estimate of global flux per unit: {np.average(batch_track_lengths['global'] / total_volume)}")

    data = []
    index = []
    for region_id in regions.keys():
        collision_flux = np.average(batch_collisions[region_id]) / (xs[region_id]['Sigma_t'] * volumes[region_id])
        tracklength_flux = np.average(batch_track_lengths[region_id]) / volumes[region_id]
        flux = [collision_flux, tracklength_flux]
        data += [flux]
        index += [region_id]
    columns = ['collision', 'tracklength']

    df = pd.DataFrame(data=data, index=index, columns=columns)
    df.index.name = 'region_id'
    df.to_csv('results.csv')


def get_vol_avg_sigma(xs, volumes):
    numerator = []
    denominator = []
    vol_frac = {}
    total_volume = np.sum(list(volumes.values()))
    for r, vol in volumes.items():
        vol_frac[r] = vol / total_volume
    partial_xs = []
    for r in xs.keys():
        partial_xs += [xs[r]['Sigma_t'] * vol_frac[r]]

    return np.sum(partial_xs), total_volume

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
    volumes = {}
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
        volumes[region_id] = values['xmax'] - values['xmin']

        # TODO: check that xmin < xmax
        xmins += [values['xmin']]
        xmaxes += [values['xmax']]

    # TODO: check that regions do not overlap and that
    # there is no void space
    global_bounds = [np.min(xmins), np.max(xmaxes)]

    return xs, xs_bins, regions, volumes, global_bounds, batches, inactive_batches, particles_per_batch

def simulate_particle(regions, global_bounds, xs, xs_bins, r_x=None, tally=True):
    """Simulate particle from birth til death

    Parameters
    ----------



    Returns
    -------
    track_length : dict
        Dictionary mapping region IDs to region-specific
        tallies of track length. Also includes a `global`
        key.
    collisions : dict
        Dictionary mapping region IDs to region-specific
        tallies of collisions. Also includes a `global`
        key.
    N_fission_neutrons : int
        Number of fission neutrons produced by the particle.
    r_fission_x : float
        Location of fission in the x-coordinate.

    """
    p = Particle(regions, global_bounds, r_x=r_x)
    track_length = {'global': 0}
    N_fission_neutrons = 0
    r_fission_x = None
    collisions = {'global': 0}
    for region_id in regions.keys():
        track_length[region_id] = 0
        collisions[region_id] = 0

    while p.alive:
        d_b = p.distance_to_boundary(regions)
        d_c = p.distance_to_collision(xs)

        if d_b <= d_c:
           p.translate(d_b)
           if tally:
               track_length['global'] += d_b
               track_length[p.region_id] += d_b
           p.get_current_region(regions, global_bounds)
        else:
           p.translate(d_c)
           if tally:
               track_length['global'] += d_c
               track_length[p.region_id] += d_c
               collisions['global'] += 1
               collisions[p.region_id] += 1
           rxn = p.sample_reaction(xs, xs_bins)
           if rxn == 'Sigma_s':
               p.sample_direction()
           elif rxn == 'Sigma_a':
               r_fission_x, N_fission_neutrons = p.sample_fission_neutrons(xs)

    return track_length, collisions, N_fission_neutrons, r_fission_x

#cProfile.run('run()')
run()
