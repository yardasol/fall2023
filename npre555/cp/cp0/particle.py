import random
import math

class Particle:
    """Class representing a particle

    Attributes
    ----------
    omega : list of float
        Particle direction

    r : list of float
        Particle position

    region_id : int
        Region ID

    alive : bool
        Flag to indicate if the particle is alive or not.

    """

    def __init__(self, regions, global_bounds, r=None):
        """Initalize new particle"""
        self.sample_direction()
        self.sample_position(global_bounds, r=r)
        self.get_current_region(regions, global_bounds)

        self.alive = True

    def sample_direction(self):
        """Sample a binary direction for the particle"""
        # Sample -1 (-x) or 1 (+x) direction
        omega_x = random.randrange(-1, 2, 2)
        self.omega = omega_x

    def sample_position(self, global_bounds, r=None):
        """Sample a position for the particle

        Parameters
        ----------
        global_bounds : list of float
            Global geometry boundaries

        r : list of float
            Pre-selected position to sample
        """

        if r is None:
            xmin, xmax = global_bounds
            r = random.uniform(xmin, xmax)

        self.r = r

    def get_current_region(self, regions, global_bounds):
        """Get ID of current region based on
        particle position

        Parameters
        ----------
        regions : dict of int to tuple
            Dictionary mapping region ID to region bounds
        """
        # Check if particle is leaving defined geometry
        if ((self.r == global_bounds[0] and self.omega == -1) or
                self.r == global_bounds[1] and self.omega == 1):
            self.kill()
            return

        # TODO: add machinery to address situation where a particle is born
        # on a region boundary
        # Naive algorithm, but okay for now
        for region_id, bounds in regions.items():
            if bounds[0] < self.r < bounds[1]:
                self.region_id = region_id
            elif bounds[1] == self.r and self.omega == -1:
                self.region_id = region_id
            elif bounds[0] == self.r and self.omega == 1:
                self.region_id = region_id


    def sample_reaction(self, xs, xs_bins):
        """Determine reaction at collision"""
        sample = random.random()
        for rxn_name, rxn_bin in xs_bins[self.region_id].items():
            if rxn_bin[0] <= sample < rxn_bin[1]:
                rxn = rxn_name
        return rxn

    def distance_to_collision(self, xs):
        xi = random.random()
        return math.log(1 - xi)/-xs[self.region_id]['Sigma_t']

    def distance_to_boundary(self, regions):
        if self.omega == -1:
            idx = 0
        else:
            idx = 1
        return abs(self.r - regions[self.region_id][idx])

    def translate(self, r):
        """Translate the particle by r in direction omega"""
        self.r += r * self.omega

    def kill(self):
        """Kill the particle"""
        self.alive = False

    def sample_fission_neutrons(self, xs):
        self.kill()
        fission_neutrons = xs[self.region_id]['nuSigma_f']/xs[self.region_id]['Sigma_a']
        I = int(fission_neutrons)
        R = fission_neutrons - I
        # Sample if we round the remainder up or down
        xi = random.random()
        if xi <= R:
            fission_neutrons = I + 1
        else:
            fission_neutrons = I
        return self.r, fission_neutrons



