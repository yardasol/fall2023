import random
import math

class Particle:
    """Class representing a particle

    Attributes
    ----------
    mu, mu_x : float
        cosine and sine of particle angle w.r.t
        the z-axis

    r_x, r_y : float
        Particle position in the x and y coordinates

    region_id : int
        Region ID

    alive : bool
        Flag to indicate if the particle is alive or not.

    """

    def __init__(self, regions, global_bounds, r_x=None):
        """Initalize new particle"""
        self.sample_direction()
        self.sample_position(global_bounds, r_x=r_x)
        self.get_current_region(regions, global_bounds)

        self.alive = True

    def sample_direction(self):
        """Sample a direction for the particle"""
        # Sample between 0 and 2\pi direction
        sign_x = random.randrange(-1, 2, 2)
        mu = random.uniform(-1, 1)
        self.mu_x = sign_x * math.sin(math.acos(mu))
        self.mu = mu

    def sample_position(self, global_bounds, r_x=None):
        """Sample a position for the particle

        Parameters
        ----------
        global_bounds : list of float
            Global geometry boundaries

        r_x : float
            Pre-selected position to sample. We only sample
            the x-coordinate and assume all y-coordinates start
            at 0
        """

        if r_x is None:
            xmin, xmax = global_bounds
            r_x = random.uniform(xmin, xmax)

        self.r_x = r_x
        self.r_y = 0

    def get_current_region(self, regions, global_bounds):
        """Get ID of current region based on
        particle position

        Parameters
        ----------
        regions : dict of int to tuple
            Dictionary mapping region ID to region bounds
        """
        # Check if particle is leaving defined geometry
        if ((self.r_x <= global_bounds[0] and self.mu_x < 0) or
                self.r_x >= global_bounds[1] and self.mu_x > 0):
            self.kill()
            return

        # TODO: add machinery to address situation where a particle is born
        # on a region boundary
        # Naive algorithm, but okay for now
        for region_id, bounds in regions.items():
            if bounds[0] < self.r_x < bounds[1]:
                self.region_id = region_id
            elif bounds[1] == self.r_x and self.mu_x < 0:
                self.region_id = region_id
            elif bounds[0] == self.r_x and self.mu_x > 0:
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
        if self.mu_x < 0:
            idx = 0
        elif self.mu_x > 0:
            idx = 1
        else:
            return math.inf
        d_x = abs(self.r_x - regions[self.region_id][idx])
        d = d_x / math.cos((math.pi/2) - math.acos(self.mu))

        return d

    def translate(self, d):
        """Translate the particle by d in the current particle direction"""
        d_x = d * self.mu_x
        d_y = d * self.mu
        self.r_x += d_x
        self.r_y += d_y

    def kill(self):
        """Kill the particle"""
        self.alive = False

    def sample_fission_neutrons(self, xs):
        self.kill()
        N_fission_neutrons = xs[self.region_id]['nuSigma_f']/xs[self.region_id]['Sigma_a']

        # Integer component
        I = int(N_fission_neutrons)
        # Remainder
        R = N_fission_neutrons - I

        # Sample if we round the remainder up or down
        xi = random.random()
        if xi <= R:
            N_fission_neutrons = I + 1
        else:
            N_fission_neutrons = I
        return self.r_x, N_fission_neutrons



