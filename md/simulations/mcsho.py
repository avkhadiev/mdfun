# -*- coding: utf-8 -*-

import numpy as np

class MCMetropolis(object):
    """A class to calculate the SHO velocity autocorrelation function using MC"""

    # creates an instance of a Metropolis MC sampler class
    def __init__(self, name, config):
        np.random.seed(160)
        self.name       = name
        self.config     = config
        if 'npoints'    not in config.keys(): config['npoints'] = 100
        if 'nrlxsteps'  not in config.keys(): config['nrlxsteps'] = 50
        if 'stepP'      not in config.keys(): config['stepP']   = 0.01
        if 'stepX'      not in config.keys(): config['stepX']   = 0.01
        if 'L'          not in config.keys(): config['L']       = 1.0
        if 'omegasq'    not in config.keys(): config['omegasq'] = 1.0
        if 'enMax'      not in config.keys(): config['enMax']   = 100
        if 'tInt'       not in config.keys(): config['tInt']    = 0.01
        if 'ntimes'     not in config.keys(): config['ntimes']  = 200
        if (config['debug'] == 'True'):
            print "SHO Monte-Carlo may have modified the configuration."
            print "The current configuration is:"
            print( config )
            print "\n"
        self.npoints    = int(config['npoints'])    # sample how many?
        self.nrlxsteps  = int(config['nrlxsteps'])  # num steps before sampling
        self.stepP      = float(config['stepP'])    # relative to ...
        self.stepX      = float(config['stepX'])    # relative to ...
        self.L          = float(config['L'])        # size of the box
        self.omegasq    = float(config['omegasq'])  # to calculate period T
        self.enMax      = float(config['enMax'])    # max energy / kT
        self.tInt       = float(config['tInt'])     # calculate for every  t/T
        self.ntimes     = int(config['ntimes'])     # tot time = tInt * ntimes
        self.points     = np.array([])              # array of sampled points
        self.trials     = 0                         # steps attempted
        self.newAcc     = 0                         # unique points accepted

    # uses an RNG to sample the first pair (p, x)
    # allowed by the input temperature PE_max / kT, omegasq, and size of box L
    def sampleFirst(self):
        pass

    # clears all sampled points and counters
    def reset(self):
        pass

    # calculates current ratio of unique accepts to number of attempts
    def calculateRatio(self):
        pass

    # calculates the normalization factor based on the current sampled points
    def normalization(self):
        pass

    # returns an array of npoints
    def sampledPoints(self):
        return self.points

    # makes nsteps new Monte-Carlo steps, saves sampled points in points array
    def evolve(self, nsteps):
        pass





