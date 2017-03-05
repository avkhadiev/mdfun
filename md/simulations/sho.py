# -*- coding: utf-8 -*-

from mdsimulation import NaiveSimulation
import numpy as np

class SHO(NaiveSimulation):
    """A class for an MD simulation of a simple harmonic oscillator"""

    # creates an instance of a simulation class
    def __init__(self, config):
        # general configs
        name = "SHO"
        if 'N'          not in config.keys(): config['N']          = 1
        if 'ndim'       not in config.keys(): config['ndim']       = 1
        if 'L'          not in config.keys(): config['L']          = 1.0
        if 'iniTemp'    not in config.keys(): config['iniTemp']    = 0.0
        if 'dt'         not in config.keys(): config['dt']         = 0.01
        if 'sampletime' not in config.keys(): config['sampletime'] = 0.01
        if 'relaxtime'  not in config.keys(): config['relaxtime']  = 10.0
        # SHO-specific configs
        if 'omegasq'    not in config.keys():
            self.omega = 1.0
        else:
            self.omega = float(config['omegasq'])
        self.kinEnArray = np.array([])
        self.potEnArray = np.array([])
        super(SHO, self).__init__(name, config)

    # calculates force from the attractor at L = 0
    def force_sho(self):
        return - 1 * self.omegasq * self.pos

    def force(self):
        return self.force_sho(self)


    # calculates potential energy of the system
    def potEn(self):
        return 0.5 * (self.pos * self.pos)

    # calculates and stores required quantities
    def recordObservables(self):
        self.kinEnArray = np.append(self.kinEnArray, self.kinEn())
        self.potEnArray = np.append(self.potEnArray, self.potEn())
        self.posArray   = np.append(self.posArray, self.pos)
        self.velArray   = np.append(self.velArray, self.vel)

    # resets required quantities
    def resetObservables(self):
        self.steps      = 0
        self.tempAcc    = 0.0
        self.tempSqAcc  = 0.0
        self.posArray   = np.array([])
        self.kinEnArray = np.array([])
        self.potEnArray = np.array([])
        self.posArray   = np.array([])
        self.velArray   = np.array([])

    # TODO
    # Plotting



