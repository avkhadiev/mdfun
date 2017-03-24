# -*- coding: utf-8 -*-

import numpy as np
import math

class MCMetropolis(object):
    """A class to calculate the SHO velocity autocorrelation function using MC"""

    # creates an instance of a Metropolis MC sampler class
    def __init__(self, name, config):
        np.random.seed()                            # give arg to remove randomness
        self.name       = name
        self.config     = config
        if 'L'          not in config.keys(): config['L']       = 1.0
        if 'omegasq'    not in config.keys(): config['omegasq'] = 1.0
        # temperature is specified as T / Tmax, where Tmax = 1/4 omegasq (L/2)^2
        # (m = kB = 1)
        if 'iniTemp'    not in config.keys(): config['iniTemp'] = 0.25
        if 'npoints'    not in config.keys(): config['npoints'] = 100
        # stepX is put in relative to L, stepP relative to L/T
        # this will multiplied by reference values automatically
        if 'stepP'      not in config.keys(): config['stepP']   = 0.001
        if 'stepX'      not in config.keys(): config['stepX']   = 0.001
        if 'tInt'       not in config.keys(): config['tInt']    = 0.01
        if 'ntimes'     not in config.keys(): config['ntimes']  = 200
        if 'stepcorr'   not in config.keys(): config['stepcorr']= 1.0
        if (config['debug'] == 'True'):
            print "SHO Monte-Carlo may have modified the configuration."
            print "The current configuration is:"
            print( config )
            print "\n"
        # PHYSICS
        self.L          = float(config['L'])        # size of the box
        self.omegasq    = float(config['omegasq'])  # to calculate period T
        self.omega      = math.sqrt( self.omegasq )
        self.tempMax    = 2 * 0.25 * self.omegasq * (self.L / 2)**2
        self.iniTemp    = float(config['iniTemp'])  # in units of T / Tmax
        self.kbT        = self.iniTemp * self.tempMax
        self.T          = 2 * math.pi / self.omega
        self.pRef       = self.L / self.T           # reference momentum
        # SIMULATION
        self.npoints    = int(config['npoints'])    # sample how many?
        self.stepP      = float(config['stepP']) * self.pRef
        self.stepX      = float(config['stepX']) * self.L
        self.tInt       = float(config['tInt'])     # calculate for every  t/T
        self.ntimes     = int(config['ntimes'])     # tot time = tInt * ntimes
        self.stepcorr   = float(config['stepcorr']) # step correction (for tests)
        # DATA
        self.point      = (0, 0)                    # current (x, p)
        self.points     = []                        # array of sampled points
        self.freq       = []                        # how many times is ith point sampled?
        self.trials     = 0                         # steps attempted
        self.newAcc     = 0                         # unique points accepted

    # uses an RNG to sample the first pair (p, x)
    # allowed by the input temperature PE_max / kT, omegasq, and size of box L
    def sampleFirst(self):
        # random position
        corr = math.sqrt( self.iniTemp )
        x  = (np.random.random() - 0.5) * self.L * corr
        # assign momentum compatible with temperature and position
        totEn = 0.5 * self.omegasq * (self.L / 2)**2 * self.iniTemp
        potEn = 0.5 * self.omegasq * x**2
        p  = math.sqrt( 2*(totEn - potEn) )
        # update points
        self.point  = (x, p)
        self.points.append( self.point )            # updates list with x0, p0
        self.freq.append( 1 )                       # x0, p0 occures 1 times
        self.newAcc += 1                            # +1 unique accepts
        self.trials += 1

    # clears all sampled points and counters
    def reset(self):
        self.points   = []
        self.npoints  = 0

    # calculates current ratio of unique accepts to number of attempts
    def calculateRatio(self):
        return self.newAcc / self.trials

    def adjustRatio(self):
        ratio = self.calculateRatio()
        if ratio < 0.5:
        # acceptance too low, decrease the step
            correction = 1 - self.stepcorr * (0.5 - ratio)
            self.stepX = self.stepX * correction
            self.stepP = self.stepP * correction
        elif ratio > 0.5:
        # acceptance too high, inscrease the step
            correction = 1 + self.stepcorr * (ratio - 0.5)
            self.stepX = self.stepX * 1.1
            self.stepP = self.stepP * 1.1

    # returns an array of npoints
    def sampledPoints(self):
        return self.points

    # updates the list of sampled points with the current self.point
    def updatePoints(self):
        if self.point == self.points[-1]:
            self.freq[-1] += 1
        else:
            self.points.append( self.point )
            self.freq.append( 1 )

    # makes a random step in phase space
    # updates respective counters
    def makeStep(self):
        self.trials += 1
        if (np.random.random() > 0.5):
            xSign = 1.
        else:
            xSign = -1.
        if (np.random.random() > 0.5):
            pSign = 1.
        else:
            pSign = -1.
        (xOld, pOld) = self.point
        xNew = xOld + xSign * self.stepX
        pNew = pOld + pSign * self.stepP
        return (xNew, pNew)

    # compute weight exp( - delta V)
    def computePba(self, xNew, pNew):
        (xOld, pOld) = self.point
        oldPE        = 0.5 * self.omegasq * (xOld * xOld)
        oldKE        = 0.5 * self.omegasq * (pOld * pOld)
        newPE        = 0.5 * self.omegasq * (xNew * pNew)
        newKE        = 0.5 * self.omegasq * (pNew * pNew)
        delPE        = newPE - oldPE
        delKE        = newKE - oldKE
        weight       = (delPE + delKE) / self.kbT
        pba          = math.exp( - weight)
        return pba

    # accepts a new trial point and updates correspoding counters / lists
    def accept(self, xTry, pTry):
        self.newAcc += 1
        self.point   = (xTry, pTry)
        self.updatePoints()

    # rejects a new trial point and updates corresponding counters / lists
    def reject(self):
        # append previous point to the list again
        self.updatePoints()

    # makes nsteps new Monte-Carlo steps, saves sampled points in points array
    def evolve(self, nsteps):
        # make nsteps
        for istep in range(nsteps):
            (xTry, pTry) = self.makeStep()
            # ACCEPT or REJECT
            pba = self.computePba(xTry, pTry)
            if pba >= 1:
                self.accept(xTry, pTry)
            else:
                if np.random.random() >= pba:
                    self.accept(xTry, pTry)
                else:
                    self.reject()

    # calculates the normalization factor based on the current sampled points
    def normalization(self):
        # convert to np.arrays
        self.points = np.array( self.points )       # (x, p)
        self.freq   = np.array( self.freq )
        # extract velocities
        v0 = self.points[:, 1]                      # second row (p)
        # integrate
        C = (v0 * self.freq).sum() / self.freq.sum()
        return C






