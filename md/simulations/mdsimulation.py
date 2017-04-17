# -*- coding: utf-8 -*-

import numpy as np
import pylab as pb
import math

class NaiveSimulation(object):
    """A base class for an MD simulation"""

    # creates an instance of a simulation class
    def __init__(self, name, config):
        # configuration
        np.random.seed()                          # RNG for velocities/positions
        self.name        = name
        self.config      = config
        self.integration = config['integration']
        self.batchmode   = config['batchmode']    # don't show
        # environment
        self.ndim        = int(config['ndim'])
        self.L           = float(config['L'])
        # time
        self.steps       = 0
        self.samplesteps = 0
        self.t           = 0.0
        self.dt          = float(config['dt'])
        self.dtsq        = self.dt ** 2
        self.sampleint   = float(config['sampleint'])
        self.runtime     = float(config['runtime'])
        self.relaxtime   = float(config['relaxtime'])
        self.tArray      = np.array([self.t])
        self.sampleTArray= np.array([])
        #objects
        self.N           = int(config['N'])
        self.pos         = np.zeros([self.ndim, self.N])
        self.vel         = np.zeros([self.ndim, self.N])
        self.acc         = np.zeros([self.ndim, self.N])
        self.velArray    = np.array([])
        self.posArray    = np.array([])
        # observables
        self.iniTemp     = float(config['iniTemp'])
        self.tempArray   = np.array([])
        self.tempAcc     = 0.0
        self.tempsqAcc   = 0.0
        self.enArray     = np.array([])

    # randomizes positions as (rand() - 1/2) * L
    # may produce overlap for large N!
    def randomPos(self):
        self.pos = (np.random.random([self.ndim, self.N]) - 0.5) * self.L

    # randomizes velocities as (rand() - 1/2) * L
    def randomVel(self):
        self.vel  = (np.random.random([self.ndim, self.N]) - 0.5) * self.L
        # T = self.temp()
        # self.vel *= math.sqrt( self.iniTemp / T )

    # zeros total momentum along each dimension
    def zeroMomentum(self):
        for dim in range(self.ndim):
            self.vel[dim][:] -= self.vel[dim][:].mean()

    # assigns initial velocities and positions, ensures zero momentum
    def initialize(self):
        self.randomPos()
        self.randomVel()
        self.zeroMomentum()

    # calculates sum of kinetic energies of each particle
    def kinEn(self):
        return 0.5 * (self.vel * self.vel).sum()

    # calculates potential energy of the system
    def potEn(self):
        pass

    # calculates ``kinetic'' temperature of the system
    def temp(self):
        return 2 * self.kinEn() / ( self.ndim * self.N )

    # TODO calculate pressure via virial theorem
    def pressure(self):
       pass

    def en(self):
        return self.kinEn() + self.potEn()

    # calculates forces on each particle
    def force(self):
        pass

    # Velocity Verlet
    # v(t + 0.5 dt) = v(t) + 1/2 dt a(t)
    # x(t + dt)     = x(t) + dt v(t + 0.5 dt)
    # a(t + dt)     = 1/m F( x(t + dt))
    # v(t)          = v(t + 0.5 dt) + 1/2 dt a(t + dt)
    def velocityVerlet(self):
        self.vel += 0.5 * self.dt * self.acc
        self.pos += self.dt * self.vel
        self.acc  = self.force()
        self.vel += 0.5 * self.dt * self.acc

    # TODO
    # Position Verlet
    # x(t + 0.5 dt) = x(t) + 1/2 dt v(t)
    # v(t + dt)     = v(t) + dt a(t + 0.5 dt)
    # x(t)          = x(t + 0.5 dt) + 1/2 dt v(t + dt)
    def positionVerlet(self):
        pass

    # TODO
    def verlet(self):
        pass

    # TODO
    def gearPC(self, iterations=0.0):
        if iterations == 0.0: iterations = self.config['iterations']
        pass

    # integrates EOMs
    def integrate(self):
        if   self.integration == "velocityVerlet":
            self.velocityVerlet()
        elif self.integration == "positionVerlet":
            self.positionVerlet()
        elif self.integration == "verlet":
            self.verlet()
        elif self.integration == "gearPC":
            self.gearPC()

    # calculates and stores requires quantities
    def recordObservables(self):
        pass

    # evolves the system a specified amount of time
    def evolve(self, time):
        steps = int(abs(time/self.dt))
        self.acc = self.force()
        for istep in xrange(steps):
            self.integrate()
            self.t      += self.dt
            self.tArray  = np.append(self.tArray, self.t)
            if (istep % self.sampleint == 0):
                self.recordObservables()
                self.sampleTArray = np.append(self.sampleTArray, self.t)
                self.samplesteps += 1
                T = self.temp()
                self.tempArray  = np.append(self.tempArray, T)
                self.tempAcc   += T
                self.tempsqAcc += T**2
            self.steps     += 1

    def reverseTime(self):
        self.dt = -self.dt


    # STATISTICS
    def resetObservales(self):
        pass

    def resetStatistics(self):
        self.steps       = 0
        self.samplesteps = 0
        self.sampleTArray= np.array([])
        self.posArray    = np.array([])
        self.velArray    = np.array([])
        self.enArray     = np.array([])
        self.tempArray   = np.array([])
        self.t           = 0.0
        self.tempAcc     = 0.0
        self.tempsqAcc   = 0.0
        self.resetObservables()

    def meanTemp(self):
        return self.tempAcc / self.samplesteps

    def meanTempsq(self):
        return self.tempscAcc / self.samplesteps

    def meanEn(self):
        return self.enArray.mean()

    def stdEn(self):
        return self.enArray.std()

    # OUTPUT

    def showPlots(self):
        pb.show()

    # writeout observables to file
    def writeObservables(self):
        pass

    def __str__(self):
        return ("%s MD simulation" % self.name)

