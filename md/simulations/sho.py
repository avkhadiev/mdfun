# -*- coding: utf-8 -*-

from mdsimulation import NaiveSimulation
import numpy as np
import pylab as pb
import math

class SHO(NaiveSimulation):
    """A class for an MD simulation of a simple harmonic oscillator"""

    # creates an instance of a simulation class
    def __init__(self, config):
        # general configs
        name = "SHO"
        if 'N'          not in config.keys(): config['N']          = 1
        if 'ndim'       not in config.keys(): config['ndim']       = 1
        if 'L'          not in config.keys(): config['L']          = 1.0
        if 'iniTemp'    not in config.keys(): config['iniTemp']    = 1.0
        if 'dt'         not in config.keys(): config['dt']         = 0.001
        if 'sampleint'  not in config.keys(): config['sampleint']  = 100
        if 'relaxtime'  not in config.keys(): config['relaxtime']  = 10.0
        if 'runtime'    not in config.keys(): config['runtime']    = 20.0
        if 'integration' not in config.keys():
            config['integration']='velocityVerlet'
        # SHO-specific configs
        if 'omegasq'    not in config.keys():
            self.omegasq = 1.0
        else:
            self.omegasq = float(config['omegasq'])
        self.omega = math.sqrt( self.omegasq )
        self.kinEnArray = np.array([])
        self.potEnArray = np.array([])
        if (config['debug'] == 'True'):
            print "SHO may have modified the configuration."
            print "The current configuration is:"
            print( config )
            print "\n"
        super(SHO, self).__init__(name, config)
        # to store information for exact solution
        self.iniEn       = 0.0
        self.iniPos = np.zeros([self.ndim, self.N])
        self.iniVel = np.zeros([self.ndim, self.N])
        self.exPosArray   = np.array([])
        self.exVelArray   = np.array([])
        self.exKinEnArray = np.array([])
        self.exPotEnArray = np.array([])
        self.xmax         = 0.0
        self.pmax         = 0.0
        # x(t) = A exp( i omega t) + B exp( - i omega t)
        self.A = 0.0 + 0.0j
        self.B = 0.0 + 0.0j
        self.T = 2 * math.pi / self.omega
        # figure counter for plotting
        self.figures = 0

    # assigns initial velocities and positions
    def initialize(self):
        self.randomPos()
        self.randomVel()

    # save current conditions as initial for the exact solution
    def initializeExact(self):
        self.iniEn = self.en()
        # maximum displacement / momentum (mass = 1)
        self.pmax  = math.sqrt( 2 * self.iniEn )
        self.xmax  = math.sqrt( 2 * self.iniEn / self.omegasq )
        # initial position / velocity
        self.iniPos = self.pos
        self.iniVel = self.vel
        # coefficients in the exact solution
        self.A = self.iniPos / 2. - 1.j * self.iniVel / (2. * self.omega)
        self.B = self.iniPos / 2. + 1.j * self.iniVel / (2. * self.omega)

    # calculates force from the attractor at L = 0
    def force_sho(self):
        return - 1 * self.omegasq * self.pos

    def force(self):
        return self.force_sho()

    # calculates potential energy of the system
    def potEn(self):
        return 0.5 * self.omegasq * (self.pos * self.pos).sum()

    # calculates and stores required quantities
    def recordObservables(self):
        self.kinEnArray = np.append(self.kinEnArray, self.kinEn())
        self.potEnArray = np.append(self.potEnArray, self.potEn())
        self.enArray    = np.append(self.enArray, self.en())
        self.posArray   = np.append(self.posArray, self.pos)
        self.velArray   = np.append(self.velArray, self.vel)

    # NEW STATISTICS METHODS

    # calculates exact solutions based on current sampleTArray
    # (call after state has been evolved)
    def calculateExact(self):
        power      = 1.j * self.omega * self.sampleTArray
        expA       = np.exp(  power )
        expB       = np.exp( -power )
        self.exPosArray = self.A * expA + self.B * expB
        self.exVelArray = 1.j * self.omega * self.A * expA - 1.j * self.omega * self.B * expB
        self.exKinEnArray = 0.5 * (self.exVelArray * self.exVelArray)
        self.exPotEnArray = 0.5 * self.omegasq * (self.exPosArray * self.exPosArray)

    # resets required accumulators
    def resetObservables(self):
        self.kinEnArray   = np.array([])
        self.potEnArray   = np.array([])

    def meanKinEn(self):
        return self.kinEnArray.mean()

    def meanPotEn(self):
        return self.potEnArray.mean()

    # TODO
    # Plotting

    # scatter plot of coordinates dim1 vs dim2 at one time
    def plotPositions(self, dim1 = 0, dim2 = 1):
        self.figures += 1
        pb.figure( self.figures )
        pb.scatter(self.pos[dim1][:], self.pos[dim2][:],
                s=5.0, marker='o', alpha=1.0)
        pb.xlabel("dim %d" % dim1)
        pb.ylabel("dim %d" % dim2)

    def plotTrajectories(self, num = 0, dim = 0):
        self.figures += 1
        f, (ax1, ax2) = pb.subplots(2, sharex=True, sharey=True)
        ax1.set_title("Trajectory, simulation (red) vs. exact (blue)")
        for ax in (ax1, ax2):
            ax.set_xlabel("t/T")
            ax.set_ylabel("x/xmax")
        x = np.reshape(self.posArray, [self.ndim, self.N, self.samplesteps])
        ax1.plot(self.sampleTArray / self.T, x[dim][num][:] / self.xmax, 'r.-')
        # plot exact solution
        exX = np.reshape(self.exPosArray, [self.ndim, self.N, self.samplesteps])
        ax2.plot(self.sampleTArray / self.T, exX[dim][num][:] / self.xmax, 'b')
        # fine-ture figure, make subplots close to each other and hide x
        # ticks for all but bottom plot.
        f.subplots_adjust(hspace=0)
        pb.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)


    # plot phase space of a given number of particles for a given dimension
    def plotPhaseSpace(self, dim=0):
        self.figures += 1
        f, (ax1, ax2) = pb.subplots(2, sharex=True, sharey=True)
        # title is above the top subplot
        ax1.set_title("Phase space, simulation (red) vs. exact (blue)")
        for ax in (ax1, ax2):
            ax.set_xlabel("x/xmax")
            ax.set_ylabel("p/pmax")
            ax.axis('equal')
        x = np.reshape(self.posArray, [self.ndim, self.N, self.samplesteps])
        p = np.reshape(self.velArray, [self.ndim, self.N, self.samplesteps])
        exX = np.reshape(self.exPosArray, [self.ndim, self.N, self.samplesteps])
        exP = np.reshape(self.exVelArray, [self.ndim, self.N, self.samplesteps])
        ax1.plot(x[dim][0][:]   / self.xmax,   p[dim][0][:] / self.pmax, "r.-")
        ax2.plot(exX[dim][0][:] / self.xmax, exP[dim][0][:] / self.pmax, "b.")
        # fine-ture figure, make subplots close to each other and hide x
        # ticks for all but bottom plot.
        f.subplots_adjust(hspace=0)
        pb.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)


    # plot temperature as a function of time
    def plotTemp(self):
        pb.figure(3)
        pb.plot(self.tArray, self.tempArray)
        pb.xlabel("time")
        pb.ylabel("temperature")

    # plot energy as a function fo time
    def plotPotEn(self):
        self.figures += 1
        f, (ax1, ax2) = pb.subplots(2, sharex=True, sharey=True)
        ax1.set_title("Potential energy, simulation (red) vs. exact (blue)")
        for ax in (ax1, ax2):
            ax.set_xlabel("t/T")
            ax.set_ylabel("PE/PE_0")
        ax1.plot(self.sampleTArray / self.T, self.potEnArray / self.iniEn, "r.-")
        ax2.plot(self.sampleTArray / self.T, self.exPotEnArray.transpose() / self.iniEn, "b")

    def printResults(self):
        print "\nRESULTS \n"
        print "time = %.3f\n"          % self.t
        print "total energy = %.3f\n"  % self.en()
        print "temperature = %.3f\n"   % self.temp()
        if (self.steps > 0):
            print "mean energy %.3f\n" % self.meanEn()
            print "std energy %.3f\n"  % self.stdEn()

    def run(self, runtime = 0.0, relaxtime = 0.0):
        if runtime   == 0.0: runtime   = self.runtime
        if relaxtime == 0.0: relaxtime = self.relaxtime
        self.initialize()
        self.evolve( relaxtime )
        self.resetStatistics()
        if self.config['debug'] == 'True':
            print "Number of sample steps after reset: ", self.samplesteps
        self.initializeExact()
        self.evolve( runtime )
        self.calculateExact()
        self.printResults()
        self.plotPotEn()
        self.plotTrajectories()
        self.plotPhaseSpace()
        self.showPlots()

