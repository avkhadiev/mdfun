# -*- coding: utf-8 -*-

from mdsimulation import NaiveSimulation
import numpy as np
import pylab as pb
import math
import os

class SHO(NaiveSimulation):
    """A class for an MD simulation of a simple harmonic oscillator"""

    # creates an instance of a simulation class
    def __init__(self, name, config):
        # general configs
        if 'N'          not in config.keys(): config['N']          = 1
        if 'ndim'       not in config.keys(): config['ndim']       = 1
        if 'L'          not in config.keys(): config['L']          = 1.0
        # temperature is specified as T / Tmax, where Tmax = 1/4 omegasq (L/2)^2
        # (m = kB = 1)
        if 'iniTemp'    not in config.keys(): config['iniTemp']    = 0.01
        if 'batchmode'  not in config.keys(): config['batchmode']  = 'False'
        if 'writeout'   not in config.keys(): config['writeout']   = 'False'
        if 'calculateExact' not in config.keys():
            config['calculateExact'] = 'True'
        if 'setInitial' not in config.keys():
            config['setInitial'] = 'False'
        else:
            if ('x0' not in config.keys()) or ('p0' not in config.keys()):
                config['x0'] = float(config['L']) / 4
                config['p0'] = 0.0
        if 'integration' not in config.keys():
            config['integration']='velocityVerlet'
        if 'outdir' not in config.keys():
            config['outdir'] = 'sho_output'
        self.outdir = config['outdir']
        # SHO-specific configs
        if 'omegasq'    not in config.keys():
            self.omegasq = 1.0
        else:
            self.omegasq = float(config['omegasq'])
        self.omega = math.sqrt( self.omegasq )
        self.kinEnArray = np.array([])
        self.potEnArray = np.array([])
        # x(t) = A exp( i omega t) + B exp( - i omega t)
        self.A = 0.0 + 0.0j
        self.B = 0.0 + 0.0j
        self.T = 2 * math.pi / self.omega
        # intervals are specified in units of the period, T
        if 'relaxtime'  not in config.keys():
            config['relaxtime'] = 1 * self.T
        else:
            config['relaxtime'] = float(config['relaxtime']) * self.T
        if 'runtime'    not in config.keys():
            config['runtime']   = 4 * self.T
        else:
            config['runtime'] = float(config['runtime']) * self.T
        if 'dt'         not in config.keys():
            config['dt']        = 0.005 * self.T
        else:
            config['dt'] = float(config['dt']) * self.T
        # sample observables every n steps
        if 'sampleint'  not in config.keys():
            config['sampleint'] = 10
        else:
            config['sampleint'] = int(config['sampleint'])
        self.relaxtime = float( config['relaxtime'] )
        self.runtime = float( config['runtime'] )
        if (config['debug'] == 'True'):
            print "SHO may have modified the configuration."
            print "The current configuration is:"
            print( config )
            print "\n"
        super(SHO, self).__init__(name, config)
        # maximum INSTANTENEOUS temperature
        self.tempMax      = 2 * 0.25  * self.omegasq * (self.L / 2)**2
        self.xmax         = 0.0
        self.pmax         = 0.0
        # to store information for exact solution
        self.iniEn       = 0.0
        self.iniPos = np.zeros([self.ndim, self.N])
        self.iniVel = np.zeros([self.ndim, self.N])
        self.exPosArray   = np.array([])
        self.exVelArray   = np.array([])
        self.exKinEnArray = np.array([])
        self.exPotEnArray = np.array([])
        self.relaxtime = float( config['relaxtime'] )
        self.runtime = float( config['runtime'] )
        # figure counter for plotting
        self.figures = 0

    def randomPos(self):
        # T/Tmax = iniTemp, where Tmax = 2 * 2 * 1/4 omegasq (L/2)^2
        corr = math.sqrt( self.iniTemp )
        if self.config['debug'] == 'True':
            print "temperature scaling correction ", corr
        self.pos = (np.random.random([self.ndim, self.N]) - 0.5) * self.L * corr

    # assigns initial velocities and positions
    def initialize(self):
        self.randomPos()
        # initial velocity is determined by total energy + initial position
        totEn = 0.5 * self.omegasq * (self.L / 2)**2 * self.iniTemp
        self.vel[0][:] = np.array( [(1 - 2 * (np.random.random() < 0.5)) * math.sqrt(2*(totEn - self.potEn()))] )
        if self.config['debug'] == 'True':
            print 'position after initialization', self.pos
            print 'velocity after initialization', self.vel
            print 'potential energy after initialization', self.potEn()
            print 'kinetic energy after initialization',   self.kinEn()
            print 'total energy from temperature scaling ', totEn
            print 'instanteneous temperature after initialization ', self.temp()

    # assign initial position and velocity to the particle
    def setInitial(self, pos, vel):
        # check conservation of energy
        energy    = 0.5 * self.omegasq * pos ** 2 + 0.5 * vel ** 2
        maxEnergy = 0.5 * self.omegasq * (self.L / 2) ** 2
        if energy > maxEnergy:
            return 1            # can't initialize these coordinates
        else:
            self.pos[0][:] = np.array([ pos ])
            self.vel[0][:] = np.array([ vel ])
            if self.config['debug'] == 'True':
                print 'position after initialization', self.pos
                print 'velocity after initialization', self.vel
                print 'potential energy after initialization',  self.potEn()
                print 'kinetic energy after initialization',    self.kinEn()
                print 'total energy from temperature scaling ', self.en()
                print 'instanteneous temperature after initialization ', self.temp()
            return 0            # initializaton went okay

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

    # Plotting

    # scatter plot of coordinates dim1 vs dim2 at one time
    def plotPositions(self, dim1 = 0, dim2 = 1):
        self.figures += 1
        pb.figure( self.figures )
        pb.scatter(self.pos[dim1][:], self.pos[dim2][:],
                s=5.0, marker='o', alpha=1.0)
        pb.xlabel("dim %d" % dim1)
        pb.ylabel("dim %d" % dim2)

    def plotTrajectories(self, overlayExact = 0, num = 0, dim = 0):
        self.figures += 1
        pb.figure( self.figures )
        pb.title("Trajectory, simulation (red) vs. exact (blue)")
        pb.xlabel("t/T")
        pb.ylabel("x/xmax")
        pb.grid(True)
        x = np.reshape(self.posArray, [self.ndim, self.N, self.samplesteps])
        pb.plot(self.sampleTArray / self.T, x[dim][num][:] / self.xmax,
                linestyle = 'None', marker = 'o', color = 'r')
        # plot exact solution
        if overlayExact:
            exX = np.reshape(self.exPosArray, [self.ndim, self.N, self.samplesteps])
            pb.plot(self.sampleTArray / self.T, exX[dim][num][:] / self.xmax,
                    linestyle = '-', color = 'b')
        # pb.savefig('trajectories.png')

    # plot phase space of a given number of particles for a given dimension
    def plotPhaseSpace(self, overlayExact = 0, dim=0):
        self.figures += 1
        pb.figure( self.figures )
        pb.title("Phase space, simulation (red) vs. exact (blue)")
        pb.xlabel("x/xmax")
        pb.ylabel("p/pmax")
        pb.grid(True)
        pb.axis('equal')
        x = np.reshape(self.posArray, [self.ndim, self.N, self.samplesteps])
        p = np.reshape(self.velArray, [self.ndim, self.N, self.samplesteps])
        pb.plot(x[dim][0][:]   / self.xmax,   p[dim][0][:] / self.pmax,
                    linestyle = 'None', marker = 'o', color = 'r')
        if overlayExact:
            exX = np.reshape(self.exPosArray, [self.ndim, self.N, self.samplesteps])
            exP = np.reshape(self.exVelArray, [self.ndim, self.N, self.samplesteps])
            pb.plot(exX[dim][0][:] / self.xmax, exP[dim][0][:] / self.pmax,
                    linestyle = '-', color = 'b')
        # pb.savefig('phaseSpace.png')

    # plot temperature as a function of time
    def plotTemp(self):
        self.figures += 1
        pb.figure( self.figures )
        pb.plot(self.sampleTArray / self.T , self.tempArray / self.tempMax,
                linestyle = '-', marker = 'o', color = 'r')
        pb.axhline(y = self.meanTemp() / self.tempMax,
                linestyle = '--',  color='b')
        pb.axhline(y = self.iniTemp,
                linestyle = '--',  color='k')
        pb.title("Instanteneous temperature (red) and Mean Temperature (blue)")
        pb.xlabel("t/T")
        pb.ylabel("Inst Temp / Max Temp")
        pb.grid(True)

    # plot energy as a function fo time
    def plotPotEn(self, overlayExact = 0):
        self.figures += 1
        pb.figure( self.figures )
        pb.title("Potential energy, simulation (red) vs. exact (blue)")
        pb.xlabel("t/T")
        pb.ylabel("PE/PE_0")
        pb.grid(True)
        pb.plot(self.sampleTArray / self.T, self.potEnArray / self.iniEn,
                linestyle = 'None', marker = 'o', color = 'r')
        if overlayExact:
            pb.plot(self.sampleTArray / self.T, self.exPotEnArray.transpose() / self.iniEn,
                    linestyle = '-', color = 'b')
        # TODO fix to include absolute / customizable directories
        # pb.savefig('potEn.png')

    def printResults(self):
        print "\nRESULTS \n"
        print "time = %.4f\n"          % self.t
        print "maximum possible temperature %.4f\n" % self.tempMax
        if (self.steps > 0):
            print "mean energy %.4f\n" % self.meanEn()
            print "std energy %.4f\n"  % self.stdEn()
            print "mean temperature = %.4f\n"   % self.meanTemp()

    # writeout observables to file
    def writeObservables(self):
        # go to output directory and save output in a subfolder
        fname = self.name
        os.chdir( self.outdir )
        subfold = os.path.join(self.outdir, fname)
        os.mkdir(subfold)
        os.chdir(subfold)
        np.savetxt('%s_%s.txt' % (fname, 'vel'),  self.velArray    )
        np.savetxt('%s_%s.txt' % (fname, 'samT'), self.sampleTArray)
