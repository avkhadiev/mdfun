# -*- coding: utf-8 -*-

from mdsimulation import NaiveSimulation
import numpy as np
import pylab as pb
import math
import os

class Barbell(NaiveSimulation):
    """A class that describes a (freely) rotating rod (around the CM)"""

    # creates an instance of a simulation class with the given name
    # sets up the variables for the exact solution
    def __init__(self, name, config):
        ############### INPUT PARAMETERS ####################
        # specify omeganotsq (in the presence of a gravitational field)
        # specify initial position and velocity
        # calculate "natural" period T
        # specify mass m, length l
        # put in some initial temperature in the config (it won't be used)
        #####################################################
        if 'omeganotsq'    not in config.keys():
            # \omega_0^2 = \frac{g}{I}
            # by default, there is no gravitational field
            self.omeganotsq = 0.0
        else:
            self.omeganotsq = float(config['omeganotsq'])
        self.omeganot = math.sqrt( self.omeganotsq )
        if self.omeganot == 0:
            # the pendulum has no period, assign a number to get a reasonable
            # number of steps
            self.T = 1.0
        else:
            self.T        = 2 * math.pi / self.omeganot
        # thetanot, in units of pi
        if 'thetanot'       not in config.keys():
            self.thetanot = 0.01
        else:
            self.thetanot = float(config['thetanot'])
        # thetanotdot, in units of pi / sec
        if 'thetanotdot'    not in config.keys():
            if self.omeganot == 0:
                # if there is no gravitational field, add some velocity
                # to have dynamics
                self.thetanotdot = 0.2
            else:
                # if there is a gravitational field, don't add angular velocity
                self.thetanotdot = 0.0
        else:
            self.thetanotdot = float(config['thetanotdot'])
        # mass of one end of the barbell
        if 'm'             not in config.keys():
            self.m = 1.0
        else:
            self.m = float(config['m'])
        # length of the barbell
        if 'l'             not in config.keys():
            self.l = 1.0
        else:
            self.l = float(config['l'])
        # initial temperature won't be used and is only required for the
        # base class initialization method
        if 'iniTemp'    not in config.keys(): config['iniTemp']    = 0.0
        ############### SIMULATION SETTINGS ################
        # specify dt in units of natural period (2pi / omeganot)
        # put in some relaxation time (it won't be used)
        # ndim, boxsize L (won't be used), number of "particles" N = 1
        # batchmode, exact solution, writeout, integration method
        # initialize with base class initialization method
        #####################################################
        # general configurations
        if 'N'          not in config.keys(): config['N']          = 1
        # the rotation is only around one axis => 1 dof, 1 dimension
        if 'ndim'       not in config.keys(): config['ndim']       = 1
        if 'L'          not in config.keys(): config['L']          = 1.0
        if 'iniTemp'    not in config.keys(): config['iniTemp']    = 0.0
        if 'batchmode'  not in config.keys(): config['batchmode']  = 'False'
        if 'writeout'   not in config.keys(): config['writeout']   = 'True'
        if 'calculateExact' not in config.keys():
            config['calculateExact'] = 'False'
        if 'integration'    not in config.keys():
            config['integration']='velocityVerlet'
        if 'runtime'        not in config.keys():
            config['runtime']   = 4 * self.T
        else:
            config['runtime'] = float(config['runtime']) * self.T
        if 'dt'             not in config.keys():
            config['dt']        = 0.005 * self.T
        else:
            config['dt'] = float(config['dt']) * self.T
        # sample the observables every sampleint steps
        if 'sampleint'      not in config.keys():
            # sample no more than a 100 points per period
            config['sampleint'] = max(1, int((self.T/50)/float(config['dt'])))
        else:
            config['sampleint'] = int(config['sampleint'])
        # the following won't be used;
        # only needed for the base class initialization and running via script
        config['setInitial'] = 'False'
        config['relaxtime'] =  0.0
        # figure counter for plotting
        self.figures = 0
        # print out modified configuration
        if (config['debug'] == 'True'):
            print "SHO may have modified the configuration."
            print "The current configuration is:"
            print( config )
            print "\n"
        super(Barbell, self).__init__(name, config)
        ###############  EXACT SOLUTION  ###################
        # prepare storage for exact solution:
        # initial position and velocity, initial energy
        # position array, velocity array
        # kinetic energy array, potential energy array
        ####################################################
        self.exEn         = 0.0
        self.exPos        = np.zeros([self.ndim, self.N])
        self.exVel        = np.zeros([self.ndim, self.N])
        self.exPosArray   = np.array([])
        self.exVelArray   = np.array([])
        self.exKinEnArray = np.array([])
        self.exPotEnArray = np.array([])

    # returns the angular acceleration for the numerical integration step
    # in the absence of a gravitational field
    def alpha_nograv(self):
        return 0.0

    # returns the angular acceleration for the numerical integration step
    # in the presence of a gravitational field
    def alpha_grav(self):
        # acc(t) = 2*omeganotsq*sin(pos(t)) - 2*(vel(t-0.5*dt)*cos(pos(t))
        # calculate alpha from torque
        alpha =  2 * self.omeganotsq * math.sin( self.pos * math.pi)
        # subtract the component along the barbell
        alpha -= 2 * (self.vel * math.cos( self.pos * math.pi ))
        return -alpha

    # standard force method; returns some acceleration (in this case, angular)
    def force(self):
        if self.omeganotsq == 0.0:
            # no gravitational field
            return self.alpha_nograv()
        else:
            # gravitational field is present
            return self.alpha_grav()

    # once initial coordinates in phase space have been specified,
    # saves initialized velocities and positions for the exact solution
    # to be called after relaxation
    def initializeExact(self):
        # FIXME
        self.iniEn = self.en()

    # assigns initial velocity and position
    # calculates the initial acceleration
    # alpha(0) = 2 * mg * (l/2) * sin( theta(0) ) / (2 * m * (l/2)^2)
    #          = 2 g/l * sin( theta(0) )
    def initialize(self):
        ############################################
        # assign position and velocity arrays
        # calculate acceleration at time 0
        ############################################
        self.pos = np.array([ self.thetanot    ])
        self.vel = np.array([ self.thetanotdot ])
        self.acc = - 2. * self.omeganotsq * math.sin(self.thetanot * math.pi)
        if self.config['debug'] == 'True':
            print 'position after initialization', self.pos
            print 'velocity after initialization', self.vel
            print 'potential energy after initialization', self.potEn()
            print 'kinetic energy after initialization',   self.kinEn()
            print 'instanteneous temperature after initialization ', self.temp()

    # redefine kinetic energy for the rotating rod
    def kinEn(self):
        # m l^2 / 4 * vel^2
        return self.m * (self.l**2) / 4 * (self.vel * self.vel).sum()

    # no potential energy for a freely rotating rod
    def potEn(self):
        # 0, if omeganotsq = 0
        if self.omeganotsq == 0:
            # there is no gravitational field
            return 0
        else:
            # position is a one-dimensional array
            return self.omeganotsq * (1 -  math.cos( self.pos * math.pi ))

    # calculates and stores required quantities
    def recordObservables(self):
        self.kinEnArray = np.append(self.kinEnArray, self.kinEn())
        self.potEnArray = np.append(self.potEnArray, self.potEn())
        self.enArray    = np.append(self.enArray, self.en())
        self.posArray   = np.append(self.posArray, self.pos)
        self.velArray   = np.append(self.velArray, self.vel)

    # resets required accumulators
    def resetObservables(self):
        self.kinEnArray   = np.array([])
        self.potEnArray   = np.array([])

    ############## STATISTICS ################

    # calculates exact solutions based on initial conditions
    def calculateExact(self):
        pass

    # returns mean kinetic energy from the array of all sampled kin energies
    def meanKinEn(self):
        return self.kinEnArray.mean()

    # returns mean potential energy from the array of all sampled pot energies
    def meanPotEn(self):
        return self.potEnArray.mean()

    ############## OUTPUTTING RESULTS ################

    def printResults(self):
        print "\nRESULTS \n"
        print "time = %.4f\n"                   % self.t
        if (self.steps > 0):
            print "mean energy %.4f\n"          % self.meanEn()
            print "std energy %.4f\n"           % self.stdEn()
            print "mean temperature = %.4f\n"   % self.meanTemp()

    # writes all required observables to a text file in a specified folder
    def writeObservables(self):
        fname = self.name
        # make a directory and go there
        os.makedirs(fname)
        os.chdir(os.path.join( os.getcwd(), fname ))
        # save arrays as text files
        np.savetxt('%s_%s.txt' % (fname, 'ken'),  self.kinEnArray  )
        np.savetxt('%s_%s.txt' % (fname, 'pen'),  self.potEnArray  )
        np.savetxt('%s_%s.txt' % (fname, 'en'),   self.enArray     )
        np.savetxt('%s_%s.txt' % (fname, 'pos'),  self.posArray    )
        np.savetxt('%s_%s.txt' % (fname, 'vel'),  self.velArray    )
        np.savetxt('%s_%s.txt' % (fname, 'temp'), self.tempArray   )
        np.savetxt('%s_%s.txt' % (fname, 'samT'), self.sampleTArray)

    def plotTrajectories(self, overlayExact = 0, num = 0, dim = 0):
        self.figures += 1
        pb.figure( self.figures )
        pb.title("Trajectory, simulation (red) vs. exact (blue)")
        pb.xlabel("t/T_0")
        pb.ylabel("theta / pi")
        pb.grid(True)
        x = np.reshape(self.posArray, [self.ndim, self.N, self.samplesteps])
        pb.plot(self.sampleTArray / self.T, x[dim][num][:],
                linestyle = 'None', marker = 'o', color = 'r')
        # plot exact solution
        # if overlayExact:
        #    exX = np.reshape(self.exPosArray, [self.ndim, self.N, self.samplesteps])
        #    pb.plot(self.sampleTArray / self.T, exX[dim][num][:] / self.xmax,
        #            linestyle = '-', color = 'b')
        # pb.savefig('trajectories.png')


    def plotPhaseSpace(self, overlayExact = 0, dim=0):
        self.figures += 1
        pb.figure( self.figures )
        pb.title("Phase space, simulation (red) vs. exact (blue)")
        pb.xlabel("theta / pi"   )
        pb.ylabel("thetadot / pi")
        pb.grid(True)
        pb.axis('equal')
        x = np.reshape(self.posArray, [self.ndim, self.N, self.samplesteps])
        p = np.reshape(self.velArray, [self.ndim, self.N, self.samplesteps])
        pb.plot(x[dim][0][:] % 1,   p[dim][0][:],
                    linestyle = 'None', marker = 'o', color = 'r')
        # if overlayExact:
        #    exX = np.reshape(self.exPosArray, [self.ndim, self.N, self.samplesteps])
        #    exP = np.reshape(self.exVelArray, [self.ndim, self.N, self.samplesteps])
        #    pb.plot(exX[dim][0][:] / self.xmax, exP[dim][0][:] / self.pmax,
        #            linestyle = '-', color = 'b')
        # pb.savefig('phaseSpace.png')


    # plot temperature as a function of time
    def plotTemp(self):
        self.figures += 1
        pb.figure( self.figures )
        # FIXME
        tempMax = 2 * self.kinEnArray.max() / (self.ndim * self.N)
        pb.plot(self.sampleTArray / self.T , self.tempArray / tempMax,
                linestyle = '-', marker = 'o', color = 'r')
        pb.axhline(y = self.meanTemp() / tempMax,
                linestyle = '--',  color='b')
        pb.title("Instanteneous temperature (red) and Mean Temperature (blue)")
        pb.xlabel("t/T_0")
        pb.ylabel("Inst Temp / Max Temp")
        pb.grid(True)

    # plot potential energy as a function fo time
    def plotPotEn(self, overlayExact = 0):
        self.figures += 1
        pb.figure( self.figures )
        pb.title("Potential energy, simulation (red) vs. exact (blue)")
        # FIXME
        pb.xlabel("t/T_0")
        pb.ylabel("PE/iniEn")
        pb.grid(True)
        pb.plot(self.sampleTArray / self.T, self.potEnArray / self.iniEn,
                linestyle = 'None', marker = 'o', color = 'r')
        # if overlayExact:
        #    pb.plot(self.sampleTArray / self.T, self.exPotEnArray.transpose() / self.iniEn,
        #            linestyle = '-', color = 'b')
        # TODO fix to include absolute / customizable directories
        # pb.savefig('potEn.png')
