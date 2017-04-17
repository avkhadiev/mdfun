# -*- coding: utf-8 -*-

import numpy as np
import pylab as pb
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
        # temperature is specified as T / Tmax, where Tmax = 2 * 1/4 omegasq (L/2)^2
        # (m = kB = 1)
        if 'iniTemp'    not in config.keys(): config['iniTemp'] = 0.10
        if 'npoints'    not in config.keys(): config['npoints'] = 10000
        # stepX is put in relative to L, stepP relative to L/T
        # this will multiplied by reference values automatically
        if 'stepP'      not in config.keys(): config['stepP']   = 2.0
        if 'stepX'      not in config.keys(): config['stepX']   = 2.0
        if 'adjustStep' not in config.keys(): config['adjustStep'] = 'False'
        if 'nadjst'     not in config.keys(): config['nadjst']  = 20
        if 'batchmode'  not in config.keys(): config['batchmode'] = 'False'
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
        self.pRef       = self.omega * (self.L / 2) # reference momentum
        # SIMULATION
        self.npoints    = int(config['npoints'])    # sample how many?
        self.stepP      = float(config['stepP']) * self.pRef
        self.stepX      = float(config['stepX']) * (self.L / 2)
        self.figures    = 0                         # number of figures plotted
        self.naccpt     = 0                         # for step adjustment
        self.nadjst     = float(config['nadjst'])   # for step adjustment
        # DATA
        self.point      = (0., 0.)                  # current (x, p)
        self.points     = np.array([])              # array of sampled points
        self.freq       = np.array([])              # ith point sampled # times
        self.trials     = 0                         # steps attempted
        self.newAcc     = 0                         # unique points accepted
        self.totAcc     = 0                         # total points accepted
        self.ratios     = np.array([])              # ratios accepts/trials
        self.stepsX     = np.array([])              # for testing purposes
        self.stepsP     = np.array([])              # for testing purposes
        self.iniEnergy  = 0.                        # initial energy
        self.enArray    = np.array([])              # monitoring total energy
        self.samplesteps= np.array([])              # steps w sampled ratios

    # uses an RNG to sample the first pair (p, x)
    # allowed by the input temperature PE_max / kT, omegasq, and size of box L
    def sampleFirst(self):
        # random position
        corr = math.sqrt( self.iniTemp )
        x    = (np.random.random() - 0.5) * self.L * corr
        # assign momentum compatible with temperature and position
        totEn = 0.5 * self.omegasq * (self.L / 2)**2 * self.iniTemp
        potEn = 0.5 * self.omegasq * x**2
        p  = (1 - 2 * (np.random.random() < 0.5)) * math.sqrt( 2*(totEn - potEn) )
        # update points
        self.point     = (x, p)
        self.iniEnergy = self.energy()
        self.points    = np.append( self.points, self.point )
        self.freq      = np.append( self.freq, 1)      # x0, p0 occures 1 times
        self.newAcc   += 1                             # +1 unique accepts
        self.totAcc   += 1                             # +1 total accepts
        self.naccpt   += 1
        self.trials   += 1

    # clears all sampled points and counters
    def reset(self):
        # use to find an optimal step size
        # only current point is left over
        self.npoints     = 1
        self.trials      = 1
        self.newAcc      = 1
        self.totAcc      = 1
        self.naccpt      = 1
        self.freq        = np.array( [1] )  # 1 points remains sampled
        self.ratios      = np.array([])     # ignore ratios before calibration
        self.stepsX      = np.array([])
        self.stepsP      = np.array([])
        self.samplesteps = np.array([])
        # only (x, p) for the last sampled point remain in the array
        self.points      = np.array([ self.points[-2], self.points[-1] ])

    # calculates total energy
    def energy(self):
        (x, p) = self.point
        return 0.5 * ( p * p ) + 0.5 * self.omegasq * ( x * x )

    # calculates mean energy from enArray
    def meanEnergy(self):
        return self.enArray.mean()

    # calculates current ratio of unique accepts to number of attempts
    def calculateRatio(self):
        return float(self.newAcc) / float(self.trials)

    # adjust the ratio and does the bookkeeping
    def adjustRatio(self):
        ratio = float(self.naccpt) / float(self.nadjst)
        if ratio < 0.5:
        # acceptance too low, decrease the step
            correction = 0.95
            self.stepX = self.stepX * correction
            self.stepP = self.stepP * correction
        elif ratio > 0.5:
        # acceptance too high, inscrease the step
            correction = 1.05
            self.stepX = self.stepX * correction
            self.stepP = self.stepP * correction
        # zero the counter
        self.naccpt = 0

    # adds ratio and steps to the lists for plotting during the test
    def addRatio(self):
        ratio = self.calculateRatio()
        self.samplesteps = np.append( self.samplesteps, self.trials )
        self.ratios      = np.append(self.ratios, ratio)

    # returns an array of npoints
    def sampledPoints(self):
        return self.points

    # makes a random step in phase space
    # updates respective counters
    def makeStep(self):
        self.trials += 1
        # generate a vector in (x, p) space
        (xOld, pOld) = self.point
        # make a step in phase space
        xNew = xOld + (2 * np.random.random() - 1) * self.stepX
        pNew = pOld + (2 * np.random.random() - 1) * self.stepP
        return (xNew, pNew)

    # compute weight exp( - delta V)
    def computePba(self, xNew, pNew):
        (xOld, pOld) = self.point
        oldPE        = 0.5 * self.omegasq * (xOld * xOld)
        oldKE        = 0.5 * (pOld * pOld)
        newPE        = 0.5 * self.omegasq * (xNew * xNew)
        newKE        = 0.5 * (pNew * pNew)
        delPE        = newPE - oldPE
        delKE        = newKE - oldKE
        weight       = delPE + delKE
        if weight <= 0:
            pba = 1
        else:
            pba = math.exp( - weight / self.kbT )
        return pba

    # accepts a new trial point and updates correspoding counters / lists
    def accept(self, xTry, pTry):
        self.newAcc += 1            # for tracking unique accepts (not zeroed)
        self.totAcc += 1            # for tracking total accepts (not zeroed)
        self.naccpt += 1            # for adjusting steps (zeroed periodically)
        self.point   = (xTry, pTry)
        self.points = np.append( self.points,  self.point )
        self.freq   = np.append( self.freq,    1          )

    # rejects a new trial point and updates corresponding counters / lists
    def reject(self):
        # append previous point to the list again
        self.totAcc   += 1
        self.freq[-1] += 1

    # makes nsteps new Monte-Carlo steps, saves sampled points in points array
    def evolve(self, nsteps):
        # make nsteps
        for istep in range(nsteps):
            (xTry, pTry) = self.makeStep()
            # ACCEPT or REJECT
            pba = self.computePba(xTry, pTry)
            if pba == 1:
                self.accept(xTry, pTry)
            elif pba < 1:
                if np.random.random() <= pba:
                    self.accept(xTry, pTry)
                else:
                    self.reject()
            else:
                print "Error, probability > 1"
            # adjust step if neccessary
            if (self.trials % self.nadjst == 0):
                if (self.config['adjustStep'] == 'True'):
                    self.adjustRatio()
                self.addRatio()
                self.stepsX      = np.append(self.stepsX, self.stepX)
                self.stepsP      = np.append(self.stepsP, self.stepP)
                self.enArray     = np.append( self.enArray, self.energy() )

    # calculates the normalization factor based on the current sampled points
    def normalization(self):
        # convert to np.arrays
        self.points = np.array( self.points )       # (x, p)
        self.freq   = np.array( self.freq )
        # extract velocities
        v0 = self.points[1::2]                      # only odd indices
        # integrate
        C = (self.freq * v0 ** 2).sum() / self.totAcc
        return C

    # plots the initial (x0, p0) distributions
    def plotPhaseSpace(self):
        self.figures += 1
        pb.figure( self.figures )
        pb.title("(x0, p0) distribution")
        pb.xlabel("x0 / (L/2)"  )
        pb.ylabel("p0 / w (L/2)")
        pb.grid(True)
        #pb.axis('equal')
        x      = self.points[::2]
        p      = self.points[1::2]
        corr   = math.sqrt( self.iniTemp )
        margin = 0.25
        xmax   = abs( x.max() ) / (self.L / 2)
        xmin   = abs( x.min() ) / (self.L / 2)
        xbound = max( xmax, xmin, 1.0 * corr)
        pmax   = abs( p.max() ) / (self.pRef)
        pmin   = abs( p.min() ) / (self.pRef)
        pbound = max( pmax, pmin, 1.0 * corr )
        pb.xlim(-xbound - margin, xbound + margin)
        pb.ylim(-pbound - margin, pbound + margin)
        if self.config['debug'] == 'True':
            print "pbound ", pbound
            print "xbound ", xbound
        # PHASE SPACE PLOT
        pb.plot(x / (self.L / 2),   p / self.pRef,
            linestyle = 'None', marker = 'o', color = 'r')
        # Maximum possible energy given the box size
        pb.axhline(1.0,
                linestyle = '--', color = 'k')
        pb.axhline(-1.0,
                linestyle = '--', color = 'k')
        pb.axvline(1.0,
                linestyle = '--', color = 'k')
        pb.axvline(-1.0,
                linestyle = '--', color = 'k')
        # Total energy given the initial temperature
        pb.axhline(1.0 * corr,
                linestyle = '--', color = 'b')
        pb.axhline(-1.0 * corr,
                linestyle = '--', color = 'b')
        pb.axvline(1.0 * corr,
                linestyle = '--', color = 'b')
        pb.axvline(-1.0 * corr,
                linestyle = '--', color = 'b')
        np.savetxt("%s_%s.txt" % (self.name, "x0"), x)
        np.savetxt("%s_%s.txt" % (self.name, "v0"), p)
        # Trajectories
        self.figures += 1
        f, axarr = pb.subplots(2, sharex=True)
        axarr[0].set_title("x0, p0  distribution")
        axarr[0].set_ylabel("x0 / L")
        axarr[0].plot( x / self.L ,
            linestyle = '-', color = 'r')
        axarr[1].set_xlabel("sample trial")
        axarr[1].set_ylabel("v0 / (L/T)")
        axarr[1].plot( p / self.pRef,
            linestyle = '-', color = 'b')

    # plots total energy monitored throughout the simulation
    def plotEnergy(self):
        self.figures += 1
        pb.figure( self.figures )
        pb.title("MC mean energy (blue) and kbT (black)")
        pb.ylabel("E/E_0")
        pb.xlabel("sample trial")
        pb.ylim(0, 2.0)
        pb.grid(True)
        pb.axhline(self.meanEnergy() / self.iniEnergy,
                linestyle = '--', color = 'r')
        pb.axhline(self.kbT / self.iniEnergy,
                linestyle = '--', color = 'k')
        if self.config['debug'] == 'True':
            print "MC mean energy ", self.meanEnergy()
            print "kT ",            self.kbT           # ndof * kbT

    # plots the steps that are adjusted during testing
    def plotSteps(self):
        self.figures += 1
        f, axarr = pb.subplots(2, sharex=True)
        axarr[0].set_title("Steps in x (red), p (blue) during testing")
        axarr[0].set_ylabel("stepX/ (L/2)")
        axarr[0].plot( self.samplesteps / self.trials, self.stepsX / (self.L / 2),
                linestyle = '-', marker = 'o', markersize = 3, color = 'r')
        axarr[0].axhline(self.stepsX.mean() / (self.L / 2),
                linestyle = '--', color = 'k')
        axarr[0].grid(True)
        axarr[1].set_xlabel("sample trial")
        axarr[1].set_ylabel("stepP / w (L/2)")
        axarr[1].plot( self.samplesteps / self.trials , self.stepsP / self.pRef,
            linestyle = '-', marker = 'd', markersize = 3, color = 'b')
        axarr[1].axhline(self.stepsP.mean() / self.pRef,
                linestyle = '--', color = 'k')
        axarr[1].grid(True)

    # plots the ratios adjusted during testing
    def plotRatios(self):
        self.figures += 1
        pb.figure( self.figures )
        pb.title(" unique accepts / trials")
        pb.xlabel("trial #")
        pb.ylabel("ratio")
        pb.grid(True)
        pb.plot( self.samplesteps / self.trials,  self.ratios,
            linestyle = '-',  color = 'r')

    def showPlots(self):
        pb.show()

