# -*- coding: utf-8 -*-

from mdsimulation import NaiveSimulation
# import numpy as np

class RotRod(NaiveSimulation):
    """A class that describes a (freely) rotating rod (around the CM)"""

    # creates an instance of a simulation class
    # sets up the variables for the exact solution
    def __init__(self, config):
        pass

    # assigns a random angle on a circle
    def randomPos(self):
        pass

    # assigns initial velocity and position
    # saves initialized velocities and positions for the exact solution
    def initialize(self):
        pass

    # redefine kinetic energy for the rotating rod
    def kinEn(self):
        pass

    # no potential energy for a freely rotating rod
    def potEn(self):
        return 0

    # calculates and stores required quantities
    def recordObservables(self):
        pass

    # resets required quantities
    def resetObservables(self):
        pass

    ############## STATISTICS ################

    # calculates exact solutions based on initial conditions
    def calculateExact(self):
        pass

    # for testing purposes, a standalone method to return mean kinetic energy
    # (it's actually constant)
    def meanKinEn(self):
        pass

    ############ PLOTTING ##################

    def printResults(self):
        pass

    def writeObservables(self):
        pass
