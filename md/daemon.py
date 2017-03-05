# -*- coding: utf-8 -*-

import simulations
import universe

def parse(config_string):
    config = dict(item.split("=") for item in config_string.split(":"))
    return config

class Configurer(object):

    """An configurer of various of parts of the simulation"""

    def __init__(self, config):
        self.config = config

    def simulation(self, name):
        config = self.config
        if name == "SHO":
            assert('timestep' in config.keys())    # ensure key components
            assert('duration' in config.keys())    # are in the config
            return simulations.sho.SHO( self )
        else:
            return "%s is not a valid simulation" % name

    def universe(self, name):
        if name == "line":
            return universe.box( name, self )
        else:
            return "%s is not a valid universe" % name

