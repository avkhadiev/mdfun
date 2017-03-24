# -*- coding: utf-8 -*-

#
# 1D Simple Harmonic Oscillator MD simulation
# launches a single simulation based on the provided configuration
# TODO calculates and saves required observables in the file
#

from optparse import OptionParser
parser = OptionParser()
#parser.add_option("-n", "--name",
#        action="store", type="string", dest="sim",
#        help="specify simulation type")
parser.add_option("-s", "--simulation",
        action="store", type="string", dest="sim",
        help="specify simulation type")
parser.add_option("-c", "--config",
        action="store", type="string", dest="config",
        help="specify configuration via ':'-separated string")
parser.add_option("-d", "--debug",
        action="store_true", dest="debug",
        help="turn on debugging mode")
parser.add_option("-v",
        action="store_true", dest="verbose",
        help ="print status messages to standartd output")
parser.add_option("-q",
        action="store_false", dest="verbose",
        help ="don't print status messages to standartd output")
parser.add_option("-t", "--test",
        action="store_true", dest="test",
        help ="run tests")
parser.set_defaults(name="newTrial",
            simulation="sho",
            config  = ":".join([
                            "dt=0.001",
                            "integration=velocityVerlet"]),
            debug   = False,
            verbose = True,
            test    = False)

(options, args) = parser.parse_args()

###############################################################################
from md import simulate

name       = options.name
sim        = options.simulation
verbose    = options.verbose
debug      = options.debug
test       = options.test

if (verbose): print "Verbose mode: on\n"
if (debug):   print "Debug mode: on\n"
if (test):    print "Test mode: on\n"

config = simulate.parse(options.config)

if (test):
    simulate.time_reversal( sim, config, verbose, debug )
else:
    simulate.run( sim, config, verbose, debug )

