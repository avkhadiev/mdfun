# -*- coding: utf-8 -*-

from simulations import sho, mcsho
import numpy as np

def parse(config_string):
    config = dict(item.split("=") for item in config_string.split(":"))
    return config

def run( sim, name, config, verbose, debug ):
    if (verbose):
        print "Using the following configuration:\n"
        print( config )
        print "\n"
    config['verbose'] = str(verbose)
    config['debug']   = str(debug)
    if sim == "sho":
        simulation = sho.SHO(name, config)
        simulation.run()
    if sim == "mcsho":
        simulation = mcsho.MCMetropolis(name)
        MCNormTest( simulation )
    else:
        print "%s is an unknown simulation" % sim

# test the MC sho
# graph sampled x0, p0 points
# graph the step size changes
# calculate the normalization and print it
def MCNormTest( simulation ):
    simulation.sampleFirst()
    while (simulation.newAcc < 100) or (simulation.trials < 100000):
        simulation.evolve(30)
        simulation.adjustRatio()
    print "MC SHO test: normalization C_0 is ", simulation.normalization
    np.savetxt('%s_%s.txt' % (simulation.name, "x0v0"), simulation.points)


def time_reversal( sim, config, verbose, debug ):
    if (verbose):
        print "Using the following configuration:\n"
        print( config )
        print "\n"
    config['verbose'] = str(verbose)
    config['debug']   = str(debug)
    if 'reverse_time' not in config.keys():
        reverse_time = 5.0
    else:
        reverse_time = float(config['reverse_time'])
    if (verbose):
        print "Using the following configuration:\n"
        print( config )
        print "\n"
    if sim == "sho":
        simulation = sho.SHO(config)
        simulation.initialize()
        simulation.evolve( 10.0 )               # relax_time
        simulation.resetStatistics()
        simulation.evolve( reverse_time )
        simulation.plotTrajectories()
        simulation.reverseTime()
        simulation.evolve( reverse_time )
        simulation.plotTrajectories()
        simulation.showPlots()
    else:
        print "%s is an unknown simulation" % sim

if __name__ == '__main__':
    import sys
    assert len(sys.argv) > 2
    sim           = sys.argv[1];
    config_string = sys.argv[2];
    config        = parse( config_string )
    config['verbose'] = 'True'
    config['debug']   = 'False'
    run( sim, config )
