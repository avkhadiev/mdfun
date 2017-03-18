# -*- coding: utf-8 -*-

from simulations import sho

def parse(config_string):
    config = dict(item.split("=") for item in config_string.split(":"))
    return config


def run( sim, config, verbose, debug ):
    if (verbose):
        print "Using the following configuration:\n"
        print( config )
        print "\n"
    config['verbose'] = str(verbose)
    config['debug']   = str(debug)
    if sim == "sho":
        simulation = sho.SHO(config)
        simulation.run()
    else:
        print "%s is an unknown simulation" % sim

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
