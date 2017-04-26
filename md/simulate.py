# -*- coding: utf-8 -*-

from simulations import sho, mcsho, rotating
import numpy as np
import pylab as pb

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
        if MDRun( simulation ):
            # if energy too high
            print "specified initial coordinates are not allowed: maximum energy exceeded"
    elif sim == "barbell":
        simulation = rotating.Barbell( name, config )
        MDRun( simulation )
    elif sim == "mcsho_test_basic":
        simulation = mcsho.MCMetropolis( name, config )
        MCBasicTest( simulation )
    elif sim == "mcsho_test_norm":
        MCNormTest()
    elif sim == "mcsho_train":
        simulation = mcsho.MCMetropolis( name, config )
        (stepX, stepP) = MCTrain( simulation )
        print "stepX: ", stepX
        print "stepP: ", stepP
    else:
        print "%s is an unknown simulation" % sim

# test the MC sho
# graph sampled x0, p0 points
# graph the step size changes
# calculate the normalization and return it
def MCBasicTest( simulation ):
    npoints = int(simulation.config['npoints'])
    (stepX, stepP) = MCTrain( simulation )
    print "Training done:"                      # training to get 50% acc ratio
    print "stepX: ", stepX
    print "stepP: ", stepP
    simulation.reset()
    simulation.config['adjustStep'] = 'False'   # no step adjustment
    simulation.stepX = stepX                    # assign steps from training
    simulation.stepP = stepP
    simulation.sampleFirst()
    simulation.evolve(1000)
    simulation.reset()
    while (simulation.totAcc < npoints):
        simulation.evolve(1)
    if simulation.config['batchmode'] == 'False':
        print "MC SHO test: normalization C_0 is ", simulation.normalization()
        print "Points sampled: ", simulation.totAcc
        print "Frequency sum: ",  simulation.freq.sum()
        simulation.plotPhaseSpace()
        simulation.plotRatios()
        simulation.plotSteps()
        simulation.plotEnergy()
        simulation.showPlots()
    if simulation.config['writeout'] == 'True':
        simulation.writeout()
    return simulation.normalization()

# prints normalization as a function of input temperature
def MCNormTest():
    iniTempArr = np.linspace(0.1, 1, 10)
    norms = np.array([])
    for iniTemp in iniTempArr:
        name    = "mcshoNormTest" + "_iniTemp_" + str(iniTemp)
        config  = parse("stepP=0.1:stepX=0.1:npoints=100000:batchmode=True:debug=False")
        config['iniTemp'] = iniTemp
        simulation = mcsho.MCMetropolis( name, config )
        norm = MCBasicTest( simulation )
        norms = np.append( norms, norm )
    tmax = simulation.tempMax
    pb.figure(1)
    pb.title("C_0 as a function of T, simulation (red) vs. kbT (black)")
    pb.xlabel("T/TempMax")
    pb.ylabel("C_0")
    pb.grid(True)
    pb.plot( iniTempArr, norms, linestyle = 'None', marker = 'o', color = 'r' )
    pb.plot( iniTempArr, iniTempArr * tmax, linestyle = '-', color = 'k' )
    pb.show()

# test the MC sho
# graph sampled x0, p0 points
# graph the step size changes
# calculate the normalization and print it
def MCTrain( simulation ):
    npoints = int(simulation.config['npoints'])
    simulation.config['adjustStep'] = 'True'
    simulation.sampleFirst()
    # relax the system first
    simulation.evolve(100)
    simulation.reset()
    while (simulation.totAcc < npoints):
        simulation.evolve(1)
    # return steps that are good to obtain a ratio of ~ 0.5
    stepX = simulation.stepsX.mean() / (simulation.L / 2)
    stepP = simulation.stepsP.mean() / (simulation.pRef )
    return (stepX, stepP)

# runs an MD simulation
# according to the config
def MDRun( simulation, runtime = -1, relaxtime = -1 ):
    if runtime   <  0.0: runtime   = simulation.runtime
    if relaxtime <  0.0: relaxtime = simulation.relaxtime
    if simulation.config['setInitial'] == 'True':
        x0 = float(simulation.config['x0'])
        p0 = float(simulation.config['p0'])
        if simulation.setInitial( x0, p0 ):
            # if energy is too high
            return 1
    else:
        simulation.initialize()
        simulation.evolve( relaxtime )
        simulation.resetStatistics()
    simulation.initializeExact()
    simulation.evolve( runtime )                        # in units of T
    if simulation.config['calculateExact'] == 'True':
        simulation.calculateExact()
    if simulation.config['writeout'] == 'True':
        simulation.writeObservables()
    if simulation.config['batchmode'] == 'False':
        simulation.printResults()
        if simulation.config['calculateExact'] == 'True':
            overlayExact = 1
        else:
            overlayExact = 0
        simulation.plotPotEn( overlayExact )
        simulation.plotTrajectories( overlayExact )
        simulation.plotPhaseSpace( overlayExact )
        simulation.showPlots()
    return 0

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
