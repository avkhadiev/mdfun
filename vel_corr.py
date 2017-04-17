# -*- coding: utf-8 -*-

from md.simulate import parse
from simulations import sho

# given a phase space point (x0, p0), a sample interval in units of period and
# a runtime in units of T,
# returns a list of points [x(t1), p(t1), ..., x(tn), p(tn)]
# where tn = runtime and t_(i+1) - t_(i) = sampleint
def launchMD(name, x0, p0, iniTemp,  sampleint=0.1, runtime=5):
    iniTemp_str   = "iniTemp="   + str(iniTemp)
    sampleint_str = "sampleint=" + str(sampleint)
    runtime_str   = "runtime="   + str(runtime)
    stdout_str    = "setInitial=True:batchmode=True:debug=False:calculateExact=False"
    config_str    = runtime_str + ":" + sampleint_str + ":" + stdout_str
    config        = parse( config_str )
    simulation = sho.SHO( name, config )
    if simulation.setInitial( x0, p0 ):
        # if energy is too high
        return 1



