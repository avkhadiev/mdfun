from md.simulate import parse
import numpy as np

if __name__ == '__main__':
    import sys, os
    assert len(sys.argv) > 3
    script        = sys.argv[1] # absolute path to the script that runs a sim
    sim           = sys.argv[2] # argument to the script --sim
    config_string = sys.argv[3] # argument to the script --config
    config        = parse( config_string )      # parse *current* config_string
    outdir        = config['outdir']
    os.chdir(outdir)                            # go to output directory

    # build the command string with configuration parameters specified
    # python script (command)
    cmd_string    = "python " + script
    # simulation type
    sim_string    = " --sim " + sim

    # configuration string
    batchmode       = ':batchmode=True'
    writeout        = ':writeout=True'
    calculateExact  = ':calculateExact=False'
    iniTemp         = ':iniTemp=0.01'
    dt              = ':dt=0.001'
    runtime         = ':runtime=3'
    relaxtime       = ':relaxtime=0'            # important! start from x0, v0
    sampleint       = ':sampleint=10'           # sample every 10 time steps
    setInitial      = ':setInitial=True'
    # make additions to the configuration string
    config_string   = (" --config "     +
                        config_string   +
                        batchmode       +
                        writeout        +
                        calculateExact  +
                        iniTemp         +
                        dt              +
                        runtime         +
                        relaxtime       +
                        sampleint       +
                        setInitial)

    # load arrays with initial points and their weights
    x0 = np.loadtxt("%s_%s.txt" % (sim, "x0"))
    p0 = np.loadtxt("%s_%s.txt" % (sim, "v0"))
    fr = np.loadtxt("%s_%s.txt" % (sim, "fr"))

    # TESTING JOB LAUNCH
    # launch jobs
    for i in range(p0.size):
        new_name    = "sho_pt%d_fr%d" % (i, fr[i])
        x0_str      = ":x0=" + str(x0[i])
        p0_str      = ":p0=" + str(p0[i])
        name_string = " -n " + new_name
        full_cmd    = (cmd_string       +
                        name_string     +
                        sim_string      +
                        config_string   +
                        x0_str          +
                        p0_str)
        print full_cmd
