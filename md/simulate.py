# -*- coding: utf-8 -*-

from . import daemon

if __name__ == '__main__':
    import sys
    assert len(sys.argv) > 2
    sim           = sys.argv[1];
    config_string = sys.argv[2];
    config        = daemon.parse( config_string )
    configure     = daemon.Configure( config )
    configure.simulation( sim )
