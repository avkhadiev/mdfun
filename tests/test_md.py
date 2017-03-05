# -*- coding: utf-8 -*-

from .context import md
import unittest
import numpy as np

class BasicMDTestSuite(unittest.TestCase):

    """Basic test cases for md package"""

    def test_simulate_setup(self):
        config     = ""
        # self.assertEqual(md.simulate.setup("SHO", config),
        #        "Run SHO simulation")
        self.assertEqual(md.simulate.setup("wrong", config),
                "wrong is not a valid simulation")

    # creates a 1d box of length size for testing purposes
    def setup_1dbox(self, size):
        size      = np.array([size]);
        return md.universe.box.Box( size )

    # given a box creates a particle in the box
    def setup_particle(self, box, pos, vel, memory=3):
        assert( pos.size == vel.size )
        assert( np.all( np.greater_equal( box, pos )) )
        particle = md.particle.point.PointParticle( box, pos, vel, memory )
        return particle

    def test_1dbox(self):
        ### SETUP ###
        line      = self.setup_1dbox( 100. )
        pos_in    = np.array([ 50.  ])
        pos_out_1 = np.array([ 101. ])
        pos_out_2 = np.array([ -1.  ])
        ### TEST ###
        self.assertTrue(line.is_inside( pos_in ))
        self.assertFalse(line.is_inside( pos_out_1 ))
        self.assertFalse(line.is_inside( pos_out_2 ))

    def test_particle(self):
        ### SETUP ###
        line     = self.setup_1dbox( 100. )                 # 1d line size 100
        pos      = np.array([ 50. ])
        vel      = np.array([ 3.  ])
        particle = self.setup_particle( line, pos, vel )
        ### TEST ###
        self.assertEqual( pos, particle.pos[0][:] )
        self.assertEqual( vel, particle.vel[0][:] )

if __name__ == '__main__':
    unittest.main()
