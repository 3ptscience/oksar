from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
import vectormath as vmath
import properties
import unittest
import oksar


class TestOksar(unittest.TestCase):

    def test_los(self):

        interferogram_refx = 741140
        interferogram_refy = 4230327
        interferogram_ref_incidence = 23
        local_earth_radius = 6386232
        satellite_altitude = 788792
        satellite_azimuth = 192
        locationUTMzone = 35
        locationE = 706216.0606
        locationN = 4269238.9999

        utmLoc = vmath.Vector3([locationE], [locationN], [0])
        refPoint = vmath.Vector3(interferogram_refx, interferogram_refy, 0)

        LOS = oksar.getLOSvector(
            utmLoc,
            locationUTMzone,
            refPoint,
            satellite_altitude,
            satellite_azimuth,
            interferogram_ref_incidence,
            local_earth_radius
        )

        true = vmath.Vector3([0.427051, -0.090772, 0.899660])
        assert (LOS - true).length < 1e-5


if __name__ == '__main__':
    unittest.main()

