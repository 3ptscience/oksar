"""
    oksar3

   Program to calcuate forward models of interferograms, strain tensor, etc.
   from Okada subroutine.

   Heritage: originally fringes.c written by Barry Parsons
             updated to oksar                                 tjw
             oksar_strain: added strain tensor calculation    tjw
             oksar3:       added new line of sight calculator tjw feb 2003
             Modified into Python by RowanCockett, 3point Science Aug 2014
"""

import numpy as np
import vectormath as vmath
from Los import getLOSvector


class Forward(object):

    def __init__(self, model):
        self.model = model

    def getLos(self, eq, utmLoc):
        refPoint = vmath.Vector3(eq.interferogram_refx, eq.interferogram_refy, 0)
        return getLOSvector(
            utmLoc,
            eq.locationUTMzone,
            refPoint,
            eq.satellite_altitude,
            eq.satellite_azimuth,
            eq.interferogram_ref_incidence,
            eq.local_earth_radius
        )

    def getDir(self, x, y):
        model = self.model
        DEG2RAD = 0.017453292519943
        alpha = (model.beta + model.mu) / (model.beta + 2.0 * model.mu)

        #  Here we could loop over models

        flt_x = model.center.x
        flt_y = model.center.y
        strike = model.strike
        dip = model.dip
        rake = model.rake
        slip = model.slip
        length = model.length
        hmin = model.depthT
        hmax = model.depthB

        rrake = (rake+90.0)*DEG2RAD
        sindip = np.sin(dip*DEG2RAD)
        w = (hmax-hmin)/sindip
        ud = slip*np.cos(rrake)
        us = -slip*np.sin(rrake)
        halflen = length/2.0
        al2 = halflen
        al1 = -al2
        aw1 = hmin/sindip
        aw2 = hmax/sindip

        if(hmin < 0.0):
            raise Exception('ERROR: Fault top above ground surface')

        if(hmin == 0.0):
            hmin = 0.00001
        sstrike = (strike+90.0)*DEG2RAD
        ct = np.cos(sstrike)
        st = np.sin(sstrike)

        X = ct*(-flt_x+x)-st*(-flt_y+y)
        Y = ct*(-flt_y+y)+st*(-flt_x+x)
        u = self.dc3d3(alpha, X, Y, -dip, al1, al2, aw1, aw2, us, ud)

        UX = ct*u.x + st*u.y
        UY = -st*u.x + ct*u.y
        UZ = u.z

        return vmath.Vector3(UX, UY, UZ)

    def dc3d3(self, alpha, X, Y, dip, al1, al2, aw1, aw2, disl1, disl2):
        F0 = 0.0
        F1 = 1.0
        F2 = 2.0
        PI2 = 6.283185307179586
        EPS = 1.0E-6

        u = vmath.Vector3(F0, F0, F0)
        dub = vmath.Vector3(F0, F0, F0)

        #  %%dccon0 subroutine
        #  Calculates medium and fault dip constants
        c0_alp3 = (F1-alpha)/alpha
        #  PI2/360
        pl8 = 0.017453292519943
        c0_sd = np.sin(dip*pl8)
        c0_cd = np.cos(dip*pl8)

        if(np.abs(c0_cd) < EPS):
            c0_cd = F0
            if(c0_sd > F0):
                c0_sd = F1

            if(c0_sd < F0):
                c0_sd = -F1

        c0_cdcd = c0_cd * c0_cd
        c0_sdcd = c0_sd * c0_cd

        #  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        p = Y * c0_cd
        q = Y * c0_sd

        jxi = ((X - al1) * (X - al2)) <= F0  # BOOLEAN
        jet = ((p - aw1) * (p - aw2)) <= F0  # BOOLEAN

        for k in [1., 2.]:
            et = 0.0
            if(k == 1):
                et = p-aw1
            else:
                et = p-aw2

            for j in [1., 2.]:
                xi = 0.0
                if(j == 1):
                    xi = X-al1
                else:
                    xi = X-al2

    #          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #          %%dccon2 subroutine
    #          % calculates station geometry constants for finite source

                dc_max = np.max(np.abs(np.c_[xi, et, q]))

                # dc_max = max(np.abs(xi),max(np.abs(et),np.abs(q)))

                xi[(np.abs(xi/dc_max) < EPS) | (np.abs(xi) < EPS)] = F0

                et[(np.abs(et/dc_max) < EPS) | (np.abs(et) < EPS)] = F0

                q[(np.abs(q/dc_max) < EPS) | (np.abs(q) < EPS)] = F0

                dc_xi = xi
                dc_et = et
                dc_q = q
                c2_r = np.sqrt(dc_xi*dc_xi + dc_et*dc_et + dc_q*dc_q)

                if np.any(c2_r == F0):
                    raise Exception('singularity error ???')

                c2_y = dc_et * c0_cd + dc_q * c0_sd
                c2_d = dc_et * c0_sd - dc_q * c0_cd
                c2_tt = np.arctan(dc_xi * dc_et / (dc_q * c2_r))
                c2_tt[dc_q == F0] = F0

                rxi = c2_r + dc_xi
                c2_x11 = F1/(c2_r*rxi)
                c2_x11[(dc_xi < F0) & (dc_q == F0) & (dc_et == F0)] = F0

                ret = c2_r + dc_et
                if np.any(ret < 1e-14):
                    raise Exception('dccon2 b %f %f %f %f' % (
                        ret, c2_r, dc_et, dc_q, dc_xi
                    ))

                c2_ale = np.log(ret)
                c2_y11 = F1/(c2_r*ret)

                ind = (dc_et < F0) & (dc_q == F0) & (dc_xi == F0)

                # if((c2_r-dc_et) < 1e-14):
                #     raise Exception('dccon2 a %f %f %f %f %f' % (
                #         c2_3-dc_et, c2_r, dc_et, dc_q, dc_xi)
                #     )

                c2_ale[ind] = -np.log(c2_r[ind]-dc_et[ind])
                c2_y11[ind] = F0

                # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if np.any(
                    (
                        (q == F0) &
                        (
                            ((jxi) & (et == F0)) |
                            ((jet) & (xi == F0))
                        )
                    ) | (c2_r == F0)
                ):
                    raise Exception('singular problems: 2')

                # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                # ub subroutine
                # part B of displacement and strain at depth due to buried
                # faults in semi-infinite medium

                rd = c2_r + c2_d
                if np.any(rd < 1e-14):
                    raise Exception('ub %f %f %f %f %f %f' % (
                        rd, c2_r, c2_d, xi, et, q
                    ))

                ai3 = 0.0
                ai4 = 0.0
                if(c0_cd != F0):
                    # xx replaces x in original subroutine
                    xx = np.sqrt(xi*xi+q*q)
                    ai4 = F1/c0_cdcd * (xi/rd*c0_sdcd + F2*np.arctan(
                        (
                            et*(xx+q*c0_cd) +
                            xx*(c2_r+xx)*c0_sd
                        ) / (xi*(c2_r+xx)*c0_cd)
                    ))
                    ai4[xi == F0] = F0

                    ai3 = (c2_y*c0_cd/rd - c2_ale + c0_sd*np.log(rd)) / c0_cdcd
                else:
                    rd2 = rd*rd
                    ai3 = (et/rd + c2_y*q/rd2 - c2_ale) / F2
                    ai4 = xi*c2_y/rd2/F2

                ai1 = -xi/rd*c0_cd - ai4*c0_sd
                ai2 = np.log(rd) + ai3*c0_sd
                qx = q*c2_x11
                qy = q*c2_y11

                # strike-slip contribution
                if(disl1 != 0.0):
                    du2x = - xi*qy - c2_tt - c0_alp3 * ai1 * c0_sd
                    du2y = - q/c2_r + c0_alp3*c2_y/rd*c0_sd
                    du2z = q*qy - c0_alp3*ai2*c0_sd
                    du2 = vmath.Vector3(du2x, du2y, du2z)
                    dub = du2 * (disl1 / PI2)
                else:
                    dub = vmath.Vector3()

                # dip-slip contribution
                if(disl2 != F0):
                    du2x = - q/c2_r + c0_alp3 * ai3 * c0_sdcd
                    du2y = - et*qx - c2_tt - c0_alp3 * xi / rd * c0_sdcd
                    du2z = q*qx + c0_alp3 * ai4 * c0_sdcd
                    du2 = vmath.Vector3(du2x, du2y, du2z)
                    dub = dub + (du2 * (disl2 / PI2))

                # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                dux = dub.x
                duy = dub.y*c0_cd - dub.z*c0_sd
                duz = dub.y*c0_sd + dub.z*c0_cd
                du = vmath.Vector3(dux, duy, duz)
                if((j+k) != 3):
                    u = + du + u
                else:
                    u = - du + u

        return u


class InSarModel(object):
    def __init__(self):
        pass

    beta = 3E10
    mu = 3E10

    sx = 432
    sy = 232
    strike = 329.6
    dip = 50
    rake = 90
    slip = 0.5
    length = 11578.907244622129
    center = vmath.Vector3(773728.2977967655, 4223586.816611591, 0)
    depthT = 0
    depthB = 15000
    wv = 0.028333
    rigidity = 30000000000

    LOS = vmath.Vector3(0.3825, 0.0780, 0.9205)


if __name__ == '__main__':
    f = Forward(InSarModel)
    n = 300
    O = vmath.Vector3(706216.0606, 4187318.9999, 0)
    U = vmath.Vector3(81920, 0, 0)
    V = vmath.Vector3(0, 81920, 0)

    vec, shape = vmath.ouv2vec(O, U, V, n)

    print vec.shape
    # strike=254
    # dip=40
    # rake=130
    # slip=0.9209931290545376
    # length=16384
    # center=747176.0606,4228278.9999
    # depthT=0.1
    # depthB=15000
    # wv=0.028333333
    # rigidity=30000000000
    # wrap=false
    # climMin=-0.621606826782
    # climMax=0.0770393759012
    # uid=dinar
    DIR = f.getDir(vec.x, vec.y)


