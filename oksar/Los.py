from __future__ import division
import utm
import numpy as np
from vectormath import Vector3


def getLOSvector(
        utmLoc,
        utmZone,
        refPoint,
        satAltitude,
        satAzimuth,
        satIncidence,
        earthRadius):
    """
        calculate beta - the angle at earth center between reference point and
        satellite nadir
    """

    DEG2RAD = np.pi / 180.
    alpha = satIncidence * DEG2RAD
    beta = (earthRadius / (satAltitude + earthRadius)) * np.sin(alpha)
    beta = alpha - np.arcsin(beta)
    beta = beta / DEG2RAD

    # calculate angular separation of (x,y) from satellite track passing
    # through (origx,origy) with azimuth satAzimuth

    # Long lat **NOT** lat long
    origy, origx = utm.to_latlon(
        refPoint.x, refPoint.y, np.abs(utmZone), northern=utmZone > 0
    )

    xy = np.array([
        utm.to_latlon(u[0], u[1], np.abs(utmZone), northern=utmZone > 0)
        for u in utmLoc
    ])
    y = xy[:, 0]
    x = xy[:, 1]

    angdist = ang_to_gc(x, y, origx, origy, satAzimuth)

    # calculate beta2, the angle at earth center between roaming point and
    # satellite nadir track, assuming right-looking satellite

    beta2 = beta-angdist
    beta2 = beta2*DEG2RAD

    # calculate alpha2, the new incidence angle

    alpha2 = np.sin(beta2) / (
        np.cos(beta2) - (earthRadius / (earthRadius + satAltitude))
    )
    alpha2 = np.arctan(alpha2)
    alpha2 = alpha2 / DEG2RAD

    # calculate pointing vector

    satIncidence = 90 - alpha2
    satAzimuth = 360 - satAzimuth

    los_x = -np.cos(satAzimuth * DEG2RAD) * np.cos(satIncidence * DEG2RAD)
    los_y = -np.sin(satAzimuth * DEG2RAD) * np.cos(satIncidence * DEG2RAD)
    los_z = np.sin(satIncidence * DEG2RAD)

    return Vector3(los_x, los_y, los_z)


def ang_to_gc(x, y, origx, origy, satAzimuth):
    """
        Calculate angular distance to great circle passing through given point
    """

    Ngc = np.zeros(3)
    cartxy = np.zeros((len(x), 3))
    satAzimuth = np.deg2rad(satAzimuth)
    origx = np.deg2rad(origx)
    origy = np.deg2rad(origy)

    x = np.deg2rad(x)
    y = np.deg2rad(y)

    # 1. calculate geocentric norm vec to great circle, Ngc = Rz*Ry*Rx*[0;1;0]
    #    where Rz = rotation of origx about geocentric z-axis
    #    where Ry = rotation of origy about geocentric y-axis
    #    where Rx = rotation of satAzimuth about geocentric x-axis
    #    and [0;1;0] is geocentric norm vec to N-S Great Circle through 0 0

    Ngc[0] = (
                (np.sin(satAzimuth) * np.sin(origy) * np.cos(origx)) -
                (np.cos(satAzimuth) * np.sin(origx))
             )
    Ngc[1] = (
                (np.sin(satAzimuth) * np.sin(origy) * np.sin(origx)) +
                (np.cos(satAzimuth) * np.cos(origx))
             )
    Ngc[2] = -np.sin(satAzimuth) * np.cos(origy)

    # 2. calculate unit vector geocentric coordinates for lon/lat
    #    position (x,y)

    cartxy[:, 0] = np.cos(x) * np.cos(y)
    cartxy[:, 1] = np.sin(x) * np.cos(y)
    cartxy[:, 2] = np.sin(y)

    # 3. Dot product between Ngc and cartxy gives angle 90 degrees
    #    bigger than what we want

    angdist = Ngc[0]*cartxy[:, 0] + Ngc[1]*cartxy[:, 1] + Ngc[2]*cartxy[:, 2]
    angdist = np.rad2deg(np.arccos(angdist)) - 90

    return angdist


if __name__ == '__main__':

    interferogram_refx = 741140
    interferogram_refy = 4230327
    interferogram_ref_incidence = 23
    local_earth_radius = 6386232
    satellite_altitude = 788792
    satellite_azimuth = 192
    locationUTMzone = 35
    locationE = 706216.0606
    locationN = 4269238.9999

    utmLoc = Vector3([locationE], [locationN], [0])
    refPoint = Vector3(interferogram_refx, interferogram_refy, 0)

    LOS = getLOSvector(
        utmLoc,
        locationUTMzone,
        refPoint,
        satellite_altitude,
        satellite_azimuth,
        interferogram_ref_incidence,
        local_earth_radius
    )

    true = Vector3([0.427051, -0.090772, 0.899660])
    print true, LOS
    assert (LOS - true).length < 1e-5
