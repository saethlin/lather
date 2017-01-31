#include <iostream>
#include <vector>
#include <math.h>
#include <stdbool.h>
#include "profile.hpp"
#include "star.hpp"
#include "spot.hpp"


const double pi = M_PI;


Spot::Spot(Star star, double latitude, double longitude, double size, bool plage, int spotResolution) {
    latitude *= pi/180.;
    longitude *= pi/180.;
    this -> star = star;
    this -> latitude = latitude;
    this -> longitude = longitude;
    this -> size = size;
    this -> spotResolution = spotResolution;
    this -> plage = plage;

    if (!plage) {
        this -> spotTemp = star.temperature - star.spotTempDiff;
        this -> intensity = planck(5293.4115e-10, spotTemp) / star.intensity;
    }

    //Compute initial location
    std::vector<std::vector<double> > centeredCoordinates(spotResolution, std::vector<double> (3));
    initialCoordinates = std::vector<std::vector<double> > (spotResolution, std::vector<double>(3));


    for (int i = 0; i < initialCoordinates.size(); i++) {
        double rho = -pi + i*2*pi/(spotResolution-1);     // For this case, phase goes from -pi to pi (no idea why)
        centeredCoordinates[i][0] = sqrt(1 - size*size);
        centeredCoordinates[i][1] = size * cos(rho);
        centeredCoordinates[i][2] = size * sin(rho);
    }

    // matrix R in the original code
    double rot_incl = pi/2. - star.inclination;
    double matrixInit[3][3] = {{cos(rot_incl)*cos(longitude)*cos(-latitude) - sin(-latitude)*sin(rot_incl),
                                -sin(longitude)*cos(rot_incl),
                                cos(rot_incl)*cos(longitude)*sin(-latitude) + sin(rot_incl)*cos(-latitude)},
                               {sin(longitude)*cos(-latitude),
                                cos(longitude),
                                sin(longitude)*sin(-latitude)},
                               {-sin(rot_incl)*cos(longitude)*cos(-latitude) - cos(rot_incl)*sin(-latitude),
                                sin(rot_incl)*sin(longitude),
                                -sin(rot_incl)*cos(longitude)*sin(-latitude) + cos(rot_incl)*cos(-latitude)}};

    // matrix R2 in the original code
    double b2 = -(pi/2 - star.inclination);
    matrixSpot[0][0] = cos(latitude) * cos(-longitude) * cos(b2) - sin(b2) * sin(latitude);
    matrixSpot[0][1] = -sin(-longitude) * cos(latitude);
    matrixSpot[0][2] = cos(latitude) * cos(-longitude) * sin(b2) + sin(latitude) * cos(b2);
    matrixSpot[1][0] = sin(-longitude) * cos(b2);
    matrixSpot[1][1] = cos(-longitude);
    matrixSpot[1][2] = sin(-longitude) * sin(b2);
    matrixSpot[2][0] = -sin(latitude) * cos(-longitude) * cos(b2) - cos(latitude) * sin(b2);
    matrixSpot[2][1] = sin(latitude) * sin(-longitude);
    matrixSpot[2][2] = -sin(latitude) * cos(-longitude) * sin(b2) + cos(latitude) * cos(b2);

    // Rotate centeredCoordinates to the correct location on the star's surface
    for (int i = 0; i < initialCoordinates.size(); i++) {
        initialCoordinates[i][0] = matrixInit[0][0]*centeredCoordinates[i][0] + matrixInit[0][1]*centeredCoordinates[i][1] + matrixInit[0][2]*centeredCoordinates[i][2];
        initialCoordinates[i][1] = matrixInit[1][0]*centeredCoordinates[i][0] + matrixInit[1][1]*centeredCoordinates[i][1] + matrixInit[1][2]*centeredCoordinates[i][2];
        initialCoordinates[i][2] = matrixInit[2][0]*centeredCoordinates[i][0] + matrixInit[2][1]*centeredCoordinates[i][1] + matrixInit[2][2]*centeredCoordinates[i][2];
    }


}


bool Spot::isVisible(double phase) {

    double axis[3] = {cos(star.inclination), 0, sin(star.inclination)};
    double rotation[3][3] = {{(1 - cos(phase)) * axis[0] * axis[0] + cos(phase),
                              (1 - cos(phase)) * axis[0] * axis[1] - sin(phase) * axis[2],
                              (1 - cos(phase)) * axis[0] * axis[2] + sin(phase) * axis[1]},
                             {(1 - cos(phase)) * axis[1] * axis[0] + sin(phase) * axis[2],
                              (1 - cos(phase)) * axis[1] * axis[1] + cos(phase),
                              (1 - cos(phase)) * axis[1] * axis[2] - sin(phase) * axis[0]},
                             {(1 - cos(phase)) * axis[2] * axis[0] - sin(phase) * axis[1],
                              (1 - cos(phase)) * axis[2] * axis[1] + sin(phase) * axis[0],
                              (1 - cos(phase)) * axis[2] * axis[2] + cos(phase)}};

    double newx, newy, newz;
    int countOn = 0, countOff = 0;
    double miny = 1;
    double maxy = -1;
    double minz = 1;
    double maxz = -1;

    // Search for bounds by applying the rotation matrix to the initial coordinates
    // x coordinate is depth, so proceed only if the spot is on the front of the star
    for (unsigned int i = 0; i < initialCoordinates.size(); i++) {
        newx = rotation[0][0]*initialCoordinates[i][0] + rotation[0][1]*initialCoordinates[i][1] + rotation[0][2]*initialCoordinates[i][2];
        if (newx > 0) {
            countOn += 1;
            newy = rotation[1][0]*initialCoordinates[i][0] + rotation[1][1]*initialCoordinates[i][1] + rotation[1][2]*initialCoordinates[i][2];
            if (newy < miny) miny = newy;
            if (newy > maxy) maxy = newy;

            newz = rotation[2][0]*initialCoordinates[i][0] + rotation[2][1]*initialCoordinates[i][1] + rotation[2][2]*initialCoordinates[i][2];
            if (newz < minz) minz = newz;
            if (newz > maxz) maxz = newz;
        }
        else {
            countOff += 1;
        }
    }

    if ((countOn > 0) && (countOff > 0)) { // There are both visible and invisible points
                                   // --> active region is on the edge
        // In this situation there are cases where the yz-area define above is
        // actually smaller than the real area of the active region on the stellar disk.
        // The minima/maxima are over/under-estimated if the active region is on one of the
        // axis (y or z). Because if on the y axis, the minimum (or maximum) won t be on the circumference of the active region. Same for z axis
        if (miny*maxy < 0) {       //active region on the z-axis because one point is on the positive side of z, and the other on the negative side of z
            if (minz < 0) {
                minz = -1;   //active region on the bottom-z axis (z<0)
            }
            else {
                maxz = 1;   //active region on the top-z axis (z>=0)
            }
        }
        if (minz*maxz < 0) {       //active region on the y-axis because one point is on the positive side of y, and the other on the negative side of z
            if (miny < 0) {
                miny = -1;   //active region on the left hand-y axis (y<0)
            }
            else {
                maxy = 1;   //active region on the right hand-y axis (y>=0)
            }
        }
    }

    // Indices of miny, minz,... on the grid
    double gridStep = 2./star.gridSize;
    iminy = (int)floor((1.+miny)/gridStep);
    iminz = (int)floor((1.+minz)/gridStep);
    imaxy = (int)ceil((1.+maxy)/gridStep);
    imaxz = (int)ceil((1.+maxz)/gridStep);

    return (countOn > 0);
}


void Spot::scan(double phase, double& flux, std::vector<double>& profile, bool observeRV) {
    unsigned int iy, iz, i;
    double y, z, xsquared;
    double v_shift;
    double spot_temp, limbSum;
    double coordinatesReal[3];
    double depth, r_cos, rSquared;
    std::vector<double> ccfShifted;
    std::vector<double> ccfActiveShifted;

    double inv_phase = phase - 2*pi;
    double inclination = star.inclination;
    double npts;

    // matrix R from spot_inverse_rotation
    double matrixPhase[3][3] = {{(1 - cos(inv_phase)) * cos(inclination) * cos(inclination) + cos(inv_phase), sin(inv_phase) * sin(inclination), (1 - cos(inv_phase)) * cos(inclination) * sin(inclination)},
                                {-sin(inv_phase) * sin(inclination), cos(inv_phase), sin(inv_phase) * cos(inclination)},
                                {(1 - cos(inv_phase)) * sin(inclination) * cos(inclination), -sin(inv_phase) * cos(inclination), (1 - cos(inv_phase)) * sin(inclination) * sin(inclination) + cos(inv_phase)}};

    /*
    double matrixScanner[3][3];
    for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 3; c++) {
            matrixScanner[r][c] = 0.;
            for (int k = 0; k < 3; k++) {
                matrixScanner[r][c] += matrixPhase[r][k] * matrixSpot[k][c];
            }
        }
    }
    */

    // TODO Why doesn't the above code do the same thing
    double matrixScanner[3][3] = {{matrixSpot[0][0]*matrixPhase[0][0]+matrixSpot[0][1]*matrixPhase[1][0]+matrixSpot[0][2]*matrixPhase[2][0],
                                   matrixSpot[0][0]*matrixPhase[0][1]+matrixSpot[0][1]*matrixPhase[1][1]+matrixSpot[0][2]*matrixPhase[2][1],
                                   matrixSpot[0][0]*matrixPhase[0][2]+matrixSpot[0][1]*matrixPhase[1][2]+matrixSpot[0][2]*matrixPhase[2][2]},
                                  {matrixSpot[1][0]*matrixPhase[0][0]+matrixSpot[1][1]*matrixPhase[1][0]+matrixSpot[1][2]*matrixPhase[2][0],
                                   matrixSpot[1][0]*matrixPhase[0][1]+matrixSpot[1][1]*matrixPhase[1][1]+matrixSpot[1][2]*matrixPhase[2][1],
                                   matrixSpot[1][0]*matrixPhase[0][2]+matrixSpot[1][1]*matrixPhase[1][2]+matrixSpot[1][2]*matrixPhase[2][2]},
                                  {matrixSpot[2][0]*matrixPhase[0][0]+matrixSpot[2][1]*matrixPhase[1][0]+matrixSpot[2][2]*matrixPhase[2][0],
                                   matrixSpot[2][0]*matrixPhase[0][1]+matrixSpot[2][1]*matrixPhase[1][1]+matrixSpot[2][2]*matrixPhase[2][1],
                                   matrixSpot[2][0]*matrixPhase[0][2]+matrixSpot[2][1]*matrixPhase[1][2]+matrixSpot[2][2]*matrixPhase[2][2]}};

    for (iy = iminy; iy < imaxy; iy++) {
        y = -1.0+iy*2.0/star.gridSize;
        coordinatesReal[1] = y;

        if (observeRV) {
            v_shift = y * star.vrot * sin(star.inclination);
            ccfShifted = star.profileQuiet.shift(v_shift);
            ccfActiveShifted = star.profileActive.shift(v_shift);
        }

        limbSum = 0;
        npts = 0;
        double intensitySum = 0;

        for (iz = iminz; iz < imaxz; iz++) {
            z = -1.0+iz*2.0/star.gridSize;
            xsquared = y*y + z*z;
            if (xsquared < 1.) { // If on the disk
                coordinatesReal[0] = sqrt(1.-xsquared);
                coordinatesReal[2] = z;

                // Rotate spot to the center of the star and check the x-coordinate, depth.
                // This is a nifty way to check if the coordinate is on the spot but assumes
                depth = matrixScanner[0][0]*coordinatesReal[0] + matrixScanner[0][1]*coordinatesReal[1] + matrixScanner[0][2]*coordinatesReal[2];

                if (depth*depth >= (1 - size*size)) { // If actually on the spot
                    rSquared = (y*y + z*z);

                    if (plage) {
                        r_cos = sqrt(1. - rSquared);
                        spot_temp = star.temperature + 250.9 - 407.7 * r_cos + 190.9 * r_cos*r_cos;
                        intensitySum += planck(5000e-10, spot_temp) / star.intensity;
                        limbSum += 1 - star.limbLinear * (1 - r_cos) - star.limbQuadratic * (1 - r_cos)*(1 - r_cos);
                        npts += 1;
                    }
                    else {
                        r_cos = sqrt(1. - rSquared);
                        limbSum += 1 - star.limbLinear * (1 - r_cos) - star.limbQuadratic * (1 - r_cos)*(1 - r_cos);
                    }
                }
            }
        }
        for (i = 0; i < ccfShifted.size(); i++) {
            profile[i] += (ccfShifted[i] - intensity * ccfActiveShifted[i]) * limbSum;
        }

        flux += (1-intensity) * limbSum;
    }
}
