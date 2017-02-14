#include "spot.hpp"


const double pi = M_PI;


Spot::Spot(Star star, double latitude, double longitude, double fillfactor, bool plage, size_t spotResolution) {
    latitude *= pi/180.;
    longitude *= pi/180.;
    size = sqrt(2*fillfactor);
    this -> star = star;
    this -> plage = plage;
    this -> temperature = star.temperature - star.spotTempDiff;
    this->latitude = latitude;
    this->longitude = longitude;
    opening_angle = acos(sqrt(1-size*size));
    double height = 1-sqrt(1-size*size);
    radius = 2*height*(height-1);
    auto h = 1-sqrt(1-size*size);
    radius_squared = 2*h - h*h + h*h;

    //Compute initial location
    std::vector<Point> centeredCoordinates(spotResolution);
    initialCoordinates = std::vector<Point>(spotResolution);


    for (auto i = 0; i < initialCoordinates.size(); i++) {
        auto rho = -pi + i*2*pi/(spotResolution-1);     // For this case, phase goes from -pi to pi
        centeredCoordinates[i].x = sqrt(1 - size*size);
        centeredCoordinates[i].y = size * cos(rho);
        centeredCoordinates[i].z = size * sin(rho);
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
        initialCoordinates[i].x = matrixInit[0][0]*centeredCoordinates[i].x + matrixInit[0][1]*centeredCoordinates[i].y + matrixInit[0][2]*centeredCoordinates[i].z;
        initialCoordinates[i].y = matrixInit[1][0]*centeredCoordinates[i].x + matrixInit[1][1]*centeredCoordinates[i].y + matrixInit[1][2]*centeredCoordinates[i].z;
        initialCoordinates[i].z = matrixInit[2][0]*centeredCoordinates[i].x + matrixInit[2][1]*centeredCoordinates[i].y + matrixInit[2][2]*centeredCoordinates[i].z;
    }
}


spot_bounds Spot::get_bounds_at(double phase) {

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
    miny = 1;
    maxy = -1;
    minz = 1;
    maxz = -1;

    // Search for bounds by applying the rotation matrix to the initial coordinates
    // x coordinate is depth, so proceed only if the spot is on the front of the star
    for (auto i = 0; i < initialCoordinates.size(); i++) {
        newx = rotation[0][0]*initialCoordinates[i].x + rotation[0][1]*initialCoordinates[i].y + rotation[0][2]*initialCoordinates[i].z;
        if (newx > 0) {
            countOn += 1;
            newy = rotation[1][0]*initialCoordinates[i].x + rotation[1][1]*initialCoordinates[i].y + rotation[1][2]*initialCoordinates[i].z;
            miny = std::min(miny, newy);
            maxy = std::max(maxy, newy);

            newz = rotation[2][0]*initialCoordinates[i].x + rotation[2][1]*initialCoordinates[i].y + rotation[2][2]*initialCoordinates[i].z;
            minz = std::min(minz, newz);
            maxz = std::max(maxz, newz);
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

    //std::cout << minz << " " << maxz << " " << miny << " " << maxy << std::endl;

    miny = floor(miny/star.grid_interval)*star.grid_interval;
    maxy = ceil(maxy/star.grid_interval)*star.grid_interval;
    minz = floor(minz/star.grid_interval)*star.grid_interval;
    maxz = ceil(maxz/star.grid_interval)*star.grid_interval;


    // Indices of miny, minz,... on the grid, eventually comment this out
    double gridStep = star.grid_interval;
    iminy = (int)floor((1.+miny)/gridStep);
    iminz = (int)floor((1.+minz)/gridStep);
    imaxy = (int)ceil((1.+maxy)/gridStep);
    imaxz = (int)ceil((1.+maxz)/gridStep);

    //return (countOn > 0);

    double theta = phase + longitude;
    double phi = M_PI_2 - latitude;
    center = Point(sin(phi)*cos(theta),
                   sin(phi)*sin(theta),
                   cos(phi));
    center = center.rotate_y(star.inclination - M_PI_2);

    //std::cout << center.x << " " << center.y << " " << center.z << std::endl;

    return {countOn > 0, miny, maxy, minz, maxz};
}


std::vector<double> Spot::get_ccf(double phase, double wavelength) {
    std::vector<double> profile(star.profileActive.size());
    auto bounds = get_bounds_at(phase);
    if (not bounds.visible) {
        return profile;
    }

    for (auto y = miny; y < maxy; y += star.grid_interval) {

        auto v_shift = y * star.vrot * sin(star.inclination);
        auto& ccfShifted = star.profileQuiet.shift(v_shift);
        auto& ccfActiveShifted = star.profileActive.shift(v_shift);

        double limbSum = 0;

        for (auto z = minz; z < maxz; z += star.grid_interval) {
            auto x = sqrt(1 - z*z - y*y);
            if (is_on_spot(x, y, z)) {
                limbSum += star.limb_brightness(x, y, z);
            }
        }

        for (auto i = 0; i < ccfShifted.size(); i++) {
            profile[i] += (ccfShifted[i] - intensity * ccfActiveShifted[i]) * limbSum;
        }
    }
    return profile;
}


double Spot::get_flux(double phase, double wavelength) {
    auto bounds = get_bounds_at(phase);
    if (not bounds.visible) {
        return 0.0;
    }

    double limb_integral = 0.0;
    for (auto y = bounds.miny; y < bounds.maxy; y += star.grid_interval) {
        for (auto z = bounds.minz; z < bounds.maxz; z += star.grid_interval) {
            auto x = sqrt(1 - y*y - z*z);
            if (is_on_spot(x, y, z)) {
                limb_integral += star.limb_brightness(x, y, z);
            }
        }
    }
    auto intensity = planck(wavelength, temperature) / star.intensity;
    return (1-intensity)*limb_integral;
}


bool Spot::is_on_spot(double x, double y, double z) {
    if ((y*y + z*z) >= 1.0) { // Check if on the star disk
        return false;
    }
    else {
        double distance_squared = (x-center.x)*(x-center.x) + (y-center.y)*(y-center.y) + (z-center.z)*(z-center.z);
        return distance_squared < radius_squared;
    }
}


bool Spot::isVisible2(double phase) {
    // Location in spherical coordinates
    double theta = phase + longitude;
    double phi = M_PI_2 - latitude;

    // Convert to cartesian
    auto center = Point(sin(phi)*cos(theta),
                        sin(phi)*sin(theta),
                        cos(phi));

    // Apply inclination
    center = center.rotate_y(star.inclination - M_PI_2);

    Point top_point = Point(sin(phi-opening_angle)*cos(theta),
                            sin(phi-opening_angle)*sin(theta),
                            cos(phi-opening_angle));

    double bottom = Point(sin(phi+opening_angle)*cos(theta),
                          sin(phi+opening_angle)*sin(theta),
                          cos(phi+opening_angle)).z;

    bool visible = center.x > -sqrt(2*size);
    return visible;
}
