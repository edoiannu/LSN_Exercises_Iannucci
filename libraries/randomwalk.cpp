#include <iostream>
#include <fstream>
#include "randomwalk.h"

using namespace std;

// Constructor that sets the starting point of the random walk at the origin of the coordinate axes
RandomWalk::RandomWalk() {
    x = 0;
    y = 0;
    z = 0;
}

// Method to set the coordinates of the random walk to the specified values (x_0, y_0, z_0)
void RandomWalk::SetCoords(double x_0, double y_0, double z_0) {
    x = x_0;
    y = y_0;
    z = z_0;
}

// Method to get the x-coordinate of the random walk
double RandomWalk::GetX() {
    return x; // Return the x-coordinate
}

// Method to get the y-coordinate of the random walk
double RandomWalk::GetY() {
    return y; // Return the y-coordinate
}

// Method to get the z-coordinate of the random walk
double RandomWalk::GetZ() {
    return z; // Return the z-coordinate
}


// Method to evolve the random walk in the x-direction on a lattice
void RandomWalk::Evolvex_lattice(int x_n, double a) {
    x += 2 * a * (x_n - 0.5);
}

// Method to evolve the random walk in the y-direction on a lattice
void RandomWalk::Evolvey_lattice(int y_n, double a) {
    y += 2 * a * (y_n - 0.5);
}

// Method to evolve the random walk in the z-direction on a lattice
void RandomWalk::Evolvez_lattice(int z_n, double a) {
    z += 2 * a * (z_n - 0.5);
}

// Method to evolve the random walk in polar coordinates (phi, theta) with a step size of 'a'
void RandomWalk::Evolve_polar(double phi, double theta, double a) {
    x += a * sin(theta) * cos(phi);
    y += a * sin(theta) * sin(phi);
    z += a * cos(theta);
}

// Method to set the coordinates of the random walk to the specified values (x_new, y_new, z_new)
void RandomWalk::Evolve_coords(double x_new, double y_new, double z_new) {
    x += x_new;
    y += y_new;
    z += z_new;
}

// Method to calculate and return the squared distance from the origin of the random walk
double RandomWalk::Getr2() {
    return x * x + y * y + z * z;
}
