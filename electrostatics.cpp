#include <iostream>
#include <cmath>
#include <vector>
#include "electrostatics.h"

double heavisideFunc(double position, double shift){

    if(position >= shift){
        return 1;
    }
    else if(position < shift){
        return 0;
    }
    return 0;
}

double Field(double bias, double region, double x, double td, double rho) {

    double epsilon0 = 8.85e-12 * 100; // F/cm
    double epsilon_r = 8.9;
    double epsi_mat = epsilon_r * epsilon0; // F/cm
    double bandgap = 3.4 * 1.602e-19; //Joules
    double q = 1.602e-19;  // Fundamental Charge in [C]

    double V_bi = bandgap/q; // V

    //Define field profile

    double field = (rho/epsi_mat) * (heavisideFunc(x, td) + (td/region - 1)) + (bias-V_bi)/region;

    return field;
}