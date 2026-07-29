#include <iostream>
#include <vector>
#include <cmath>
#include <fstream> 
#include "GaussianElimination.h"
#include "interpolation.h"
#include "electrostatics.h"

using namespace std;

int main(){

    const int bins = 25;
    const double region = 0.5e-4; //cm

    double t_d = region/3; // cm

    int columns = bins + 1;

    vector<double> deltaJmatrix((bins + 1) * bins);

    double photon_energy = 9.613e-19; //Photon energy in Joules
    double q = 1.602e-19; //Fundamental Charge in C 

    double temp = 300;  //K

    double GaN_frac = 0; // Al Mole Fraction 

    //Boundary Conditions

    double jp0 = 0; // A/m^2

    double jn0 = 50 * 1000; // A/m^2

    double bias = -90; //V

    double step = region/bins;

    double sigma = q * 10e12; // C/cm^2

    double C = 1; //ARBITRARY CONSTANT IN UNITS OF INTENSITY [eV / cm^2]

    double GaNgamma = gamma(0, photon_energy/(1.602e-19)); // 1/cm, absorption interpolation takes energy values of [eV]
    //double GaNgamma = 0;

    //CALCULATING FIELD VALUES
    vector<double> field_vals(bins + 1);

    for(int i = 0; i < bins+1; i++) {
        field_vals[i] = -1 * Field(bias, region, i * step, t_d, sigma); // [V/cm]
        //field_vals[i] = 2.5e6; If you want a constant E field [V/cm]
    }

    //CALCULATING AVALANCHE COEFF VALUES

    vector<double> alpha_vals(bins+1);

    for(int i = 0; i < bins+1; i++) {
        alpha_vals[i] = alpha(GaN_frac, abs(field_vals[i]), temp);
    }

    vector<double> beta_vals(bins + 1);

    for(int i = 0; i < bins + 1; i++) {
        beta_vals[i] = beta(GaN_frac, abs(field_vals[i]), temp);
    }

    vector<double> avSums(bins + 1);

    for(int i = 0; i < bins + 1; i ++) {
        avSums[i] = alpha_vals[i] + beta_vals[i];
    }

    vector<double> avDiffs(bins + 1);

    for(int i = 0; i < bins + 1; i++) {
        avDiffs[i] = alpha_vals[i] - beta_vals[i];
    }

    /*
    BUILDING MATRIX!!!!!!!!!!!!!!!!!!!!!!!!
    */

    //BUILDING FIRST ROW

    deltaJmatrix[0 * columns + 0] = 1;

    deltaJmatrix[0 * columns + (bins - 1)] = (1 + (avSums[0]/avDiffs[0]) * (1 - exp(avDiffs[0] * step))) * exp(-avDiffs[0] * step);

    deltaJmatrix[0 * columns + bins] = (2 * q * C * GaNgamma/photon_energy) * (1/(avDiffs[0] + GaNgamma)) * exp(-GaNgamma * region)
            * (exp((avDiffs[0] + GaNgamma) * step) - 1) *  exp(-avDiffs[0] * step) + 2 * jp0 * exp(-avDiffs[0] * step);

    //deltaJmatrix[0 * columns + bins] = 2 * jp0 * exp(-(avDiffs[0] * step));

    //BUILDING LAST ROW 

    deltaJmatrix[(bins - 1) * columns + (bins - 2)] = 1;

    deltaJmatrix[(bins - 1) * columns + (bins - 1)] = -1 * (exp(avDiffs[bins - 1] * step) + (avSums[bins - 1]/avDiffs[bins - 1])
            * (1 - exp(avDiffs[bins - 1] * step)));
    
    deltaJmatrix[(bins - 1) * columns + bins] = (-1 * 2 * q * GaNgamma * C/photon_energy) * (1/(avDiffs[bins - 1] + GaNgamma))
            * (exp((avDiffs[bins - 1] + GaNgamma) * step) - 1) + -1 * 2 * jn0 * exp(avDiffs[bins - 1] * step);

    //BUILDING MATRIX

    for(int i = 1; i < bins - 1; i ++) {
        deltaJmatrix[i * columns + (i - 1)] = 1;
        deltaJmatrix[i * columns + i] = -1 * exp((avDiffs[i]) * step);
        deltaJmatrix[i * columns + (bins - 1)] = -1 * (avSums[i]/avDiffs[i]) * (1 - exp(avDiffs[i] * step));
        deltaJmatrix[i * columns + bins] = (-1 * 2 * GaNgamma * C * q/photon_energy) * (1/(avDiffs[i] + GaNgamma))
                * exp(-1 * GaNgamma * (region - i * step)) * (exp((avDiffs[i] + GaNgamma) * step) - 1);
    }

    /*for (auto deltaCheck : deltaJmatrix) {
     cout << deltaCheck << ", ";
    }*/

    cout << "\n";

    rowEchelon(deltaJmatrix, bins);

    vector<double> deltaJsol = solVector(deltaJmatrix, bins);

    double const j_tot = deltaJsol[bins - 1];

    //Once the solution vector is found, turn it into Jp and Jn matrix along each partition

    vector<double> Jp(bins + 1);

    Jp[0] = jp0;

    Jp[bins] = j_tot - jn0;

    vector<double> Jn(bins + 1);

    Jn[0] = j_tot - jp0;

    Jn[bins] = jn0;

    for(int i=1; i < bins; i++) {
        Jp[i] = 0.5 * (j_tot + deltaJsol[i-1]);
    }

    for(int i=1; i < bins; i++) {
        Jn[i] = 0.5 * (j_tot - deltaJsol[i-1]);
    }
    
    //VECTOR CHECK

    cout << "JP CHECK" << "\n";
    for (auto JpCheck : Jp) {
     cout << JpCheck << ", ";
    }
     cout << "\n";

    cout << "JN CHECK" << "\n";

    for (auto JnCheck : Jn) {
     cout << JnCheck << ", ";
    }

    //cout << "\n";

    /*cout << "FIELD CHECK" << "\n";

    for (auto fieldCheck : field_vals) {
     cout << fieldCheck << ", ";
    }*/

    //cout << "\n";

    /*cout << "AVSUM CHECK" << "\n";

    for (auto avsumCheck : avSums) {
     cout << avsumCheck << ", ";
    }

    cout << "AVDIFF CHECK" << "\n";

    for (auto avdiffCheck : avDiffs) {
     cout << avdiffCheck << ", ";
    }*/
}