#include <iostream>
#include <cmath>
#include <vector> 

#include "interpolation.cpp"
#include "avaProb.cpp"

double ReflectProb(double RefInd);

std::vector<double> AbsorptionProbDist(double Width, int Steps, std::vector<double> AlFracProfile, double PhotonEnergy);

std::vector<double> AvaPairProb(double Width, int Steps, double P_h_0, double Accuracy);

double SPDE(double RefInd, double Width, int Steps, std::vector<double> AlFracProfile, double PhotonEnergy, double P_h_0, double Accuracy);