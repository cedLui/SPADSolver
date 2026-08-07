#include <iostream>
#include <cmath>
#include <vector> 

#include "avaProb.h"

double ReflectProb(double RefInd);

std::vector<double> AbsorptionProbDist(double Width, int Steps, std::vector<double> AlFracProfile, double PhotonEnergy);

std::vector<double> AvaPairProb(double Width, int Steps, double Accuracy, double Bias, double td, double rho, std::vector<double> AlFracProfile);

double SPDE(double RefInd, double Width, int Steps, std::vector<double> AlFracProfile, double PhotonEnergy, double Accuracy, double Bias, double td, double rho);