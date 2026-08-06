#pragma once

#include <iostream>
#include <cmath>
#include <vector>

#include "interpolation.h"
#include "electrostatics.h"


double alphaX(double EField, double T);

double betaX(double EField, double T);

std::vector<double> guess100(double Width, int Steps, double EField, std::vector<double> AlFracProf);

std::vector<std::vector<double>> avaProb(double Width, int Steps, double Accuracy, double EField, std::vector<double>AlFracProf);