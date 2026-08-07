#pragma once

#include <iostream>
#include <cmath>
#include <vector>

#include "interpolation.h"
#include "electrostatics.h"

std::vector<double> guess100(double Width, int Steps, double Bias, double td, double rho, std::vector<double> AlFracProf);

std::vector<std::vector<double>> avaProb(double Width, int Steps, double Accuracy, double Bias, double td, double rho, std::vector<double>AlFracProf);