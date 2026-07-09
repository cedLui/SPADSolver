#pragma once

#include <iostream>
#include <cmath>
#include <vector>

double alphaX(double Pos);

double betaX(double Pos);

std::vector<std::vector<double>> avaProb(double Width, double StepSize, double P_h_0);