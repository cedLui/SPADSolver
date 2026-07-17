#pragma once

#include <iostream>
#include <cmath>
#include <vector>

double alphaX(double Pos);

double betaX(double Pos);

std::vector<std::vector<double>> avaProb(double Width, int Steps, double P_h_0, double Accuracy);