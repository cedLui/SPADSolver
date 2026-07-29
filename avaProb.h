#pragma once

#include <iostream>
#include <cmath>
#include <vector>

double alphaX(double DeltaV, double Pos, double Width);

double betaX(double DeltaV, double Pos, double Width);

std::vector<double> guess100(double Width, int Steps, double DeltaV);

std::vector<std::vector<double>> avaProb(double Width, int Steps, double P_h_0, double Accuracy, double DeltaV);