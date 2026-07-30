#pragma once

#include <iostream>
#include <cmath>
#include <vector>

double alphaX(double EField, double T);

double betaX(double EField, double T);

std::vector<double> guess100(double Width, int Steps, double EField);

std::vector<std::vector<double>> avaProb(double Width, int Steps, double Accuracy, double EField);