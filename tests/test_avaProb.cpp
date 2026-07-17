#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>

#include "avaProb.cpp"
#include "avaProb.h"

double analyiticalAvaProbTest(double Width, int Steps, double P_h_0, double Accuracy, double alpha, double beta){
    //When alpha and beta are constant, we have an analytical solution
    //We need to solve for the constant k using our Guess from the numerical solution, and check and see if the equations match values close enough

    //First, get the numerical answers for P_e and P_h
    std::vector<std::vector<double>> AvaProbVec = avaProb(Width, Steps, P_h_0, Accuracy);
    std::vector<double> P_e = AvaProbVec[0];
    std::vector<double> P_h = AvaProbVec[1];

    //Plug P_h[0] in for k
    double k = (std::pow((1.0 - P_h[0]), alpha/beta - 1.0) - std::exp((beta - alpha) * Width))/(1.0 - std::exp((beta - alpha) * Width));

    double P_e_analytic;
    double P_h_analytic;

    //Check that P_e(W) and P_h(W) are consistent
    //Unfinished
    return 0.0;
}

