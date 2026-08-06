#include <iostream>
#include <vector>

#include "interpolation.h"
#include "interpolation.cpp"

//Returns all gamma values for mole fractions between 0 and 1 in increments of stepSize
std::vector<double> plot_gamma(double stepSize, double Energy){
    int numSteps = 1 / stepSize;
    std::vector<double> gammaValues(numSteps);
    for(int i = 0; i < numSteps; i++){
        double AlFrac = i * stepSize;
        gammaValues[i] = gamma(AlFrac, Energy);
    }
    return gammaValues;
}

int main(){
    std::vector<double> gammaValues = plot_gamma(0.005, 6);
    for(int i = 0; i < gammaValues.size(); i++){
        std::cout << gammaValues[i] << std::endl;
    }
    return 0;
}


