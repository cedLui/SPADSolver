#include <iostream>
#include <cmath>
#include <vector> 

#include "interpolation.cpp"
#include "avaProb.cpp"
#include "SPDE.h"

double ReflectProb(double RefInd){
    //Use Fresnel's equation for reflectivity
    double RefIndAir = 1;
    return std::pow(((RefInd - RefIndAir)/(RefInd + RefIndAir)),2);
}

std::vector<double> AbsorptionProbDist(double Width, int Steps, std::vector<double> AlFracProfile, double PhotonEnergy){
    //Given the composition of the material, determine the Probability distribution for absorption
    //The size of AlFracProfile should be Steps + 1
    
    double StepSize = Width/Steps;
    double Gamma;
    double IntGamma;
    std::vector<double> AbsProbDist;

    //For every step in the Al mole fraction profile
    for (int i = 0; i < StepSize + 1; i++){
        Gamma = gamma(AlFracProfile[i], PhotonEnergy); //Calculate absorption coefficient

        if (i == 0){
            IntGamma = 0;
        } else{
            IntGamma += Gamma * StepSize; //Do a right hand riemman sum to approximate the integral of Gamma over distance
        }

        AbsProbDist.push_back(Gamma * std::exp(-1*IntGamma));
    }

    return AbsProbDist;
}

std::vector<double> AvaPairProb(double Width, int Steps, double P_h_0, double Accuracy){

    std::vector<std::vector<double>> PVec = avaProb(Width, Steps, P_h_0, Accuracy); //Calculate P_e and P_h
    std::vector<double> P_e = PVec[0];
    std::vector<double> P_h = PVec[1];
    std::vector<double> P_p;

    for (int i = 0; i < Steps + 1; i++){ //Calculate the probability that a carrier pair at any x-coor in the device can avalanche
        P_p.push_back(P_e[i] + P_h[i] - P_e[i]*P_h[i]);
    }

    return P_p;
}

double SPDE(double RefInd, double Width, int Steps, std::vector<double> AlFracProfile, double PhotonEnergy, double P_h_0, double Accuracy){
    double refProb = ReflectProb(RefInd); //Calculate probability of reflection
    std::vector<double> absProbDist = AbsorptionProbDist(Width, Steps, AlFracProfile, PhotonEnergy); //Calculate absorption probability distribution
    std::vector<double> avaPairProb = AvaPairProb(Width, Steps, P_h_0, Accuracy); //Calculate probability of an avalanche at points through the device

    double StepSize = Width/Steps; //Calculate the size of a step

    double IntAbsAva;
    
    for (int i = 0; i < Steps + 1; i++){ //Calculate IntAbsAva using a left handed riemann sum
        IntAbsAva += absProbDist[i] * avaPairProb[i] * StepSize;
    }

    return (1 - refProb) * IntAbsAva;

}