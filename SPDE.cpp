#include <iostream>
#include <cmath>
#include <vector> 

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
    double prevGamma = gamma(AlFracProfile[Steps], PhotonEnergy);
    double Gamma;
    double IntGamma = 0;
    std::vector<double> AbsProbDist;

    for (int k = 0; k < Steps + 1; k++){ //Initialize the probability distribution vectors
        AbsProbDist.push_back(0);
    }

    AbsProbDist[Steps] = gamma(AlFracProfile[Steps], PhotonEnergy); //Initialize the rightmost value

    //For every step in the Al mole fraction profile
    for (int i = Steps - 1; i > -1; i--){
        Gamma = gamma(AlFracProfile[i], PhotonEnergy); //Calculate absorption coefficient
        
        IntGamma += (prevGamma + Gamma) * StepSize/2.0; //Do a trapazoidal sum

        prevGamma = Gamma;

        AbsProbDist[i] = (Gamma * std::exp(-1*IntGamma));
    }

    return AbsProbDist;
}

std::vector<double> AvaPairProb(double Width, int Steps, double Accuracy, double EField, std::vector<double> AlFracProf){

    std::vector<std::vector<double>> PVec = avaProb(Width, Steps, Accuracy, EField, AlFracProf); //Calculate P_e and P_h
    std::vector<double> P_e = PVec[0];
    std::vector<double> P_h = PVec[1];
    std::vector<double> P_p;

    for (int i = 0; i < Steps + 1; i++){ //Calculate the probability that a carrier pair at any x-coor in the device can avalanche
        P_p.push_back(P_e[i] + P_h[i] - P_e[i]*P_h[i]);
    }

    return P_p;
}

double SPDE(double RefInd, double Width, int Steps, std::vector<double> AlFracProfile, double PhotonEnergy, double Accuracy, double EField){
    double refProb = ReflectProb(RefInd); //Calculate probability of reflection
    std::vector<double> absProbDist = AbsorptionProbDist(Width, Steps, AlFracProfile, PhotonEnergy); //Calculate absorption probability distribution
    std::vector<double> avaPairProb = AvaPairProb(Width, Steps, Accuracy, EField, AlFracProfile); //Calculate probability of an avalanche at points through the device

    double StepSize = Width/Steps; //Calculate the size of a step

    double IntAbsAva = 0.0;

    double prev = absProbDist[0] * avaPairProb[0]; //For trapazoidal sum
    double curr;
    
    for (int i = 1; i < Steps+1; i++){ //Calculate IntAbsAva using a trapazoidal sum 
        curr = absProbDist[i] * avaPairProb[i];
        IntAbsAva += (prev+curr) * StepSize/2.0;
        prev = curr;
    }

    return (1 - refProb) * IntAbsAva;
}

int main(){
    int Steps = 1000000;
    double spde;
    std::vector<double> AlFracProf;
    for (int i = 0; i < Steps + 1; i++){
        AlFracProf.push_back(0);
    }

    for (double j = 0.000034; j <= 0.00009; j += 0.000001){
        spde = SPDE(1.0, j, Steps, AlFracProf, 6.0, 0.000001, 3000000);
        std::cout << spde << std::endl;
    }
    
    
    return 0;
}