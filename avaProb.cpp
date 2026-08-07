#include <iostream>
#include <cmath>
#include <vector> 

#include "avaProb.h"

//These are the impact ionization coefficients for the material. They are taken from the Maeda paper
double alphaX(double EField, double T){ //In cm^-1. EField is in V/cm, T is in Kelvin
    double a_n =  2.69e7; // cm^-1
    double b_n = 2.00e-3; // K^-1
    double c_n = 2.27e7; // V * cm^-1
    double d_n = 5.00e-4; // K^-1

    return a_n*(1 + b_n*(T - 298)) * std::exp(-c_n*(1+d_n*(T - 298))/EField);

    /*This was to replicate Oldham's findings
    double V_break = 27.08;
    double E = (V_break + DeltaV)/Width; //Avg electric field (testing the constant case)
    return 3.8 * 1000000 * std::exp(-1.75*1000000/E);*/
}

double betaX(double EField, double T){ //In cm^-1. EField is in V/cm, T is in Kelvin
    double a_p =  4.32e6; // cm^-1
    double b_p = 2.00e-3; // K^-1
    double c_p = 1.31e7; // V * cm^-1
    double d_p = 9.00e-4; // K^-1

    return a_p*(1 + b_p*(T - 298)) * std::exp(-c_p*(1+d_p*(T - 298))/EField);

    /*This was to replicate Oldham's findings
    double V_break = 27.08;
    double E = (V_break + DeltaV)/Width; //Avg electric field profile
    return 2.25 * 10000000 * std::exp(-3.26 * 1000000/E);*/
}

std::vector<double> guess100(double Width, int Steps, double Bias, double td, double rho, std::vector<double> AlFracProf){
    double StepSize = Width/Steps;
    double EField = 0; //Initialize Electric Field

    std::vector<double> Guesses;
    std::vector<double> F; //Holds P_h.at(Steps)

    std::vector<double> P_e;
    std::vector<double> P_h;

    for (int i = 0; i < 100; i++){ //Initializes a guess for every 0.01
        Guesses.push_back(i/100.0);
    }

    for (int j = 0; j < Steps + 1; j++){ //Initialize P_e and P_h
        P_e.push_back(0);
        P_h.push_back(0);
    }

    for (double Guess : Guesses){ //For each guess, run the solver
        P_e.at(0) = 0.0; //Reinit P_e(0) = 0
        P_h.at(0) = Guess; //Init guess

        for (int k = 0 ; k<Steps; k++){
            EField = Field(Bias, Width, k*StepSize, td, rho);//Calculate Efield for each point in the APD

            P_e.at(k+1) = P_e.at(k) + (1- P_e.at(k)) * alpha(AlFracProf[k], EField, 298) * (P_e.at(k)+P_h.at(k)-P_e.at(k)*P_h.at(k)) * StepSize;
            P_h.at(k+1) = P_h.at(k) - (1- P_h.at(k)) * beta(AlFracProf[k], EField, 298) * (P_h.at(k)+P_e.at(k)-P_e.at(k)*P_h.at(k)) * StepSize;
        }
        F.push_back(P_h.at(Steps)); //Creates a vector of all P_h(W), which should be close to zero
    }

    for (int i = 2; i<Steps + 1; i++){
        if (std::abs(F[i-1] + F[i]) != std::abs(F[i-1]) + std::abs(F[i])){
            std::cout << i-1 << std::endl;
            std::cout << i << std::endl;
        }
    }

    return F;
}

//This function will return 2 vectors: 1st a descretized version of P_e, and the 2nd a descretized version of P_h
std::vector<std::vector<double>> avaProb(double Width, int Steps, double Accuracy, double Bias, double td, double rho, std::vector<double> AlFracProf){
    int LoopCount = 0; //Counts how many iterations
    double StepSize = Width / Steps;
    double Guess = 1; //Guess for P_h(0)
    double prevGuess = 1; //Previous guess for P_h(0)

    double EField = 0; //Initialize EField

    double left; //For binary search
    double right;

    double PhW = 0; //Current and previous P_h(W)
    double prevPhW = 0;

    bool firstIter = true;
    bool binarySearch = false;


    std::vector<std::vector<double>> Pair; //Vector of vectors to return

    //Equations:
    //P_e(i+1) = P_e(i) + (1- P_e(i))*alphaX(x)*(P_e(i)+P_h(i)-P_e(i)*P_h(i)) * StepSize
    //P_h(i+1) = P_h(i) - (1- P_h(i))*betaX(x)*(P_h(i)+P_e(i)-P_e(i)*P_h(i)) * StepSize

    //P_e(0) = 0
    //P_h(Width) = 0
    //P_h(0) = P_h_0

    //Define Lists for P_e and P_h
    std::vector<double> P_e;
    std::vector<double> P_h;

    //Initialize P_e and P_h, creating a list of zeros of size Steps + 1
    for (int i = 0; i <= Steps; i++){
        P_e.push_back(0);
        P_h.push_back(0);
    }

    //Start shooting method while loop
    while (Guess > 0){ //All values should be between 0 and 1, since they are probabilities
        LoopCount++; //Update the loop counter
        
        P_h.at(0) = Guess; //Update guess for P_h(0)
        P_e.at(0) = 0; //Reinitialize guess for P_e(0)
        
        for (int i = 0 ; i<Steps; i++){
            EField = Field(Bias, Width, i*StepSize, td, rho); //Calculate the electric field at this point

            P_e.at(i+1) = P_e.at(i) + (1- P_e.at(i)) * alpha(AlFracProf[i], EField, 298) * (P_e.at(i)+P_h.at(i)-P_e.at(i)*P_h.at(i)) * StepSize;
            P_h.at(i+1) = P_h.at(i) - (1- P_h.at(i)) * beta(AlFracProf[i], EField, 298) * (P_h.at(i)+P_e.at(i)-P_e.at(i)*P_h.at(i)) * StepSize;
        }

        prevPhW = PhW; //Update PhWs
        PhW = P_h.at(Steps);

        if (std::abs(P_h.at(Steps)) < Accuracy){ //If P_h(Width) is close enough to 0, return the lists
            Pair.push_back(P_e);
            Pair.push_back(P_h);
            return Pair;
        }


        if (firstIter){ 
            firstIter = false;
            prevGuess = Guess;
            Guess = Guess - 0.01;
        } else{
            if (!binarySearch){ //If the binary search hasn't been activated
                if (std::abs(PhW + prevPhW) != std::abs(PhW) + std::abs(prevPhW)){ //If there was a sign change, we know the real P_h_0 is between the current and previous guesses
                    binarySearch = true;
                    if (PhW < 0){
                        left = Guess;
                        right = prevGuess;
                    } else{
                        right = Guess;
                        left = prevGuess;
                    }
                    Guess = (prevGuess + Guess)/2.0;
                } else{
                    prevGuess = Guess;
                    Guess = Guess - 0.01;
                }
            } else{ //The binary search has been activated
                    if (PhW < 0){
                        left = Guess;
                        Guess = (left + right)/2.0;
                    } else{
                        right = Guess;
                        Guess = (left + right)/2.0;
                    }
            }
        }
    }
    //std::cout << "No solution found" << std::endl;
    //Pair.push_back(P_e);
    //Pair.push_back(P_h);
    return Pair; //If the while loop ends, return the lists
}

int main(){
    std::vector<double> AlFracProf;
    for (int x = 0; x < 10000; x++){
        AlFracProf.push_back(0);
    }
    
    std::vector<std::vector<double>> Pair = avaProb(0.00005, 10000, 0.00000001, 150, 0.00001, -182e-9, AlFracProf);
    for (int i = 0; i <= 10000; i = i + 200){
        std::cout << Pair[0][i] << std::endl;
    }
}

