#include <iostream>
#include <cmath>
#include <vector> 
#include "avaProb.h"

//These are the impact ionization coefficients for the material. Currently they are constant, but they can be changed to be functions of position
double alphaX(double Pos){
    return 0.2;
}

double betaX(double Pos){
    return 0.2;
}

//This function will return 2 vectors: 1st a descretized version of P_e, and the 2nd a descretized version of P_h
std::vector<std::vector<double>> avaProb(double Width, double Steps, double P_h_0, double Accuracy){
    double StepSize = Width / Steps;
    double Guess = P_h_0; //Guess for P_h(0)
    double F = 0; //Value of P_h(Width)

    double previousGuess = 0; //Previous guess & previous value of P_h(Width) for secant method
    double previousF = 0;

    double tempGuess = 0;

    bool firstIter = true;

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
    while (Guess >= 0 && Guess <= 1 && Guess != previousGuess){ //All values should be between 0 and 1, since they are probabilities
        P_h.at(0) = Guess; //Update guess for P_h(0)
        P_e.at(0) = 0; //Reinitialize guess for P_e(0)
        
        for (int i = 0 ; i<Steps; i++){
            P_e.at(i+1) = P_e.at(i) + (1- P_e.at(i)) * alphaX(i*StepSize) * (P_e.at(i)+P_h.at(i)-P_e.at(i)*P_h.at(i)) * StepSize;
            P_h.at(i+1) = P_h.at(i) - (1- P_h.at(i)) * betaX(i*StepSize) * (P_h.at(i)+P_e.at(i)-P_e.at(i)*P_h.at(i)) * StepSize;
        }

        if (std::abs(P_h.at(Steps)) < Accuracy){ //If P_h(Width) is close enough to 0, return the lists
            Pair.push_back(P_e);
            Pair.push_back(P_h);
            return Pair;
        }

        F = P_h.at(Steps); //Updates F to the current value of P_h(Width)

        //If we need to fix the guess
        if (firstIter){ //On the first iteration, do not use secant method
            firstIter = false;
            previousGuess = Guess;
            previousF = P_h.at(Steps);
            if (Guess < 0.9){
                Guess += 0.1;
            }
            else{
                Guess -= 0.1;
            }
        }
        //Use secant method otherwise
        else{
            tempGuess = Guess;
            Guess = Guess - (Guess - previousGuess)/(F - previousF) * F;
            previousGuess = tempGuess;
            previousF = F;
        }
    }
    std::cout << "No solution found" << std::endl;
    Pair.push_back(P_e);
    Pair.push_back(P_h);
    return Pair; //If the while loop ends, return the lists
}

int main(){
    std::vector<std::vector<double>> Pair = avaProb(5, 10000, 0.5, 0.00001);
    for (std::vector<double> P_x : Pair){
        std::cout << "List: " << std::endl;
        for (double Step : P_x){
            std::cout << Step << std::endl;
        }
    }
    return 0;
}