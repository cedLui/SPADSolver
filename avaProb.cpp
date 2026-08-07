#include <iostream>
#include <cmath>
#include <vector> 

#include "avaProb.h"


std::vector<double> guess100(double Width, int Steps, double Bias, double td, double rho, std::vector<double> AlFracProf){ //This function was a debugging function
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
    double Guess = 1; //Guess for P_h(0), intialized at 1
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
        
        for (int i = 0 ; i<Steps; i++){ //Use FDM to solve equations
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


        if (firstIter){ //We need to check at least two guesses to even have a chance at starting binary search
            firstIter = false;
            prevGuess = Guess;
            Guess = Guess - 0.01; //Move on to the next guess
        } else{
            if (!binarySearch){ //If the binary search hasn't been activated, check to see if it should be
                if (std::abs(PhW + prevPhW) != std::abs(PhW) + std::abs(prevPhW)){ //If there was a sign change, we know the real P_h_0 is between the current and previous guesses
                    binarySearch = true;
                    if (PhW < 0){ //Depending on whether the change in P_h(W) was from - to + or vice versa, the binary search bounds would change
                        left = Guess;
                        right = prevGuess;
                    } else{
                        right = Guess;
                        left = prevGuess;
                    }
                    Guess = (prevGuess + Guess)/2.0;
                } else{
                    prevGuess = Guess;
                    Guess = Guess - 0.01; //If binary search wasn't activated, continue sweeping
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

/*int main(){
    std::vector<double> AlFracProf;
    for (int x = 0; x < 10000; x++){
        AlFracProf.push_back(0);
    }
    
    std::vector<std::vector<double>> Pair = avaProb(0.00005, 10000, 0.00000001, 150, 0.00001, 0, AlFracProf);
    for (int i = 0; i <= 10000; i = i + 200){
        std::cout << Pair[0][i] << std::endl;
    }
}*/

