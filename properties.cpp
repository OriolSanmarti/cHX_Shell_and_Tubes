#include <iostream>
#include <string>
#include <cmath>
#include "properties.h" 

double calcRho(double T, std::string fluid){
        double value = 0;
        if(fluid == "solarsalt"){
            value = 2090 - 0.636*T;
        }else if(fluid == "yaramost"){
            value = 2085 - 0.74*T;
        }
        else if(fluid == "hitecxl"){
            value = 2240 - 0.8266*T;
        }
        else if(fluid == "ternary"){
            value = 2190 - 0.9*T;
        }
        else if(fluid == "quaternary"){
            value = 2170 - 0.688*T;
        }
        return value;
}

double calcCp(double T, std::string fluid){
        double value = 0;
        if(fluid == "solarsalt"){
            value = 1443 + 0.172*T;
        }else if(fluid == "yaramost"){
            value = 1549 - 0.15*T;
        }
        else if(fluid == "hitecxl"){
            value = 1536 - 0.2624*T;
        }
        else if(fluid == "ternary"){
            value = 1600 - 0.122*T;
        }
        else if(fluid == "quaternary"){
            value = 1660 - 0.719*T;
        }
        return value;
}

double calcLambda(double T, std::string fluid){
        double value = 0;
        if(fluid == "solarsalt"){
            value = 0.443 + 1.93E-4*T;
        }else if(fluid == "yaramost"){
            value = 0.697 - 4.61E-4*T;
        }
        else if(fluid == "hitecxl"){
            value = 0.519;
        }
        else if(fluid == "ternary"){
            value = 0.5;
        }
        else if(fluid == "quaternary"){
            value = 0.5;
        }
        return value;
}

double calcMu(double T, std::string fluid){
        double value = 0;
        if(fluid == "solarsalt"){
            value = 22.714E-3 - 0.120E-3*T + 2.281E-7*pow(T, 2.0) -1.474E-10*pow(T, 3.0);
        }else if(fluid == "yaramost"){
            value = (31.59 - 0.1948*T + 0.000425*pow(T,2.0) - 0.0000002122*pow(T,3.0) )/ 1000.0;
        }
        else if(fluid == "hitecxl"){
            value = 1372000*pow(T,-3.364);
        }
        else if(fluid == "ternary"){
            value = (15 - 0.0557*T + 0.0000575*pow(T,2.0)) / 1000.0;
        }
        else if(fluid == "quaternary"){
            value = (54 - 0.311*T + 0.000496*pow(T,2.0)) / 1000.0;
        }
        return value;
}


double calcH(double T, std::string fluid){
        double value = 0;
        double Tref = 0;
        if(fluid == "solarsalt"){
            double cpmid = (T*calcCp(T, fluid) - Tref*calcCp(Tref, fluid)) / (T-Tref);
            value =  cpmid * (T-Tref);
        }else if(fluid == "yaramost"){
            value = 0;
        }
        else if(fluid == "hitecxl"){
            value = 0;
        }
        else if(fluid == "ternary"){
            value = 0;
        }
        else if(fluid == "quaternary"){
            value = 0;
        }
        return value;
}

double calcS(double T, std::string fluid){
        double value = 0;
        double Tref = 0;
        if(fluid == "solarsalt"){
            double cpmid = (T*calcCp(T, fluid) - Tref*calcCp(Tref, fluid)) / (T-Tref);
            value = cpmid * log(T-Tref);
        }else if(fluid == "yaramost"){
            value = 0;
        }
        else if(fluid == "hitecxl"){
            value = 0;
        }
        else if(fluid == "ternary"){
            value = 0;
        }
        else if(fluid == "quaternary"){
            value = 0;
        }
        return value;
}

