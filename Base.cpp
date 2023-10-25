#include <iostream>
#include <vector>
#include <cmath>
#include "properties.h"
#include "cHXshelltube.h"
#include <fstream>

using namespace std;

int main(){
    cout << "La lambda del fluid(quaternary)a 250 ºC es: " << calcLambda(250, "quaternary") << endl;
    cout << "La mu del fuid(quaternary) a 250 ºC es: "<< calcMu(250, "quaternary") << endl;
    cout << "La mu del fuid(quaternary) a 350 ºC es: "<< calcMu(350, "quaternary") << endl;
    cout << "La mu del fuid(quaternary) a 450 ºC es: "<< calcMu(450, "quaternary") << endl;
    cout << "La cp min es: " << calcCp(180, "quaternary") << endl;
    cout << "La cp max es: " << calcCp(350, "quaternary") << endl;
    double T1in = 180.0;
    double T1o = 350.0;
        
    T1in = 20;
    T1o = 50.0;
    double T2in = 350;
    double T2o = T2in;
    string fluid1 = "air";
    string fluid2 = "quaternary";
    cHXshelltube(T1in, T1o, T2in, T2o, fluid1, fluid2);
}
