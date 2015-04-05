//
//  main.cpp
//  Varmeplader
//
//  Created by Anders Launer Bæk on 24/03/15.
//  Copyright (c) 2015 Anders Launer Bæk. All rights reserved.
//

#include <iostream>
#include <iomanip>      // std::setw
#include "nr3.h"
#include "quadrature.h"

// DATA
Doub T1 = 1000;
Doub T2 = 500;
Doub epsilon1 = 0.80;
Doub epsilon2 = 0.60;
Doub sigma = 1.712*10e-09;
Doub d = 1.00;
Doub w = 1.00;
VecDoub x(5);
VecDoub y(5);
VecDoub interval(2);










Doub eksfunc(Doub x){
    return sqrt(x)*cos(x*x)*exp(-x);
}

struct Funcd {
    Doub operator() (const Doub x) {
        return eksfunc(x);
    }
};








int main() {
    x[0]=0; x[1]=0.25; x[2]=-0.25; x[3]=0.5; x[4]=-0.5;
    y[0]=0; y[1]=0.25; y[2]=-0.25; y[3]=0.5; y[4]=-0.5;
    interval[0] = -0.5*w;
    interval[1] = 0.5*w;
    
    // a)
    
    
    //Doub precision = 1e-6 ; // Set precision
    int maxiterations = 20;
    Funcd fx;
//    Trapzd< Funcd > Tra(fx, interval[0], interval[1]);
    Trapzd< Funcd > Tra(fx, 0, 1);
    Tra.setww = 14;
    
    
    std::cout<< std::setprecision(5)<<"k"<<std::setw(Tra.setww)<<"S(hk)"<<std::setw(Tra.setww)<<"ROE"<<std::setw(Tra.setww)<<"REE" <<std::endl;
    for (int i=0; i<maxiterations; i++) {
        Tra.next();
    }
    
    
    
    
    // b)
    // c)
    
    
    
    
    return 0;
}
