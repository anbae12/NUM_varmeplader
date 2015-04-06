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
#include "svd.h"

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


VecDoub b(2);
MatDoub A(10,2);








Doub Func(Doub xxx,Doub yyy){
    return 0.5*(pow(d, 2)/(2*pow(pow(d, 2)+pow(xxx-yyy,2),(3/2))));
}

Doub eksfunc(Doub x){
    return sqrt(x)*cos(x*x)*exp(-x);
}

struct Funcd {
    Doub operator() (const Doub x) {
        return eksfunc(x);
    }
};








int main() {
    x[2]=0; x[3]=0.25; x[1]=-0.25; x[4]=0.5; x[0]=-0.5;
    y[2]=0; y[3]=0.25; y[1]=-0.25; y[4]=0.5; y[0]=-0.5;
    
    
    // ANDERS
    
    interval[0] = -0.5*w;
    interval[1] = 0.5*w;
    int N = 4;
    Doub h =(interval[1]-interval[0])/N;
    
    
    A[0][0] = 0;
    A[1][0] = 0;
    A[2][0] = 0;
    A[3][0] = 0;
    A[4][0] = 0;
    A[5][0] = 0.5*Func(x[2], interval[0]);
    A[6][0] = Func(x[2], interval[0]+h);
    A[7][0] = Func(x[2], interval[0]+h+h);
    A[8][0] = Func(x[2], interval[0]+h+h+h);
    A[9][0] = 0.5*Func(x[2], interval[0]+h+h+h+h);
    
    A[0][1] = 0.5*Func(interval[0],y[2]);
    A[1][1] = Func(interval[0]+h,y[2]);
    A[2][1] = Func(interval[0]+h+h,y[2]);
    A[3][1] = Func(interval[0]+h+h+h,y[2]);
    A[4][1] = 0.5*Func(interval[0]+h+h+h+h,y[2]);
    A[5][1] = 0;
    A[6][1] = 0;
    A[7][1] = 0;
    A[8][1] = 0;
    A[9][1] = 0;
    
    
    std::cout<<A<<std::endl;
    
    b[0] = (-epsilon1*sigma*pow(T1,4))/((1-epsilon1)*h);
    b[1] = (-epsilon2*sigma*pow(T2,4))/((1-epsilon2)*h);
    
    std::cout<<"\n"<<b<<std::endl;
    
    
    
    
    
    SVD obj(A);
//    VecDoub W=obj.w;
//    MatDoub V=obj.v;
//    MatDoub U=obj.u;
//    
//    std::cout<<"W\n"<<W<<std::endl;
//    std::cout<<"V\n"<<V<<std::endl;
//    std::cout<<"U\n"<<U<<std::endl;
//    
    
    /* Estimate the parameters q = (x0, y0, a, b) and state your results.
     State also the residual error ∥Aq − z∥.
     */
//    VecDoub q(dataColumn);
    VecDoub z(10);
    obj.solve(z, b);

    std::cout<<"\n"<<z<<std::endl;
  
    
    // MIKKEL
    VecDoub I(5);
    I[0] = (h*(0.5*Func(interval[0], interval[0])*z[5])+0.5*Func(interval[0], interval[1])*z[9]+Func(interval[0], interval[0]+h)*z[6]+Func(interval[0], interval[0]+h+h)*z[7]+Func(interval[0], interval[0]+h+h+h))*z[8];
    I[1] = (h*(0.5*Func(interval[0]+h, interval[0])*z[5])+0.5*Func(interval[0]+h, interval[1])*z[9]+Func(interval[0]+h, interval[0]+h)*z[6]+Func(interval[0]+h, interval[0]+h+h)*z[7]+Func(interval[0]+h, interval[0]+h+h+h))*z[8];
    I[2] = (h*(0.5*Func(interval[0]+h+h, interval[0])*z[5])+0.5*Func(interval[0]+h+h, interval[1])*z[9]+Func(interval[0]+h+h, interval[0]+h)*z[6]+Func(interval[0]+h+h, interval[0]+h+h)*z[7]+Func(interval[0]+h+h, interval[0]+h+h+h))*z[8];
    I[3] = (h*(0.5*Func(interval[0]+h+h+h, interval[0])*z[5])+0.5*Func(interval[0]+h+h+h, interval[1])*z[9]+Func(interval[0]+h+h+h, interval[0]+h)*z[6]+Func(interval[0]+h+h+h, interval[0]+h+h)*z[7]+Func(interval[0]+h+h+h, interval[0]+h+h+h))*z[8];
    I[4] = (h*(0.5*Func(interval[0]+h+h+h+h, interval[0])*z[5])+0.5*Func(interval[0]+h+h+h+h, interval[1])*z[9]+Func(interval[0]+h+h+h+h, interval[0]+h)*z[6]+Func(interval[0]+h+h+h+h, interval[0]+h+h)*z[7]+Func(interval[0]+h+h+h+h, interval[0]+h+h+h))*z[8];
    
    
    Doub Q1 = h*(0.5*(z[5]-I[0])+0.5*(z[9]-I[4])+(z[6]-I[1])+(z[7]-I[2])+(z[8]-I[3]));
    
    std::cout<<Q1<<std::endl;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
//    
//    
//    // a)
//    
//    
//    //Doub precision = 1e-6 ; // Set precision
//    int maxiterations = 20;
//    Funcd fx;
//    //    Trapzd< Funcd > Tra(fx, interval[0], interval[1]);
//    Trapzd< Funcd > Tra(fx, 0, 1);
//    Tra.setww = 14;
//    
//    
//    std::cout<< std::setprecision(5)<<"k"<<std::setw(Tra.setww)<<"S(hk)"<<std::setw(Tra.setww)<<"ROE"<<std::setw(Tra.setww)<<"REE" <<std::endl;
//    for (int i=0; i<maxiterations; i++) {
//        Tra.next();
//    }
//    
//
//    
//    
    // b)
    // c)
    
    
    
    
    return 0;
}
