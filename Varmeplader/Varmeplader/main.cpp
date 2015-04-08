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

Doub make_h(int iteratations, Doub llim, Doub ulim){
    return (ulim-llim)/iteratations;
}
MatDoub make_A(int xD,int yD, Doub XXX, Doub YYY, Doub h){
    int rows = xD;
    int cols = yD;
    int tempH1 = 0;
    int tempH2 = 0;

    MatDoub A(cols,rows);
    for (int i = 0; i<cols; i++) {
        for (int j = 0; j<rows; j++) {
            A[j][i]=0;
        }
    }
    
    int i = 0;
    for (int j = 0; j<rows; j++) {
        A[j][i]=1;
        i++;
    }
    
    for (int i = cols/2; i<cols; i++) {
        for (int j = 0; j<rows/2; j++) {
            if (j <= rows/2 && i == cols/2) {
                A[j][i] = 0.5*Func(XXX,YYY);
            }
            if (j <= rows/2 && i < cols-1 && i > cols/2) {
                ++tempH1;
                A[j][i] = Func(XXX,YYY+(h*tempH1));
            }
            if (j <= rows/2 && i == cols-1) {
                ++tempH1;
                A[j][i] = 0.5*Func(XXX,YYY+(h*tempH1));
            }
        }
    }
    
    for (int i = 0; i<cols/2; i++) {
        for (int j = rows/2; j<rows; j++) {
            if (j >= rows/2 && i == 0) {
                A[j][i] = 0.5*Func(XXX,YYY);
            }
            if (j >= rows/2 && i < cols/2 && i > 0) {
                ++tempH2;
                A[j][i] = Func(XXX,YYY+(h*tempH2));
            }
            if (j >= rows/2 && i == cols/2-1) {
                ++tempH2;
                A[j][i] = 0.5*Func(XXX,YYY+(h*tempH2));
            }
        }
    }
    return A;
}

VecDoub make_b(int N, Doub h){
    VecDoub b(N);
    for (int i = 0; i<N; i++) {
        if (i<N/2) {
            b[i] = (-epsilon1*sigma*pow(T1,4))/((1-epsilon1)*h);
        }else{
            b[i] = (-epsilon2*sigma*pow(T2,4))/((1-epsilon2)*h);
        }
    }
   return b;
}


// MIKKEL Trapetz
template<class T>

Doub trapez_integral(T &func, Doub a, Doub b,int N){
    Doub interval = b-a;
    Doub h = interval /N;
    Doub sum = 0;
    for(int i = 1; i<N;i++){
        sum+= func(a+i*h);
    }
    return h*(0.5*func(a)+0.5*func(b)+sum);
}


Doub f(Doub c){
    return c*c;
}

int main() {
    // TEST!!!
    std::cout<<"Result: "<<trapez_integral(f,0,10,10000)<<std::endl;
    
    // DATA POINTS
    x[2]=0; x[3]=0.25; x[1]=-0.25; x[4]=0.5; x[0]=-0.5;
    y[2]=0; y[3]=0.25; y[1]=-0.25; y[4]=0.5; y[0]=-0.5;
    interval[0] = -0.5*w;
    interval[1] = 0.5*w;
    int N = 4;
    
    
    // CREATE A and b
    Doub h = make_h(N, interval[0], interval[1]);
    MatDoub AA;
    VecDoub bb;
    AA = make_A(N,N, interval[0], interval[1], h);
    bb = make_b(N,h);


    std::cout<<AA<<std::endl;
    std::cout<<bb<<std::endl;
    
    //
    //
    //
    //
    //
    //    SVD obj(AA);
    //    VecDoub W=obj.w;
    //    MatDoub V=obj.v;
    //    MatDoub U=obj.u;
    //
    //    std::cout<<"W\n"<<W<<std::endl;
    //    std::cout<<"V\n"<<V<<std::endl;
    //    std::cout<<"U\n"<<U<<std::endl;
    //
    
    //    /* Estimate the parameters q = (x0, y0, a, b) and state your results.
    //     State also the residual error ∥Aq − z∥.
    //     */
    //    //    VecDoub q(dataColumn);
    //    VecDoub z(10);
    //    obj.solve(z, b);
    //
    //    std::cout<<"\n"<<z<<std::endl;
    //
    //
    //    // MIKKEL
    //    VecDoub I(5);
    //    I[0] = (h*(0.5*Func(interval[0], interval[0])*z[5])+0.5*Func(interval[0], interval[1])*z[9]+Func(interval[0], interval[0]+h)*z[6]+Func(interval[0], interval[0]+h+h)*z[7]+Func(interval[0], interval[0]+h+h+h))*z[8];
    //    I[1] = (h*(0.5*Func(interval[0]+h, interval[0])*z[5])+0.5*Func(interval[0]+h, interval[1])*z[9]+Func(interval[0]+h, interval[0]+h)*z[6]+Func(interval[0]+h, interval[0]+h+h)*z[7]+Func(interval[0]+h, interval[0]+h+h+h))*z[8];
    //    I[2] = (h*(0.5*Func(interval[0]+h+h, interval[0])*z[5])+0.5*Func(interval[0]+h+h, interval[1])*z[9]+Func(interval[0]+h+h, interval[0]+h)*z[6]+Func(interval[0]+h+h, interval[0]+h+h)*z[7]+Func(interval[0]+h+h, interval[0]+h+h+h))*z[8];
    //    I[3] = (h*(0.5*Func(interval[0]+h+h+h, interval[0])*z[5])+0.5*Func(interval[0]+h+h+h, interval[1])*z[9]+Func(interval[0]+h+h+h, interval[0]+h)*z[6]+Func(interval[0]+h+h+h, interval[0]+h+h)*z[7]+Func(interval[0]+h+h+h, interval[0]+h+h+h))*z[8];
    //    I[4] = (h*(0.5*Func(interval[0]+h+h+h+h, interval[0])*z[5])+0.5*Func(interval[0]+h+h+h+h, interval[1])*z[9]+Func(interval[0]+h+h+h+h, interval[0]+h)*z[6]+Func(interval[0]+h+h+h+h, interval[0]+h+h)*z[7]+Func(interval[0]+h+h+h+h, interval[0]+h+h+h))*z[8];
    //
    //
    //    Doub Q1 = h*(0.5*(z[5]-I[0])+0.5*(z[9]-I[4])+(z[6]-I[1])+(z[7]-I[2])+(z[8]-I[3]));
    //
    //    std::cout<<Q1<<std::endl;
    //
    //
    //
    //
    //
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
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
