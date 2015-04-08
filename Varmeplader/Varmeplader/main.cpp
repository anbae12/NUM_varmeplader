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
Doub sigma = 1.712e-09;
Doub d = 1.00;
Doub w = 1.00;
VecDoub x(5);
VecDoub y(5);
VecDoub interval(2);

Doub Func(Doub xxx,Doub yyy){
    return 0.5*pow(d, 2)/(pow(pow(d,2)+pow(xxx-yyy,2),1.5));
}

Doub make_h(int iteratations, Doub llim, Doub ulim){
    return (ulim-llim)/iteratations;
}

MatDoub make_A(VecDoub beta, int N, Doub l_bound, Doub h_bound, Doub h){
    int rows = 2*(N+1);
    int cols = 2*(N+1);
    
    
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
                A[j][i] = -beta[0]*0.5*Func(l_bound+(h*j),l_bound);
            }
            if (j <= rows/2 && i < cols-1 && i > cols/2) {
                A[j][i] = -beta[0]*Func(l_bound+(h*j),l_bound+(h*(i-(cols/2))));
            }
            if (j <= rows/2 && i == cols-1) {
                A[j][i] = -beta[0]*0.5*Func(l_bound+(h*j),h_bound);
            }
        }
    }
    
    for (int i = 0; i<cols/2; i++) {
        for (int j = rows/2; j<rows; j++) {
            if (j >= rows/2 && i == 0) {
                A[j][i] = -beta[1]*0.5*Func(l_bound,l_bound+(h*(j-(rows/2))));
            }
            if (j >= rows/2 && i < cols/2 && i > 0) {
                A[j][i] = -beta[1]*Func(l_bound+(h*i),l_bound+(h*(j-(rows/2))));
            }
            if (j >= rows/2 && i == cols/2-1) {
                A[j][i] = -beta[1]*0.5*Func(h_bound,l_bound+(h*(j-(rows/2))));
            }
        }
    }
    return A;
}

VecDoub make_beta( Doub h){
    VecDoub beta(2);
            beta[0] = (1-epsilon1)*h;
            beta[1] = (1-epsilon2)*h;
   return beta;
}

VecDoub make_b(int N_1, Doub h){
    int N = 2*(N_1+1);
    VecDoub b(N);
    for (int i = 0; i<N; i++) {
        if (i<N/2) {
            b[i] = epsilon1*sigma*pow(T1,4);
        }else{
            b[i] = epsilon2*sigma*pow(T2,4);
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
    int NN = N*2+2;
    
    
    // CREATE A and b
    Doub h = make_h(N, interval[0], interval[1]);
    std:cout<<"h: "<<h<<std::endl;
    MatDoub AA;
    VecDoub bb;
    AA = make_A(make_beta(h),N, interval[0], interval[1], h);
    bb = make_b(N,h);


    std::cout<<AA<<std::endl;
    std::cout<<bb<<std::endl;
    
    //
       SVD obj(AA);
       VecDoub z(NN);
       obj.solve(bb, z);
    //
        std::cout<<"\n"<<z<<std::endl;

    return 0;
}
