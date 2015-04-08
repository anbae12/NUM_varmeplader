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


// Anders
Doub make_h(int iteratations, Doub llim, Doub ulim){
    return (ulim-llim)/iteratations;
}
MatDoub make_A(int x, Doub XXX, Doub YYY, Doub h){
      int rows = (iteratations+1)*2;
        int cols = 2;
        int tempH1 = 0;
        int tempH2 = 0;
    //
    //
    //    MatDoub A(rows,cols);
    //    for (int i=0; i<cols; i++) {
    //        for (int j=0; j<rows; j++) {
    //            if (i==0) {
    //                if (j<iteratations+1) {
    //                    A[j][i] = 0;
    //                }
    //                if (j==iteratations+1) {
    //
    //                    A[j][i] = 0.5*Func(XXX,YYY);
    //                }
    //                if (j>iteratations+1 && j<rows-1) {
    //                    ++tempH1;
    //                    A[j][i] = Func(XXX,YYY+(h*tempH1));
    //                }
    //                if (j==rows-1) {
    //                    ++tempH1;
    //                    A[j][i] = 0.5*Func(XXX,YYY+(h*tempH1));
    //                }
    //            }
    //            if (i==1) {
    //                if (j>iteratations) {
    //                    A[j][i] = 0;
    //                }
    //                if (j==0) {
    //                    A[j][i] = 0.5*Func(XXX, YYY);
    //                }
    //                if (j>0 && j<iteratations) {
    //                    ++tempH2;
    //                    A[j][i] = Func(XXX,YYY+(h*tempH2));
    //                }
    //                if (j==iteratations) {
    //                    ++tempH2;
    //                    A[j][i] = 0.5*Func(XXX,YYY+(h*tempH2));
    //                }
    //            }
    //        }
    //    }
    //    return A;
    int cols = 12;
    int rows = 12;
    
    
    MatDoub A(cols,rows);
    //    A[0][0] = 0; A[0][1] = 0; A[0][2] = 0; A[0][3] = 0; A[0][4] = 0; A[0][5] = 0; A[0][6] = 0; A[0][7] = 0; A[0][8] = 0; A[0][9] = 0;
    //    A[1][0] = 0; A[1][1] = 0; A[1][2] = 0; A[1][3] = 0; A[1][4] = 0; A[1][5] = 0; A[1][6] = 0; A[1][7] = 0; A[1][8] = 0; A[1][9] = 0;
    //    A[2][0] = 0; A[2][1] = 0; A[2][2] = 0; A[2][3] = 0; A[2][4] = 0; A[2][5] = 0; A[2][6] = 0; A[2][7] = 0; A[2][8] = 0; A[2][9] = 0;
    //    A[3][0] = 0; A[3][1] = 0; A[3][2] = 0; A[3][3] = 0; A[3][4] = 0; A[3][5] = 0; A[3][6] = 0; A[3][7] = 0; A[3][8] = 0; A[3][9] = 0;
    //    A[4][0] = 0; A[4][1] = 0; A[4][2] = 0; A[4][3] = 0; A[4][4] = 0; A[4][5] = 0; A[4][6] = 0; A[4][7] = 0; A[4][8] = 0; A[4][9] = 0;
    //    A[5][0] = 0; A[5][1] = 0; A[5][2] = 0; A[5][3] = 0; A[5][4] = 0; A[5][5] = 0; A[5][6] = 0; A[5][7] = 0; A[5][8] = 0; A[5][9] = 0;
    //    A[6][0] = 0; A[6][1] = 0; A[6][2] = 0; A[6][3] = 0; A[6][4] = 0; A[6][5] = 0; A[6][6] = 0; A[6][7] = 0; A[6][8] = 0; A[6][9] = 0;
    //    A[7][0] = 0; A[7][1] = 0; A[7][2] = 0; A[7][3] = 0; A[7][4] = 0; A[7][5] = 0; A[7][6] = 0; A[7][7] = 0; A[7][8] = 0; A[7][9] = 0;
    //    A[8][0] = 0; A[8][1] = 0; A[8][2] = 0; A[8][3] = 0; A[8][4] = 0; A[8][5] = 0; A[8][6] = 0; A[8][7] = 0; A[8][8] = 0; A[8][9] = 0;
    //    A[9][0] = 0; A[9][1] = 0; A[9][2] = 0; A[9][3] = 0; A[9][4] = 0; A[9][5] = 0; A[9][6] = 0; A[9][7] = 0; A[9][8] = 0; A[9][9] = 0;
    //
    
    
    int i = 0;
    for (int j = 0; j<rows; j++) {
        A[j][i]=1;
        i++;
    }
    
    
    for (int i = cols/2; i<cols; i++) {
        for (int j = 0; j<rows/2; j++) {
            A[j][i]=9;
            
        }
    }
    for (int i = 0; i<cols/2; i++) {
        for (int j = rows/2; j<rows; j++) {
            A[j][i]=8;
            
        }
    }
    
    
    return A;
}











int main() {
    // DATA POINTS
    x[2]=0; x[3]=0.25; x[1]=-0.25; x[4]=0.5; x[0]=-0.5;
    y[2]=0; y[3]=0.25; y[1]=-0.25; y[4]=0.5; y[0]=-0.5;
    interval[0] = -0.5*w;
    interval[1] = 0.5*w;
    int N = 4;
    
    
    // CREATE A and B
    Doub h = make_h(N, interval[0], interval[1]);
    MatDoub AA;
    AA = make_A(N, interval[0], interval[1], h);
    
    std::cout<<AA<<std::endl;
    
    //    b[0] = (-epsilon1*sigma*pow(T1,4))/((1-epsilon1)*h);
    //    b[1] = (-epsilon2*sigma*pow(T2,4))/((1-epsilon2)*h);
    //
    //    std::cout<<b<<std::endl;
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
