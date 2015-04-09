//
//  main.cpp
//  Varmeplader
//
//  Created by Anders Launer Baek on 24/03/15.
//  Copyright (c) 2015 Anders Launer Baek. All rights reserved.
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
Doub const1_1 = epsilon1*sigma*pow(T1,4);
Doub const1_2 = 1-epsilon1;
Doub const2_1 = epsilon2*sigma*pow(T2,4);
Doub const2_2 = 1-epsilon2;

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

VecDoub make_b(int N_1){
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

Doub Q(VecDoub list, Doub const_1, Doub const_2, Doub h){
    int N = list.size();
    Doub sum = 0;
    for(int i = 0; i<N; i++){
        if(i==0){
            sum += 0.5*h* (list[i]- (list[i]-const_1)/const_2);
        }
        else if(i==N-1){
            sum += 0.5*h* (list[i]- (list[i]-const_1)/const_2);
        }
        else{
            sum += h*(list[i]- (list[i]-const_1)/const_2);
        }
    }
    return sum;
}

int main() {
    // DATA POINTS
    interval[0] = -0.5*w;
    interval[1] = 0.5*w;
    
    Doub Q1, Q2, Q1_last, Q2_last, Q1_lastlast, Q2_lastlast, alpha_k1, alpha_k2, alpha_k, error1, error2;
    for(int N = 4; N<1040; N = N*2){
        int NN = N*2+2;
        
        // CREATE A and b
        Doub h = make_h(N, interval[0], interval[1]);
        MatDoub AA;
        VecDoub bb;
        AA = make_A(make_beta(h),N, interval[0], interval[1], h);
        bb = make_b(N);
        
        //Solve with SVD
        SVD obj(AA);
        VecDoub z(NN);
        obj.solve(bb, z);
        
        //Calculate Q1 and Q2
        VecDoub z_u(NN/2),z_v(NN/2);
        for(int i = 0; i<NN; i++){
            if(i<NN/2){
                z_u[i] = z[i];
            }
            else{
                z_v[i-NN/2] = z[i];
            }
        }
        Q1_lastlast = Q1_last;
        Q2_lastlast = Q2_last;
        Q1_last = Q1;
        Q2_last = Q2;
        
        //Calculate Q1 and Q2
        Q1 = Q(z_u,const1_1,const1_2,h);
        Q2 = Q(z_v,const2_1,const2_2,h);
        
        //Richardson alpha_k estimate
        alpha_k1 = (Q1_lastlast - Q1_last)/(Q1_last-Q1);
        alpha_k2 = (Q2_lastlast - Q2_last)/(Q2_last-Q2);
        
        
        // It is seen that both alpha_k converges to 4!
        alpha_k = 4;
        
        //Richardson error estimate
        error1 = (Q1 - Q1_last)/(alpha_k-1);
        error2 = (Q2 - Q2_last)/(alpha_k-1);
        
        int NNN=NN-1;
        int NNNN=NNN/2;
        
        std::cout<<scientific<<setprecision(16)
        << "N = "       << N                        << "\n"
        << "u(-1/2) = " << z[0]                     << "\n"
        << "u(-1/4) = " << z[NNNN-(((NNNN/2)-1))]   << "\n"
        << "u(0) = "    << z[((NNNN/2)/2)]          << "\n"
        << "u(1/4) = "  << z[NNNN-((NNNN/2)/2)]     << "\n"
        << "u(1/2) = "  << z[NNNN]                  << "\n"
        << "v(-1/2) ="  << z[NNNN+1]                << "\n"
        << "v(-1/4) ="  << z[NNN-(((NNN/2)-1))]     << "\n"
        << "v(0) = "    << z[NNN-((NNN/2)/2)]       << "\n"
        << "v(1/4) ="   << z[((NNN/2)/2)+NNNN]      << "\n"
        << "v(1/2) ="   << z[NNN]                   << "\n"
        << "Q1 = "      << Q1                       << "\n"
        << "Q2 = "      << Q2                       << "\n"
        << "Alp_k1 = "  << alpha_k1                 << "\n"
        << "Alp_k2 = "  << alpha_k2                 << "\n"
        << "Error1 = "  << error1                   << "\n"
        << "Error2 = "  << error2                   << "\n\n";
//        // Print to table
//        std::cout<<scientific<<setprecision(8)
//        <<N<<"&"<<z[0]<<"&"<<z[NNNN-(((NNNN/2)-1))]<<"&"<<z[((NNNN/2)/2)]<< "&"<< z[NNNN-((NNNN/2)/2)]<<"&"<< z[NNNN] <<"&"<<z[NNNN+1]<<"&"<<z[NNN-(((NNN/2)-1))]<<"&"<<z[NNN-((NNN/2)/2)]<<"&"<<z[((NNN/2)/2)+NNNN]<<"&"<<z[NNN]<< "\n"
//        << N<<"&"<<Q1<<"&"<<alpha_k1<<"&"<<error1<<"&"<<Q2<<"&"<<alpha_k2<<"&"<<error2<< "\n\n";
    }
    
    return 0;
}
