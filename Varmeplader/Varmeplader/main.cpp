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
    /*
    int N = 4;
    int NN = N*2+2;

    // CREATE A and b
    Doub h = make_h(N, interval[0], interval[1]);
    MatDoub AA;
    VecDoub bb;
    AA = make_A(make_beta(h),N, interval[0], interval[1], h);
    bb = make_b(N);


    std::cout<<AA<<std::endl;
    std::cout<<bb<<std::endl;
    SVD obj(AA);
    VecDoub z(NN);
    obj.solve(bb, z);
    std::cout<<"\n"<<z<<std::endl;

    //Calculate Q1 and Q2
    VecDoub z_u(NN/2),z_v(NN/2);
    Doub Q1, Q2;
    for(int i = 0; i<NN; i++){
        if(i<NN/2){
            z_u[i] = z[i];
        }
        else{
            z_v[i-NN/2] = z[i];
        }
    }
    z_v.print();

    Q1 = Q(z_u,const1_1,const1_2,h);
    Q2 = Q(z_v,const2_1,const2_2,h);
    std::cout<<"Q1: "<<Q1<<endl;
    std::cout<<"Q2:"<<Q2<<endl;
    */


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

        //std::cout<<"Alpha_k1 = "<<alpha_k1<<std::endl;
        //std::cout<<"Alpha_k2 = "<<alpha_k2<<std::endl;

        // It is seen that both alpha_k converges to 4!
        alpha_k = 4;

        //Richardson error estimate
        error1 = (Q1 - Q1_last)/(alpha_k-1);
        error2 = (Q2 - Q2_last)/(alpha_k-1);

        //OUTPUT
        std::cout<<scientific<<setprecision(16)
                << "N = "       << N                << setw(15)
                << "x(-1/4) = " << z[NN/8]          << setw(15)
                << "x(0) = "    << z[NN/4]          << setw(15)
                << "x(1/4) = "  << z[NN/4+NN/8]     << setw(15)
                << "y(-1/4) ="  << z[NN/2+NN/8]     << setw(15)
                << "y(0) = "    << z[NN/2+NN/4]     << setw(15)
                << "y(1/4) ="   << z[NN/2+NN/4+NN/8]<< "\n"
                << "Q1 = "      << Q1               << setw(15)
                << "Q2 = "      << Q2               << setw(15)
                << "Alp_k1 = "  << alpha_k1         << setw(15)
                << "Alp_k2 = "  << alpha_k2         << setw(15)
                << "Error1 = "  << error1           << setw(15)
                << "Error2 = "  << error2
                << std::endl<<std::endl<<std::endl;






    }





    //Part c


    return 0;
}
