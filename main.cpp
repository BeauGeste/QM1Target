//
//  main.cpp
//  1-Target MT

//  Small updates by Nicholas Harmon on 6/16/20.
//  Small updates by Richard Gerst on 8/--/19.
//  Created by Stephen McMillan on 5/23/17.
//  Copyright © 2017 Stephen McMillan. All rights reserved.

//  Yu averaged rates are hard-coded in (if uncommented) which will give Yu results. When actual hopping matrix elements are used, Yu is not correct.

#include <iostream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <gsl/gsl_rng.h>	// GSL random number generators
#include <gsl/gsl_randist.h>	// GSL random distributions
#include <fstream>
#include <sstream>
//#include <bsumtree.h> // not being use in code
//#include "/Users/smcmillan/Documents/Physics/Flatté Research/Code/MH Spin Relaxation/bsumtree.h"

using namespace std;

#define Pi 3.1415926535897932384626433832795028
#define HBar (6.582119514 * pow(10,-7))  //units of eV*ns
#define GammaSqr (1.07 * pow(10,-3))  //Alq3 electron polaron value (Yu) (1.07 * pow(10,-3))
#define Lambda 	 sqrt(1.07 * pow(10,-3)/2) //Lambda = Sqrt(2*GammaSqr) sqrt(1.07 * pow(10,-3)/2)
#define V	0.001 //0.001 // units of eV (Yu) (0.1)

gsl_rng * rng_ptr;
//pointer to random number generator(rng)

double          random_exp(double a)
{
    return gsl_ran_exponential(rng_ptr, a);
}

double
random_gaussian(double sigma)
{
    return gsl_ran_gaussian(rng_ptr, sigma);
}

double
random_uniform()
{
    return gsl_ran_flat(rng_ptr, 0, 1);
}


//All V matrix elements are calculated with the following approx: vijx=vijy=vijz=V0

double
ModSqrVpp(double itheta, double iphi, double ftheta, double fphi)
{
    
    double MSvpp;
    
    //PRB
    // MSvpp = pow(V,2.0)*(pow(cos(itheta)*cos(ftheta)+cos(iphi-fphi)*sin(itheta)*sin(ftheta),2.0)+pow(Lambda,  2.0)*(1+sin(2*itheta))*pow(sin(ftheta)*sin(iphi-fphi),2.0));
    
    //MINE/SHE
    MSvpp = pow(V,2.0)*(pow(cos(itheta)*cos(ftheta),2.0)+2*cos(itheta)*cos(ftheta)*sin(itheta)*sin(ftheta)*cos(iphi-fphi)+pow(sin(itheta)*sin(ftheta)*cos(iphi-fphi),2.0)+pow(2*Lambda*sin(itheta)*sin(ftheta)*sin(iphi-fphi),2.0));
    
    return MSvpp;
}

double
ModSqrVmm(double itheta, double iphi, double ftheta, double fphi)
{
    double MSvmm;
    
    //PRB
    //MSvmm = pow(V,2.0)*(pow(cos(itheta)*cos(ftheta)+cos(iphi-fphi)*sin(itheta)*sin(ftheta),2.0)+pow(Lambda,2.0)*(1+sin(2*itheta))*pow(sin(ftheta)*sin(iphi-fphi),2.0));
    
    //MINE/SHE
    MSvmm = pow(V,2.0)*(pow(cos(itheta)*cos(ftheta),2.0)+2*cos(itheta)*cos(ftheta)*sin(itheta)*sin(ftheta)*cos(iphi-fphi)+pow(sin(itheta)*sin(ftheta)*cos(iphi-fphi),2.0)+pow(2*Lambda*sin(itheta)*sin(ftheta)*sin(iphi-fphi),2.0));
    
    return MSvmm;
}

double
ModSqrVpm(double itheta, double iphi, double ftheta, double fphi)
{
    double MSvpm;
    
    
    // MSvpm = 4*pow(V*Lambda*cos(ftheta)*sin(itheta)*sin(iphi-fphi),2.0);
    
    //MINE/SHE
    MSvpm = pow(V*Lambda*2,2.0)*(pow(cos(ftheta)*sin(itheta),2.0)-2*cos(ftheta)*sin(itheta)*cos(itheta)*sin(ftheta)*cos(iphi-fphi)+pow(cos(itheta)*sin(ftheta),2.0));
    
    return MSvpm;
}

double
ModSqrVmp(double itheta, double iphi, double ftheta, double fphi)
{
    
    double MSvmp;
    
    //MSvmp = 4*pow(V*Lambda*cos(ftheta)*sin(itheta)*sin(iphi-fphi),2.0);
    
    //MINE/SHE
    MSvmp = pow(V*Lambda*2,2.0)*(pow(cos(ftheta)*sin(itheta),2.0)-2*cos(ftheta)*sin(itheta)*cos(itheta)*sin(ftheta)*cos(iphi-fphi)+pow(cos(itheta)*sin(ftheta),2.0));
    
    
    return MSvmp;
}

double
rate_for_E(double E, double k0)
{
    if (E < 0) {
        return k0;
    } else {
        return k0 * exp(- E);
    }
}


int
main(int argc, char *argv[])
{
    clock_t t1,t2;
    t1=clock();
    
    
    rng_ptr = gsl_rng_alloc(gsl_rng_taus2);
    //allocate the rng
    gsl_rng_set(rng_ptr, 1099);
    //seed the rng
    
    //string input;
    double SIGMA;
    int NDIS;
    int max_hops;
    int count = 0;
    int totalhop = 0;
    int flip_num = 0;
    int con_num = 0;
    
    
    cout << "Enter the disorder (in units of T), sigma: ";
    cin >> SIGMA;
    cout << "Enter the number of disorder configurations to average over: ";
    cin >> NDIS;
    cout << "Enter the maximum number of hops: ";
    cin >> max_hops;
    
    
    std::ostringstream strs2;
    strs2 << SIGMA;
    std::string input2 = strs2.str();
    std::ostringstream strs3;
    strs3 << NDIS;
    std::string input3 = strs3.str();
    std::ostringstream strs4;
    strs4 << max_hops;
    std::string input4 = strs4.str();
    
    std::ofstream out_file;
    out_file.open(("1MT_s"+input2+".dat").c_str());
    cout << "The output file is " << "1MT_s"+input2+".dat" << endl;
    
    std::ofstream angleWTout_file;
    angleWTout_file.open(("1MT_s"+input2+".dat").c_str());
    cout << "The total rate output file is " << "out1MT_s"+input2+".dat" << endl;
    
    int HOPS;
    double time, tOld, tNew;
//    double time_max = 0.015; //
//    double dt = 0.0001; //
    double time_max = 100.0;
    double dt = 5;
    
    for (int id = 0; id < NDIS; id++)
    {
        HOPS = 0;
        
        int state = 0;
        double sim_time = 0;
        double phi_i = random_uniform()*2*Pi;
        double theta_i = acos(2*random_uniform()-1);
        double Sz = 1.; //Start in spin-up state
        
        while (HOPS < max_hops + 1)
        {
            double rates[2];
            double Vrate[2];
            double total_rate = 0;
           // double vs = 0.00907915;
            //double vc = 30.18649;

            double phi_f = random_uniform()*2*Pi;
            double theta_f = acos(2*random_uniform()-1);//random_uniform()*Pi;
      //      double phi_f2 = random_uniform()*2*Pi;
      //      double theta_f2 = acos(2*random_uniform()-1);//random_uniform()*Pi;
      //      double phi_i2 = random_uniform()*2*Pi;
      //      double theta_i2 = acos(2*random_uniform()-1);
            double E = (SIGMA != 0) ? -random_exp(SIGMA) : 0; //as SIGMA --> 0, random_exp --> dirac delta distribution; i.e., Prob(E=0) = 1. R Gerst, 15 Aug 2019
            double MTrate_for_E = exp(E);
            
            if (state == 0)
            {
                Vrate[0] = ModSqrVmp(theta_i, phi_i, theta_f, phi_f);
                Vrate[1] = ModSqrVpp(theta_i, phi_i, theta_f, phi_f);
                
           //     Vrate[0] = vs;
           //     Vrate[1] = vc;
//                rates[0] = 0.00907915;//2*Pi/HBar * Vrate[0];//= 2*Pi/HBar * Vrate[0];0.00907915;//supposedly coming from Yu's averaging
//                rates[1] = 3.18649;//2*Pi/HBar * Vrate[1];//3.18649; // these hard code values give relaxation rate of Yu using V = 0.001
                rates[0] = 2*Pi/HBar * Vrate[0] * MTrate_for_E;
                rates[1] = 2*Pi/HBar * Vrate[1] * MTrate_for_E;
               // std::cout << "old" << rates[0]  << std::endl;

                
              //  rates[0] = vs* MTrate_for_E;
              //  rates[1] = vc * MTrate_for_E;
              //  std::cout << "new" << rates[0]  << std::endl;
            }
            else if (state == 1)
            {
                Vrate[0] = ModSqrVpm(theta_i, phi_i, theta_f, phi_f);
                Vrate[1] = ModSqrVmm(theta_i, phi_i, theta_f, phi_f);


//               rates[0] = 0.00907915;//2*Pi/HBar * Vrate[0]; //= 0.00907915;
//                rates[1] = 3.18649;// = 2*Pi/HBar * Vrate[1];
                rates[0] = 2*Pi/HBar * Vrate[0] * MTrate_for_E;
                rates[1] = 2*Pi/HBar * Vrate[1] * MTrate_for_E;
              //  std::cout << "old" << rates[0]  << std::endl;
                
              //  Vrate[0] = vs;
              //  Vrate[1] = vc;
                
             //   rates[0] = vs* MTrate_for_E;
             //   rates[1] = vc * MTrate_for_E;
         //       std::cout << "new" << rates[0]  << std::endl;
                
            }
            
            total_rate = rates[0]+rates[1];  // NICK HAS THIS INSIDE STATE if STATEMENTS
            
            //choose hop:
            double          rand_rate = random_uniform() * total_rate;
            int             ihop = 0;
            
            
            
            if (rates[0] > 0.0) //Adding this 'if' was done on 21Sept2016
            {
                if (rand_rate <= rates[0])
                {
                    ihop = 0; //flip
                    flip_num++;
                }
                
                else
                {
                    ihop = 1; //no flip
                    con_num++;
                }
            }
            else
            {
                ihop = 1; //no flip
                con_num++;
            }
            totalhop++;
            
            time = random_exp(1.0 / total_rate );  //was 1/total_rate
            tNew = sim_time + time;
            tOld = sim_time;
            sim_time = tNew;
//            time_max = 200.0;//0.015; //100.;.05
//            dt = 5.0;//0.0005; //3.;.001
            
            //  if (ihop == 1)
            // {
           // angleWTout_file << theta_f << "\t" << rates[1] << std::endl;
            // }
            
            

            for (double t = 0.; t <= time_max; t+= dt) {
                
                if (t >= tOld && t < tNew) {
                    
                    out_file << t << "\t" << Sz << std::endl;

                }
            }
            
            if (state == 0){
                if (ihop == 1){
                    state = 0;
                }
                
                else if (ihop == 0){
                    state = 1;
                }
            }
            
            else if (state == 1){
                if (ihop == 1){
                    state = 1;
                }
                
                else if (ihop == 0){
                    state = 0;
                }
            }
            
            //Calculate Spin along z (Yu PRB 2.10 and 2.11)
            
            if (state == 0){
                
                Sz = 1 - GammaSqr * pow(cos(theta_f),2.0);
            }
            
            else if (state == 1){
                
                Sz = -(1 - GammaSqr * pow(cos(theta_f),2.0));
            }
            
            HOPS++;
            
            phi_i = phi_f;
            theta_i = theta_f;
            
            
            //For Uncorrelated Hopping
            /* if (ihop == 0)
             {
             theta_i = theta_f;
             phi_i = phi_f;
             }
             else
             {
             theta_i = theta_f2;
             phi_i = phi_f2;
             }
             */
        }
        
        
        
    }
    
    out_file.close();
    angleWTout_file.close();
    
    std::cout << "disorder configs = " << "\t" << NDIS << std::endl;
    std::cout << "sigma = " << "\t" << SIGMA << std::endl;
    std::cout << "Times flip > con = " << "\t" << count <<"/"<<totalhop <<" = "<< count/totalhop * 100<<"%"<< std::endl;
    std::cout << "flip/con = " << "\t" << flip_num <<"/"<< con_num <<std::endl;
    
#ifdef WIN32
    char            c;
    std::cin >> c;
#endif
    
    gsl_rng_free(rng_ptr);
    //free the random number generator
    
    t2=clock();
    float diff ((float)t2-(float)t1);
    float seconds = diff / CLOCKS_PER_SEC;
    std::cout<< "timer:\t" << seconds<< " seconds"<<std::endl;
    
    return 0;
    
}



