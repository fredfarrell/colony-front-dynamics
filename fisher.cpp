//2D Fisher equation solver

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <vector>
#include"string.h"
using namespace std;


#define SIZEX 500
#define SIZEY 500

double **phi,**dphi, **p, **dp; //phi is the total density, p the proportion of fitter mutants (phi=p+q, q not considered explicitly)
double dt=0.001;
double dx=0.1;
double v=1.0;
double D=1.0;
double alpha_1=12.0;
double alpha_2=10.0;

int totalshift=0;

double f(int x) {

    if(x>100 && x<400) return 5*exp( -(x-SIZEX/2)*(x-SIZEX/2)/500.0 );
	else return 0;

}

void init() {

	phi = new double*[SIZEX]; dphi = new double*[SIZEX]; p = new double*[SIZEX]; dp = new double*[SIZEX];

	for(int i=0;i<SIZEX; i++) {
		phi[i]=new double[SIZEY]; dphi[i]=new double[SIZEY]; p[i]=new double[SIZEY]; dp[i]=new double[SIZEY];
	}


	for(int i=0;i<SIZEX;i++) {
		for(int j=0;j<SIZEY;j++) {

		//phi[i][j]=tanh( (j-SIZE/2)/10.0 - 5*sin(i*2*M_PI/((double)SIZE)))  ;
        phi[i][j]=0.5*tanh( (SIZEY/4-j)/10.0) + 0.5 ;
        p[i][j]=(0.5*tanh( (SIZEY/4-j)/10.0) + 0.5) * exp (-(i-SIZEX/2)*(i-SIZEX/2)/100.0) ;

	}}

}

int up(int x) { //functions for PBCs

	if(x==SIZEX-1) return 0;
	else return x+1;
}

int dwn(int x) {

	if(x==0) return SIZEX-1;
	else return x-1;

}


void timestep() {



	for(int i=0;i<SIZEX;i++) {
		for(int j=1;j<SIZEY-1;j++) {

            double d2phidx2 = (phi[up(i)][j]+phi[dwn(i)][j]-2*phi[i][j])/(dx*dx);  
            double d2phidy2 = (phi[i][j+1]+phi[i][j-1]-2*phi[i][j])/(dx*dx);
            double d2pdx2 = (p[up(i)][j]+p[dwn(i)][j]-2*p[i][j])/(dx*dx);  
            double d2pdy2 = (p[i][j+1]+p[i][j-1]-2*p[i][j])/(dx*dx);

			dphi[i][j] = dt*alpha_2*phi[i][j]*(1-phi[i][j]) + dt*(alpha_1-alpha_2)*p[i][j]*(1-phi[i][j]) + D*dt*(d2phidx2+d2phidy2);

            dp[i][j] = dt*alpha_1*p[i][j]*(1-phi[i][j]) + D*dt*(d2pdx2+d2pdy2);

	}}

	for(int i=0;i<SIZEX;i++) {
		for(int j=1;j<SIZEY-1;j++) {
	
			phi[i][j]+=dphi[i][j];	
            p[i][j]+=dp[i][j];	

	}}


	for(int i=0;i<SIZEX;i++) {
	
            //zero flux BC
			phi[i][0]=phi[i][1];
            phi[i][SIZEY-1]=phi[i][SIZEY-2];		
			p[i][0]=p[i][1];
            p[i][SIZEY-1]=p[i][SIZEY-2];	

	}

}

void printgrid(int t) {

	char str[30];
	sprintf(str,"out%i.dat",t);
	ofstream outp(str);

	for(int i=1;i<SIZEX-1;i++) {
		for(int j=1;j<SIZEY-1;j++) {

			outp<<i<<" "<<j<<" "<<phi[i][j]<<" "<<p[i][j]<<endl;

	}outp<<endl;}

}

void shiftEverythingDown() { //shift everything to keep the interface roughly centred in the simulation

    double tolerance = 0.01;
    int ymin=SIZEY;
    int buffer=100; //number of lattice sites between the bottom of the contour and the bottom of the sim
    
    for(int i=0;i<SIZEX;i++) {
        for(int j=0;j<SIZEY;j++) {
            
            if(fabs(phi[i][j])<tolerance && j<ymin) ymin=j;
            
        }
    }
    
    int shift = ymin-buffer;
    
    /*if(shift<0) {
	printgrid(3);
        cout<<"Error: hitting bottom of sim box!" <<ymin<<endl;
        exit(0);
    }*/

    if(shift>0) {            
	    for(int i=0;i<SIZEX;i++) {
		for(int j=0;j<SIZEY;j++) {
		    
		    int x=j+shift;
		    if(x>=SIZEY) {dphi[i][j]=0; dp[i][j]=0;}
		    else {dphi[i][j]=phi[i][x]; dp[i][j]=p[i][x];}
	    
		}
	    }
     }

    else if(shift<0) {

    	for(int i=0;i<SIZEX;i++) {
		for(int j=0;j<SIZEY;j++) {
		    
		    int x=j+shift;
		    if(x<0) {dphi[i][j]=1; dp[i][j]=1;}
		    else {dphi[i][j]=phi[i][x]; dp[i][j]=p[i][x];}
	    
		}
	    }

    }

    else return;

    for(int i=0;i<SIZEX;i++) {
	    for(int j=0;j<SIZEY;j++) {
		    
                phi[i][j]=dphi[i][j]; //using dphi as a temp array to update grid simultaneously
                p[i][j]=dp[i][j];
	    
	    }
	}


    
}

int main() {

	init();
	printgrid(1);

	for(int i=0;i<20001;i++) {
		
		cout<<i<<endl;
		timestep();
		shiftEverythingDown();

        if(i%1000==0) printgrid(i);
	
	}


	return 0;

}


