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


#define SIZE 500

double **phi,**dphi;
double dt=0.001;
double dx=0.1;
double v=1.0;
double D=1.0;
double alpha=10.0;

int totalshift=0;

double f(int x) {

    if(x>100 && x<400) return 5*exp( -(x-250)*(x-250)/500.0 );
	else return 0;

}

void init() {

	phi = new double*[SIZE]; dphi = new double*[SIZE]; 

	for(int i=0;i<SIZE; i++) {
		phi[i]=new double[SIZE]; dphi[i]=new double[SIZE];
	}


	for(int i=0;i<SIZE;i++) {
		for(int j=0;j<SIZE;j++) {

		//phi[i][j]=tanh( (j-SIZE/2)/10.0 - 5*sin(i*2*M_PI/((double)SIZE)))  ;
        phi[i][j]=0.5*tanh( f(i) + (SIZE/4-j)/10.0) + 0.5 ;

	}}

}

int up(int x) { //functions for PBCs

	if(x==SIZE-1) return 0;
	else return x+1;
}

int dwn(int x) {

	if(x==0) return SIZE-1;
	else return x-1;

}


void timestep() {

	for(int i=0;i<SIZE;i++) {
		for(int j=1;j<SIZE-1;j++) {

            double d2phidx2 = (phi[up(i)][j]+phi[dwn(i)][j]-2*phi[i][j])/(dx*dx);  
            double d2phidy2 = (phi[i][j+1]+phi[i][j-1]-2*phi[i][j])/(dx*dx);

			dphi[i][j] = dt*alpha*phi[i][j]*(1-phi[i][j]) + D*dt*(d2phidx2+d2phidy2);
	
	}}

	for(int i=0;i<SIZE;i++) {
		for(int j=1;j<SIZE-1;j++) {
	
			phi[i][j]+=dphi[i][j];		

	}}


	for(int i=0;i<SIZE;i++) {
	
            //zero flux BC
			phi[i][0]=phi[i][1];
            phi[i][SIZE-1]=phi[i][SIZE-2];		

	}
	
}

void printgrid(int t) {

	char str[30];
	sprintf(str,"out%i.dat",t);
	ofstream outp(str);

	for(int i=1;i<SIZE-1;i++) {
		for(int j=1;j<SIZE-1;j++) {

			outp<<i<<" "<<j<<" "<<phi[i][j]<<endl;

	}outp<<endl;}

}

void shiftEverythingDown() { //shift everything to keep the interface roughly centred in the simulation

    double tolerance = 0.01;
    int ymin=SIZE;
    int buffer=100; //number of lattice sites between the bottom of the contour and the bottom of the sim
    
    for(int i=0;i<SIZE;i++) {
        for(int j=0;j<SIZE;j++) {
            
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
	    for(int i=0;i<SIZE;i++) {
		for(int j=0;j<SIZE;j++) {
		    
		    int x=j+shift;
		    if(x>=SIZE) dphi[i][j]=0;
		    else dphi[i][j]=phi[i][x];
	    
		}
	    }
     }

    else if(shift<0) {

    	for(int i=0;i<SIZE;i++) {
		for(int j=0;j<SIZE;j++) {
		    
		    int x=j+shift;
		    if(x<0) dphi[i][j]=1;
		    else dphi[i][j]=phi[i][x];
	    
		}
	    }

    }

    else return;

    for(int i=0;i<SIZE;i++) {
	    for(int j=0;j<SIZE;j++) {
		    
                phi[i][j]=dphi[i][j]; //using dphi as a temp array to update grid simultaneously
	    
	    }
	}


    
}

int main() {

	init();
	printgrid(1);

	for(int i=0;i<21001;i++) {
		
		cout<<i<<endl;
		timestep();
		shiftEverythingDown();

        if(i%1000==0) printgrid(i);
	
	}


	return 0;

}


