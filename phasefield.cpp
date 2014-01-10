//code to solve the equation d\phi/dt = v |grad \phi| in 2D

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
double dt=0.01;
double dx=0.1;
double v=1.0;

int totalshift=0;

void init() {

	phi = new double*[SIZE]; dphi = new double*[SIZE]; 

	for(int i=0;i<SIZE; i++) {
		phi[i]=new double[SIZE]; dphi[i]=new double[SIZE];
	}


	for(int i=0;i<SIZE;i++) {
		for(int j=0;j<SIZE;j++) {

		phi[i][j]=tanh( (j-SIZE/2)/10.0 - 5*sin(i*2*M_PI/((double)SIZE)))  ;

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
	
			double dphidx = (phi[up(i)][j]-phi[dwn(i)][j])/(2*dx); //use FTCS scheme
			double dphidy = (phi[i][j+1]-phi[i][j-1])/(2*dx);

			dphi[i][j] = -dt*v*sqrt(dphidx*dphidx+dphidy*dphidy);
	
	}}

	for(int i=0;i<SIZE;i++) {
		for(int j=1;j<SIZE-1;j++) {
	
			phi[i][j]+=dphi[i][j];		

	}}
	
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
    
    if(shift<0) {
        cout<<"Error: hitting bottom of sim box!"<<endl;
        exit(0);
    }
            
    for(int i=0;i<SIZE;i++) {
        for(int j=0;j<SIZE;j++) {
            
            int x=j+shift;
            if(x>=SIZE) phi[i][j]=1;
            else phi[i][j]=phi[i][j+shift];
    
        }
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

int main() {

	init();
	printgrid(1);

	for(int i=0;i<900;i++) {
		
		cout<<i<<endl;
		timestep();
		shiftEverythingDown();
	
	}

	printgrid(2);

	return 0;

}


