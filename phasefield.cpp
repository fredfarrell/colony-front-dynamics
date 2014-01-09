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

void init() {

	phi = new double*[SIZE]; dphi = new double*[SIZE]; 

	for(int i=0;i<SIZE; i++) {
		phi[i]=new double[SIZE]; dphi[i]=new double[SIZE];
	}


	for(int i=0;i<SIZE;i++) {
		for(int j=0;j<SIZE;j++) {

		phi[i][j]=tanh((j-SIZE/2)/10.0);

	}}

}

int up(int x) {

	if(x==SIZE-1) return 0;
	else return x+1;
}

int dwn(int x) {

	if(x==0) return SIZE-1;
	else return x-1;

}


void timestep() {

	for(int i=1;i<SIZE-1;i++) {
		for(int j=1;j<SIZE-1;j++) {
	
			double dphidx = (phi[up(i)][j]-phi[dwn(i)][j])/(2*dx);
			double dphidy = (phi[i][up(j)]-phi[i][dwn(j)])/(2*dx);

			dphi[i][j] = -dt*v*sqrt(dphidx*dphidx+dphidy*dphidy);
	
	}}

	for(int i=1;i<SIZE-1;i++) {
		for(int j=1;j<SIZE-1;j++) {
	
			phi[i][j]+=dphi[i][j];		

	}}
	
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

	for(int i=0;i<100;i++) {
		
		timestep();
	
	}

	printgrid(2);

	return 0;

}


