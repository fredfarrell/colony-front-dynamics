#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <vector>
#include"string.h"
using namespace std;

#define SIZE 1000
#define SGN(x) (double)((x > 0) - (x < 0))


vector<double> rx,ry,deltarx,deltary;
double dt=0.01;
double dx=0.1;
double v=1.0;
double ymin,ymax;



void init() {

	
	for(int i=0;i<SIZE;i++) {

		rx.push_back(i);
		ry.push_back(100*exp(-(i-500)*(i-500)/1000.0));
		//ry.push_back(0);

		deltarx.push_back(0);
		deltary.push_back(0);		

		ymin=0;
		ymax=100;

	}

}

void print(int j) {

	char str[30];

	sprintf(str,"out%i.dat",j);
	ofstream outp(str);

	for(int i=0;i<rx.size();i++) {
	
		outp<<i<<" "<<rx[i]<<" "<<ry[i]<<endl;

	}

	outp.close();

}



void timestep(int tm) {


	for(int i=0;i<rx.size();i++) {


		int iup = (i==rx.size()-1)?0:i+1;
		int idwn = (i==0)?rx.size()-1:i-1;

		double drx = rx[i]-rx[idwn]; //v is +ve for delta-rx
		if (fabs(drx)>SIZE*0.5) drx = drx - SIZE*SGN(drx); //minimum image for PBCs

		double dry = ry[iup]-ry[i]; //v is -ve for delta-ry
		double norm = sqrt(drx*drx+dry*dry);

		//cout<<i<<" "<<tm<<" "<<drx<<" "<<dry<<" "<<norm<<endl;

		drx/=norm; dry/=norm;

		deltarx[i] =  - dt*v*dry/dx; //update with upwind scheme
		deltary[i] =    dt*v*drx/dx;


	}


	for(int i=0;i<rx.size();i++) { //move the points

		int iup = (i==rx.size()-1)?0:i+1;
		int idwn = (i==0)?rx.size()-1:i-1; 

		rx[i]+=deltarx[i]; //update with upwind scheme
		ry[i]+=deltary[i];
	}


	for(int i=0;i<rx.size();i++) { //PBCs
		if (rx[i]>SIZE) rx[i]-=SIZE;
		if (rx[i]<0) rx[i]+=SIZE;
	}


}



int main() {

	init();
	print(1000);

	for (int i=0;i<50000;i++) {
		
		timestep(i);
		print(i);
		cout<<i<<" "<<rx.size()<<endl;
	}

	print(2000);


	return 0;

}
