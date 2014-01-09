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
double ds=1.0;
double v0=1.0;
double D=0.0;
double ymin,ymax;



double v(double x) {

	if(x>300 && x<700) return v0 + 0.05*exp( -(x-500)*(x-500)/200 );
	else return v0;

}

void init() {

	
	for(int i=0;i<SIZE;i++) {

		rx.push_back(i);
		//ry.push_back(100*exp(-(i-500)*(i-500)/1000.0));
		ry.push_back(0);

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


//functions to calculate i+1 and i-1 with PBCs
int up(int x) {

	if(x==rx.size()-1) return 0;
	else return x+1;
}

int dwn(int x) {

	if(x==0) return rx.size()-1;
	else return x-1;

}


void timestep(int tm) {


	for(int i=0;i<rx.size();i++) {


		//calculate positions of points 1 or 2 up or down from current cell using minimum image
		//sees if the difference between the current cell and (e.g.) i+1 is > SIZE/2, returns a shifted point if so
		double rxup,rxdwn,rx2up,rx2dwn;
		rxup=(rx[up(i)]-rx[i]<-SIZE*0.5)?rx[up(i)]+SIZE:rx[up(i)]; 
		rxdwn=(rx[dwn(i)]-rx[i]>SIZE*0.5)?rx[dwn(i)]-SIZE:rx[dwn(i)]; 
		//rx2up=(rx[up(up(i))]-rx[i]<-SIZE*0.5)?rx[up(up(i))]+SIZE:rx[up(up(i))];  not actually used
		//rx2dwn=(rx[dwn(dwn(i))]-rx[i]>SIZE*0.5)?rx[dwn(dwn(i))]-SIZE:rx[dwn(dwn(i))]; ditto

		double drxds = (rxup - rxdwn)/(2.0*ds); //centred difference 
		double dryds = (ry[up(i)] - ry[dwn(i)])/(2.0*ds); 

		double d2rxds2 = (rxup+rxdwn-2*rx[i])/(ds*ds);
		double d2ryds2 = (ry[up(i)]+ry[dwn(i)]-2*ry[i])/(ds*ds);

		double norm = sqrt(drxds*drxds+dryds*dryds);

		drxds/=norm; dryds/=norm;

		deltarx[i] =  - dt*v(i)*dryds + D*dt*d2rxds2;
		deltary[i] =    dt*v(i)*drxds + D*dt*d2ryds2;

	}


	for(int i=0;i<rx.size();i++) { //move the points

		rx[i]+=deltarx[i]; 
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
		if(i%1000==0) print(i);
	}


	print(2000);


	return 0;

}
