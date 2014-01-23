#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <vector>
#include"string.h"
using namespace std;

int n[1000];
double psi[1000];
double dt=0.01;
double dx=1;
double D = 1;
    
double pmin=0.0001;

void print() {
    for(int i=0;i<1000;i++) cout<<i<<" "<<n[i]<<" "<<psi[i]<<endl;
}


int main() {

    double R = D*dt/(dx*dx);
    
    for(int i=0;i<1000;i++) {
        n[i]= (exp(-(i-500)*(i-500)/100.0))/pmin;
        psi[i]=0;
    }
    
    for(int t=0;t<10000;t++) {
        
        for(int i=0;i<1000;i++) {
        
        
        }
    
    }
    
    print();
    
   // for(int t=0;t<1000;t++)


}