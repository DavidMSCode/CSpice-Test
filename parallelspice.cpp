/**
 * @file parallelspice.cpp
 * @author your name (you@domain.com)
 * @brief Testing CSpice thread collisions
 * @version 0.1
 * @date 2022-04-14
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "parallelspice.h"
#include "SpiceUsr.h"

#include <omp.h>
#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <cstdlib>

//using std namespace to avoid std:: before a lot of functions
using namespace std;
//using vec3dDouble as an alias for a 3d vector containing doubles

void main(int argc, char* argv[]){
    //default values
    int N = 1000;
    int num_runs = 2;
    int max_threads = 2;

    //Parse input to replace default values
    switch(argc) {
        case 2:
            //Just specify number of reads per run
            N=atoi(argv[1]);
            break;
        case 3:
            //Specify number of reads per run and num of parallel runs
            N = atoi(argv[1]);
            num_runs = atoi(argv[2]);
            break;
        case 4:
            //Specify number of reads per run and num of parallel runs
            N = atoi(argv[1]);
            num_runs = atoi(argv[2]);
            max_threads = atoi(argv[3]);
            break;
        default:
            //use default values
            break;
    }
    cout << "Running with " + to_string(N) + " SPICE reads per run, " + to_string(num_runs) + " parallel runs and " + to_string(max_threads) + " CPU threads\n";


    serialRun(int(N/2), num_runs);          //dividing N by 2 because the read function reads for moon and sun

    parallelRun(int(N/2), num_runs, max_threads);
}

void serialRun(int N, int num_runs){
    //runs in series a bunch of ephemeris data read runs in series to confirm that CSPICE is working correctly
    cout << "Starting serial run\n";
    //Load Kernels
    furnsh_c("de440.bsp");
    furnsh_c("naif0012.tls");
    for(int i=0;i<num_runs;i++){
        //an arbitrary time interval
        double t0 = 0.0;
        double tf = 1500;
        //read data
        vec3dDouble theData = readALotOfData(t0,tf, N);
    }
    //Unload Kernels
    kclear_c();
    cout << "Finished serial run\n";
}

void parallelRun(int N, int num_runs, int max_threads){
    //runs a bunch of parralel threads to encourage thread collisions 
    cout << "Starting parallel runs\n";
    omp_set_num_threads(max_threads);
    #pragma omp parallel for shared(N)
        for(int i=0; i<num_runs; i++){
            furnsh_c("de440.bsp");
            furnsh_c("naif0012.tls");
            //an arbitrary time interval
            double t0 = 0.0;
            double tf = 1500;
            //read data
            vec3dDouble theData = readALotOfData(t0,tf,N);
            kclear_c();
        } 
    cout << "Ending parallel runs\n";
}


//Function where majority of work occurs
vec3dDouble readALotOfData(double t0, double tf, int N){
    //too lazy to add real arrays so this function returns a vector of vectors of vectors of doubles 2xNx3
    //reads a lot of positions for the sun and moon relative to the earth to simulate CSPICE usage
    int threadnum = omp_get_thread_num();
    //cout << "Starting data read on thread " + to_string(threadnum) + "\n";
    
    //create vector of times
    vector<double> times = linspace(t0,tf,N);
    //init empty vector of vectors
    vector<vector<double> > sunList(N, vector<double>(3,0.0));
    vector<vector<double> > moonList(N, vector<double>(3,0.0));
    for(int i=0;i<N;i++){
        //read sun position at time i
        double sunpos[3];
        double lts;
        spkpos_c("SUN",times[i],"J2000","LT+S","EARTH",sunpos, &lts);
        vector<double> sunvec(begin(sunpos), end(sunpos));
        //store in vector
        sunList[i] = sunvec;
        //read moon position at time i
        double moonpos[3];
        double ltm;
        spkpos_c("MOOn",times[i],"J2000","LT+S","EARTH",moonpos, &ltm);
        vector<double> moonvec(begin(moonpos), end(moonpos));
        //store in vector
        moonList[i] = moonvec;
    }
    vec3dDouble output = {sunList,moonList};
    //cout << "Finished data read on thread " + to_string(threadnum) + "\n";
    return output;
}

//Utility functions
vector<double> linspace(double a, double b, int N){
    //returns a vector of N evenly spaced doubles over the interval [a,b]
    vector<double> out(N);
    double step = (b-a)/(N-1);
    out[0] = a;
    for (int i=1;i<N;i++){
        out[i] = out[i-1]+step;
    }
    return out;
}
