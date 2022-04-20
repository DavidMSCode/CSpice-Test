/**
 * @file parallelspice.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2022-04-14
 * 
 * @copyright Copyright (c) 2022
 * 
 */


#ifndef __PSPICE__
#define __PSPICE__
#include <vector>
using vec3dDouble = std::vector<std::vector<std::vector<double> > >;
/**
 * @brief Runs some spice functions on a single thread
 * 
 */
void serialRun(int N);

vec3dDouble readALotOfData(double t0, double tf, int N);

std::vector<double> linspace(double a, double b, int N);

#endif