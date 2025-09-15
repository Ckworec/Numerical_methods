#pragma once
#include <stdio.h>
#include <math.h>
#include <numeric>
#include <algorithm>
#include <cmath> 

#include <thread>
#include <omp.h>
#include <typeinfo>
#include <chrono>

#include <iostream>
#include <fstream>

#include <vector>
#include <algorithm>
#include <functional>

using namespace std;

class Vector;
class Matrix;

double l1_norm(Vector& x);
double l2_norm(Vector& x);
double l2_norm_square(Vector& x);
double l_inf_norm(Vector& x);

double f1(Vector x);

double g(Vector x);