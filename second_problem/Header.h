#pragma once
#include <stdio.h>
#include <math.h>
#include <numeric>
#include <algorithm>
#include <cmath> 
#include <cstdlib>

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

ostream& operator<<(ostream& cout, const vector<double>& b);

class Matrix;
class Fredholm;

double l2_norm_square(vector<double>& x);

Fredholm Fredholm1();
Fredholm Fredholm2();
Fredholm Fredholm3();
Fredholm Fredholm5();
Fredholm Fredholm7();
Fredholm Fredholm8();
Fredholm Fredholm9();