#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <functional>
#include <limits>
#include <tuple>
#include <algorithm>
#include <fstream>
#include <sstream>

// Определим тип для точки в n-мерном пространстве
using Point = std::vector<double>;
// Определим тип для функции f(x)
using Function = std::function<double(const Point&)>;

#define eps 1.e-4
#define tao 0.1 * sqrt(eps)
#define r (3.0 - sqrt(5.0)) / 2.0
#define max_iter 10000

double norm(const Point& v);
double dot(const Point& a, const Point& b);
Point operator-(const Point& a, const Point& b);
Point operator+(const Point& a, const Point& b);
Point operator*(double scalar, const Point& v);

Point project_to_box(const Point& x, double a, double b);

Point numerical_gradient(const Function& f, const Point& x);

double line_search_backtracking(const Function& f, const Point& x_k, const Point& h_k, double a, double b, double alpha_init, double beta, double sigma);
double line_search_golden_section(const Function& f, const Point& x_k, const Point& h_k, double a, double b, double a_alpha, double b_alpha, double tol);
std::tuple<Point, double, int> minimize(const Function& f, const Point& x0, double a, double b);

Function create_function_1();
Function create_function_2();
Function create_function_3();
Function create_function_4();
Function create_function_5();
Function create_function_6();
Function create_function_7();
Function create_function_8();
Function create_function_9();
Function create_function_10();
Function create_function_11();
Function create_function_12();
Function create_function_13();
Function create_function_14();
Function create_function_15();

std::string format_result(int func_num, const Point& start_point, const Point& min_point, double min_value, int iterations, double a, double b);