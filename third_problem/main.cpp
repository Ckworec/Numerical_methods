#include <iostream>
#include <cmath>

using namespace std;

// Правая часть дифференциального уравнения y' = f(x,y)
double f(double x, double y)
{
    return y;
}

double RK_variant1(double x, double y, double h)
{
    double k1 = h * f(x, y);

    double k2 = h * f(x + h/2.0,
                      y + k1/2.0);

    double k3 = h * f(x + h,
                      y - k1 + 2.0*k2);

    return y + (k1 + 4.0*k2 + k3)/6.0;
}

double RK_variant2(double x, double y, double h)
{
    double k1 = h * f(x, y);

    double k2 = h * f(x + h/2.0,
                      y + k1/2.0);

    double k3 = h * f(x + h/2.0,
                      y + k2/2.0);

    double k4 = h * f(x + h,
                      y + k3);

    return y + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
}

double rungeKuttaControlStep(double x, double y, double h, double &error)
{
    double k1 = h * f(x, y);

    double k2 = h * f(x + h/3.0,
                      y + k1/3.0);

    double k3 = h * f(x + h/3.0,
                      y + (k1 + k2)/6.0);

    double k4 = h * f(x + h/2.0,
                      y + k1/8.0 + 3.0*k3/8.0);

    double k5 = h * f(x + h,
                      y + k1/2.0 - 3.0*k3/2.0 + 2.0*k4);

    double y_next = y + (k1 + 4.0*k4 + k5)/6.0;

    error = (2.0*k1 - 9.0*k3 + 8.0*k4 - k5) / 30.0;

    return y_next;
}

int main()
{
    double a = 0.0;
    double b = 1.0;
    double y = 1.0;

    double h = 0.2;
    double eps = 0.001;

    double x = a;

    cout << "x\t y\t h\n";

    while (x < b)
    {
        if (x + h > b)
            h = b - x;

        double err;
        double y1 = RK_variant1(x, y, h);
        // double y1 = RK_variant2(x, y, h);
        // double y1 = rungeKuttaControlStep(x, y, h, err);

        double rho = fabs(err);
        double limit = eps * h / (b - a);

        if (rho > limit)
        {
            h /= 2.0;
            continue;
        }

        x += h;
        y = y1;

        cout << x << "\t" << y << "\t" << h << endl;
    }

    return 0;
}
