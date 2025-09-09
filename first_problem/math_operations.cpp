#include "func.hpp"

// --- Вспомогательные функции для работы с векторами ---
double norm(const Point& v) {
    double sum = 0.0;

    for (double elem : v) {
        sum += elem * elem;
    }

    return std::sqrt(sum);
}

double dot(const Point& a, const Point& b) {
    double sum = 0.0;

    for (size_t i = 0; i < a.size(); ++i) {
        sum += a[i] * b[i];
    }

    return sum;
}

Point operator-(const Point& a, const Point& b) {
    Point result(a.size());

    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] - b[i];
    }

    return result;
}

Point operator+(const Point& a, const Point& b) {
    Point result(a.size());
    
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] + b[i];
    }

    return result;
}

Point operator*(double scalar, const Point& v) {
    Point result(v.size());

    for (size_t i = 0; i < v.size(); ++i) {
        result[i] = scalar * v[i];
    }

    return result;
}

// --- Проекция точки на ограничения ---
// Проектирует первую координату на отрезок [a, b]
Point project_to_box(const Point& x, double a, double b) {
    Point x_proj = x;
    x_proj[0] = std::max(a, std::min(b, x[0]));
    return x_proj;
}

// --- Численное дифференцирование ---
Point numerical_gradient(const Function& f, const Point& x) {
    size_t n = x.size();
    Point grad(n);
    Point x_plus_h = x;
    Point x_minus_h = x;

    for (size_t i = 0; i < n; ++i) {
        x_plus_h[i] = x[i] + tao;
        x_minus_h[i] = x[i] - tao;
        double f_plus = f(x_plus_h);
        double f_minus = f(x_minus_h);
        grad[i] = (f_plus - f_minus) / (2 * tao);

        // Восстанавливаем x_plus_h
        x_plus_h[i] = x[i];
        x_minus_h[i] = x[i];
    }

    return grad;
}