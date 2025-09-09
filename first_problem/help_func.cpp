#include "func.hpp"

// --- Функция для формирования строки с результатами ---
std::string format_result(int func_num, const Point& start_point, const Point& min_point, double min_value, int iterations, double a, double b) {
    std::ostringstream oss;
    oss << "Функция " << func_num << "\n";
    oss << "Ограничение на x[0]: [" << a << ", " << b << "]\n";
    oss << "Стартовая точка: ";
    for (size_t i = 0; i < start_point.size(); ++i) {
        oss << std::fixed << std::setprecision(6) << start_point[i];
        if (i < start_point.size() - 1) oss << ", ";
    }
    oss << "\n";
    oss << "Найденная точка минимума: ";
    for (size_t i = 0; i < min_point.size(); ++i) {
        oss << std::fixed << std::setprecision(6) << min_point[i];
        if (i < min_point.size() - 1) oss << ", ";
    }
    oss << "\n";
    oss << "Значение функции в минимуме: " << std::fixed << std::setprecision(10) << min_value << "\n";
    oss << "Количество итераций: " << iterations << "\n";
    oss << "----------------------------------------\n";
    return oss.str();
}