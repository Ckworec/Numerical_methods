#include "func.hpp"

// --- Реализация тестовых функций ---
Function create_function_1() {
    return [](const Point& p) -> double {
        if (p.size() < 2) return std::numeric_limits<double>::infinity();
        return p[0] * p[0] - p[1] * p[1]; // x1^2 - x2^2
    };
}

Function create_function_2() {
    return [](const Point& p) -> double {
        if (p.size() < 2) return std::numeric_limits<double>::infinity();
        return p[0] * p[0] + p[1] * p[1]; // x1^2 + x2^2
    };
}

Function create_function_3() {
    return [](const Point& p) -> double {
        if (p.size() < 2) return std::numeric_limits<double>::infinity();
        return (p[0] - 1.0) * (p[0] - 1.0) + (p[1] + 1.0) * (p[1] + 1.0); // (x1-1)^2 + (x2+1)^2
    };
}

Function create_function_4() {
    return [](const Point& p) -> double {
        if (p.size() < 2) return std::numeric_limits<double>::infinity();
        return p[0] * p[0] * p[0] + p[1] * p[1] * p[1] - 3.0 * p[0] * p[1]; // x1^3 + x2^3 - 3*x1*x2
    };
}

Function create_function_5() {
    return [](const Point& p) -> double {
        if (p.size() < 2) return std::numeric_limits<double>::infinity();
        return p[0] * p[0] - p[0] * p[1] + p[1] * p[1] - 2.0 * p[0] + p[1]; // x1^2 - x1*x2 + x2^2 - 2*x1 + x2
    };
}

Function create_function_6() {
    return [](const Point& p) -> double {
        if (p.size() < 2) return std::numeric_limits<double>::infinity();
        return p[0] * p[0] - p[1] * p[1] - 4.0 * p[0] + 6.0 * p[1]; // x1^2 - x2^2 - 4*x1 + 6*x2
    };
}

Function create_function_7() {
    return [](const Point& p) -> double {
        if (p.size() < 2) return std::numeric_limits<double>::infinity();
        return 2.0 * p[0] * p[0] + p[0] * p[1] + p[1] * p[1]; // 2*x1^2 + x1*x2 + x2^2
    };
}

Function create_function_8() {
    return [](const Point& p) -> double {
        if (p.size() < 2) return std::numeric_limits<double>::infinity();
        return (1.0 - p[0]) * (1.0 - p[0]) + 10.0 * (p[1] - p[0] * p[0]) * (p[1] - p[0] * p[0]); // (1-x1)^2 + 10*(x2-x1^2)^2
    };
}

Function create_function_9() {
    return [](const Point& p) -> double {
        if (p.size() < 2) return std::numeric_limits<double>::infinity();
        return (p[1] - p[0] * p[0]) * (p[1] - p[0] * p[0]) + (1.0 - p[1]) * (1.0 - p[1]); // (x2-x1^2)^2 + (1-x2)^2
    };
}

Function create_function_10() {
    return [](const Point& p) -> double {
        if (p.size() < 2) return std::numeric_limits<double>::infinity();
        return 3.0 * p[0] * p[1] - p[0] * p[0] * p[1] - p[1] * p[1] * p[0]; // 3*x1*x2 - x1^2*x2 - x1*x2^2
    };
}

Function create_function_11() {
    return [](const Point& p) -> double {
        if (p.size() < 2) return std::numeric_limits<double>::infinity();
        return 3.0 * p[0] * p[0] + 4.0 * p[0] * p[1] + p[1] * p[1] - 8.0 * p[0] - 12.0 * p[1]; // 3*x1^2 + 4*x1*x2 + x2^2 - 8*x1 - 12*x2
    };
}

Function create_function_12() {
    return [](const Point& p) -> double {
        if (p.size() < 3) return std::numeric_limits<double>::infinity();
        return p[0] * p[0] + 5.0 * p[1] * p[1] + 3.0 * p[2] * p[2] + 4.0 * p[0] * p[1] - 2.0 * p[1] * p[2] - 2.0 * p[0] * p[2]; // x1^2 + 5*x2^2 + 3*x3^2 + 4*x1*x2 - 2*x2*x3 - 2*x1*x3
    };
}

Function create_function_13() {
    return [](const Point& p) -> double {
        if (p.size() < 2) return std::numeric_limits<double>::infinity();
        return p[0] * p[0] - p[0] * p[1] + p[1] * p[1] - 2.0 * p[0] + 3.0 * p[1] - 4.0; // x1^2 - x1*x2 + x2^2 - 2*x1 + 3*x2 - 4
    };
}

Function create_function_14() {
    return [](const Point& p) -> double {
        if (p.size() < 3) return std::numeric_limits<double>::infinity();
        return -p[0] * p[0] - p[1] * p[1] - p[2] * p[2] - p[0] + p[0] * p[1] + 2.0 * p[1] * p[2]; // -x1^2 - x2^2 - x3^2 - x1 + x1*x2 + 2*x2*x3
    };
}

Function create_function_15() {
    return [](const Point& p) -> double {
        if (p.size() < 3) return std::numeric_limits<double>::infinity();
        return p[0] * p[0] * p[0] + p[1] * p[1] * p[1] + p[2] * p[2] * p[2] + p[1] * p[2] - 3.0 * p[0] + 6.0 * p[1] + 2.0; // x1^3 + x2^3 + x3^3 + x2*x3 - 3*x1 + 6*x2 + 2
    };
}