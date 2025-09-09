#include "func.hpp"

int main() {
    std::ofstream output_file("results.txt");

    output_file << "Результаты тестирования алгоритма минимизации с ограничением на x[0]\n";
    output_file << "=====================================================================\n";

    // Список функций и их стартовых точек (далеко от минимума)
    // Добавим границы отрезка [a, b] для x[0]
    // func, point, boundary
    std::vector<std::tuple<Function, Point, double, double>> test_cases = {
        {create_function_1(), {-10.0, 5.0}, -15.0, 15.0}, // f1: Минимум вдоль прямой x1=-x2
        {create_function_2(), {10.0, -15.0}, -20.0, 20.0}, // f2: Минимум в (0,0)
        {create_function_3(), {-5.0, 10.0}, -10.0, 10.0}, // f3: Минимум в (1,-1)
        {create_function_4(), {-2.0, 3.0}, -5.0, 5.0}, // f4: Минимум в (1,1) или (-1,-1)
        {create_function_5(), {10.0, -5.0}, 0.0, 20.0}, // f5: Минимум в (1,1)
        {create_function_6(), {5.0, -10.0}, -5.0, 15.0}, // f6: Минимум в (2,3)
        {create_function_7(), {-3.0, 2.0}, -10.0, 10.0}, // f7: Минимум в (0,0)
        {create_function_8(), {-1.0, 5.0}, -2.0, 2.0}, // f8: Минимум в (1,1)
        {create_function_9(), {3.0, -2.0}, 0.0, 2.0}, // f9: Минимум в (1,1)
        {create_function_10(), {5.0, -3.0}, 0.0, 10.0}, // f10: Минимум в (1,1) или (0,0)
        {create_function_11(), {10.0, -15.0}, -5.0, 15.0}, // f11: Минимум в (2,3)
        {create_function_12(), {-5.0, 10.0, -8.0}, -10.0, 10.0}, // f12: Минимум в (0,0,0)
        {create_function_13(), {5.0, -10.0}, 0.0, 5.0}, // f13: Минимум в (2,1)
        {create_function_14(), {5.0, -10.0, 8.0}, -10.0, 10.0}, // f14: Минимум в около (0.5, 0.5, 0.5)
        {create_function_15(), {-5.0, 10.0, -8.0}, -10.0, 10.0} // f15: Минимум в около (1, -1, -1)
    };

    for (size_t i = 0; i < test_cases.size(); ++i) {
        auto [func, start_point_data, a, b] = test_cases[i];
        Point start_point = start_point_data; // Копируем, так как minimize может модифицировать
        auto [min_point, min_value, iterations] = minimize(func, start_point, a, b);

        std::string result_str = format_result(i + 1, start_point, min_point, min_value, iterations, a, b);

        output_file << result_str;
    }

    output_file.close();

    std::cout << "Тестирование завершено. Результаты сохранены в файле 'results.txt'." << std::endl;

    return 0;
}