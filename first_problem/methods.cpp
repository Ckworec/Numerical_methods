#include "func.hpp"

// --- Методы одномерной минимизации ---

// 1. Метод дробления шага с проекцией
double line_search_backtracking(const Function& f, const Point& x_k, const Point& h_k, double a, double b, double alpha_init = 1.0, double beta = 0.5, double sigma = 1e-4) {
    double alpha = alpha_init;
    Point x_k_proj = project_to_box(x_k, a, b);
    double fxk = f(x_k_proj);
    
    // Используем скалярное произведение градиента и направления для критерия убывания
    Point grad_fk = numerical_gradient(f, x_k_proj);
    double directional_derivative = dot(grad_fk, h_k); // Производная по направлению h_k

    while (alpha > 1e-16) {
        Point x_new = x_k + alpha * h_k; // Пробный шаг без проекции
        Point x_new_proj = project_to_box(x_new, a, b); // Проектируем на допустимую область
        double fx_new = f(x_new_proj);
        
        // Критерий Армихо: f(x_new_proj) <= f(x_k_proj) + sigma * alpha * (grad_f * h_k)
        if (fx_new <= fxk + sigma * alpha * directional_derivative) {
            return alpha;
        }
        alpha *= beta;
    }

    return alpha;
}

// 2. Метод золотого сечения с проекцией
double line_search_golden_section(const Function& f, const Point& x_k, const Point& h_k, double a, double b, double a_alpha, double b_alpha, double tol) {
    const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
    const double resphi = 2.0 - phi;

    double x1 = a_alpha + resphi * (b_alpha - a_alpha);
    double x2 = b_alpha - resphi * (b_alpha - a_alpha);
    
    Point p1_trial = x_k + x1 * h_k;
    Point p1 = project_to_box(p1_trial, a, b); // Проектируем
    
    Point p2_trial = x_k + x2 * h_k;
    Point p2 = project_to_box(p2_trial, a, b); // Проектируем
    
    double f1 = f(p1);
    double f2 = f(p2);

    while (std::abs(b_alpha - a_alpha) > tol) {
        if (f1 < f2) {
            b_alpha = x2;
            x2 = x1;
            f2 = f1;
            x1 = a_alpha + resphi * (b_alpha - a_alpha);
            
            Point p1_trial = x_k + x1 * h_k;
            p1 = project_to_box(p1_trial, a, b);
            f1 = f(p1);
        } else {
            a_alpha = x1;
            x1 = x2;
            f1 = f2;
            x2 = b_alpha - resphi * (b_alpha - a_alpha);
            
            Point p2_trial = x_k + x2 * h_k;
            p2 = project_to_box(p2_trial, a, b);
            f2 = f(p2);
        }
    }

    return (a_alpha + b_alpha) / 2.0;
}


// --- Основная функция минимизации с ограничением на первую координату ---
std::tuple<Point, double, int> minimize(const Function& f, const Point& x0, double a, double b) {
    Point x_k = x0;
    // Убеждаемся, что начальная точка удовлетворяет ограничению
    x_k[0] = std::max(a, std::min(b, x_k[0]));
    
    int total_iterations = 0;
    const double sqrt_eps = std::sqrt(eps);
    const size_t n = x0.size(); // Определяем размерность

    // --- Этап 1: Метод дробления шага ---
    double current_eps = sqrt_eps;
    int stage1_iter = 0;

    while (true) {
        Point grad = numerical_gradient(f, x_k);
        double grad_norm = norm(grad);
        Point h_k = -1.0 * grad; // Направление антиградиента

        if (grad_norm < current_eps) {
            // std::cout << "Этап 1 завершен. Норма градиента < sqrt(eps)." << std::endl;
            break;
        }

        // Одномерная минимизация методом дробления шага с проекцией
        double alpha_k = line_search_backtracking(f, x_k, h_k, a, b);

        Point x_k1 = x_k + alpha_k * h_k;
        x_k1 = project_to_box(x_k1, a, b); // Проектируем результат шага
        total_iterations++;
        stage1_iter++;

        x_k = x_k1;

        if (total_iterations > max_iter) {
            // std::cout << "Превышено максимальное количество итераций на этапе 1." << std::endl;
            break;
        }
    }

    // --- Этап 2: Метод золотого сечения ---
    current_eps = eps;
    int stage2_iter = 0;
    double last_alpha_guess = 1.0; // Сохраняем для эвристики

    while (true) {
        Point grad = numerical_gradient(f, x_k);
        double grad_norm = norm(grad);
        Point h_k = -1.0 * grad; // Направление антиградиента

        if (grad_norm < current_eps) {
            // std::cout << "Этап 2 завершен. Норма градиента < eps." << std::endl;
            break;
        }

        // Одномерная минимизация методом золотого сечения с проекцией
        double alpha_guess = last_alpha_guess;
        if (stage2_iter > 0) {
             alpha_guess = 2.0 * last_alpha_guess; // Эвристика
        }
        // Определяем отрезок поиска [a_alpha, b_alpha] для alpha
        double a_alpha = 0.0;
        // Убедимся, что b_alpha > a_alpha
        double b_alpha = (alpha_guess > 0) ? (3.0 * alpha_guess) : 1.0;
        if (b_alpha <= a_alpha) b_alpha = a_alpha + 1.0;

        // Точность поиска alpha
        double alpha_tol = eps / (1.0 + std::abs(dot(grad, h_k))); // Адаптивная точность
        alpha_tol = std::max(alpha_tol, 1e-12); // Минимальная точность

        double alpha_k = line_search_golden_section(f, x_k, h_k, a, b, a_alpha, b_alpha, alpha_tol);
        last_alpha_guess = alpha_k; // Сохраняем для следующей итерации

        Point x_k1 = x_k + alpha_k * h_k;
        x_k1 = project_to_box(x_k1, a, b); // Проектируем результат шага
        total_iterations++;
        stage2_iter++;

        // Проверка на сходимость по аргументу (в проекции)
        Point diff = x_k1 - x_k;
        Point diff_proj = project_to_box(diff, a-a, b-b);
        // Более простая проверка: норма разности
        if (norm(x_k1 - x_k) < eps) {
             // std::cout << "Сходимость по аргументу достигнута." << std::endl;
             x_k = x_k1;
             break;
        }

        x_k = x_k1;

        if (total_iterations > max_iter) {
            // std::cout << "Превышено максимальное количество итераций на этапе 2." << std::endl;
            break;
        }
    }

    return std::make_tuple(x_k, f(x_k), total_iterations);
}