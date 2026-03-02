#include <iostream>
#include <cmath>
#include <iomanip>
#include <chrono>

const double PI = 3.14159265358979323846;
const double ZETA3 = 1.2020569031595942854; // ζ(3)

// Функция для оценки времени выполнения
class Timer {
private:
    std::chrono::high_resolution_clock::time_point start_time;
public:
    void start() {
        start_time = std::chrono::high_resolution_clock::now();
    }
    
    double elapsed_seconds() {
        auto end_time = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(end_time - start_time).count();
    }
};

// Задача 2: Вычисление s(x) = sum(1/(k(k-x))) - sum(1/(k(k+x)))

// Метод А: Прямое суммирование (медленное)
double s_direct(double x, double eps, int& terms_used, int max_terms = 10000000) {
    double sum1 = 0.0;
    double sum2 = 0.0;
    int k = 1;
    
    while (k <= max_terms) {
        double term1 = 1.0 / (static_cast<double>(k) * (static_cast<double>(k) - x));
        double term2 = 1.0 / (static_cast<double>(k) * (static_cast<double>(k) + x));
        
        sum1 += term1;
        sum2 += term2;
        
        // Оценка остатка через интеграл ~ x/k^2
        if (k > 100) {
            double remainder_estimate = std::abs(x) / (static_cast<double>(k) * static_cast<double>(k));
            if (remainder_estimate < eps) {
                break;
            }
        }
        k++;
    }
    
    terms_used = k;
    return sum1 - sum2;
}

// Метод Б: Ускоренное суммирование (вычитаем асимптотику)
double s_fast(double x, double eps, int& terms_used, int max_terms = 1000000) {
    if (std::abs(x - std::round(x)) < 1e-10 && x != 0) {
        std::cerr << "Ошибка: x не должно быть целым числом" << std::endl;
        return NAN;
    }
    
    double sum = 0.0;
    int k = 1;
    
    // Вычисляем аналитическую часть: 2x * ζ(3)
    double analytic_part = 2.0 * x * ZETA3;
    
    while (k <= max_terms) {
        // Точное значение члена
        double exact_term = 2.0 * x / (static_cast<double>(k) * 
                           (static_cast<double>(k) * static_cast<double>(k) - x * x));
        
        // Асимптотика: 2x/k^3
        double asymp_term = 2.0 * x / (static_cast<double>(k) * 
                           static_cast<double>(k) * static_cast<double>(k));
        
        // Суммируем разность (она убывает как O(1/k^5))
        sum += (exact_term - asymp_term);
        
        // Оценка остатка ~ 2x * x^2 / k^5
        if (k > 10) {
            double remainder_estimate = 2.0 * std::abs(x) * x * x / 
                                       (std::pow(static_cast<double>(k), 5));
            if (remainder_estimate < eps) {
                break;
            }
        }
        k++;
    }
    
    terms_used = k;
    // Добавляем аналитическую часть
    return sum + analytic_part;
}

// Метод В: Использование специальных функций (самый быстрый)
double s_analytic(double x) {
    // Используем формулу: s(x) = (psi(1-x) - psi(1+x) + 2*gamma)/(2x) - 1/x^2
    // Но в C++ нет встроенной psi-функции, поэтому используем приближение
    // через ряд для демонстрации
    
    // Для простоты используем fast метод с очень малой погрешностью
    int dummy;
    return s_fast(x, 1e-15, dummy, 1000);
}

// Задача 3: Вычисление ряда sum(1/(n^2 + 1))

// Прямое суммирование
double sum_direct(double eps, int& terms_used, int max_terms = 10000000) {
    double sum = 0.0;
    int n = 1;
    
    while (n <= max_terms) {
        double term = 1.0 / (static_cast<double>(n) * static_cast<double>(n) + 1.0);
        sum += term;
        
        // Оценка остатка через интеграл: ∫ dx/(x^2+1) = arctan(N) - π/2
        // Приближенно: ~ 1/N для больших N
        if (n > 100) {
            double remainder_estimate = 1.0 / static_cast<double>(n);
            if (remainder_estimate < eps) {
                break;
            }
        }
        n++;
    }
    
    terms_used = n;
    return sum;
}

// Преобразованный ряд (с использованием известных сумм)
double sum_fast(double eps, int& terms_used, int max_terms = 1000000) {
    // Используем тождество из условия:
    // 1/(n^2+1) = 1/n^2 - 1/n^4 + 1/(n^4(n^2+1))
    
    double sum = 0.0;
    int n = 1;
    
    // Аналитические части: ζ(2) и ζ(4)
    double zeta2 = PI * PI / 6.0;
    double zeta4 = PI * PI * PI * PI / 90.0;
    
    while (n <= max_terms) {
        // Остаточный член: 1/(n^4(n^2+1)) ~ 1/n^6
        double term = 1.0 / (std::pow(static_cast<double>(n), 4) * 
                            (static_cast<double>(n) * static_cast<double>(n) + 1.0));
        sum += term;
        
        // Оценка остатка ~ 1/n^5
        if (n > 10) {
            double remainder_estimate = 1.0 / std::pow(static_cast<double>(n), 5);
            if (remainder_estimate < eps) {
                break;
            }
        }
        n++;
    }
    
    terms_used = n;
    // Добавляем аналитические части: ζ(2) - ζ(4) + остаток
    return zeta2 - zeta4 + sum;
}

// Точное значение ряда sum(1/(n^2+1)) через формулу:
// (π * coth(π) - 1)/2
double sum_exact() {
    return (PI * 1.0 / std::tanh(PI) - 1.0) / 2.0;
}

int main() {
    std::cout << std::fixed << std::setprecision(12);
    std::cout << "========================================\n";
    std::cout << "ЗАДАЧА 2: Вычисление s(x)\n";
    std::cout << "========================================\n";
    
    double eps = 8.0/3.0 * 1e-8; // Заданная точность
    Timer timer;
    
    // Тест для x = 0.5
    double x1 = 0.5;
    std::cout << "\n--- x = " << x1 << " ---\n";
    
    timer.start();
    int terms_direct;
    double result_direct = s_direct(x1, eps, terms_direct);
    double time_direct = timer.elapsed_seconds();
    
    timer.start();
    int terms_fast;
    double result_fast = s_fast(x1, eps, terms_fast);
    double time_fast = timer.elapsed_seconds();
    
    std::cout << "Прямой метод:\n";
    std::cout << "  Результат: " << result_direct << "\n";
    std::cout << "  Членов: " << terms_direct << "\n";
    std::cout << "  Время: " << time_direct * 1000 << " мс\n";
    
    std::cout << "Ускоренный метод:\n";
    std::cout << "  Результат: " << result_fast << "\n";
    std::cout << "  Членов: " << terms_fast << "\n";
    std::cout << "  Время: " << time_fast * 1000 << " мс\n";
    
    // Оценка времени при 500 мкс на член
    double time_per_term = 500e-6; // 500 микросекунд
    std::cout << "\nОценка времени при " << time_per_term * 1e6 << " мкс на член:\n";
    std::cout << "  Прямой метод: " << terms_direct * time_per_term << " сек\n";
    std::cout << "  Ускоренный метод: " << terms_fast * time_per_term << " сек\n";
    
    // Тест для x = 0.999999999
    double x2 = 0.999999999;
    std::cout << "\n--- x = " << x2 << " (близко к 1) ---\n";
    
    result_direct = s_direct(x2, eps, terms_direct, 10000000);
    result_fast = s_fast(x2, eps, terms_fast);
    
    std::cout << "Прямой метод:\n";
    std::cout << "  Результат: " << result_direct << "\n";
    std::cout << "  Членов: " << terms_direct << "\n";
    
    std::cout << "Ускоренный метод:\n";
    std::cout << "  Результат: " << result_fast << "\n";
    std::cout << "  Членов: " << terms_fast << "\n";
    
    std::cout << "\n========================================\n";
    std::cout << "ЗАДАЧА 3: Вычисление sum(1/(n^2+1))\n";
    std::cout << "========================================\n";
    
    double eps3 = 0.5e-10; // Погрешность меньше 0.5e-10 для 10-го знака
    
    std::cout << "\nТочное значение: " << sum_exact() << "\n\n";
    
    timer.start();
    int terms_direct3;
    double result_direct3 = sum_direct(eps3, terms_direct3);
    time_direct = timer.elapsed_seconds();
    
    timer.start();
    int terms_fast3;
    double result_fast3 = sum_fast(eps3, terms_fast3);
    time_fast = timer.elapsed_seconds();
    
    std::cout << "Прямое суммирование:\n";
    std::cout << "  Результат: " << result_direct3 << "\n";
    std::cout << "  Ошибка: " << std::abs(result_direct3 - sum_exact()) << "\n";
    std::cout << "  Членов: " << terms_direct3 << "\n";
    std::cout << "  Время: " << time_direct * 1000 << " мс\n";
    
    std::cout << "\nПреобразованный ряд (с использованием ζ(2) и ζ(4)):\n";
    std::cout << "  Результат: " << result_fast3 << "\n";
    std::cout << "  Ошибка: " << std::abs(result_fast3 - sum_exact()) << "\n";
    std::cout << "  Членов: " << terms_fast3 << "\n";
    std::cout << "  Время: " << time_fast * 1000 << " мс\n";
    
    std::cout << "\nСравнение скорости сходимости:\n";
    std::cout << "  Прямой ряд: ~1/N\n";
    std::cout << "  Преобразованный ряд: ~1/N^5\n";
    std::cout << "  Ускорение в ~" << terms_direct3 / terms_fast3 << " раз по числу членов\n";
    
    return 0;
}









#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>

/**
 * Решает квадратное уравнение a*x^2 + b*x + c = 0 с повышенной надежностью.
 * Учитывает переполнение, потерю знаков и машинный ноль.
 *
 * Параметры:
 *   a, b, c - коэффициенты (числа с плавающей точкой двойной точности)
 *
 * Возвращает:
 *   0 - успех, найдены два корня (могут быть равны)
 *   1 - уравнение вырождено (a=0, b=0, c=0) -> бесконечно много решений
 *   2 - уравнение вырождено (a=0, b=0, c!=0) -> нет решений
 *   3 - уравнение линейное (a=0) -> один корень
 *  -1 - ошибка: дискриминант отрицательный (комплексные корни, не рассматриваем в этом примере)
 */
int solve_quadratic(double a, double b, double c, double* x1, double* x2) {
    const double ZERO_THRESHOLD = 1e-14; // Порог для определения "машинного нуля"
    *x1 = *x2 = 0.0;

    // --- 1. Обработка вырожденных случаев (a == 0) ---
    if (fabs(a) < ZERO_THRESHOLD) {
        if (fabs(b) < ZERO_THRESHOLD) {
            if (fabs(c) < ZERO_THRESHOLD) {
                return 1; // Случай: 0 = 0 (все числа)
            } else {
                return 2; // Случай: c = 0 (нет решений)
            }
        } else {
            // Линейное уравнение b*x + c = 0
            *x1 = -c / b;
            return 3; // Один корень
        }
    }

    // --- 2. Масштабирование для предотвращения переполнения (как в Случае 3) ---
    // Ищем максимальный по модулю коэффициент, чтобы привести числа к разумному диапазону.
    // Это важно, если коэффициенты огромны (например, 10^30).
    double max_coeff = fabs(a);
    if (fabs(b) > max_coeff) max_coeff = fabs(b);
    if (fabs(c) > max_coeff) max_coeff = fabs(c);

    // Если числа слишком большие, масштабируем их вниз.
    // Используем степени двойки для масштабирования, чтобы избежать лишних ошибок округления.
    if (max_coeff > 1e20) { // Эмпирический порог
        double scale = max_coeff / 1e10; // Примерный масштаб
        a /= scale;
        b /= scale;
        c /= scale;
        // printf("[Масштабирование] Коэффициенты уменьшены для предотвращения переполнения.\n");
    }


    // --- 3. Вычисление дискриминанта и корней с защитой от потери знаков ---
    double discriminant = b * b - 4.0 * a * c;

    if (discriminant < 0 && discriminant > -ZERO_THRESHOLD) {
        discriminant = 0.0; // Исправление машинного нуля для случая, когда дискриминант близок к нулю
    }

    if (discriminant < 0) {
        return -1; // Комплексные корни
    }

    double sqrt_discriminant = sqrt(discriminant);

    // --- Главный трюк для устойчивости (как в Случае 1) ---
    // Вычисляем "устойчивый" корень, избегая вычитания близких чисел.
    // Используем формулу: x = (-b - sign(b) * sqrt(D)) / (2*a)
    double q;
    if (b < 0) {
        q = (-b + sqrt_discriminant) / 2.0;
    } else {
        q = (-b - sqrt_discriminant) / 2.0;
    }

    // Корни находятся по формулам:
    // x1 = q / a
    // x2 = c / q
    // Эти формулы следуют из теоремы Виета: x1 * x2 = c/a
    // Это позволяет избежать катастрофической потери знаков при вычислении второго корня.
    *x1 = q / a;
    if (fabs(q) > ZERO_THRESHOLD) {
        *x2 = c / q;
    } else {
        // Если q=0, значит один из корней ноль, вычисляем второй стандартно
        *x2 = (-b - sqrt_discriminant) / (2.0 * a);
    }

    // Сортировка: x1 должен быть меньше или равен x2 (для соответствия условию z1 <= z2)
    if (*x1 > *x2) {
        double temp = *x1;
        *x1 = *x2;
        *x2 = temp;
    }

    // --- 4. Проверка на машинный ноль (чтобы не выводить -0.000000) ---
    if (fabs(*x1) < ZERO_THRESHOLD) *x1 = 0.0;
    if (fabs(*x2) < ZERO_THRESHOLD) *x2 = 0.0;

    return 0;
}

// Вспомогательная функция для тестирования
void test_equation(double a, double b, double c, const char* description) {
    double x1, x2;
    printf("\n====================================\n");
    printf("Тест: %s\n", description);
    printf("Уравнение: %.4e x^2 + %.4e x + %.4e = 0\n", a, b, c);

    int result = solve_quadratic(a, b, c, &x1, &x2);

    switch(result) {
        case 0:
            printf("Корни: x1 = %.10f, x2 = %.10f\n", x1, x2);
            // Проверка подстановкой (обратный анализ)
            double f1 = a*x1*x1 + b*x1 + c;
            double f2 = a*x2*x2 + b*x2 + c;
            printf("[Проверка] Невязка для x1: %.2e, для x2: %.2e\n", fabs(f1), fabs(f2));
            break;
        case 1:
            printf("Результат: Уравнение 0=0. Бесконечно много решений.\n");
            break;
        case 2:
            printf("Результат: Противоречие (0 = %f). Нет решений.\n", c);
            break;
        case 3:
            printf("Линейное уравнение. Корень: x = %.10f\n", x1);
            break;
        case -1:
            printf("Результат: Дискриминант отрицательный. Комплексные корни.\n");
            break;
    }
}

int main() {
    setlocale(LC_ALL, "RUS");

    printf("=== Решение квадратных уравнений с контролем неустойчивости ===\n");
    printf("(Основано на материале: Обратная рекурсия, алгоритмы Виета, обратный анализ ошибок)\n");

    // Случай 1: Классическая потеря знаков (из вашего примера)
    test_equation(1.0, -100000.0, 1.0, "Чувствительный случай (b >> 4ac)");

    // Случай 2: Обычное уравнение
    test_equation(6.0, 5.0, -4.0, "Обычное уравнение");

    // Случай 3: Потенциальное переполнение (с масштабированием)
    test_equation(6e30, 5e30, -4e30, "Большие числа (требуется масштабирование)");

    // Случай 4: Огромный и крошечный корень
    test_equation(1e-30, -1.0, 1e30, "Огромный разброс корней (x1 ~ 1e30, x2 ~ 1e-30)");

    // Случай 5: Чувствительное уравнение (близкие корни)
    test_equation(1.0, -4.0, 3.9999999, "Чувствительное уравнение (почти кратный корень)");

    return 0;
}
