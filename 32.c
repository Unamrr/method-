#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

// Метод прогонки для решения СЛАУ с трёхдиагональной матрицей
bool tridiagonal_solve(const double alpha[], const double beta[], const double gamma[],
    const double b[], double x[], int n) {

    double* c = new double[n + 2];  // прогоночные коэффициенты
    double* d = new double[n + 2];

    // ----- ПРЯМАЯ ПРОГОНКА -----
    // Начальные значения c[2] и d[2] (в индексации методички)
    // В коде: c[2] соответствует i=2, используем сдвиг индексов
    if (fabs(beta[0]) < 1e-12) {
        cout << "Ошибка: beta[1] = 0!" << endl;
        delete[] c; delete[] d;
        return false;
    }

    c[2] = -gamma[0] / beta[0];   // c2 = -γ1/β1 с2 и d2 это отправные точки со всей прогонки
    d[2] = b[0] / beta[0];         // d2 = b1/β1

    // Вычисляем c[i+1] и d[i+1] для i = 2, 3, 4 (до n-1)
    for (int i = 2; i <= n - 1; i++) {
        double znamenatel = alpha[i - 2] * c[i] + beta[i - 1];  // αi*ci + βi  (исправлено: alpha[i-2])
        if (fabs(znamenatel) < 1e-12) {
            cout << "Ошибка: деление на ноль!" << endl;
            delete[] c; delete[] d;
            return false;
        }
        c[i + 1] = -gamma[i - 1] / znamenatel;                    // c_{i+1}
        d[i + 1] = (b[i - 1] - alpha[i - 2] * d[i]) / znamenatel;   // d_{i+1} (исправлено: alpha[i-2])
    }
    //xi = ci+1 *xi+1 + di+1 мы находим с и d тк они все связаны между собой предполагаем, идем с конца типа тк чтоб найти х1 надо знать х2 и так далее
    //находим чтобы обратная прогонка могла начаться с самого правогох и дойти до самого левого х
    // ----- ОБРАТНАЯ ПРОГОНКА -----
    // Находим xn
    double denom_last = alpha[n - 2] * c[n] + beta[n - 1];  // исправлено: alpha[n-2]
    if (fabs(denom_last) < 1e-12) {
        cout << "Ошибка: деление на ноль!" << endl;
        delete[] c; delete[] d;
        return false;
    }
    x[n - 1] = (b[n - 1] - alpha[n - 2] * d[n]) / denom_last;  // xn (исправлено: alpha[n-2])

    // Находим остальные xi (с конца)
    for (int i = n - 1; i >= 1; i--) {
        x[i - 1] = c[i + 1] * x[i] + d[i + 1];  // пример x3 = c4*x4 + d4 и так до х1
    }

    delete[] c;
    delete[] d;
    return true;
}

int main() {
    setlocale(LC_ALL, "ru");
    cout << fixed << setprecision(4);

    // Система из методички (индексация: α2, α3, α4, α5; β1..β5; γ1..γ4)
    // α2=-2, α3=2, α4=1, α5=3
    // β1=1, β2=4, β3=-2, β4=1, β5=-1  
    // γ1=3, γ2=-1, γ3=1, γ4=1
    // b1=5, b2=1, b3=3, b4=-2, b5=-1

    const int n = 5;

    // ВАЖНО: индексы массивов в C++ начинаются с 0
    // alpha[i] соответствует α_{i+2} (потому что α1 не существует)
    double alpha[] = { -2, 2, 1, 3 };   // α2, α3, α4, α5
    double beta[] = { 1, 4, -2, 1, -1 };   // β1, β2, β3, β4, β5
    double gamma[] = { 3, -1, 1, 1 };   // γ1, γ2, γ3, γ4
    double b[] = { 5, 1, 3, -2, -1 };   // b1, b2, b3, b4, b5

    double x[n];

    cout << "===== МЕТОД ПРОГОНКИ =====" << endl;
    cout << "\nИсходная система:" << endl;
    cout << "  x1 + 3x2 = 5" << endl;
    cout << "-2x1 + 4x2 -  x3 = 1" << endl;
    cout << "      2x2 - 2x3 +  x4 = 3" << endl;
    cout << "            x3 +  x4 + x5 = -2" << endl;
    cout << "                 3x4 - x5 = -1" << endl;

    if (tridiagonal_solve(alpha, beta, gamma, b, x, n)) {
        cout << "\nРешение системы:" << endl;
        for (int i = 0; i < n; i++) {
            cout << "x" << i + 1 << " = " << x[i] << endl;
        }
    }
    else {
        cout << "Ошибка!" << endl;
    }

    return 0;
}
