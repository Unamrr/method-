#include <iostream>   // для ввода/вывода
#include <cmath>      // для fabs, sqrt, acos
#include <iomanip>    // для форматирования вывода

using namespace std;

const double PI = acos(-1.0);  // число π

//  функция erf из ЛР1
double erf_series(double x) {
    double sum = 0.0;      // накопленная сумма ряда
    double term = x;       // первый член ряда (n=0)
    double prev_sum;       // предыдущая сумма для проверки
    int n = 0;

    do {
        prev_sum = sum;
        sum += term;
        n++;
        // рекуррентная формула для следующего члена
        term *= -x * x * (2 * n - 1) / (n * (2.0 * n + 1));
    } while (sum != prev_sum);  // останов, когда член слишком мал

    return (2.0 / sqrt(PI)) * sum;
}


// Функция для решения системы линейных уравнений методом Гаусса
// Возвращает true, если решение найдено, иначе false (вырожденная матрица)
bool gauss_solve(double** A, double* b, double* x, int n) {
    // Создаём копию матрицы и правой части, чтобы не портить оригинал
    double** M = new double* [n];
    for (int i = 0; i < n; i++) {
        M[i] = new double[n + 1];  // расширенная матрица (A|b)
        for (int j = 0; j < n; j++)
            M[i][j] = A[i][j];
        M[i][n] = b[i];
    }

    //метод Гаусса Жордана
    for (int i = 0; i < n; i++) { // Поиск главного элемента (для устойчивости).
        // ищем в строке i строку с самым большим по модулю элементом
        int max_row = i; 
        for (int k = i + 1; k < n; k++) //k эт оиндекс строки которую мы проверяем в данный млмент
            if (fabs(M[k][i]) > fabs(M[max_row][i]))
                max_row = k; 
        // Меняем строки местами как только нашли строку с бОльшим элементом
        if (max_row != i) {
            for (int j = i; j <= n; j++)
                swap(M[i][j], M[max_row][j]);
        }
        // проверяем,  если гл элемент почти ноль — матрица вырождена
        if (fabs(M[i][i]) < 1e-12) {
            for (int k = 0; k < n; k++) delete[] M[k]; // очищаем память и ретурн фолс тк определитель 0 решения нет или бесконечно много
            delete[] M;
            return false;
        }
        //делим всю строку i на диагональный элемент и как итог на диагонали M[i][i]ровно  1 
        double div = M[i][i];
        for (int j = i; j <= n; j++)
            M[i][j] /= div;
        // Обнуляем остальные строки в столбце i 
        // //для каждой строки k, кроме текущей i, берёт коэффициент factor = M[k][i] и вычитает из строки k строку i, умноженную на этот коэффициент
        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = M[k][i];
                for (int j = i; j <= n; j++)
                    M[k][j] -= factor * M[i][j];
            }
        }
    }//как итог вов сех остальных строках в столбце i становится 0
    // копируем реешние из расшир матрицы в массив x
    for (int i = 0; i < n; i++)
        x[i] = M[i][n];

    // Освобождаем память
    for (int i = 0; i < n; i++) delete[] M[i];
    delete[] M;
    return true;
}

// Вычисление нормы матрицы (максимум суммы модулей по столбцам)
double matrix_norm(double** A, int n) {
    double max_col_sum = 0.0;
    for (int j = 0; j < n; j++) {
        double col_sum = 0.0;
        for (int i = 0; i < n; i++)
            col_sum += fabs(A[i][j]);// складываем модули элементов в столбце j
        if (col_sum > max_col_sum)
            max_col_sum = col_sum;// запоминаем максимум по столбцам
    }
    return max_col_sum;
}

// Вычисление обратной матрицы (для обусловленности) методом Гаусса-Жордана
// Возвращает true, если обратная найдена
bool inverse_matrix(double** A, double** inv, int n) {
    // Создаём расширенную матрицу [A | I]
    double** M = new double* [n];
    for (int i = 0; i < n; i++) {
        M[i] = new double[2 * n]; //выделяем память под расш матрицу
        for (int j = 0; j < n; j++)
            M[i][j] = A[i][j]; // копируем левую часть
        for (int j = n; j < 2 * n; j++) // заполняем правуюч асть
            M[i][j] = (j - n == i) ? 1.0 : 0.0;
    }
    // Прямой ход
    for (int i = 0; i < n; i++) {
        // Поиск главного элемента
        int max_row = i;
        for (int k = i + 1; k < n; k++) //избегаем деления на ноль или очень маленькое число
            if (fabs(M[k][i]) > fabs(M[max_row][i]))
                max_row = k;
        if (max_row != i) {//Если нашли строку с большим элементом не на месте i — меняем строки местами
            for (int j = i; j < 2 * n; j++)
                swap(M[i][j], M[max_row][j]);
        }
        if (fabs(M[i][i]) < 1e-12) { // проверяем не равен ли диаг эл нулю(машинный ноль)
            for (int k = 0; k < n; k++) delete[] M[k]; //Если почти ноль → матрица вырождена, обратной не существует
            delete[] M;
            return false;
        }
        double div = M[i][i];
        for (int j = i; j < 2 * n; j++)
            M[i][j] /= div;// на диагонали м [i] [i] становится 1
        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = M[k][i];
                for (int j = i; j < 2 * n; j++)
                    M[k][j] -= factor * M[i][j]; //итог в столбце i везде 0 кроме позиции гле 1
            }
        }
    }
    // Извлекаем обратную матрицу
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inv[i][j] = M[i][j + n];/*j + n — сдвиг на n столбцов вправо Теперь inv содержит обратную матрицу A⁻¹*/
    //очищаем память и возвращаем результат
    for (int i = 0; i < n; i++) delete[] M[i];
    delete[] M;
    return true;
}

// Вычисление невязки: max |Ax - b|
double nev(double** A, double* x, double* b, int n) {// n размерность 3
    double max_diff = 0.0;//для хранения максим разности, начинаем с 0
    for (int i = 0; i < n; i++) { /*Ax_i — это скалярное произведение i-й строки матрицы A на вектор X

*/
        double Ax_i = 0.0;
        for (int j = 0; j < n; j++)
            Ax_i += A[i][j] * x[j];
        double diff = fabs(Ax_i - b[i]);
        if (diff > max_diff) max_diff = diff; //если текущая разность больше чем пред то запоминаем
    }
    return max_diff; // возврат результата
}

int main() {
    setlocale(LC_ALL, "ru");
    cout << fixed << setprecision(10);
    const int n = 3;

    // -------- Пункт 1: система из ЛР3.1 --------
    // Матрица A
    double** A = new double* [n];
    for (int i = 0; i < n; i++)
        A[i] = new double[n];

    A[0][0] = 1.00; A[0][1] = 0.80; A[0][2] = 0.64;
    A[1][0] = 1.00; A[1][1] = 0.90; A[1][2] = 0.81;
    A[2][0] = 1.00; A[2][1] = 1.10; A[2][2] = 1.21;

    // Правая часть b — используем  функцию erf_series
    double b[3] = { erf_series(0.80), erf_series(0.90), erf_series(1.10) };

    // Решаем систему
    double x[3] = { 0,0,0 };
    bool solved = gauss_solve(A, b, x, n);

    cout << "===== Пункт 1 =====" << endl;
    if (!solved) {
        cout << "Ошибка: матрица вырождена!" << endl;
    }
    else {
        cout << "а) Решение системы:" << endl;
        cout << "   x1 = " << x[0] << endl;
        cout << "   x2 = " << x[1] << endl;
        cout << "   x3 = " << x[2] << endl;

        // Вычисляем обусловленность cond(A) = ||A|| * ||A^{-1}||
        double norm_A = matrix_norm(A, n);
        double** invA = new double* [n];
        for (int i = 0; i < n; i++) invA[i] = new double[n];
        bool inv_ok = inverse_matrix(A, invA, n);
        if (inv_ok) {
            double norm_invA = matrix_norm(invA, n);
            double cond_A = norm_A * norm_invA;
            cout << "   Обусловленность cond(A) = " << cond_A << endl;
            cout << "   (норма A = " << norm_A << ", норма A^{-1} = " << norm_invA << ")" << endl;
        }
        else {
            cout << "   Не удалось вычислить обратную матрицу" << endl;
        }

        // Невязка
        double res = nev(A, x, b, n);
        cout << "б) Невязка |Ax - b| = " << res << endl;

        // Сумма решений и сравнение с erf(1.0)
        double sum_x = x[0] + x[1] + x[2];
        double erf_1 = erf_series(1.0);
        cout << "в) x1 + x2 + x3 = " << sum_x << endl;
        cout << "   erf(1.0) = " << erf_1 << endl;
        cout << "   Разность = " << fabs(sum_x - erf_1) << endl;//сравниваем
       

        // Освобождаем память для обратной матрицы
        for (int i = 0; i < n; i++) delete[] invA[i];
        delete[] invA;
    }

   

    // Освобождаем память
    for (int i = 0; i < n; i++) delete[] A[i];
    delete[] A;

    return 0;
}

/*[1.00  0.80  0.64 | 1.0  0.0  0.0]
[1.00  0.90  0.81 | 0.0  1.0  0.0]
[1.00  1.10  1.21 | 0.0  0.0  1.0]*/

/*[1.00  0.80  0.64 | 1.0  0.0  0.0]
[0.00  0.10  0.17 | -1.0 1.0  0.0]
[0.00  0.30  0.57 | -1.0 0.0  1.0]*/

/*[1.00  0.00  -0.72| 9.0  -8.0  0.0]
[0.00  1.00   1.7 |-10.0  10.0  0.0]
[0.00  0.00   0.06| 2.0  -3.0  1.0]*/
/*[1.00  0.00  0.00 | 6.25  -12.5  6.25]
[0.00  1.00  0.00 |-17.5   35.0 -17.5]
[0.00  0.00  1.00 | 12.5  -25.0  12.5] правая половина это обратная матрица*/
/*Невязка — это насколько сильно 
A * X отличаетсч от B
Невязка=∣A⋅X−B∣
измеряет Насколько точно найденное X удовлетворяет системе*/

/*Формула:
(A* X)i = Ej = 1 до n Aij * xj
/* пример Ax_0 = A[0][0]*x[0] + A[0][1]*x[1] + A[0][2]*x[2];
Ax_0 = 1.00 * x1 + 0.80 * x2 + 0.64 * x3*/;
