#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace std;

double norm(const vector<vector<double>>& A,
    const vector<double>& x,
    const vector<double>& b) {
    int n = x.size();
    double sum = 0;
    for (int i = 0; i < n; ++i) {
        double Ax = 0;
        for (int j = 0; j < n; ++j) {
            Ax += A[i][j] * x[j];
        }
        double r = Ax - b[i];
        sum += r * r;
    }
    return sqrt(sum);
}


void jacobi(const vector<vector<double>>& A,
    const vector<double>& b,
    vector<double>& x,
    double eps,
    int max_iter,
    vector<double>& residuals) {

    int n = x.size();
    vector<double> x_new(n);
    residuals.clear();

    for (int iter = 0; iter < max_iter; ++iter) {
        for (int i = 0; i < n; ++i) {
            double sum = 0;
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    sum += A[i][j] * x[j];
                }
            }
            x_new[i] = (b[i] - sum) / A[i][i];
        }

        double res = norm(A, x_new, b);
        residuals.push_back(res);

        if (res < eps) {
            x = x_new;
            return;
        }

        x = x_new;
    }
}


void seidel(const vector<vector<double>>& A,
    const vector<double>& b,
    vector<double>& x,
    double eps,
    int max_iter,
    vector<double>& residuals) {

    int n = x.size();
    residuals.clear();

    for (int iter = 0; iter < max_iter; ++iter) {
        for (int i = 0; i < n; ++i) {
            double sum = 0;
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    sum += A[i][j] * x[j];
                }
            }
            x[i] = (b[i] - sum) / A[i][i];
        }

        double res = norm(A, x, b);
        residuals.push_back(res);

        if (res < eps) {
            return;
        }
    }
}


string doubleToStr(double val) {
    stringstream ss;
    ss << fixed << scientific << val;
    string s = ss.str();
    // Заменяем точку на запятую
    for (char& c : s) {
        if (c == '.') c = ',';
    }
    return s;
}


void saveToCSV(const string& filename,
    const vector<double>& res_jac,
    const vector<double>& res_sei,
    const vector<double>& start,
    double eps) {
    ofstream file(filename);

    // Заголовок
    file << "Номер итерации;Невязка Якоби;Невязка Зейделя\n";

  
    size_t max_len = max(res_jac.size(), res_sei.size());

    for (size_t i = 0; i < max_len; ++i) {
        file << i << ";";

        if (i < res_jac.size()) {
            file << doubleToStr(res_jac[i]);
        }
        file << ";";

        if (i < res_sei.size()) {
            file << doubleToStr(res_sei[i]);
        }

        file << "\n";
    }

    file.close();
    cout << "   Сохранено в файл: " << filename << endl;
}



void printResults(const string& name, const vector<double>& x,
    const vector<double>& residuals, double eps) {
    cout << "\n--- " << name << " ---\n";
    cout << "Решение: ";
    for (double val : x) {
        cout << fixed << setprecision(8) << val << " ";
    }
    cout << "\nКоличество итераций: " << residuals.size() << endl;
    cout << "Достигнутая точность: " << scientific << residuals.back()
        << " (задана " << eps << ")" << endl;
}


int main() {
    setlocale(LC_ALL, "ru");
    // Матрица A (твоя система из 4 уравнений)
    vector<vector<double>> A = {
        {12.14,  1.32, -0.78, -2.75},
        {-0.89, 16.75,  1.88, -1.55},
        {2.65,  -1.27, -15.64, -0.64},
        {2.44,   1.52,  1.93, -11.43}
    };

    // Вектор правых частей b
    vector<double> b = { 14.78, -12.14, -11.65, 4.26 };

    double eps = 1e-6;
    int max_iter = 1000;

    // Разные начальные приближения
    vector<vector<double>> starts = {
        {0.0, 0.0, 0.0, 0.0},
        {10.0, 10.0, 10.0, 10.0},
        {100.0, 0.0, -100.0, 50.0}
    };

    vector<string> start_names = {
        "start_0_0_0_0",
        "start_10_10_10_10",
        "start_100_0_-100_50"
    };

    for (size_t k = 0; k < starts.size(); ++k) {
        cout << "Начальное приближение: [";
        for (size_t i = 0; i < starts[k].size(); ++i) {
            cout << starts[k][i] << (i < starts[k].size() - 1 ? ", " : "");
        }
        cout << "]" << endl;

        vector<double> res_jac, res_sei;

        // Метод Якоби
        vector<double> x_jac = starts[k];
        jacobi(A, b, x_jac, eps, max_iter, res_jac);
        printResults("Метод Якоби", x_jac, res_jac, eps);

        // Метод Зейделя
        vector<double> x_sei = starts[k];
        seidel(A, b, x_sei, eps, max_iter, res_sei);
        printResults("Метод Зейделя", x_sei, res_sei, eps);

        // Сохраняем в CSV файл (с запятыми для Excel)
        string filename = "residuals_" + start_names[k] + ".csv";
        saveToCSV(filename, res_jac, res_sei, starts[k], eps);
    }

  

    return 0;
}
//чем больше итераций тем меньше невязка приближается к 0
//Зейдель использует новые значения сразу, а Якоби все считает последовательно по старым
//также Зейдель ближе к истие
