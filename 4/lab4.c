import numpy as np
import matplotlib.pyplot as plt

# Точность (0.000001)
# Метод Якоби x_new[i] = (b[i] - сумма_остальных) / диагональный_элемент
def jacobi(A, b, x0, eps=1e-6, max_iter=1000):
    n = len(b)
    x = x0.copy()
    residuals = []
    
    for it in range(max_iter):
        x_new = x.copy()#копируем начальное приближение
        for i in range(n):#Считаем сумму всех слагаемых, кроме диагонального.
            sum_ = 0
            for j in range(n):
                if j != i:
                    sum_ += A[i][j] * x[j]
            x_new[i] = (b[i] - sum_) / A[i][i]
        
        residual = np.linalg.norm(A @ x_new - b)
        residuals.append(residual)
        
        if residual < eps:
            return x_new, residuals
        x = x_new #в якоби мы писалив  отдельный массив потом обновляли все сразу а в Зейделе обновляем сразу без отдельного массива
    
    return x, residuals

# ------------------------------------------------------------
# Метод Зейделя 
def seidel(A, b, x0, eps=1e-6, max_iter=1000):
    n = len(b)
    x = x0.copy()
    residuals = []
    
    for it in range(max_iter):
        for i in range(n):
            # Сумма всех слагаемых, кроме диагонального
            sum_ = 0
            for j in range(n):
                if j != i:
                    sum_ += A[i][j] * x[j]  # используем уже обновлённые x[j]!
            x[i] = (b[i] - sum_) / A[i][i]
        
        residual = np.linalg.norm(A @ x - b)
        residuals.append(residual)
        
        if residual < eps: #если достигли точность то возвращаем решение и историю невязок
            return x, residuals
    
    return x, residuals

# ------------------------------------------------------------
def print_results(method_name, x, residuals, eps):
    print(f"\n{method_name}")
    # Формируем строку для любого количества переменных
    x_str = ", ".join([f"{val:.8f}" for val in x])
    print(f"Решение: [{x_str}]")
    print(f"Количество итераций: {len(residuals)}")
    print(f"Достигнутая точность: {residuals[-1]:.2e} (задана {eps})")
# ------------------------------------------------------------
# Основная программа
def main():
   A = np.array([
    [12.14, 1.32, -0.78, -2.75],
    [-0.89, 16.75, 1.88, -1.55],
    [2.65, -1.27, -15.64, -0.64],
    [2.44, 1.52, 1.93, -11.43]
], dtype=float)

b = np.array([14.78, -12.14, -11.65, 4.26], dtype=float)
    eps = 1e-6
    
   starts = [
    np.array([0.0, 0.0, 0.0, 0.0]),        
    np.array([10.0, 10.0, 10.0, 10.0]),   
    np.array([100.0, 0.0, -100.0, 50.0])   
]
    
    for start in starts:
        print("\n" + "="*60)
        print(f"Начальное приближение: [{start[0]}, {start[1]}, {start[2]}, {start[3]}]")
        print("="*60)
        
        # Метод Якоби
        x_jac, res_jac = jacobi(A, b, start, eps)
        print_results("Метод Якоби", x_jac, res_jac, eps)
        
        # Метод Зейделя
        x_sei, res_sei = seidel(A, b, start, eps)
        print_results("Метод Зейделя", x_sei, res_sei, eps)
        
        # Построение графика
        plt.figure(figsize=(10, 5))
        plt.semilogy(res_jac, label='Якоби', marker='o', markersize=3)
        plt.semilogy(res_sei, label='Зейдель', marker='s', markersize=3)
        plt.xlabel('Номер итерации')
        plt.ylabel('Норма невязки (лог. шкала)')
        plt.title(f'Сравнение методов (начальное приближение: {start})')
        plt.grid(True, which="both", ls="--", alpha=0.5)
        plt.axhline(y=eps, color='r', linestyle='--', label=f'Точность {eps}')
        plt.legend()
        plt.grid(True)
        plt.show()

# ------------------------------------------------------------
# ЗАПУСК
if __name__ == "__main__":
    main()

#||r|| = √(r₁² + r₂² + r₃²)
#A @ x_new — умножаем матрицу на вектор (получаем левую часть)
#A @ x_new - b — вычитаем правую часть (получаем невязку)
#np.linalg.norm(...) — вычисляем длину вектора (√(суммы квадратов))
