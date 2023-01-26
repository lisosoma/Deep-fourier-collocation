import numpy as np
import matplotlib.pyplot as plt

class solver:
    def __init__(self, nu=1, L=35, N=100, t0=0, tN=300, dt=0.05):
        self.nu = nu    # параметр, имеющий физический смысл в уравнении
        self.L = L      # длина отрезка интегрирования по х
        self.N = N      # количество узлов расчетной сетки
        self.t0 = t0    # начальный момент времени
        self.tN = tN    # конечный момент времени
        self.dt = dt    # шаг по времени
        self.k = np.arange(-self.N/2, self.N/2, 1) # коэффициенты k для вычисления операторов дифференцирования
        self.FL = ((((2 * np.pi) / self.L) * self.k) ** 2 
                   - self.nu * (((2 * np.pi) / L) * self.k) ** 4 )  # оператор дифференцирования для линейной части
        self.FN = -(1 / 2) * (1j) * ((2 * np.pi) / self.L) * self.k # оператор дифференцирования для нелинейной части

    def solve(self):
        # количество шагов по времени
        nt = int((self.tN - self.t0) /self.dt) 

        # сетка
        self.t = np.linspace(start=self.t0, stop=self.tN, num=nt)
        self.x = np.linspace(start=0, stop=self.L, num=self.N) 

        # начальные условия
        self.u0 = -np.cos((2 * np.pi * self.x) / self.L) + np.sin((4 * np.pi * self.x) / self.L)

        # коэффициенты пространства Фурье
        u0_hat = (1 / self.N) * np.fft.fftshift(np.fft.fft(self.u0))
        # отдельно для нелинейной части
        u0_hat2 = (1 / self.N) * np.fft.fftshift(np.fft.fft(self.u0**2))

        # массивы с решениями
        u = np.zeros((self.N, nt))
        u_hat = np.zeros((self.N, nt), dtype=complex)
        u_hat2 = np.zeros((self.N, nt), dtype=complex)

        # инициализация массивов
        u[:,0] = self.u0
        u_hat[:,0] = u0_hat
        u_hat2[:,0] = u0_hat2

        for i in range(0, nt-1):
            uhat_curr = u_hat[:,i]
            uhat_curr2 = u_hat2[:,i]
            if i == 0:
                uhat_prev2 = u_hat2[:,0]
            else:
                uhat_prev2 = u_hat2[:,i-1]
            # схема Кранка-Николсона для линейной части и Адамса-Башфорта для нелинейной части
            # таймстепинг в пространстве коэффициентов Фурье
            u_hat[:,i+1] = (1 / (1 - (self.dt / 2) * self.FL)) * ((1 + (self.dt / 2) * self.FL) * uhat_curr
                            + (((3 / 2) * self.FN) * (uhat_curr2) - ((1 / 2) * self.FN) * (uhat_prev2)) * self.dt)
            # возврат в физическое пространство
            u[:,i+1] = np.real(self.N * np.fft.ifft(np.fft.ifftshift(u_hat[:,i+1])))
            # корректировка коэффициента
            u_hat[:,i+1] = (1 / self.N) * np.fft.fftshift(np.fft.fft(u[:,i+1]))
            # вычисление нелинейной части
            u_hat2[:,i+1] = (1 / self.N) * np.fft.fftshift(np.fft.fft(u[:,i+1]**2))

        self.u = u


    def show(self):
        fig, ax = plt.subplots(figsize=(25,6))
        tt, xx = np.meshgrid(self.t, self.x)
        levels = np.arange(-6, 6, 0.01)
        cs = ax.contourf(tt, xx, self.u, cmap='plasma')
        ax.set_xlabel("t", fontsize=20)
        ax.set_ylabel("x",  fontsize=20)
        plt.rcParams.update({'font.size': 20})
        ax.set_title(f"Kuramoto-Sivashinsky: N = {self.N}, L = {self.L}, nu = {self.nu}")
