import numpy as np
import matplotlib.pyplot as plt

class solver:
    def __init__(self, nu=1, L=35, N=100, t0=0, tN=300, dt=0.05):
        self.nu = nu    
        self.L = L      # length of interval
        self.N = N      # number of grid points
        self.t0 = t0    # start time mmoment
        self.tN = tN    # end time moment
        self.dt = dt    # time step
        self.k = np.arange(-self.N/2, self.N/2, 1) 
        self.FL = ((((2 * np.pi) / self.L) * self.k) ** 2 
                   - self.nu * (((2 * np.pi) / L) * self.k) ** 4 )  # linear differectiaction operator
        self.FN = -(1 / 2) * (1j) * ((2 * np.pi) / self.L) * self.k # nonlinear differenciation operator

    def solve(self):
        nt = int((self.tN - self.t0) /self.dt) 
        self.t = np.linspace(start=self.t0, stop=self.tN, num=nt)
        self.x = np.linspace(start=0, stop=self.L, num=self.N) 

        self.u0 = -np.cos((2 * np.pi * self.x) / self.L) + np.sin((4 * np.pi * self.x) / self.L)

        u0_hat = (1 / self.N) * np.fft.fftshift(np.fft.fft(self.u0))
        u0_hat2 = (1 / self.N) * np.fft.fftshift(np.fft.fft(self.u0**2))

        u = np.zeros((self.N, nt))
        u_hat = np.zeros((self.N, nt), dtype=complex)
        u_hat2 = np.zeros((self.N, nt), dtype=complex)

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

            u_hat[:,i+1] = (1 / (1 - (self.dt / 2) * self.FL)) * ((1 + (self.dt / 2) * self.FL) * uhat_curr
                            + (((3 / 2) * self.FN) * (uhat_curr2) - ((1 / 2) * self.FN) * (uhat_prev2)) * self.dt)
            # back to real space
            u[:,i+1] = np.real(self.N * np.fft.ifft(np.fft.ifftshift(u_hat[:,i+1])))
            # corretion of coefficient
            u_hat[:,i+1] = (1 / self.N) * np.fft.fftshift(np.fft.fft(u[:,i+1]))
            # computing of nonlinear part
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
