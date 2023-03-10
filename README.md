# Deep-fourier-collocation

Under-resolved simulations of chaotic dynamics by deep spectral methods.

## Introduction

Chaotic dynamics is a popular object of research. For discribe chaotic dymanimcs by math ussualy used partial differectial equations (PDE). The Kuramoto-Sivashinsky equation will be considered as an object of chaotic dynamics research in this paper.

## Methods

Fourier Collocation method:


$$
	\partial_tu(t, x) = -\partial_{xx}u(t, x)- \partial_{xxxx}u(t, x)- u(t, x)\partial_xu(t, x), \ u(t, x_j) = \sum\limits_{k = N/2}^{N/2-1}\hat u_k(t)e^{ikx_j\frac{2}{L}},
$$

$$
	\partial_t\left(\sum\limits_{k =- N/2}^{N/2-1}\hat u_k(t)e^{ikx_j\frac{2}{L}}\right) = - \partial_{xx} \left(\sum\limits_{k =- N/2}^{N/2-1}\hat u_k(t)e^{ikx_j\frac{2}{L}}\right) - \partial_{xxxx} \left(\sum\limits_{k =- N/2}^{N/2-1}\hat u_k(t)e^{ikx_j\frac{2}{L}}\right) -
$$

$$
-\partial_{x}\left(\sum\limits_{k = -N/2}^{N/2-1}\hat u_k(t)e^{ikx_j\frac{2}{L}}\right) \cdot \sum\limits_{k = -N/2}^{N/2-1}\hat u_k(t)e^{ikx_j\frac{2}{L}} = 
$$

$$
	= \sum\limits_{k = N/2}^{N/2-1}k\left(\frac{2}{L} \right)^2\hat u_k(t)e^{ikx_j\frac{2}{L}} - \sum\limits_{k =- N/2}^{N/2-1}\left(\frac{2k}{L} \right)^4\hat u_k(t)e^{ikx_j\frac{2k}{L}} - 
$$

$$
\- i\sum\limits_{k =- N/2}^{N/2-1}\left(\frac{2k}{L} \right)\hat u_ke^{ikx_j\frac{2}{L}} \cdot \sum\limits_{k =- N/2}^{N/2-1}\hat u_k(t)e^{ikx_j\frac{2}{L}} =
$$

$$
	= \sum\limits_{k =- N/2}^{N/2-1}\left(\frac{2}{L} \right)^2\hat u_k(t)e^{ikx_j\frac{2k}{L}} - \sum\limits_{k =- N/2}^{N/2-1}\left(\frac{2}{L} \right)^4\hat u_kk(t)e^{ikx_j\frac{2k}{L}} -
$$

$$
	- i\left(\frac{2}{L} \right)\sum\limits_{p = -N/2}^{N/2}\sum\limits_{l = -N/2}^{N/2}\hat u_p(t)\hat u_l(t)e^{i(p+l)x_j\frac{2}{L}} = 
$$

$$
	= \sum\limits_{k =- N/2}^{N/2-1}\left(\frac{2k}{L} \right)^2\hat u_k(t)e^{ikx_j\frac{2}{L}} - \sum\limits_{k =- N/2}^{N/2-1}\left(\frac{2k}{L} \right)^4\hat u_k(t)e^{ikx_j\frac{2}{L}} -
$$

$$
	- i\left(\frac{2k}{L} \right) \sum\limits_{p+l=k}\hat u_p(t)\hat u_l(t)e^{i(p+l)x_j\frac{2}{L}}
$$

$$
	\partial_t \hat u_k(t) = \left(\frac{2k}{L} \right)^2\hat u_k(t) - \left(\frac{2k}{L} \right)^4 \hat u_k(t) - \frac{1}{2}\cdot \frac{2ik}{L}\sum\limits_{l+p = k}\hat u_p(t) \hat u_l(t).
$$

## Results of Under-resolved simulations

![KS_20_AB_10pi_nu1](https://user-images.githubusercontent.com/60492990/214952806-d9c472fe-3808-4206-8b09-795e59c003da.png)

![KS_50_AB_10pi_nu1](https://user-images.githubusercontent.com/60492990/214952928-b39ca89b-2a31-4d82-bb81-2030681098b0.png)

![KS_100_AB_10pi_nu1](https://user-images.githubusercontent.com/60492990/214952954-079d4c43-cbcc-4a6e-b1cf-658fbaeaa973.png)


