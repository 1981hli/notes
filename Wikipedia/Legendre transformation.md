> 2022.10.03

# Legendre transformation

- involutory function / self-inverse function: 反函数是自己的函数.

- Legendre transformation is an involutive transformation on real-valued convex functions of one real variable.

- 物理中, Legendre 变换用来将一些量 (速度, 压强, 温度) 的函数变成其共轭量 (动量, 体积, 熵) 的函数.
  常用来由 Lagrange 公式给出 Hamilton 公式; 给出热力学势; 给出多元微分方程的解.

- 光滑函数 $f$ 的 Legendre 变换 $f^*$ 满足它们的一阶导函数互为反函数.
  $$Df(\cdot)=(Df^*)^{-1}(\cdot)$$

## 1 定义

### 1.1 通过求导理解 Legendre 变换

- 函数 $f$ 的导函数 $f'$ 与其的 Legendre 变换 $f^*$ 的导函数 $f^{*'}$ 互为反函数

## 2 性质

- 凸函数的 Legendre 变换是凸函数.

- Legendre 变换是 involution, 即 $f^{**}=f$

## 3 例

## 4 Legendre 变换下的微分

## 5 应用

## 6 几何解释

## 7 高维 Legendre 变换

## 8 流形上的 Legendre 变换

- 对于矢量丛 $(M,E,\pi)$, $L:E\to\mathbb{R}$ 是丛上的光滑函数.
  对于矢量丛的特殊情况 $(\mathbb{R},\mathbb{R}\times\mathbb{R},\pi)$, $L:E\to\mathbb{R},(x,v)\mapsto\frac{1}{2}mv^2-V(x)$ 就是系统的 Lagrangian.

- $E^*$ 是 $E$ 的对偶矢量丛; $E_x$ 是 $x$ 点的 fiber; $L_{E_{x}}:E_{x}\to\mathbb{R}$ 是 $L$ 在 $E_{x}$ 的限制.

- $L$ 的 Legendre 变换是
  $$\begin{aligned}
  \mathrm{F}L:E&\to E^* \\
  v&\mapsto \mathrm{d}\left(L_{E_x}\right)(v)
  \end{aligned}$$

## 9 其它性质

## A 程序

- The Legendre transform is a way to describe a function in terms of it’s supporting hyperplanes.

- Example: Legendre transformation of  $x^2$ where  $x\in(a,b)$ 
  ```Mathematica
  f[x_]:=x^2
  {a,b}={-4,4};
  Lf[k_]:=MaxValue[{k*x-f[x],a<x<b},x]
  Plot[f[x],{x,a,b}]
  Plot[Lf[k],{k,-5,5}]
  ```

- Example: Legendre transformation of  $\mathrm{e}^x$ where  $x\in(a,b)$ , comparing with $k\ln k-k$
  ```Mathematica
  f[x_]:=Exp[x]
  {a,b}={-4,4};
  Lf[k_]:=MaxValue[{k*x-f[x],a<x<b},x]
  Plot[f[x],{x,a,b}]
  Plot[{Lf[k],k*Log[k]-k},{k,0.01,5}]
  ```


<!--stackedit_data:
eyJoaXN0b3J5IjpbLTkyODIwMTYwNywtODc2MDIwOTI1LC0xMj
M5MTcwNzkwLDEwNzkyNDkxMjMsMTcwMTcxNzczMCwtMTM0ODIw
Nzc4NSwtNzM4ODYwMjk4LC0yNTU0MjM3NzEsLTIwMjMxOTI5Nz
AsLTk4ODgxMTgzMSw5MDA5NzY1MThdfQ==
-->