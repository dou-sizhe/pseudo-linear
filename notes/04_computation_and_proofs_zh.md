# 两阶段时域延拓 iLQR：计算细节与证明链

> 2026-09-06 研究范围更新：原问题的控制取值为 $u_k\in\mathbb R^m$，不施加输入范围约束。终端集、CLF 进度条件和局部校正邻域属于证书或算法条件，不是原问题的输入限制。本文保留历史路线的职责；当前修订及证明边界以 [06_review_and_revised_method_zh.md](06_review_and_revised_method_zh.md) 为准。

## 0. 文档范围

本文只整理算法中的计算和理论证明，不讨论研究背景、文献综述和数值实验。
研究对象是固定初值、无折扣、无输入范围限制的离散时间控制仿射非线性系统。全文所有
接受状态均由原始非线性动力学重新 rollout；线性化只用于构造 iLQR 局部模型。

需要始终区分以下五类结论：

1. 终端反馈可行且尾部代价有限；
2. 第一阶段在有限步内获得终端安全前缀；
3. 固定有限时域上的近似 KKT 性；
4. 随时域增长的控制序列收敛与无限时域局部最优性；
5. 全局最优性。

本文在明确附加假设下证明前四类中的相应结论，不证明第五类结论。

---

## 1. 问题、符号与精确 rollout

考虑

$$
x_{k+1}=F(x_k,u_k):=a(x_k)+B(x_k)u_k,
\qquad x_0=\bar x,
$$

其中 $x_k\in\mathbb R^n$，$u_k\in\mathbb R^m$。原问题不设置输入幅值上下界。

阶段代价为

$$
\ell(x,u)=x^\top Qx+u^\top Ru,
\qquad Q\succ0,\quad R\succ0.
$$

目标无限时域问题为

$$
\min_{U=(u_0,u_1,\ldots)}
J_\infty(U)
:=
\sum_{k=0}^{\infty}\ell(x_k,u_k),
$$

并满足原始动力学。优化变量 $u_k$ 不受取值范围约束。

对有限控制序列

$$
U_N=(u_0,\ldots,u_{N-1}),
$$

定义精确 rollout 映射

$$
X_N(U_N)=(x_0,\ldots,x_N),
\qquad
x_{k+1}=a(x_k)+B(x_k)u_k.
$$

后文的有限时域增广目标是

$$
\widehat J_N(U_N)
=
\sum_{k=0}^{N-1}\ell(x_k,u_k)+V_f(x_N).
$$

$V_f$ 是局部终端函数，不预先假设它等于真实无限时域值函数。

### 1.1 基础正则性

计算与终端证书使用以下基础条件：

1. $a(0)=0$，$a$ 和 $B$ 在有关轨迹邻域内二阶连续可微；
2. 它们的一、二阶导数在有关紧集上局部 Lipschitz；
3. 原点线性化

$$
A_0=a_x(0),
\qquad
B_0=B(0)
$$

可稳定；
4. 算法访问的精确 rollout 留在上述光滑区域的某个紧子集内。

因为 $Q\succ0$，线性二次型中的可检测性条件自动满足；结合 $(A_0,B_0)$ 可稳定，
离散代数 Riccati 方程存在稳定化解 $P\succ0$。

---

## 2. 局部终端三元组的计算与证明

### 2.1 DARE、终端增益与闭环恒等式

求解离散代数 Riccati 方程

$$
\begin{aligned}
P={}&Q+A_0^\top PA_0\\
&-A_0^\top PB_0
(R+B_0^\top PB_0)^{-1}
B_0^\top PA_0.
\end{aligned}
$$

定义

$$
G:=R+B_0^\top PB_0\succ0,
$$

$$
K_f:=G^{-1}B_0^\top PA_0,
\qquad
A_c:=A_0-B_0K_f,
$$

以及

$$
S:=Q+K_f^\top RK_f\succ0.
$$

下面逐项验证

$$
A_c^\top PA_c-P=-S.
$$

先展开左端：

$$
\begin{aligned}
A_c^\top PA_c-P
={}&A_0^\top PA_0-P
-A_0^\top PB_0K_f
-K_f^\top B_0^\top PA_0\\
&+K_f^\top B_0^\top PB_0K_f.
\end{aligned}
$$

由 $GK_f=B_0^\top PA_0$，有

$$
A_0^\top PB_0K_f=K_f^\top GK_f,
$$

$$
K_f^\top B_0^\top PA_0=K_f^\top GK_f.
$$

DARE 又给出

$$
A_0^\top PA_0-P=-Q+K_f^\top GK_f.
$$

代回后得到

$$
\begin{aligned}
A_c^\top PA_c-P
={}&-Q-K_f^\top GK_f
+K_f^\top B_0^\top PB_0K_f\\
={}&-Q-K_f^\top RK_f=-S.
\end{aligned}
$$

### 2.2 缩放终端函数

取任意 $\gamma>1$，定义

$$
V_f(x)=\gamma x^\top Px,
\qquad
\kappa_f(x)=-K_fx.
$$

终端集合与带裕量终端集合分别为

$$
\mathcal X_f
=
\{x:V_f(x)\le\alpha_f\},
$$

$$
\mathcal X_f^\delta
=
\{x:V_f(x)\le(1-\delta)\alpha_f\},
\qquad \delta\in(0,1).
$$

$\gamma>1$ 的作用可以从后面的余项展开直接看出：它留下严格的负二次裕量
$-(\gamma-1)x^\top Sx$，用于压住非线性三阶及更高阶项。

### 2.3 非线性闭环余项

终端闭环映射为

$$
F_f(x):=a(x)+B(x)\kappa_f(x)=a(x)-B(x)K_fx.
$$

由二阶光滑性和 $a(0)=0$，在足够小的球内可写成

$$
F_f(x)=A_cx+r(x),
$$

并存在 $c_r>0$ 使

$$
\|r(x)\|\le c_r\|x\|^2.
$$

这里 $B(x)$ 的状态依赖并不会破坏二阶余项，因为

$$
[B(x)-B_0]K_fx=O(\|x\|^2).
$$

计算终端 Lyapunov 差：

$$
\begin{aligned}
V_f(F_f(x))-V_f(x)
={}&\gamma(A_cx+r)^\top P(A_cx+r)-\gamma x^\top Px\\
={}&\gamma x^\top(A_c^\top PA_c-P)x
+2\gamma x^\top A_c^\top Pr
+\gamma r^\top Pr\\
={}&-\gamma x^\top Sx
+2\gamma x^\top A_c^\top Pr
+\gamma r^\top Pr.
\end{aligned}
$$

另一方面

$$
\ell(x,\kappa_f(x))
=x^\top Qx+x^\top K_f^\top RK_fx
=x^\top Sx.
$$

因此

$$
\begin{aligned}
&V_f(F_f(x))-V_f(x)+\ell(x,\kappa_f(x))\\
&=-(\gamma-1)x^\top Sx
+2\gamma x^\top A_c^\top Pr(x)
+\gamma r(x)^\top Pr(x).
\end{aligned}
$$

令

$$
s_{\min}=\lambda_{\min}(S),
\qquad
p_{\max}=\|P\|_2,
\qquad
c_A=\|A_c^\top P\|_2.
$$

利用 Cauchy--Schwarz 不等式和余项界，得到

$$
2\gamma x^\top A_c^\top Pr(x)
\le
2\gamma c_Ac_r\|x\|^3,
$$

$$
\gamma r(x)^\top Pr(x)
\le
\gamma p_{\max}c_r^2\|x\|^4.
$$

所以

$$
\begin{aligned}
&V_f(F_f(x))-V_f(x)+\ell(x,\kappa_f(x))\\
&\le
-\Big[
(\gamma-1)s_{\min}
-2\gamma c_Ac_r\|x\|
-\gamma p_{\max}c_r^2\|x\|^2
\Big]\|x\|^2.
\end{aligned}
$$

选择 $\rho_f>0$ 足够小，使

$$
2\gamma c_Ac_r\rho_f
+\gamma p_{\max}c_r^2\rho_f^2
\le
\frac{\gamma-1}{2}s_{\min},
$$

便有

$$
V_f(F_f(x))-V_f(x)
\le
-\ell(x,\kappa_f(x))
-\frac{\gamma-1}{2}s_{\min}\|x\|^2
$$

对 $\|x\|\le\rho_f$ 成立。忽略最后的额外负项即可得到所需终端下降条件

$$
V_f(F_f(x))-V_f(x)
\le
-\ell(x,\kappa_f(x)).
$$

### 2.4 如何由半径计算 $\alpha_f$

记

$$
p_{\min}=\lambda_{\min}(P)>0.
$$

若 $x\in\mathcal X_f$，则

$$
\gamma p_{\min}\|x\|^2
\le V_f(x)\le\alpha_f,
$$

从而

$$
\|x\|
\le
\sqrt{\frac{\alpha_f}{\gamma p_{\min}}}.
$$

因此只需取

$$
\alpha_f\le\gamma p_{\min}\rho_f^2,
$$

即可保证整个 $\mathcal X_f$ 位于余项界有效的球内。

这里 $\rho_f$ 仅由光滑邻域和非线性下降余项界确定，无需输入裕量。由终端下降和 $\ell\ge0$，

$$
V_f(F_f(x))\le V_f(x)\le\alpha_f,
$$

所以 $F_f(x)\in\mathcal X_f$。至此得到：

$$
\boxed{
x\in\mathcal X_f
\Longrightarrow
\begin{cases}
F_f(x)\in\mathcal X_f,\\
V_f(F_f(x))-V_f(x)\le-\ell(x,\kappa_f(x)).
\end{cases}}
$$

理论上“存在足够小的 $\alpha_f$”已经得到证明；实际给出具体数值时，仍必须用
解析余项界、区间界或其他可靠证书计算 $c_r$ 和有效半径。有限采样只能提供数值
支持，不能替代上述全集合不等式。

---

## 3. 终端尾部的可行性、代价界和实现误差

### 3.1 可实现无限策略

若某个有限前缀满足 $x_N\in\mathcal X_f$，定义实现策略

$$
\Pi_N(U_N)
=(u_0,\ldots,u_{N-1})\oplus\kappa_f,
$$

其中 $N$ 以后每一步都把 $\kappa_f$ 作为状态反馈施加于原始非线性系统。

由终端反馈定义良好和正不变性，整个无限尾部均可行。对任意 $T>N$，逐项求和

$$
V_f(x_{k+1})-V_f(x_k)
\le-\ell(x_k,\kappa_f(x_k))
$$

得到

$$
\sum_{k=N}^{T-1}\ell(x_k,\kappa_f(x_k))
\le
V_f(x_N)-V_f(x_T)
\le V_f(x_N).
$$

令 $T\to\infty$，由非负级数的单调收敛性，

$$
\sum_{k=N}^{\infty}\ell(x_k,\kappa_f(x_k))
\le V_f(x_N).
$$

加上有限前缀代价可得

$$
\boxed{
J_\infty(\Pi_N(U_N))
\le
\widehat J_N(U_N)<\infty.}
$$

又因为 $R\succ0$，

$$
\lambda_{\min}(R)
\sum_{k=0}^{\infty}\|u_k\|^2
\le J_\infty(\Pi_N(U_N)),
$$

故实现控制序列属于 $\ell^2$。

### 3.2 终端二次型与真实尾代价之差

定义终端 Bellman 余项

$$
r_f(x)
:=
\ell(x,\kappa_f(x))
+V_f(F_f(x))-V_f(x).
$$

终端证书给出 $r_f(x)\le0$，但一般并没有 $r_f(x)=0$。沿终端闭环

$$
z_0=x_N,
\qquad
z_{t+1}=F_f(z_t),
$$

第 3.1 节已经证明终端闭环的阶段代价级数有限。由 $Q\succ0$，

$$
\lambda_{\min}(Q)
\sum_{t=0}^{\infty}\|z_t\|^2
\le
\sum_{t=0}^{\infty}\ell(z_t,\kappa_f(z_t))
<\infty,
$$

所以 $z_t\to0$，进而 $V_f(z_t)\to0$。对有限 $T$ 有精确恒等式

$$
\sum_{t=0}^{T-1}r_f(z_t)
=
\sum_{t=0}^{T-1}\ell(z_t,\kappa_f(z_t))
+V_f(z_T)-V_f(z_0).
$$

令 $T\to\infty$，所以

$$
J_{\kappa_f}(x_N)-V_f(x_N)
=
\sum_{t=0}^{\infty}r_f(z_t).
$$

因此对整条实现策略，

$$
\boxed{
J_\infty(\Pi_N(U_N))-\widehat J_N(U_N)
=
\sum_{t=0}^{\infty}r_f(F_f^t(x_N)).}
$$

若还可证明存在 $C_r>0$、$p_f\ge2$、$C_z>0$ 和 $\sigma\in(0,1)$ 使

$$
|r_f(x)|\le C_r\|x\|^{p_f},
$$

$$
\|F_f^t(x)\|\le C_z\sigma^t\|x\|,
$$

则

$$
\begin{aligned}
|J_{\kappa_f}(x_N)-V_f(x_N)|
&\le
C_r\sum_{t=0}^{\infty}
\|F_f^t(x_N)\|^{p_f}\\
&\le
\frac{C_rC_z^{p_f}}{1-\sigma^{p_f}}
\|x_N\|^{p_f}.
\end{aligned}
$$

这个误差随时域消失还需要移动终端状态满足 $x_N\to0$。仅有
$x_N\in\mathcal X_f$ 并不能推出这一点。

---

## 4. 有限代价域与有限步终端可达

定义有限代价域

$$
\mathcal D_{\rm fin}
=
\left\{
x:\exists U,\ J_\infty(x;U)<\infty
\right\}.
$$

在终端三元组已经成立时，对固定初值 $\bar x$，下列两件事等价：

1. $\bar x\in\mathcal D_{\rm fin}$；
2. 存在有限 $T$ 和可行前缀 $U_T$，使 $x_T\in\mathcal X_f^\delta$。

证明第一方向。若 $J_\infty(\bar x;U)<\infty$，则

$$
\lambda_{\min}(Q)
\sum_{k=0}^{\infty}\|x_k\|^2
\le
\sum_{k=0}^{\infty}x_k^\top Qx_k
\le J_\infty(\bar x;U)<\infty.
$$

故 $\sum_k\|x_k\|^2<\infty$，从而 $x_k\to0$。因为
$V_f(x)=\gamma x^\top Px$ 连续，$V_f(x_k)\to0$，所以必存在有限 $T$ 使

$$
V_f(x_T)\le(1-\delta)\alpha_f.
$$

证明反方向。若有限前缀进入 $\mathcal X_f^\delta\subset\mathcal X_f$，在第 $T$
步后接入 $\kappa_f$，则

$$
J_\infty(\Pi_T(U_T))
\le
\sum_{k=0}^{T-1}\ell(x_k,u_k)+V_f(x_T)<\infty.
$$

等价性只说明某条有限可达前缀存在，不说明局部 iLQR 能从任意初始化找到它。

---

## 5. 第一阶段一步问题的精确计算

### 5.1 严格凸二次结构

对固定状态 $x$，一步增广目标为

$$
\phi_x(u)
:=
\ell(x,u)+V_f(a(x)+B(x)u).
$$

把二次项完全展开：

$$
\begin{aligned}
\phi_x(u)
={}&x^\top Qx+u^\top Ru\\
&+\gamma[a(x)+B(x)u]^\top
P[a(x)+B(x)u]\\
={}&x^\top Qx+\gamma a(x)^\top Pa(x)\\
&+2\gamma a(x)^\top PB(x)u\\
&+u^\top[R+\gamma B(x)^\top PB(x)]u.
\end{aligned}
$$

定义

$$
H(x):=R+\gamma B(x)^\top PB(x),
$$

$$
h(x):=\gamma B(x)^\top Pa(x),
$$

$$
c(x):=x^\top Qx+\gamma a(x)^\top Pa(x).
$$

于是

$$
\phi_x(u)=c(x)+2h(x)^\top u+u^\top H(x)u.
$$

因为 $R\succ0$ 且 $B^\top PB\succeq0$，

$$
H(x)\succeq R\succ0.
$$

所以该一步问题是无约束严格凸二次问题，具有唯一解。

### 5.2 配方与解析解

由

$$
\phi_x(u)=c(x)-h(x)^\top H(x)^{-1}h(x)
+[u+H(x)^{-1}h(x)]^\top H(x)[u+H(x)^{-1}h(x)],
$$

以及 $H(x)\succ0$，唯一最小点为

$$
u^g(x)=-H(x)^{-1}h(x).
$$

实际计算用 Cholesky 分解求解 $H(x)u^g=-h(x)$，不显式求逆。

### 5.3 一步最优性验证

无约束一阶条件为

$$
\nabla_u\phi_x(u)=2H(x)u+2h(x)=0.
$$

对任意 $u$，配方恒等式给出

$$
\phi_x(u)-\phi_x(u^g)=(u-u^g)^\top H(x)(u-u^g)\ge0,
$$

且等号仅在 $u=u^g$ 成立。这直接验证唯一全局最优性，不需要输入边界或活动集乘子。

### 5.4 一步延拓缺陷恒等式

在当前精确终端状态 $x_N$ 上取

$$
u_N^g=u^g(x_N),
\qquad
x_{N+1}=a(x_N)+B(x_N)u_N^g.
$$

定义一步延拓缺陷

$$
\Delta_N^g
:=
\ell(x_N,u_N^g)+V_f(x_{N+1})-V_f(x_N).
$$

则

$$
\begin{aligned}
&\widehat J_{N+1}(U_N\oplus u_N^g)\\
&=\sum_{k=0}^{N-1}\ell(x_k,u_k)
+\ell(x_N,u_N^g)+V_f(x_{N+1})\\
&=\widehat J_N(U_N)+\Delta_N^g.
\end{aligned}
$$

这是精确恒等式，不是局部近似。

若 $\kappa_f(x_N)$ 可行，由 $u_N^g$ 的全局一步最优性，

$$
\Delta_N^g
\le
\ell(x_N,\kappa_f(x_N))
+V_f(F_f(x_N))-V_f(x_N).
$$

当 $x_N\in\mathcal X_f$ 时，右端不大于零；当 $x_N\notin\mathcal X_f$ 时，
这个右端未必非正。因此“一步贪心不劣于终端反馈”不能被改写成全局 Lyapunov
下降。

### 5.5 累计正缺陷与条件性有限命中

设第一阶段访问一列时域 $N_j$，所有扩展均逐步使用上述精确一步解，所有被接受的
iLQR 修正满足增广代价不增加。定义

$$
D_j^+
:=
\sum_{k\in\mathcal E_j}[\Delta_k^g]_+,
\qquad
[z]_+:=\max\{z,0\},
$$

其中 $\mathcal E_j$ 是截至 $N_j$ 已执行的所有一步延拓索引。由一步恒等式，负缺陷
只会降低目标，而 correction 也不增目标，所以

$$
\widehat J_{N_j}(U_{N_j})
\le
\widehat J_1(U_1)+D_j^+.
$$

若

$$
D_j^+=o(N_j),
$$

则

$$
\frac{\widehat J_{N_j}(U_{N_j})}{N_j}\to0.
$$

又因为 $V_f\ge0$ 且 $u_k^\top Ru_k\ge0$，

$$
\widehat J_{N_j}(U_{N_j})
\ge
\sum_{k=0}^{N_j-1}x_k^\top Qx_k
\ge
\lambda_{\min}(Q)
\sum_{k=0}^{N_j-1}\|x_k\|^2.
$$

因此

$$
\min_{0\le k<N_j}\|x_k\|^2
\le
\frac{1}{N_j}
\sum_{k=0}^{N_j-1}\|x_k\|^2
\le
\frac{\widehat J_{N_j}(U_{N_j})}
{N_j\lambda_{\min}(Q)}
\to0.
$$

$\mathcal X_f^\delta$ 包含原点的一个邻域，所以对足够大的有限 $j$，当前精确
rollout 中至少有一个状态落入 $\mathcal X_f^\delta$。若算法在每次扩展和接受修正
后扫描整条轨迹并在首次命中处截断，则

$$
D_j^+=o(N_j)
\Longrightarrow
N_{\rm hit}<\infty.
$$

这是性能通道的充分条件，不是一步 QP 自动保证的性质，也不说明末端序列
$x_{N_j}$ 本身收敛。

---

## 6. 固定时域 iLQR 的逐项计算

### 6.1 名义轨迹与精确 Jacobian

给定名义控制

$$
\bar U_N=(\bar u_0,\ldots,\bar u_{N-1}),
$$

先按原始动力学计算

$$
\bar x_0=\bar x,
\qquad
\bar x_{k+1}=a(\bar x_k)+B(\bar x_k)\bar u_k.
$$

记 $B_{\cdot j}(x)$ 为 $B(x)$ 的第 $j$ 列。因为

$$
F(x,u)=a(x)+\sum_{j=1}^mB_{\cdot j}(x)u_j,
$$

所以名义点上的精确一阶 Jacobian 为

$$
A_k
=F_x(\bar x_k,\bar u_k)
=a_x(\bar x_k)
+\sum_{j=1}^m\bar u_{k,j}[B_{\cdot j}]_x(\bar x_k),
$$

$$
B_k
=F_u(\bar x_k,\bar u_k)
=B(\bar x_k).
$$

扰动动力学为

$$
\delta x_{k+1}=A_k\delta x_k+B_k\delta u_k,
\qquad
\delta x_0=0.
$$

### 6.2 阶段代价与终端代价的导数

二次阶段代价给出

$$
q_k:=\ell_x(\bar x_k,\bar u_k)=2Q\bar x_k,
$$

$$
r_k:=\ell_u(\bar x_k,\bar u_k)=2R\bar u_k,
$$

$$
Q_k:=\ell_{xx}=2Q,
\qquad
R_k:=\ell_{uu}=2R,
\qquad
S_k:=\ell_{ux}=0.
$$

终端初始化为

$$
V_{x,N}=2\gamma P\bar x_N,
\qquad
V_{xx,N}=2\gamma P.
$$

### 6.3 从下游价值展开到局部 $Q$ 函数

假设下游价值变化近似为

$$
\delta V_{k+1}
=
V_{x,k+1}^\top\delta x_{k+1}
+\frac12\delta x_{k+1}^\top
V_{xx,k+1}\delta x_{k+1}.
$$

把 $\delta x_{k+1}=A_k\delta x_k+B_k\delta u_k$ 代入线性项：

$$
V_{x,k+1}^\top\delta x_{k+1}
=
(A_k^\top V_{x,k+1})^\top\delta x_k
+(B_k^\top V_{x,k+1})^\top\delta u_k.
$$

再展开二次项：

$$
\begin{aligned}
&\frac12(A_k\delta x_k+B_k\delta u_k)^\top
V_{xx,k+1}(A_k\delta x_k+B_k\delta u_k)\\
&=
\frac12\delta x_k^\top A_k^\top V_{xx,k+1}A_k\delta x_k\\
&\quad+
\delta u_k^\top B_k^\top V_{xx,k+1}A_k\delta x_k\\
&\quad+
\frac12\delta u_k^\top B_k^\top V_{xx,k+1}B_k\delta u_k.
\end{aligned}
$$

与阶段代价的二次展开合并，得到

$$
Q_{x,k}=q_k+A_k^\top V_{x,k+1},
$$

$$
Q_{u,k}=r_k+B_k^\top V_{x,k+1},
$$

$$
Q_{xx,k}=Q_k+A_k^\top V_{xx,k+1}A_k,
$$

$$
Q_{uu,k}=R_k+B_k^\top V_{xx,k+1}B_k,
$$

$$
Q_{ux,k}=S_k+B_k^\top V_{xx,k+1}A_k.
$$

因此局部二次模型是

$$
\begin{aligned}
\delta Q_k
={}&Q_{x,k}^\top\delta x_k
+Q_{u,k}^\top\delta u_k\\
&+\frac12\delta x_k^\top Q_{xx,k}\delta x_k
+\delta u_k^\top Q_{ux,k}\delta x_k\\
&+\frac12\delta u_k^\top Q_{uu,k}\delta u_k.
\end{aligned}
$$

iLQR 在这里保留动力学一阶项和价值二次项，舍去价值梯度与 $F_{xx}$、$F_{xu}$
的二阶张量收缩。控制仿射性只给出 $F_{uu}=0$；一般仍有 $F_{xx}\ne0$ 和
$F_{xu}\ne0$。保留这些二阶动力学项得到的是 full DDP，而不是当前 iLQR。

### 6.4 无约束正则化反向步

取 Levenberg--Marquardt 正则参数 $\lambda\ge0$，令

$$
H_k:=Q_{uu,k}+\lambda I.
$$

若 $H_k$ 非正定，则增大 $\lambda$ 并重做反向递推。在 $H_k\succ0$ 时，局部控制增量满足

$$
H_k\delta u_k+Q_{u,k}+Q_{ux,k}\delta x_k=0.
$$

因此

$$
d_k=-H_k^{-1}Q_{u,k},\qquad L_k=-H_k^{-1}Q_{ux,k},
\qquad \delta u_k=d_k+L_k\delta x_k.
$$

实现中解正定线性方程，无需活动集、控制边界或投影。

### 6.5 价值递推的完整代入

代入局部策略

$$
\delta u_k=d_k+L_k\delta x_k.
$$

常数项为

$$
\Delta V_{0,k}
=Q_{u,k}^\top d_k+\frac12d_k^\top H_kd_k.
$$

关于 $\delta x_k$ 的一次项来自四处：

$$
Q_{x,k}^\top\delta x_k,
$$

$$
d_k^\top Q_{ux,k}\delta x_k,
$$

$$
\delta x_k^\top L_k^\top Q_{u,k},
$$

$$
\delta x_k^\top L_k^\top H_kd_k.
$$

因此

$$
V_{x,k}
=
Q_{x,k}
+Q_{ux,k}^\top d_k
+L_k^\top Q_{u,k}
+L_k^\top H_kd_k.
$$

二次项合并为

$$
V_{xx,k}
=
Q_{xx,k}
+L_k^\top H_kL_k
+L_k^\top Q_{ux,k}
+Q_{ux,k}^\top L_k.
$$

数值上使用

$$
V_{xx,k}\leftarrow
\frac12(V_{xx,k}+V_{xx,k}^\top)
$$

消除浮点误差造成的非对称部分。

### 6.6 预测下降、真实下降和接受比

反向递推累计

$$
\Delta_1
=
\sum_{k=0}^{N-1}Q_{u,k}^\top d_k,
$$

$$
\Delta_2
=
\sum_{k=0}^{N-1}d_k^\top H_kd_k.
$$

对步长 $\alpha\in(0,1]$，局部模型预测下降为

$$
\Delta m_N(\alpha)
=
-\alpha\Delta_1
-\frac12\alpha^2\Delta_2.
$$

trial 控制和状态由

$$
u_k^+(\alpha)
=
\bar u_k+\alpha d_k
+L_k(x_k^+(\alpha)-\bar x_k),
$$

$$
x_{k+1}^+(\alpha)
=
a(x_k^+(\alpha))+B(x_k^+(\alpha))u_k^+(\alpha),
\qquad
x_0^+(\alpha)=\bar x
$$

逐步计算。这是原始非线性 rollout，不是线性化状态递推。

真实下降与比值为

$$
\Delta J_N(\alpha)
=
\widehat J_N(\bar U_N)-\widehat J_N(U_N^+(\alpha)),
$$

$$
\rho_N(\alpha)
=
\frac{\Delta J_N(\alpha)}{\Delta m_N(\alpha)}.
$$

只有同时满足

$$
\Delta m_N(\alpha)>0,
\qquad
\Delta J_N(\alpha)>0,
\qquad
\rho_N(\alpha)\ge\eta_\rho
$$

时才接受。进入认证阶段后还必须满足

$$
x_N^+(\alpha)\in\mathcal X_f.
$$

若所有正则化和步长均失败，保持上一条可行轨迹并报告 corrector 失败，不能把失败
解释成收敛。

---

## 7. 精确有限时域 KKT 诊断

### 7.1 伴随梯度的推导

对精确 rollout，构造拉格朗日函数

$$
\begin{aligned}
\mathcal L
={}&
\sum_{k=0}^{N-1}\ell(x_k,u_k)+V_f(x_N)\\
&+\sum_{k=0}^{N-1}
p_{k+1}^\top[F(x_k,u_k)-x_{k+1}].
\end{aligned}
$$

对 $x_N$ 求导得

$$
p_N=\nabla V_f(x_N)=2\gamma Px_N.
$$

对 $x_k$ 求导并令其为零，得到反向伴随方程

$$
p_k
=
2Qx_k+F_x(x_k,u_k)^\top p_{k+1},
$$

即

$$
p_k
=
2Qx_k
+\left[
a_x(x_k)
+\sum_{j=1}^mu_{k,j}[B_{\cdot j}]_x(x_k)
\right]^\top p_{k+1}.
$$

对控制求导得到精确约化梯度块

$$
g_{u,k}
=
2Ru_k+B(x_k)^\top p_{k+1}.
$$

堆叠后

$$
g_N=(g_{u,0},\ldots,g_{u,N-1})
=
\nabla\widehat J_N(U_N).
$$

该梯度来自原始非线性系统的一阶导数，不是冻结 LTV 问题的替代梯度。

### 7.2 原问题的真实梯度残差

控制空间为 $\mathbb R^{mN}$，故定义

$$
\mathcal G_N(U_N):=g_N(U_N)=\nabla\widehat J_N(U_N).
$$

原有限时域问题的一阶驻点条件直接是 $\mathcal G_N=0$；实际精度门检查 $\|\mathcal G_N\|$。这里无需投影梯度、控制边界乘子或输入法锥。

终端成员条件是另一个问题：若显式要求 $c_N(U_N)=V_f(x_N(U_N))-\alpha_f\le0$，其 KKT 系统仍须包含

$$
g_N+\nu_N\nabla c_N=0,\quad c_N\le0,\quad \nu_N\ge0,\quad \nu_Nc_N=0.
$$

只有终端条件不活动或其乘子项为零时，该系统才退化为原问题的梯度条件。严格内部的起点不能保证后续边界不活动，详见第 06 号修订稿。

### 7.3 Riccati 预条件残差

在 $\lambda=0$ 的未正则化局部模型上重新做一次反向递推，定义

$$
\chi_N(U_N)
:=
\max_{0\le k<N}\|d_k^{(\lambda=0)}\|_\infty.
$$

若所有所需控制 Hessian 块正定，且每个阶段无约束二次问题精确求解，则

$$
\chi_N(U_N)=0
\Longleftrightarrow
\mathcal G_N(U_N)=0.
$$

证明可用反向归纳。若所有 $d_k=0$，价值梯度递推退化为真实伴随递推，阶段 QP
在 $d_k=0$ 处的最优性条件正是

$$
g_{u,k}=0.
$$

反之，若整个控制序列满足 KKT，则从 $k=N-1$ 开始，未正则化阶段 QP 在
$d_k=0$ 处满足最优性条件；价值梯度因此与伴随一致，逐步向前归纳得到全部
$d_k=0$。

但在残差非零时，$\chi_N$ 只是 Riccati 预条件量。若没有时域一致的曲率上下界和
反向递推的统一有界性，不能自动得到

$$
c_1\|\mathcal G_N\|
\le\chi_N
\le c_2\|\mathcal G_N\|
$$

且常数与 $N$ 无关。因此 continuation gate 可以使用 $\chi_N$ 控制 corrector
工作量，但有限时域一阶驻点性必须同时用 $\mathcal G_N$ 审计。

---

## 8. 第一阶段：见证保护的终端安全前缀获取

### 8.1 可达见证假设

假设给定有限 $T_b$ 和候选控制 $U_{T_b}^b$。算法必须重新用原始非线性动力学
rollout，并逐项验证：

$$
u_k^b\in\mathbb R^m\text{ 且数值有限},
$$

$$
x_{k+1}^b=a(x_k^b)+B(x_k^b)u_k^b,
$$

状态均有限，且存在首次命中索引

$$
T_b
=
\min\{k:x_k^b\in\mathcal X_f^\delta\}<\infty.
$$

验证后在首次命中处截断并独立保存。见证只证明有限可达，不含任何最优性结论。

### 8.2 双通道计算

性能通道从 $N=1$ 的严格凸一步 QP 开始。每当时域从 $N$ 增加到 $M$ 时，对每个
新状态重复计算

$$
u_k^g
=
-H(x_k)^{-1}h(x_k),
$$

并以原始动力学得到 $x_{k+1}$。随后在当前时域最多尝试有限次 ratio-globalized
iLQR correction。任何 correction 失败都保留原轨迹；任何接受轨迹都重新扫描

$$
T_\delta(X_N)
=
\min\{0\le k\le N:x_k\in\mathcal X_f^\delta\}.
$$

证书通道始终保存已验证见证，不允许性能通道覆盖。

### 8.3 有限终止定理

**定理 1（第一阶段有限终止）** 设见证验证、时域调度和每个时域的 correction
次数都有有限上界，并且调度最终允许到达或回退到 $T_b$。则第一阶段在有限工作量
后返回某个

$$
N_{\rm hit}\le T_b
$$

及精确可行前缀，满足

$$
x_{N_{\rm hit}}\in\mathcal X_f^\delta.
$$

**证明。** 若性能通道先命中，则被返回的控制与状态数值有限，状态来自精确
rollout，终端成员关系已经显式检查。若性能通道在有限预算内没有命中，则返回
独立保存且已验证的见证前缀。见证长度有限，每个 correction 块有限，时域调度
有限，因此总工作量有限。整个证明不需要 iLQR 收敛。$\square$

由第 3 节的尾部证明，第一阶段返回后立即得到一条可行、有限无限时域代价的回退
策略。此时尚未证明有限时域 KKT 性。

### 8.4 第一阶段工作量

验证长度为 $T_b$ 的见证需要 $O(T_b)$ 次动力学计算。若性能通道访问

$$
1=N_0<N_1<\cdots<N_s\le T_b
$$

且每个时域最多做 $K_{\rm I}$ 次 Riccati sweep，则 corrector 工作量为

$$
\sum_{j=0}^sO(K_{\rm I}N_j).
$$

逐点调度 $N_j=j$ 时，

$$
\sum_{j=1}^{T_b}O(K_{\rm I}j)
=O(K_{\rm I}T_b^2).
$$

几何调度 $N_{j+1}\ge rN_j$、$r>1$ 时，几何级数给出

$$
\sum_jO(K_{\rm I}N_j)
=O(K_{\rm I}T_b).
$$

逐步一步 QP 总数不超过 $T_b$，因此这部分关于新增时刻是线性的。

---

## 9. 第二阶段：残差强迫的认证时域延拓

### 9.1 首次认证与强迫序列

选取

$$
\tau_N=c_\tau(N+1)^{-p},
\qquad c_\tau>0,\quad p>0.
$$

在时域 $N$ 上，认证门为

$$
x_N\in\mathcal X_f,
\qquad
\chi_N(U_N)\le\tau_N.
$$

第一阶段只保证

$$
x_{N_{\rm hit}}\in\mathcal X_f^\delta
\subset\operatorname{int}\mathcal X_f.
$$

第二阶段先固定 $N=N_{\rm hit}$，运行保持终端成员关系的安全 corrector，直到首次
通过残差门。该时域记为 $N_{\rm cert}$。数值上可能有
$N_{\rm cert}=N_{\rm hit}$，但逻辑上二者分别表示“已有安全尾部”和“安全尾部加
有限时域残差门”。

### 9.2 时域一致局部 corrector 假设

假设每个第二阶段时域 $N$ 上存在局部集合

$$
\mathcal B_N\subset\mathbb R^{mN}
$$

和与 $N$ 无关的 $\theta\in(0,1)$，满足：

1. 从第一阶段输出出发，有限次被接受的终端安全 correction 能进入
   $\mathcal B_{N_{\rm hit}}$；
2. 若 $U_N\in\mathcal B_N$、$x_N\in\mathcal X_f$ 且 $\chi_N(U_N)>0$，则存在被
   接受的安全 correction $U_N^+$，使

$$
U_N^+\in\mathcal B_N,
\qquad
x_N(U_N^+)\in\mathcal X_f,
$$

$$
\widehat J_N(U_N^+)<\widehat J_N(U_N),
$$

$$
\chi_N(U_N^+)\le\theta\chi_N(U_N).
$$

这是第二阶段的新假设，不由终端 Lyapunov 下降、局部光滑性或可稳定性自动推出。
它概括了统一正曲率、反向递推的统一有界性、导数 Lipschitz 性和统一吸引域等局部求解器
条件。

### 9.3 首次认证的有限步证明

**命题 2（首次认证有限完成）** 在上述 corrector 假设下，第一阶段安全前缀经过
有限次被接受 correction 后通过认证门。

**证明。** 先用假设的第一部分在有限步内进入 $\mathcal B_{N_{\rm hit}}$。设进入
该集合时的残差为 $\chi^0$。重复使用收缩不等式得到

$$
\chi^i\le\theta^i\chi^0.
$$

若 $\chi^0\le\tau_{N_{\rm hit}}$，不再需要 correction。否则取

$$
i
\ge
\left\lceil
\frac{
\log(\chi^0/\tau_{N_{\rm hit}})
}{
\log(1/\theta)
}
\right\rceil,
$$

便有 $\theta^i\chi^0\le\tau_{N_{\rm hit}}$。corrector 不变性在全过程保持
$x_N\in\mathcal X_f$，故认证门在有限步后成立。$\square$

### 9.4 精确终端反馈延拓

从安全前缀 $U_N$ 延拓到 $M>N$，定义

$$
\mathcal E_N^MU_N
=
(u_0,\ldots,u_{N-1},
\kappa_f(x_N),\ldots,\kappa_f(x_{M-1})),
$$

其中所有新增状态都按

$$
x_{k+1}=F_f(x_k)
=a(x_k)+B(x_k)\kappa_f(x_k)
$$

真实 rollout。由于 $\mathcal X_f$ 正不变，新增控制和状态均可行。

计算延拓前后增广代价之差：

$$
\begin{aligned}
&\widehat J_M(\mathcal E_N^MU_N)-\widehat J_N(U_N)\\
&=
\sum_{k=N}^{M-1}\ell(x_k,\kappa_f(x_k))
+V_f(x_M)-V_f(x_N)\\
&=
\sum_{k=N}^{M-1}
\left[
\ell(x_k,\kappa_f(x_k))
+V_f(x_{k+1})-V_f(x_k)
\right]\\
&=
\sum_{k=N}^{M-1}r_f(x_k)
\le0.
\end{aligned}
$$

因此

$$
\boxed{
\widehat J_M(\mathcal E_N^MU_N)
\le
\widehat J_N(U_N).}
$$

该式只说明安全性和增广代价非增，不说明延拓后的 $M$ 时域控制已经满足 KKT。

### 9.5 候选时域与实测工作量

在认证时域 $N$ 上设置候选集合

$$
\mathcal H(N)
=
\{N+1\}
\cup
\{\lceil r_1N\rceil,\ldots,\lceil r_sN\rceil\},
\qquad
1<r_i\le r_{\max}.
$$

对每个候选 $M$，实际构造

$$
W_M=\mathcal E_N^MU_N,
$$

重新 rollout，并计算真实的 $\chi_M(W_M)$。给定保守经验收缩率
$\widehat\theta\in(0,1)$，预测达到门槛所需 correction 数为

$$
\widehat K(M)
=
\left\lceil
\frac{
\log[\chi_M(W_M)/\tau_M]
}{
\log(1/\widehat\theta)
}
\right\rceil_+,
$$

其中当 $\chi_M(W_M)\le\tau_M$ 时定义为零。等价地，候选满足目标工作预算
$K_{\rm tar}$ 当且仅当

$$
\widehat\theta^{K_{\rm tar}}
\chi_M(W_M)
\le\tau_M.
$$

选择通过预算的最大候选；若没有几何候选通过，仍保留 $M=N+1$ 作为必选候选。
预测值只用于选择时域，不能替代 corrector 后重新检查真实认证门。

---

## 10. 延拓正则性与一致 correction 次数

### 10.1 强迫归一化延拓假设

终端 Lyapunov 下降不控制延拓后的驻点残差。为此单独假设存在与 $N,M$ 无关的
$D_E<\infty$，使对每个已认证检查点和每个 $M\in\mathcal H(N)$，

$$
W_M=\mathcal E_N^MU_N\in\mathcal B_M,
$$

$$
\chi_M(W_M)\le D_E\tau_M.
$$

第一式保证 warm start 落入新时域的统一局部吸引域；第二式只要求延拓残差与新
强迫容差的比值一致有界。

### 10.2 旧的指数加性界为何是充分条件

一个更强但不是主定理所必需的条件是

$$
\chi_M(W_M)
\le
\chi_N(U_N)+C_e\beta^N,
\qquad
C_e>0,\quad\beta\in(0,1).
$$

在认证点 $\chi_N(U_N)\le\tau_N$，所以

$$
\begin{aligned}
\frac{\chi_M(W_M)}{\tau_M}
&\le
\frac{\tau_N}{\tau_M}
+\frac{C_e\beta^N}{\tau_M}\\
&=
\left(\frac{M+1}{N+1}\right)^p
+\frac{C_e}{c_\tau}(M+1)^p\beta^N.
\end{aligned}
$$

候选满足 $M\le\lceil r_{\max}N\rceil$，故第一项一致有界。第二项是“多项式乘
指数衰减”；序列 $N^p\beta^N$ 有有限上确界。因此整个右端存在有限一致上界，
可取为 $D_E$。这证明指数加性界蕴含强迫归一化延拓界。

反方向一般不成立：归一化界只约束算法实际查询的检查点和候选，不要求给出指数
衰减率。

### 10.3 一致 correction 上界

**定理 3（一致 correction 次数）** 在时域一致 corrector 假设和强迫归一化延拓
假设下，首次认证后从任一已认证 $N$ 延拓到任一候选 $M$，最多需要

$$
K_E
:=
\left\lceil
\frac{\log D_E}{\log(1/\theta)}
\right\rceil_+
$$

次被接受 correction 即可通过新门槛。

**证明。** 延拓后

$$
\chi_M(U_M^0)\le D_E\tau_M.
$$

连续应用收缩不等式 $K$ 次，

$$
\chi_M(U_M^K)
\le
\theta^K\chi_M(U_M^0)
\le
\theta^KD_E\tau_M.
$$

$K=K_E$ 时 $\theta^{K_E}D_E\le1$，故

$$
\chi_M(U_M^{K_E})\le\tau_M.
$$

延拓位于 $\mathcal B_M$，corrector 又保持该集合和终端成员关系，所以每一步均可
继续，且新检查点安全。$\square$

每次完成一个时域块后都有

$$
N_{j+1}\ge N_j+1.
$$

因此理想的、没有软件上限 $N_{\max}$ 的算法产生严格递增整数序列，必有

$$
N_j\to\infty.
$$

达到有限的 $N_{\max}$ 只表示计算停止条件，不等于上述理想极限已经在数学上达到。

### 10.4 几何调度的线性 sweep 复杂度

再假设：

1. $\widehat\theta\ge\theta$，即预测收缩不比真实收缩更乐观；
2. 候选集中包含固定 $r>1$ 的 $M_N=\lceil rN\rceil$；
3. 目标工作量满足

$$
\widehat\theta^{K_{\rm tar}}D_E\le1.
$$

则对任一候选，

$$
\widehat\theta^{K_{\rm tar}}\chi_M(W_M)
\le
\widehat\theta^{K_{\rm tar}}D_E\tau_M
\le\tau_M,
$$

所以所有候选都通过预测预算，算法至少选择 $\lceil rN\rceil$。真实收缩满足
$\theta\le\widehat\theta$，故实际也不超过 $K_{\rm tar}$ 次 correction。

到达有限上限 $N_{\max}$ 前，认证时域至少几何增长，所以阶段数为
$O(\log N_{\max})$。每次 Riccati sweep 的时标工作量为 $O(N_j)$，从而

$$
\sum_jO(K_{\rm tar}N_j)
=
O(K_{\rm tar}N_{\max}).
$$

若始终采用 $N\mapsto N+1$，即使每个时域只需常数次 sweep，也会有

$$
\sum_{N=1}^{N_{\max}}O(N)
=O(N_{\max}^2).
$$

---

## 11. 被接受增广代价的收敛

第二阶段有两种会更新当前检查点的操作：

1. 终端反馈延拓，其增广代价不增加；
2. 被 ratio test 接受的 correction，其真实增广代价严格下降。

因此，把所有被保存的第二阶段增广代价按实际发生顺序记为 $c_0,c_1,\ldots$，有

$$
0\le c_{j+1}\le c_j.
$$

单调有下界序列必收敛，故存在 $c_\infty\ge0$ 使

$$
c_j\to c_\infty.
$$

该结论只涉及一列标量，不能单独推出控制序列收敛。控制收敛需要下一节的分支
假设。

---

## 12. 从有限时域检查点到无限时域控制极限

### 12.1 无限控制空间与嵌入

定义控制 Hilbert 空间

$$
\mathscr U
=
\ell^2(\mathbb R^m).
$$

对安全有限前缀 $U_N$，用终端反馈生成其真实无限尾部，并把得到的无限开环控制
序列记为

$$
\mathcal I_NU_N\in\mathscr U.
$$

它不是简单的零填充，而是固定初值下由 $\kappa_f$ 与原始非线性闭环实际产生的
控制尾部。第 3 节已经证明它属于 $\ell^2$。

记 $\mathcal P_T$ 为提取前 $T$ 个控制块的坐标投影，其作用与输入限幅无关。无输入范围约束时，$T_{\mathscr U}(U)=\ell^2$。

**历史证明边界。** $U\in\ell^2$ 不能自动保证原系统有有限代价，也不能保证 $J_\infty$ 在开环 $\ell^2$ 邻域可微，尤其对开环不稳定系统。下面的分支与光滑性结论只在显式附加假设成立时有效；不能把取消输入范围约束当作这些假设的证明。当前修订采用稳定反馈坐标，见第 06 号文档。

### 12.2 选定局部分支假设

为识别极限并证明局部最优，需要显式假设存在：

- 足够大的起始时域 $N_0$；
- 每个有关时域上的选定有限 KKT 点 $U_N^\star$；
- 极限 $U_\infty^\star\in\mathscr U$；
- 非降函数 $\omega:[0,\infty)\to[0,\infty)$，且
  $\omega(s)\to0$ 当 $s\downarrow0$；
- 沿访问时域的 $\varepsilon_N\to0$；
- 常数 $\mu>0$。

它们满足以下五项。

**分支 KKT：**

$$
\nabla\widehat J_N(U_N^\star)=0,
$$

且 $x_N(U_N^\star)\in\mathcal X_f$。

**残差误差界：** 对第二阶段局部集合中的任意 $U_N$，

$$
\|\mathcal I_NU_N-\mathcal I_NU_N^\star\|_{\ell^2}
\le
\omega(\chi_N(U_N)).
$$

这里允许任意趋零模函数 $\omega$，不要求线性误差界。

**分支收敛：**

$$
\|\mathcal I_NU_N^\star-U_\infty^\star\|_{\ell^2}
\le\varepsilon_N.
$$

**固定前缀梯度一致性：** $J_\infty$ 在分支附近连续 Fr\'echet 可微，且对每个
固定 $T$，

$$
\left\|
\mathcal P_T\nabla J_\infty(\mathcal I_NU_N^\star)
-
\mathcal P_T\nabla\widehat J_N(U_N^\star)
\right\|
\to0.
$$

**可行切向曲率：** 在 $U_\infty^\star$ 的一个开凸邻域 $\mathcal O$ 上，
$J_\infty$ 二阶连续 Fr\'echet 可微，并且

$$
\langle v,\nabla^2J_\infty(U)v\rangle
\ge
\mu\|v\|_{\ell^2}^2
$$

对所有 $U\in\mathcal O$ 和
$v\in T_{\mathscr U}(U_\infty^\star)$ 成立。

前两阶段的光滑性、可达见证和终端证书都不能自动推出这五项；它们是从有限局部
解分支过渡到原无限维问题所需的独立接口。

### 12.3 实现策略的强收敛

**定理 4（实现策略强收敛）** 在前述 corrector、延拓和分支假设下，认证检查点
$U_{N_j}^j$ 满足

$$
\|\mathcal I_{N_j}U_{N_j}^j-U_\infty^\star\|_{\ell^2}
\le
\omega(\tau_{N_j})+\varepsilon_{N_j}
\to0.
$$

**证明。** 插入同一时域的分支点并使用三角不等式：

$$
\begin{aligned}
&\|\mathcal I_{N_j}U_{N_j}^j-U_\infty^\star\|_{\ell^2}\\
&\le
\|\mathcal I_{N_j}U_{N_j}^j
-\mathcal I_{N_j}U_{N_j}^\star\|_{\ell^2}\\
&\quad+
\|\mathcal I_{N_j}U_{N_j}^\star
-U_\infty^\star\|_{\ell^2}\\
&\le
\omega(\chi_{N_j}(U_{N_j}^j))
+\varepsilon_{N_j}\\
&\le
\omega(\tau_{N_j})+\varepsilon_{N_j}.
\end{aligned}
$$

因为 $N_j\to\infty$，有 $\tau_{N_j}\to0$；再用 $\omega(s)\to0$ 和
$\varepsilon_{N_j}\to0$，右端趋于零。

对一次延拓到 $M$ 后的第 $i$ 个中间 correction，

$$
\chi_M(U_M^i)
\le
\theta^iD_E\tau_M
\le D_E\tau_M.
$$

同理可得

$$
\|\mathcal I_MU_M^i-U_\infty^\star\|_{\ell^2}
\le
\omega(D_E\tau_M)+\varepsilon_M
\to0.
$$

因此有限初始块之后的延拓点和中间 correction 也收敛到同一极限。$\square$

由于 $\ell^2$ 强收敛蕴含任意有限维投影收敛，对每个固定 $T$，

$$
\mathcal P_T\mathcal I_{N_j}U_{N_j}^j
\to
\mathcal P_TU_\infty^\star.
$$

这给出了固定控制前缀收敛。它比只说每个有限时域残差小更强，但仍依赖显式分支
假设。

---

## 13. 无限时域 KKT 与严格局部最优性

### 13.1 从有限前缀梯度到无限时域驻点

**定理 5（无限时域驻点）** 在第 12 节全部假设成立时，极限满足

$$
\nabla J_\infty(U_\infty^\star)=0.
$$

**证明。** 令 $Y_N=\mathcal I_NU_N^\star$。由分支收敛，$Y_N\to U_\infty^\star$ 于 $\ell^2$。固定任意 $T$；在 $N\ge T$ 时，有限时域驻点条件给出 $\mathcal P_T\nabla\widehat J_N(U_N^\star)=0$。因此

$$
\begin{aligned}
\|\mathcal P_T\nabla J_\infty(U_\infty^\star)\|
&\le\|\mathcal P_T[\nabla J_\infty(U_\infty^\star)-\nabla J_\infty(Y_N)]\|\\
&\quad+\|\mathcal P_T\nabla J_\infty(Y_N)-\mathcal P_T\nabla\widehat J_N(U_N^\star)\|\to0.
\end{aligned}
$$

第一项由梯度连续性趋零，第二项由固定前缀梯度一致性趋零。每个有限前缀梯度均为零，而无限梯度属于 $\ell^2$，故整个梯度为零。$\square$

### 13.2 二次增长与严格局部最优

**定理 6（严格局部最优性）** 若再满足可行切向曲率条件，则对所有充分接近
$U_\infty^\star$ 的可行 $U$，

$$
J_\infty(U)
\ge
J_\infty(U_\infty^\star)
+\frac\mu2
\|U-U_\infty^\star\|_{\ell^2}^2.
$$

因此 $U_\infty^\star$ 是严格局部极小点。

**证明。** 令

$$
h:=U-U_\infty^\star.
$$

由无限时域梯度为零，

$$
\langle\nabla J_\infty(U_\infty^\star),h\rangle
\ge0.
$$

$\mathscr U$ 凸且 $U$ 可行，所以 $h$ 属于
$T_{\mathscr U}(U_\infty^\star)$。当 $U$ 足够接近，使线段
$U_\infty^\star+th$ 完全位于 $\mathcal O$ 内，积分型 Taylor 公式给出

$$
\begin{aligned}
&J_\infty(U)-J_\infty(U_\infty^\star)\\
&=
\langle\nabla J_\infty(U_\infty^\star),h\rangle\\
&\quad+
\int_0^1(1-t)
\langle
h,
\nabla^2J_\infty(U_\infty^\star+th)h
\rangle\,dt.
\end{aligned}
$$

第一项非负，第二项由切向强制性满足

$$
\begin{aligned}
&\int_0^1(1-t)
\langle
h,
\nabla^2J_\infty(U_\infty^\star+th)h
\rangle\,dt\\
&\ge
\int_0^1(1-t)\mu\|h\|_{\ell^2}^2dt
=
\frac\mu2\|h\|_{\ell^2}^2.
\end{aligned}
$$

故二次增长不等式成立。$\square$

二次增长只在 $U_\infty^\star$ 的局部邻域和可行切向上成立，所以结论是严格局部
最优，不是全局最优。

---

## 14. 完整条件性定理

把前面的结果合并，可得到以下主结论。

**定理 7（两阶段算法的完整条件性保证）** 假设：

1. 第 1 节的控制仿射系统正则性、原点可稳定性和二次正定代价成立；
2. 按第 2 节选择了经过可靠余项验证的非平凡终端三元组；
3. 第 8 节的有限可达见证存在并通过精确 rollout 验证；
4. 第 9 节的时域一致局部 corrector 假设成立；
5. 第 10 节的强迫归一化延拓正则性成立；
6. 第 12 节的选定分支一致性、固定前缀梯度一致性和切向曲率成立。

则理想的两阶段算法满足：

$$
N_{\rm hit}<\infty,
\qquad
N_{\rm cert}<\infty;
$$

每个第二阶段保存点都可接入一条可行、有限代价的无限尾部；被接受的增广代价
收敛；每个认证时域的 correction 数一致有界于 $K_E$；理想认证时域满足

$$
N_j\to\infty;
$$

实现无限策略满足

$$
\mathcal I_{N_j}U_{N_j}^j
\to U_\infty^\star
\quad\text{于 }\ell^2;
$$

并且 $U_\infty^\star$ 是原无限时域问题的严格局部极小点。

若还满足几何候选和

$$
\widehat\theta^{K_{\rm tar}}D_E\le1,
$$

则首次认证后的每个时域至多需要 $K_{\rm tar}$ 次 correction，到有限软件上限
$N_{\max}$ 的 Riccati-stage 总工作量为

$$
O(K_{\rm tar}N_{\max}).
$$

**证明。** 有限 $N_{\rm hit}$ 来自定理 1；有限 $N_{\rm cert}$ 来自命题 2；可行
有限尾部来自第 3 节的望远镜求和；增广代价收敛来自第 11 节的单调有下界；一致
correction 上界和 $N_j\to\infty$ 来自定理 3；策略强收敛来自定理 4；无限 KKT
和严格局部最优分别来自定理 5、定理 6；复杂度结论来自第 10.4 节。$\square$

---

## 15. 证明边界汇总

各结论与其最少使用的证明接口如下。

| 结论 | 必需接口 | 不能顺带推出的结论 |
|---|---|---|
| 局部终端证书 | DARE、非线性余项界、有效局部半径 | 终端函数等于最优值函数 |
| 有限代价尾部 | 终端成员关系、终端下降 | 有限前缀满足 KKT |
| 第一阶段有限命中 | 经验证有限可达见证 | iLQR 全局收敛 |
| 性能通道条件性命中 | $D_j^+=o(N_j)$、全轨迹扫描 | 末端序列 $x_{N_j}\to0$ |
| 有限时域近似驻点 | $\chi_N\le\tau_N$，并用 $\mathcal G_N$ 审计 | 局部极小或全局极小 |
| 一致 correction 次数 | 统一局部收缩、归一化延拓界 | 控制序列自动收敛 |
| 实现策略强收敛 | 残差误差界、选定分支收敛 | 全局最优 |
| 无限时域 KKT | 固定前缀梯度一致性、梯度连续性 | 二阶局部最优 |
| 严格局部最优 | 无限 KKT、可行切向正曲率 | 全局最优 |

特别需要保留以下三条否定性边界：

1. $x_N\in\mathcal X_f$ 只证明尾部可行和有限代价，不证明移动终端
   $x_N\to0$；
2. $\chi_N$ 是 Riccati 预条件残差，不是原控制坐标下的真实梯度；
3. 到达有限 $N_{\max}$ 只表示计算上限，不能称为无限时域收敛。

---

## 16. 按计算顺序的最小实现检查表

对一次完整运行，计算与证明接口应按以下顺序检查：

1. 求 $P$、$K_f$，数值验证 DARE 余量

$$
\|A_c^\top PA_c-P+Q+K_f^\top RK_f\|;
$$

2. 给出可靠的 $c_r$ 和 $\rho_f$，据此选择 $\alpha_f$，验证正不变性和下降不等式；
3. 对可达见证做独立精确 rollout，并截断在首次
   $\mathcal X_f^\delta$ 命中；
4. 从 $N=1$ 求解严格凸一步 QP，后续每个新增状态重复该 QP，并记录每个
   $\Delta_N^g$、累计 $D_j^+$ 和 $\widehat J_{N_j}/N_j$；
5. 每次 iLQR correction 都重新计算真实 Jacobian、反向 Riccati 方向、预测下降、
   原始非线性 trial、真实下降和 ratio；
6. 第一阶段命中后保存回退策略，记录 $N_{\rm hit}$；
7. 第二阶段同时计算 $\chi_N$ 和精确伴随梯度 $\mathcal G_N$，通过安全与残差门后
   记录 $N_{\rm cert}$；
8. 每个候选时域都实际做终端反馈延拓并计算
   $\chi_M(\mathcal E_N^MU_N)/\tau_M$，再选择时域；
9. 每个新时域必须在有限 correction 后重新通过自己的认证门；
10. 若要支持超出有限时域近似驻点的结论，还必须单独检验或证明分支误差、固定
    前缀梯度一致性、移动终端行为和无限维二阶条件。

这套顺序把“可计算量”和“理论接口”一一对应，避免用终端安全证书替代驻点证明，
也避免用有限时域数值收敛替代无限时域局部最优性证明。
