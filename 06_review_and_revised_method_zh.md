# 时域延拓 iLQR 的审计与改进：从安全前缀到同一无限时域目标

初稿：2026-09-05；无输入约束修订：2026-09-06；完整流程整理：2026-09-07。本文核对了当前代码、`03`—`05` 笔记及论文，以 `05_clf_guided_horizon_continuation_zh.md` 为被审计的最新候选方案，并给出建议采用的修订路线；这不表示旧生产代码已经完成迁移。

本次按原问题统一取 $u_k\in\mathbb R^m$，没有输入集合、幅值上限或饱和约束。终端集、CLF 收缩和局部认证球是算法与证明条件，不是原问题的输入取值范围。

## 1. 判断与研究目的

当前思路值得保留的核心是：**先找到可以接稳定尾部的真实非线性前缀，再增加可优化的前缀长度，并用 iLQR 校正。** 主要问题不在 Riccati 计算，而在不同阶段使用的目标、残差、可行域和无限时域函数空间尚未完全对齐。

推荐将研究目的写成：

> 对固定初值的非线性无折扣最优控制问题，逐步构造可执行的有限前缀加稳定反馈尾部；在明确的局部可行性、稳定性和曲率条件下，用有限计算认证真实无限时域策略的局部精度，并证明随优化前缀增长，算法解收敛到该邻域中的严格局部最优解。

本次改进有两层。第一层修正 CLF 阶段的接受和失败逻辑，保留可达性证书。第二层把优化目标改成给定反馈尾部的**实际代价**，采用稳定反馈扰动变量，获得同一目标上的嵌套优化问题。第二层的强结论限定在可认证的局部邻域；第一阶段命中终端集不自动满足这一入口条件。

最终成果应包含：可行策略、真实残差、尾部误差和结论适用范围。不得把“达到最大时域”或“小步停止”列为最优性证书。

阅读上，第 4 节给阶段一，第 5—8 节给阶段二所需的坐标、目标、有限尾计算和 corrector，第 9—10 节给证明与最终证书；第 11 节按实际执行顺序把这些部件重新串成完整算法。

### 1.1 方法的五层职责

| 层次 | 输入 | 主要任务 | 可以得到 | 不能据此得到 |
|---|---|---|---|---|
| 终端准备 | 局部模型、$Q,R$ | 构造并验证 $P,\kappa_f,\mathcal X_f,T_f$ | 稳定有限代价尾部 | 远端初值可达、最优性 |
| 阶段一 | $x_0$、CLF 参数 | 搜索并提交真实收缩短块 | 可执行的安全前缀加反馈尾部 | 驻点、局部最优 |
| 阶段间入口 | 阶段一策略 | 验证稳定坐标、可行球和一致曲率 | 第二阶段的局部证明环境 | 全局吸引域 |
| 阶段二 | 参考反馈与 $v=0$ | 固定 $N$ 校正、几何增加自由变量 | 有限前缀真实残差和收敛检查点 | 仅凭前缀残差得到无限精度 |
| 最终认证 | 当前检查点 | 合并前缀梯度和未开放尾梯度 | 局部距离与性能差上界 | 全局最优性 |

全文始终把可行性、真实下降、固定时域驻定、无限时域收敛和局部最优性作为五种不同结论。

## 2. 对修订前方案的审计及仍需保留的结论

| 位置 | 已有内容 | 审计结论及修正 |
|---|---|---|
| `05` §2.3—2.4 | 任意 $W=x^TSx$ 定义终端集，使用 DARE 矩阵 $P$ 构造 $V_f$ | $V_f$ 下降不保证任意 $W$ 子水平集不变。主线固定 $S=P$；否则单独证明 $W$ 不变性。 |
| `05` §8—9 | 从正松弛开始探索，又要求试探块立即硬收缩 | 区分内部可行性搜索与外层提交；搜索允许正松弛，只有提交必须通过真实硬证书。 |
| `05` §11、§15.3 | 全前缀校正后声称保留所有历史块收缩 | 实际只保持当前终点预算。应对不同版本轨迹的当前终点归纳，不能把历史端点当成未变化。 |
| `03`、`05` §18 | 在固定首次命中时域上安全校正直到原目标残差小 | 可能根本不存在既满足终端门又满足原目标驻点条件的点；内部初值不能修复。 |
| `continuation.py` | 用未正则化前馈量 $\chi$ 过门 | 有限阈值下不等于真实梯度；保留作数值诊断，不能直接作原问题精度证书。 |
| `fixed_horizon.py` | 小轨迹变化或小相对代价变化即 `converged` | 大阻尼即可制造假收敛；应报告停滞，成功必须通过相应真实残差。 |
| `04` 无限时域部分 | 假设原开环控制的 $\ell^2$ 邻域内目标为 $C^2$ | 对积分器、开环不稳定系统可能不存在这样的有限代价邻域。改用稳定反馈变量或可行轨迹流形。 |
| 延拓分析 | 另设有限分支趋向无限分支 | 作为条件传递命题可以成立，但未证明分支收敛。改进版从同一目标的稠密子空间逼近直接推导。 |
| 计算优势 | 热启动、每段少量校正、线性工作量 | 热启动不证明常数校正次数。主结论改为有条件的 $O(N\log N)$，常数次数作为额外加强条件。 |

### 2.1 固定时域安全校正可能必然卡住

考虑完全线性的反例：

$$
x^+=x+u,\quad \ell=x^2+u^2,\quad x_0=2,\quad u\in\mathbb R,
$$

$$
P=\frac{1+\sqrt5}{2},\quad V_f(x)=2Px^2,\quad
\kappa_f(x)=-x/P,\quad \mathcal X_f=[-0.1,0.1].
$$

DARE 给出 $P$，闭环为 $x^+=x/P^2$，故该终端集不变，$V_f$ 满足所需 Lyapunov 下降。时域一的目标是

$$
J_1(u)=4+u^2+2P(2+u)^2,\qquad
J_1'(u)=2u+4P(2+u).
$$

无终端约束的唯一驻点满足

$$
x_1^\star=\frac{2}{1+2P}\approx0.472136>0.1.
$$

终端安全只允许 $u\in[-2.1,-1.9]$，该区间内 $J_1'$ 恒负。即使用 $u=-2$ 从终端集严格内部开始，也不可能在保留终端门的同时使原梯度趋零。安全约束问题的最优点是 $u=-1.9$；它需要非零终端约束乘子才能满足完整 KKT。

这里的区间仅由终端成员条件诱导，原问题没有输入范围限制。因此，“安全筛选 + 原目标梯度残差 + 固定时域必过门”不能一起作为普遍结论。只增加阻尼、线搜索或缩小初始终端集，均不能消除该逻辑冲突。

### 2.2 不稳定系统中的函数空间问题

对 $x_{k+1}=x_k+u_k$，在一条有限代价轨迹上仅将 $u_0$ 增加 $\epsilon\ne0$，其余开环控制保持不变。扰动控制的 $\ell^2$ 范数只有 $|\epsilon|$，但所有后续状态增加 $\epsilon$，正定状态代价使新代价为无穷。

所以不能在这类系统上直接假设原开环目标在普通控制 $\ell^2$ 开球内是有限的光滑函数。下面改用稳定反馈扰动，后续控制会随扰动后的状态反馈变化。

## 3. 统一符号与结论范围

原问题为

$$
\begin{aligned}
\min_{\{u_k\}_{k\ge0},\ u_k\in\mathbb R^m}\quad&
\sum_{k=0}^{\infty}(x_k^TQx_k+u_k^TRu_k)\\
\mathrm{s.t.}\quad&x_{k+1}=a(x_k)+B(x_k)u_k,\qquad x_0\ \text{给定}.
\end{aligned}
$$

这里 $u_k\in\mathbb R^m$ 只声明变量维数，没有任何输入幅值限制。有限总代价由 $R\succ0$ 推出 $\sum_k\|u_k\|^2<\infty$，这是有限代价的后果，不是预先施加的输入取值范围。以下沿用

$$
x_{k+1}=F(x_k,u_k)=a(x_k)+B(x_k)u_k,\qquad
J_\infty=\sum_{k=0}^\infty\ell(x_k,u_k),\qquad
\ell=x^TQx+u^TRu,\quad Q,R\succ0.
$$

$a(x)\in\mathbb R^n$，$B(x)\in\mathbb R^{n\times m}$。假设 $F(0,0)=0$。$k$ 为系统时刻，$b$ 为第一阶段已提交块的编号，$j$ 为第二阶段时域编号，$i$ 为固定时域校正编号。所有状态均由原系统 rollout。

终端构造使用 $W(x)=x^TPx$、$V_f=\gamma W$，$\gamma>1$。在 $\mathcal X_f=\{W\le\alpha_f\}$ 上验证反馈有定义、不变性及

$$
V_f(F(x,\kappa_f(x)))-V_f(x)\le-\ell(x,\kappa_f(x)).
$$

取消输入范围后，无需输入可行性检验，也无需用反馈幅值进一步缩小终端集；终端集半径只由真实非线性余项、不变性和 Lyapunov 下降决定。若平衡输入非零，应先平移平衡点，并在 $F_x$ 中保留输入加权的 $B_x$ 项。

原问题只保留动力学等式，不增加状态或输入的硬取值范围。阶段一的 CLF 收缩门与阶段二的局部认证邻域用于保证构造和证明；它们仍须显式说明。若在固定时域子问题中保留终端硬约束，其乘子仍存在，不能因取消输入约束就把该约束问题的 KKT 与原目标梯度混同。

## 4. 第一阶段：证书优先的短块扩展及其求解

### 4.1 核心子问题

当前终点记为 $z_b=x_{N_b}^{(b)}$。固定 $0<\bar\rho<1$、有限步候选 $h=1,\ldots,h_{\max}$ 和允许的中间增长系数 $c\ge1$。寻找真实控制块使

$$
z_{t+1}=F(z_t,w_t),\quad z_0=z_b,\quad w_t\in\mathbb R^m,
$$

$$
W(z_h)\le\bar\rho W(z_b),\qquad
W(z_t)\le cW(z_b)\quad(1\le t<h).
$$

一步问题中 $F(z,w)=a(z)+B(z)w$，故 $W(F(z,w))$ 为凸二次函数。记 $C=P^{1/2}B(z)$、$d=P^{1/2}a(z)$，则无输入约束时

$$
\min_{w\in\mathbb R^m}W(F(z,w))
=\min_w\|d+Cw\|^2
=\|(I-CC^\dagger)d\|^2.
$$

一个最小范数解为 $w_{\min}=-C^\dagger d$。即使 $C$ 不满列秩，最小二乘解仍存在，但一般不唯一；这里不能称为严格凸问题。若最小值不超过 $\bar\rho W(z)$，就得到一步收缩可行性。

在硬收缩条件下进一步最小化 $w^TRw+\lambda_bW(F(z,w))$，仍是凸 QCQP；其控制能量项使目标强制且严格凸，故只要收缩可行，就存在唯一性能最优解。这里唯一的不等式是算法要求的 CLF 收缩，没有输入箱约束。

若只求不带 CLF 硬条件的一步加权问题，则直接求解线性方程

$$
(R+\lambda_bB^TPB)w_{\lambda_b}=-\lambda_bB^TPa,
$$

即

$$
w_{\lambda_b}=-(R+\lambda_bB^TPB)^{-1}\lambda_bB^TPa.
$$

矩阵因 $R\succ0$ 始终正定，不需要投影或裁剪；但该加权解仍须另外检验是否满足 CLF 收缩。无输入幅值限制也不意味着欠驱动系统一定能一步收缩。

有限步流通常非凸。下面把固定块长的求解拆成“先找零违反证书，再在硬约束下改善性能”两层。这个拆分使可行性和性能职责分离，也给出了局部求解失败时应保留的回退点。

### 4.2 固定块长的第一层：最小违反量搜索

固定当前块起点 $z_b$ 和长度 $h$，记

$$
W_b=W(z_b),\qquad
U_h=\operatorname{col}(w_0,\ldots,w_{h-1}),\qquad
X_h=\operatorname{col}(z_1,\ldots,z_h).
$$

先求共同松弛问题

$$
\begin{aligned}
\inf_{X_h,U_h,s}\quad&s\\
\mathrm{s.t.}\quad
&z_0=z_b,\\
&z_{t+1}=F(z_t,w_t),\qquad t=0,\ldots,h-1,\\
&W(z_h)-\bar\rho W_b\le s,\\
&W(z_t)-cW_b\le s,\qquad t=1,\ldots,h-1,\\
&s\ge0.
\end{aligned}
$$

也可以给终点和各中间时刻使用独立松弛，再最小化它们的最大值。共同松弛形式对应真实最大违反量

$$
\Phi_h(U_h)=
\max\left\{
0,\,
W(z_h(U_h))-\bar\rho W_b,\,
\max_{1\le t<h}\bigl[W(z_t(U_h))-cW_b\bigr]
\right\}.
$$

当 $h=1$ 时不存在中间时刻约束，上式中的内部最大项直接省略。

第一层的目的不是证明全局最小值等于零，而是找到任意一个满足 $\Phi_h(U_h)=0$ 的实际控制块。无界控制域上，下确界未必取得，因此 $s_h^{\inf}=0$ 本身不是证书；反过来，局部求解器停在 $s>0$ 也不能证明该块长不可行。

不能把第一层随意改成单一加权和

$$
\min\ s+\epsilon\|U_h\|^2.
$$

即使存在 $s=0$ 的控制块，加权和仍可能偏好控制较小但 $s>0$ 的点。控制正则只能作为字典序第二目标、SQP 子问题的近端项或信赖域机制；第一优先级必须保持为找到零违反点。

### 4.3 固定块长的第二层：硬收缩下优化性能

第一层一旦找到经真实 rollout 验证的零违反块，立即把它保存为 $U_h^{\rm cert}$。随后以它为初值求

$$
\begin{aligned}
\min_{X_h,U_h}\quad&
\sum_{t=0}^{h-1}
\left(z_t^TQz_t+w_t^TRw_t\right)
+\lambda_bW(z_h)\\
\mathrm{s.t.}\quad
&z_0=z_b,\\
&z_{t+1}=F(z_t,w_t),\\
&W(z_h)\le\bar\rho W_b,\\
&W(z_t)\le cW_b,\qquad t=1,\ldots,h-1.
\end{aligned}
$$

记上式目标为 $J_b(U_h)$。$R\succ0$ 使目标对控制块具有强制性；若真实硬可行集非空且闭，则相应子水平集紧，故全局极小值存在。但有限步问题仍通常非凸，局部 SQP 只给局部结果。性能求解失败、终止在不可行点或没有取得真实下降时，直接返回已经保存的 $U_h^{\rm cert}$。

因此两层的严格优先级是

$$
\text{真实零违反证书}
\quad\succ\quad
\text{硬证书内的性能改善}.
$$

### 4.4 多重射击变量、稀疏 Jacobian 与约束梯度

推荐把

$$
z_1,\ldots,z_h,\qquad
w_0,\ldots,w_{h-1}
$$

都作为优化变量，并显式保留动力学等式。这种多重射击形式的约束 Jacobian 是块带状的，比完全消去状态后形成稠密控制 Hessian 更适合稀疏 SQP。

定义动力学残差

$$
r_{t+1}=z_{t+1}-F(z_t,w_t).
$$

在名义点 $(\bar z_t,\bar w_t)$ 上记

$$
A_t=F_x(\bar z_t,\bar w_t),\qquad
B_t=F_u(\bar z_t,\bar w_t).
$$

对于控制仿射系统，

$$
A_t=a_x(\bar z_t)
+\sum_{q=1}^m\bar w_{t,q}
\frac{\partial B_{\cdot q}}{\partial x}(\bar z_t),
\qquad
B_t=B(\bar z_t).
$$

动力学等式的非零 Jacobian 块为

$$
\frac{\partial r_{t+1}}{\partial z_{t+1}}=I,\qquad
\frac{\partial r_{t+1}}{\partial z_t}=-A_t,\qquad
\frac{\partial r_{t+1}}{\partial w_t}=-B_t.
$$

由于 $W(z)=z^TPz$，

$$
\nabla W(z)=2Pz,\qquad
\nabla^2W(z)=2P.
$$

阶段代价和块终端项的直接导数为

$$
\nabla_{z_t}\ell=2Qz_t,\qquad
\nabla_{w_t}\ell=2Rw_t,\qquad
\nabla_{z_h}\bigl[\lambda_bW(z_h)\bigr]=2\lambda_bPz_h,
$$

对应的直接 Hessian 分别为 $2Q$、$2R$ 和 $2\lambda_bP$。所以目标、终点约束和中间 CLF 约束的解析导数都可以直接提供给 SQP。

若采用只以 $U_h$ 为变量的消元形式，则控制到状态的敏感度满足

$$
S_{0,q}=0,\qquad
S_{t+1,q}=A_tS_{t,q}+B_t\mathbf 1_{\{q=t\}},
$$

其中 $S_{t,q}=\partial z_t/\partial w_q$。于是

$$
\nabla_{w_q}W(z_t)=2S_{t,q}^TPz_t.
$$

短块可以用这一消元形式快速实现；随着 $h$ 增大，多重射击通常具有更好的稀疏性和数值条件。

### 4.5 逐次凸化子问题

令当前名义控制块经过精确 rollout 得到 $\bar z_{0:h}$。线性化动力学为

$$
\delta z_0=0,\qquad
\delta z_{t+1}=A_t\delta z_t+B_t\delta w_t.
$$

因为 $W$ 本身是凸二次型，可以在线性化状态上保留完整二次表达，而不必再把 $W$ 线性化。第一层的一个凸 QCQP 子问题为

$$
\begin{aligned}
\operatorname{lexmin}_{\delta X,\delta U,s}\quad&
\left(s,\,
\frac12\|\delta X\|^2+
\frac12\|\delta U\|^2\right)\\
\mathrm{s.t.}\quad&
\delta z_{t+1}=A_t\delta z_t+B_t\delta w_t,\\
&W(\bar z_h+\delta z_h)\le\bar\rho W_b+s,\\
&W(\bar z_t+\delta z_t)\le cW_b+s,\quad 1\le t<h,\\
&s\ge0,\qquad
\|D_X\delta X\|_\infty\le\Delta_X,\quad
\|D_U\delta U\|_\infty\le\Delta_U.
\end{aligned}
$$

这里 $\operatorname{lexmin}$ 表示先最小化 $s$，再在所得最小松弛层内选择较小增量。实际软件若不支持字典序目标，可以先单独求最小 $s$，再固定 $s\le s_{\rm found}+\tau_s$ 求最小增量；不能只用未经论证的大权重代替字典序。

第二层使用相同的线性化动力学和信赖域，将目标换成

$$
\begin{aligned}
\widehat J_b(\delta X,\delta U)
=&
\sum_{t=0}^{h-1}
\left[
(\bar z_t+\delta z_t)^TQ(\bar z_t+\delta z_t)
+(\bar w_t+\delta w_t)^TR(\bar w_t+\delta w_t)
\right]\\
&+\lambda_bW(\bar z_h+\delta z_h)
+\frac{\sigma}{2}
\left(\|\delta X\|^2+\|\delta U\|^2\right),
\end{aligned}
$$

其中 $\sigma\ge0$ 是近端正则。配合硬 CLF 二次约束后仍是凸 QCQP。这里的凸性只属于当前局部子问题，不代表原非线性多步问题凸。

### 4.6 真实 rollout、接受比与信赖域更新

QCQP 给出的预测状态不能直接提交。取

$$
w_t^{\rm trial}=\bar w_t+\delta w_t
$$

后，必须从同一个 $z_b$ 出发重新计算

$$
z_{t+1}^{\rm trial}=F(z_t^{\rm trial},w_t^{\rm trial}).
$$

第一层比较预测和真实违反量下降：

$$
\operatorname{pred}_\Phi
=\Phi_h(\bar U_h)-\widehat\Phi_h(\bar U_h+\delta U_h),
$$

$$
\operatorname{ared}_\Phi
=\Phi_h(\bar U_h)-\Phi_h(U_h^{\rm trial}),
\qquad
r_\Phi=\frac{\operatorname{ared}_\Phi}
{\max\{\operatorname{pred}_\Phi,\epsilon_{\rm den}\}}.
$$

其中 $\widehat\Phi_h$ 是 QCQP 在线性化动力学上的预测最大违反量。只有真实 rollout 有定义、$\operatorname{pred}_\Phi>0$、$\operatorname{ared}_\Phi>0$ 且 $r_\Phi$ 达到预设门槛时才接受发现步；预测可行而真实不可行时缩短信赖域。只要真实约束已经带数值裕量通过，就结束第一层，不必继续把松弛目标求到机器零。

第二层只接受同时满足

$$
\Phi_h(U_h^{\rm trial})=0,
\qquad
J_b(U_h^{\rm trial})<J_b(\bar U_h)
$$

的真实试探点，并可使用相应的目标实际/预测下降比更新信赖域。任何未通过真实硬证书的性能步都被拒绝。

有限精度下应采用区间或保守裕量，例如

$$
\overline W(z_h^{\rm trial})
\le\bar\rho\,\underline W(z_b)-\delta_{\rm cert},
\qquad
\delta_{\rm cert}>0,
$$

而不是把浮点计算的 $s\le10^{-k}$ 直接等同于数学上的 $s=0$。

### 4.7 块长、初值与完整失败语义

块长外循环优先取

$$
h=1,2,\ldots,h_{\max}.
$$

若假设只保证这些整数中存在某个可行长度，就必须覆盖全部整数。几何候选可以减少尝试次数，但可能跳过唯一可行长度，除非另外假设可行长度属于候选集，否则不能据此声称搜索完备。

每个块长可以使用多组初值：

1. 全零控制块；
2. 一步加权控制的重复、衰减或平移；
3. 冻结 LTV 模型及其块敏感度矩阵给出的控制；
4. 前一块或相邻块长的解作 warm start；
5. 少量确定性扰动或随机多启动。

具体地，冻结仿射 LTV 预测可写成

$$
z_h\approx d_h+\mathcal C_hU_h,\qquad
\bar R_h=\operatorname{diag}(R,\ldots,R).
$$

忽略硬 CLF 约束的加权问题具有显式解

$$
U_{h,\lambda}
=-\lambda
\left(
\bar R_h+\lambda\mathcal C_h^TP\mathcal C_h
\right)^{-1}
\mathcal C_h^TPd_h.
$$

它可以作为初值；若连同终端和中间 CLF 条件一起处理，则得到一个凸 QCQP 预测器。冻结 LTV 解只用于预测和初始化，必须经过真实非线性 rollout。若某次局部求解停在正违反量，可以更换初值、扩大允许的控制信赖域或尝试下一块长；不能立即报告不可行。

运行状态必须区分：

- finite_step_certificate_found：已经找到并验证真实零违反块；
- local_search_failed：当前初值和局部方法没有找到证书；
- candidate_lengths_exhausted：尝试完预定块长仍未找到；
- globally_infeasible：只有完备全局方法给出不可行证明时才能使用。

主算法把第二、第三种情形汇总为 certificate_not_found，但其语义始终是“有限搜索未找到证书”，不是“系统不可达”。

### 4.8 去掉不必要的权重循环

主算法用固定有限的 $\lambda_b>0$ 即可。可行性由硬收缩保证，性能由 iLQR 改善，不需要为证明可达性再引入乘子驱动权重和 Bellman Hessian 匹配。它们可以作为独立加速消融，不进入主证明。

搜索到可行块后立即保存。后续块性能优化或全前缀 iLQR 失败，均保留该可行块；不得因为可选性能优化失败而丢失已经找到的证书。

### 4.9 校正只需保持当前终点预算

设

$$
\beta_0=W(x_0),\qquad \beta_{b+1}=\bar\rho\beta_b.
$$

扩展后的全前缀 iLQR 可运行有限次，接受条件为：同一冻结目标真实下降、真实动力学 rollout 有定义，并且

$$
W(x_{N_{b+1}}^{(b+1)})\le\beta_{b+1}.
$$

没有合格校正则保留提交块。只承诺当前版本终点满足预算，不声称改动后的历史端点仍满足原块间收缩。

这一接受门正是第 4.13 节归纳证明能够跨越可选校正步骤的原因。有限命中还需要第 4.12 节的有限步 CLF 存在性和求解器发现条件，不能从局部求解器名称推出。

### 4.10 数值误差进入接受式

可使用区间界验证 $\overline W(z_h)\le\bar\rho\,\underline W(z_b)$，或预留足够严格的相对裕量。若实现只能保证带统一加性误差的收缩，则必须采用第 4.14 节的误差地板、修正块数和切换条件；不能继续沿用精确几何收缩结论。

### 4.11 第一阶段的完整算法

第一阶段的输入是：

- 原始动力学 $F$、阶段代价 $\ell$ 和初值 $x_0$；
- 进度函数 $W(x)=x^TPx$、终端反馈 $\kappa_f$、终端集 $\mathcal X_f=\{W\le\alpha_f\}$；
- 切换阈值 $0<\alpha_{\rm sw}<\alpha_f$；
- 收缩率 $0<\bar\rho<1$、中间增长系数 $c\ge1$、最大块长 $h_{\max}$；
- 固定的块性能权重 $\lambda_b>0$、数值证书裕量和可选的全前缀校正调度。

初始化

$$
b=0,\qquad N_0=0,\qquad
z_0=x_0,\qquad
\beta_0=W(x_0),
$$

并令保存前缀只包含 $x_0$。随后执行以下循环。

**步骤一：检查是否已经到达切换集。** 若

$$
W(z_b)\le\alpha_{\rm sw},
$$

则停止短块搜索，令 $N_{\rm hit}=N_b$，把 $\kappa_f$ 接在当前前缀之后。若 $x_0$ 一开始就在切换集内，则 $N_{\rm hit}=0$。

**步骤二：从短到长搜索下一块。** 对

$$
h=1,2,\ldots,h_{\max}
$$

依次调用第 4.2—4.7 节的求解器。对每个 $h$：

1. 由全零控制、冻结 LTV 预测、相邻块 warm start 或多启动产生初值；
2. 用最小违反量 SQP/SCP 搜索 $\Phi_h(U_h)=0$；
3. 每个候选都从当前真实端点 $z_b$ 重新 rollout；
4. 一旦得到带数值裕量的真实零违反块，立即保存为 $U_h^{\rm cert}$；
5. 可选地从 $U_h^{\rm cert}$ 出发求解硬收缩性能问题；成功则使用改进块，失败则恢复 $U_h^{\rm cert}$。

默认提交第一个通过的最短块。若所有长度和预定初值均未找到真实证书，则保持原前缀不变并返回 certificate_not_found。该状态不等于不可达。

**步骤三：原子地提交控制块。** 设选中的块长为 $h_b$。把整个控制块及其真实状态一次性拼接到当前前缀，

$$
N_{b+1}=N_b+h_b,
$$

并立即保存新前缀。提交不是逐步移动截止时刻的滚动承诺；若实现必须逐时刻写入，也必须保留整块剩余控制和原终点收缩承诺，直到 $N_{b+1}$。

**步骤四：更新进度预算。** 令

$$
\beta_{b+1}=\bar\rho\beta_b.
$$

由提交块的真实收缩以及 $W(z_b)\le\beta_b$，未经全前缀校正的新终点满足

$$
W(x_{N_{b+1}})\le\bar\rho W(z_b)
\le\bar\rho\beta_b
=\beta_{b+1}.
$$

**步骤五：可选地校正整个当前前缀。** 可以冻结当前时域 $N_{b+1}$ 和阶段一性能目标，做有限次 iLQR。每个试探步必须用原动力学重新 rollout，并同时满足：

1. 冻结目标真实下降；
2. rollout 全程有定义；
3. 当前终点继续满足

   $$
   W(x_{N_{b+1}})\le\beta_{b+1}.
   $$

没有合格步就恢复步骤三保存的前缀。这里不要求校正后的历史块端点继续满足旧收缩式；第一阶段归纳只使用当前版本的当前终点预算。

**步骤六：进入下一块。** 用校正后或回退后的真实终点定义

$$
z_{b+1}=x_{N_{b+1}},
$$

令 $b\leftarrow b+1$，返回步骤一。

第一阶段最终输出

$$
\Pi^{\rm safe}
=
(u_0,\ldots,u_{N_{\rm hit}-1})
\oplus\kappa_f,
$$

即一个有限真实控制前缀加稳定反馈尾部；它是第二阶段参考策略的中心。

### 4.12 第一阶段依赖的假设

为避免把不同强度的条件混在一起，按其用途列出假设。

**A1：基本适定性。** $F$ 在算法访问区域连续，有限控制块产生唯一且数值有限的真实 rollout；$P\succ0$，$Q,R\succ0$。若使用第 4.4—4.6 节的 SQP/SCP 导数和模型误差分析，还要求相应区域内 $F$ 至少为 $C^1$，并具有足够的局部 Jacobian 连续性或 Lipschitz 界。连续性负责证书和存在性，导数正则性负责具体局部求解器，不应混成一个无差别的全局 $C^2$ 假设。

**A2：终端尾证书。** 已验证

$$
\mathcal X_f=\{x:W(x)\le\alpha_f\}
$$

在 $\kappa_f$ 下正向不变，且在 $\mathcal X_f$ 上

$$
V_f(F(x,\kappa_f(x)))-V_f(x)
\le-\ell(x,\kappa_f(x)).
$$

切换集严格包含于终端集，即 $0<\alpha_{\rm sw}<\alpha_f$。这个假设只负责进入切换集后的无限尾部可执行性和有限代价。

**A3：沿访问区域的有限步 CLF 存在性。** 对第一阶段可能访问的每个

$$
z\notin\{W\le\alpha_{\rm sw}\},
$$

存在某个 $h(z)\in\{1,\ldots,h_{\max}\}$ 和实际有限控制块，使

$$
W(z_{h(z)})\le\bar\rho W(z),\qquad
W(z_t)\le cW(z),\quad 1\le t<h(z).
$$

这是排除算法在某个访问状态因不存在下一块而阻塞的核心控制假设。它弱于区域一致一步 CLF，但强于“从初值存在某条最终可达轨迹”。

**A4：有限发现性或求解器条件。** 对 A3 保证存在的块，长度外循环和内层求解器能在有限计算内实际找到一个经真实 rollout 验证的零违反块。A3 是控制问题的存在性，A4 是算法的发现能力；普通局部 SQP、iLQR 或多启动本身不自动证明 A4。若没有 A4，仍可对每个实际找到的块给出运行时证书，但不能预先证明算法永不返回 certificate_not_found。

**A5：原子提交与保护性校正。** 控制块整体提交；块性能优化和全前缀 iLQR 都保留最后一个已验证前缀，并且只接受保持当前终点预算的真实 rollout。由于校正是可选的，拒绝全部校正步也不破坏阶段一的终止证明。

**A6：数值证书可靠性。** 理想分析采用精确 rollout 和精确 $W$。有限精度实现必须使用区间界或明确裕量，使被接受块在真实量上满足收缩。若只能保证统一加性误差

$$
W(z_{b+1})\le\bar\rho W(z_b)+e,
$$

则所有结论必须改用第 4.10 节的误差地板版本。

这些假设与结论的依赖关系为：

| 结论 | 所需条件 |
|---|---|
| 某个已提交块是真实 CLF 证书 | A1、该块通过真实验证、A6 |
| 校正不会丢失当前进度预算 | A5 |
| 第一阶段不会在切换集外阻塞 | A1、A3、A4、A5、A6 |
| 有限块数进入切换集 | 上一行再加 $0<\bar\rho<1$、$\alpha_{\rm sw}>0$ |
| 得到有限代价无限策略 | 有限进入切换集再加 A2 |

### 4.13 第一阶段有限终止与可行尾部定理

**定理（精确证书情形）。** 假设 A1—A5 成立，并且每个提交块精确满足收缩门。则第一阶段具有以下性质。

1. 对每个已提交块编号 $b$，当前保存前缀的终点满足

   $$
   W(z_b)\le\beta_b=\bar\rho^bW(x_0).
   $$

2. 算法至多经过

   $$
   b_{\rm hit}
   =
   \max\left\{
   0,\,
   \left\lceil
   \frac{\log(W(x_0)/\alpha_{\rm sw})}
   {\log(1/\bar\rho)}
   \right\rceil
   \right\}
   $$

   个已提交块进入切换集；当 $W(x_0)=0$ 时直接取 $b_{\rm hit}=0$。
3. 因为每块 $h_b\le h_{\max}$，切换时的物理时域满足

   $$
   N_{\rm hit}
   =\sum_{b=0}^{b_{\rm hit}-1}h_b
   \le h_{\max}b_{\rm hit}.
   $$

4. 将 $\kappa_f$ 接在 $N_{\rm hit}$ 后得到满足原动力学的无限策略，而且

   $$
   \sum_{k=N_{\rm hit}}^\infty
   \ell(x_k,\kappa_f(x_k))
   \le V_f(x_{N_{\rm hit}})<\infty.
   $$

   有限前缀只含有限个数值有限的状态和控制，所以整个策略总代价有限。

**证明。** 初始时 $W(z_0)=\beta_0$。假设第 $b$ 个当前终点满足 $W(z_b)\le\beta_b$。提交块的真实收缩给出

$$
W(z_{b+1}^{\rm raw})
\le\bar\rho W(z_b)
\le\bar\rho\beta_b
=\beta_{b+1}.
$$

若不做全前缀校正，直接保存该端点；若做校正，A5 的接受门仍强制校正后端点不超过 $\beta_{b+1}$。因此由归纳得到第一条。

当 $b\ge b_{\rm hit}$ 时，

$$
W(z_b)\le\bar\rho^bW(x_0)\le\alpha_{\rm sw},
$$

所以在不超过 $b_{\rm hit}$ 个提交块后触发切换。A3 保证每个非终端访问状态存在长度不超过 $h_{\max}$ 的块，A4 保证算法有限找到它，因此在达到该块数之前不会阻塞。这给出第二、第三条。

切换时 $W(x_{N_{\rm hit}})\le\alpha_{\rm sw}<\alpha_f$，故端点属于 $\mathcal X_f$。A2 的不变性使后续反馈轨迹始终留在终端区；对 Lyapunov 下降从 $N_{\rm hit}$ 到 $N_{\rm hit}+T-1$ 求和，

$$
\sum_{k=N_{\rm hit}}^{N_{\rm hit}+T-1}
\ell(x_k,\kappa_f(x_k))
\le
V_f(x_{N_{\rm hit}})
-V_f(x_{N_{\rm hit}+T})
\le V_f(x_{N_{\rm hit}}).
$$

令 $T\to\infty$ 即得第四条。

### 4.14 有限精度版本和复杂度边界

若每块只能保证

$$
W(z_{b+1})\le\bar\rho W(z_b)+e,
$$

归纳得到

$$
W(z_b)
\le
\bar\rho^bW(x_0)
+\frac{e(1-\bar\rho^b)}{1-\bar\rho}.
$$

记误差地板

$$
E_\infty=\frac{e}{1-\bar\rho}.
$$

只有当 $\alpha_{\rm sw}>E_\infty$ 时，上界才可能在有限块数内进入切换集。若 $W(x_0)>\alpha_{\rm sw}>E_\infty$，充分块数为

$$
b\ge
\left\lceil
\frac{
\log\bigl((W(x_0)-E_\infty)/
(\alpha_{\rm sw}-E_\infty)\bigr)}
{\log(1/\bar\rho)}
\right\rceil.
$$

第一阶段能够证明块数和最终物理时域上界，但不能仅由上述假设给出非凸短块搜索的统一多项式复杂度。若另外给定每个访问状态上有限步求解器的统一工作上界 $C_{\rm block}$，才可写出

$$
\operatorname{work}_{\rm I}
\le b_{\rm hit}C_{\rm block}
+\operatorname{work}_{\rm optional\ corrector}.
$$

可选全前缀 iLQR 的工作量必须单列，不能隐藏在短块证书复杂度中。

### 4.15 第一阶段究竟证明了什么

第一阶段在上述条件下证明：

- 被提交控制块来自真实非线性 rollout；
- 当前保存前缀的块端点进度按几何预算下降；
- 条件性地在有限块数、有限物理时域内进入切换集；
- 最终得到一个可执行、满足原动力学且总代价有限的“有限前缀加反馈尾部”策略；
- 任何可选性能校正失败都不会删除最后一个已认证前缀。

第一阶段没有证明：

- 短块 NLP 的局部解是全局最优；
- 未找到控制块意味着系统不可达；
- 校正后的所有历史块端点仍分别满足原来的收缩式；
- 当前前缀满足原无限时域目标的驻点或 KKT 条件；
- 所得策略是局部或全局最优；
- 第二阶段的稳定坐标、强凸性和完整梯度入口自动成立。

因此，第一阶段的准确定位是“构造并保护一个可执行的有限代价参考策略”，而不是“求解无限时域最优控制问题”。局部驻定、时域延拓和最优性从第 5 节开始另行建立。

## 5. 第二阶段：稳定反馈变量和可行邻域

### 5.1 如何保留阶段一的远端前缀

阶段一交付安全前缀 $(\bar x_k,\bar u_k)_{k<N_0}$。定义固定的参考反馈

$$
\mu_k(x)=
\begin{cases}
\bar u_k+K_k^{\rm ref}(x-\bar x_k),&k<N_0,\\
\kappa_f(x),&k\ge N_0.
\end{cases}
$$

前 $N_0$ 个反馈只需在该有限前缀附近光滑且使真实轨迹有定义；$K_k^{\rm ref}=0$ 也可以，但可能导致较差的敏感性界。不要假设局部 $\kappa_f$ 在远端初始区域也有效。

用新变量

$$
u_k=\mu_k(x_k)+v_k,\qquad
x_{k+1}=\mathcal F_k(x_k,v_k):=F(x_k,\mu_k(x_k)+v_k),\qquad
v\in \mathcal H:=\ell^2(\mathbb N;\mathbb R^m).
$$

$v=0$ 正好实现保存的安全前缀和反馈尾部。记 $\mathcal J(v)$ 为这条真实无限轨迹的总代价。

### 5.2 局部认证入口，不把命中终端集当作全部条件

在 $\mathcal B_{r_{\mathrm B}}=\{v:\|v\|_2\le r_{\mathrm B}\}$ 的开邻域内，需建立以下条件：

1. 所有轨迹由原动力学真实生成且有定义，尾部位于稳定反馈的定义域；轨迹映射 $v\mapsto(X(v),U(v))$ 取值于 $\ell^2\times\ell^2$，为 $C^2$ 且导数有界。
2. 同一个真实目标满足一致曲率界

$$
0<mI\preceq D^2\mathcal J(v)\preceq LI<\infty.
$$

3. 中心的**完整无限梯度**具有经有限计算认证的上界

$$
\eta\ge\|\nabla\mathcal J(0)\|_2,\qquad \eta<mr_{\mathrm B}/2.
$$

这些条件不假设有限解分支或无限最优解已经存在。轨迹适定性、局部稳定域及曲率界是真正的局部假设；必须用模型界、稳定敏感性和二阶分析验证。$Q,R\succ0$ 本身不证明非线性约化目标强凸。

中心证书未通过时，可以另做保留终端尾部的性能改善或更换参考反馈，再重新认证；不能宣称这些尝试一定成功。此时仍可返回阶段一的可行有限代价策略，但不附加本节的局部最优结论。

在新变量下，每个 $v_k\in\mathbb R^m$ 也没有输入幅值限制，更新 $u_k=\mu_k(x_k)+v_k$ 不做裁剪。$\mathcal B_{r_{\mathrm B}}$ 是用于建立稳定性、可微性和局部最优性的序列空间邻域；它没有被写入原最优控制问题。后文证明下降迭代及极小点严格位于该球内部，因此最终局部结论不依赖人工球边界乘子。

动力学消元后的无约束驻点条件直接是 $\nabla\mathcal J_N(v)=0$；在原输入变量中则为 $\nabla_UJ_N(U)=0$。这里不再使用控制集合投影残差或控制法锥。第 6 节的 $\mathsf P_N$ 只表示取前 $N$ 个坐标的线性投影，和输入约束无关。

### 5.3 条件的具体来源

稳定反馈的变分系统具有有界输入—状态算子，是第一条的关键。有限长度的远端前缀只改变常数；稳定尾部避免这些常数随优化时域无限增长。零输入下稳定不自动意味着对所有 $\ell^2$ 扰动稳定，仍需局部鲁棒性界。

令 $X_1=DX(v)h$、$U_1=DU(v)h$。真实 Hessian 为

$$
D^2\mathcal J(v)[h,h]
=2\langle QX_1,X_1\rangle+2\langle RU_1,U_1\rangle
+2\langle QX,D^2X[h,h]\rangle
+2\langle RU,D^2U[h,h]\rangle.
$$

若 $\|D\mu_k\|\le K$，由 $h_k=\delta u_k-D\mu_k\delta x_k$ 得

$$
\|h\|_2^2\le(1+K^2)(\|X_1\|_2^2+\|U_1\|_2^2).
$$

所以前两项至少为 $m_0\|h\|_2^2$，其中 $m_0=2\min\{\lambda_{\min}(Q),\lambda_{\min}(R)\}/(1+K^2)$。若后二阶余项的绝对值在球内不超过 $c_2\|h\|_2^2$ 且 $c_2<m_0$，可取 $m=m_0-c_2$。这是验证路线，不能以若干采样 Hessian 正定代替球内的一致界。

## 6. 同一目标上的时域延拓

### 6.1 给定反馈的真实尾价值

在终端区定义

$$
T_f(z)=\sum_{t=0}^\infty\ell(z_t,\kappa_f(z_t)),\qquad
z_{t+1}=F_f(z_t):=F(z_t,\kappa_f(z_t)).
$$

它是给定反馈的代价，不是未知的最优价值函数。Lyapunov 求和给出

$$
0\le T_f(z)\le V_f(z),\qquad
T_f(z)=\ell(z,\kappa_f(z))+T_f(F_f(z)).
$$

第一条只需可行与下降；导数公式需要额外的闭环敏感性界，见下一节。

令 $\mathcal H_N=\{v:v_k=0,\ k\ge N\}$，$N\ge N_0$，$\mathsf P_N$ 为正交投影。有限维目标取

$$
\mathcal J_N(v)=\sum_{k=0}^{N-1}\ell(x_k,\mu_k(x_k)+v_k)+T_f(x_N).
$$

对 $v\in \mathcal H_N$，尾部自然执行 $\kappa_f$，因此严格有

$$
\mathcal J_N(v)=\mathcal J(v),\qquad
\mathcal J_M(v,0)=\mathcal J_N(v),\qquad
\mathsf P_N\nabla\mathcal J_M(v,0)=\nabla\mathcal J_N(v).
$$

旧变量的目标及梯度在零填充时完全一致；新增变量的梯度可以非零。因此热启动不是已经解好大时域问题。真正的改进在于：不再需要用移动 $V_f$ 的边界缺陷修补不同优化目标之间的关系。

终端集仍负责定义安全局部尾部，但局部认证球已经保证整条策略可行，且下降子水平集与球边界有严格距离。由此不再依赖“固定 $N_{\rm hit}$ 加终端硬筛选一定能达到原驻点”的错误接口。**仅换成 $T_f$ 而保留原硬筛选，并不足以解决第 2.1 节反例；可行邻域和内部性同样必不可少。**

## 7. 无限尾部怎样用有限计算认证

设对终端区内的 $z$ 有经验证的界

$$
\|F_f^t(z)\|\le C_xq_x^t\|z\|,\qquad
\|D F_f^t(z)\|\le C_dq_d^t,\qquad 0<q_x,q_d<1,
$$

$$
\ell(z,\kappa_f(z))\le C_\ell\|z\|^2,\qquad
\|\nabla_z\ell(z,\kappa_f(z))\|\le C_r\|z\|.
$$

计算有限和 $T_{f,H}(z)=\sum_{t=0}^{H-1}\ell(z_t,\kappa_f(z_t))$。每一尾项非负，故

$$
0\le T_f(z)-T_{f,H}(z)
\le E_0(H,z):=\frac{C_\ell C_x^2q_x^{2H}}{1-q_x^2}\|z\|^2.
$$

每一项的初值导数为 $(D F_f^t(z))^T\nabla\ell(z_t,\kappa_f(z_t))$。一致几何上界允许逐项求导，并给出

$$
\|\nabla T_f(z)-\nabla T_{f,H}(z)\|
\le E_1(H,z):=\frac{C_rC_xC_d(q_xq_d)^H}{1-q_xq_d}\|z\|.
$$

这一步独立证明导数误差，未从值误差猜测梯度误差。二阶光滑性需要对应的二阶敏感性可求和界。

对前缀 $v$，记 $S_N=D_vx_N$。有限尾计算的梯度误差满足

$$
\|\widehat g_N-\nabla\mathcal J_N(v)\|_2\le
\delta_{g,N}:=\overline{\|S_N\|}E_1(H,x_N)+e_{g,\rm arith}.
$$

$S_N$ 可由线性敏感性递推计算或用稳定算子界控制；横线表示经验证的上界。$e_{g,\rm arith}$ 覆盖状态 rollout、有限和、伴随及敏感性数值计算产生的误差。函数值区间也必须覆盖完整 rollout 的舍入误差，而非只给末端截断项加一个误差。尾部的控制是固定反馈，不是新的自由变量；计算尾价值及其导数时不能对这些控制做最小化。

最终有限前缀残差门是

$$
\boxed{\|\widehat g_N\|_2+\delta_{g,N}\le\varepsilon_N.}
$$

函数值使用包含真实值的区间 $[\underline J,\overline J]$ 做下降判定，并把浮点舍入误差计入区间；普通双精度“差值为零”不等于精确证书。

## 8. 固定时域校正：iLQR 主方向，梯度后备

### 8.1 真实梯度

在当前轨迹上记

$$
\widetilde A_k=F_x+F_uD\mu_k,\qquad
\widetilde B_k=F_u,\qquad
c_k(x,v)=\ell(x,\mu_k(x)+v).
$$

则 $c_x=2Qx+2D\mu_k^TRu$，$c_v=2Ru$。精确伴随满足

$$
p_N=\nabla T_f(x_N),\qquad
p_k=c_{x,k}+\widetilde A_k^Tp_{k+1},\qquad
g_k=c_{v,k}+\widetilde B_k^Tp_{k+1}.
$$

实际用有限尾的导数替代，并加上第 7 节误差界。此 $p_k$ 是真实梯度伴随；不要和下面优化二次模型时的价值梯度混写。

### 8.2 Riccati 二次方向

在物理代价上采用 Gauss–Newton 曲率，局部矩阵为

$$
C_{xx}=2Q+2D\mu^TRD\mu,\quad C_{vv}=2R,\quad C_{vx}=2RD\mu.
$$

它保留真实一阶项，忽略造成不定性的动力学二阶项。尾部使用有限反馈 rollout 的半正定 Gauss–Newton 曲率 $P_N$；不需要把可能不定的真实 $D^2T_f$ 硬当正定模型。

初始化二次模型价值梯度 $s_N=\nabla T_{f,H}(x_N)$，终端曲率取上述 $P_N$。记下一时刻的模型价值梯度为 $s_{k+1}$，反向递推

$$
q_x=c_x+\widetilde A^Ts_{k+1},\quad q_v=c_v+\widetilde B^Ts_{k+1},
$$

$$
Q_{xx}=C_{xx}+\widetilde A^TP_{k+1}\widetilde A,\quad
Q_{vv}=C_{vv}+\widetilde B^TP_{k+1}\widetilde B,\quad
Q_{vx}=C_{vx}+\widetilde B^TP_{k+1}\widetilde A.
$$

以 $\widetilde Q_{vv}=Q_{vv}+\lambda I\succ0$ 求

$$
d_k=-\widetilde Q_{vv}^{-1}q_v,\qquad
K_k=-\widetilde Q_{vv}^{-1}Q_{vx},
$$

$$
s_k=q_x-Q_{xv}\widetilde Q_{vv}^{-1}q_v,\qquad
P_k=Q_{xx}-Q_{xv}\widetilde Q_{vv}^{-1}Q_{vx}.
$$

若正则化加入整个局部目标，上式对应正则化后的模型；预测值也应使用同一模型。正则化有上界或方向检查，不允许依靠 $\lambda\to\infty$ 制造成功残差。

从 $\delta x_0=0$ 线性前推

$$
h_k=d_k+K_k\delta x_k,\qquad
\delta x_{k+1}=\widetilde A_k\delta x_k+\widetilde B_kh_k.
$$

得到完整方向 $h\in \mathcal H_N$。试探策略明确取 $v^{\rm trial}=v+\alpha h$，再做原系统 rollout。这样线搜索路径就是参数空间直线；不需要隐含假设标准非线性反馈 trial 路径具有统一二阶界。

稳定敏感性和第 5.3 节的正下界给出约化模型矩阵

$$
b_-I\preceq M_N\preceq b_+I,\qquad b_->0,
$$

常数与 $N$ 无关。因此精确模型方向 $h=-M_N^{-1}g$ 满足

$$
-g^Th\ge\|g\|^2/b_+,\qquad \|h\|\le\|g\|/b_-.
$$

实现需定量检查 $-\widehat g^Th\ge\|\widehat g\|^2/b_+$ 和 $\|h\|\le\|\widehat g\|/b_-$，其中常数与 $N$ 无关；仅检查负内积不够。将界保守扩为 $b_-\le1\le b_+$，则 $h=-\widehat g$ 一定满足相同的方向检查，可作为模型或数值失败时的后备。这些检查是计算条件，不是用来替代真实残差的停止指标。

### 8.3 有限精度下的实际接受规则

令 $s=-\widehat g^Th>0$。增加尾计算长度，直到

$$
\delta_g\le\varepsilon_N/4,\qquad
\delta_g\|h\|\le s/4.
$$

精化尾长度或算术精度后，重新计算 $\widehat g$、方向 $h$、$s$ 和相应误差界；不能沿用精化前尚未验证的方向关系。对回溯步长 $\alpha$，检查试探点在认证球内且真实 rollout 有定义，并把两个值区间的总宽度降到至多 $\alpha s/8$。仅当

$$
\overline J(v+\alpha h)\le\underline J(v)-\alpha s/4
$$

才接受。梯度门通过则直接结束，不继续要求严格下降。

为何这个条件可以实现：真实方向满足 $g^Th\le-3s/4$；$L$ 光滑性给出

$$
\mathcal J(v+\alpha h)-\mathcal J(v)
\le-3\alpha s/4+L\alpha^2\|h\|^2/2.
$$

当 $\alpha\le s/(2L\|h\|^2)$ 时，右端不超过 $-\alpha s/2$。加上值区间总宽度 $\alpha s/8$ 后，仍可通过上述更宽松的 $-\alpha s/4$ 接受门。尾误差几何衰减，故每次所需精度用有限尾长度可达到；算术精度也需随证书要求调整。

## 9. 完整局部定理与逐步证明

本节假设第 5.2 节的局部认证入口成立，固定使用该参考反馈和认证球；所有校正由第 8 节的下降规则接受。

### 9.1 存在性和不会碰到人工边界

由二阶积分余项，对于 $r=\|v\|\le r_{\mathrm B}$，

$$
\mathcal J(v)\ge\mathcal J(0)-\eta r+\frac m2r^2.
$$

故球边界的值严格大于中心值，所有下降迭代满足

$$
\mathcal J(v)\le\mathcal J(0)\Longrightarrow
\|v\|\le r_s:=2\eta/m<r_{\mathrm B}.
$$

球是 Hilbert 空间中的弱紧集；$\mathcal J$ 在球上连续、凸，因此将其与球的指标函数相加得到弱下半连续函数，故在球上取得极小值。边界不可能取到极小值，强凸性给出唯一内部点 $v^\star$，满足 $\nabla\mathcal J(v^\star)=0$。

同样，对每个 $N\ge N_0$，$\mathcal B_{r_{\mathrm B}}\cap \mathcal H_N$ 上存在唯一内部极小点 $v_N^\star$，满足 $\mathsf P_N\nabla\mathcal J(v_N^\star)=0$。**这些解的存在和内部性由本节证明，并未预设一条收敛分支。**

### 9.2 每个固定时域都能有限达到正阈值

先看精确计算。方向界、$L$ 光滑性、梯度界 $\|g\|\le\eta+Lr_{\mathrm B}$ 及内部裕量 $r_{\mathrm B}-r_s>0$ 共同给出与 $N$ 无关的足够小步长；几何回溯因此有统一下界 $\alpha_{\min}>0$。Armijo 给出某个 $a>0$：

$$
\mathcal J_N(v^+)\le\mathcal J_N(v)-a\|g_N(v)\|^2.
$$

目标下界为零，所以只要 $\|g_N\|>\varepsilon_N>0$，便不可能无限次都减少至少 $a\varepsilon_N^2$。因此固定时域在有限次校正后过门。

有限精度时，门未通过且 $\delta_g\le\varepsilon_N/4$ 蕴含 $\|\widehat g\|>3\varepsilon_N/4$。模型谱界给 $s\ge\|\widehat g\|^2/b_+$；第 8.3 节保证真代价至少减少 $\alpha s/4$。统一步长下界和相同的下有界论证仍给有限终止。这里依赖可提高尾计算和算术精度；固定双精度只能支持达到其可认证的有限阈值。

### 9.3 真实残差控制有限维解误差

设返回 $v_N\in \mathcal H_N$ 且 $\|\mathsf P_N\nabla\mathcal J(v_N)\|\le\varepsilon_N$。由强单调性，

$$
m\|v_N-v_N^\star\|^2
\le\langle \mathsf P_N\nabla\mathcal J(v_N),v_N-v_N^\star\rangle
\le\varepsilon_N\|v_N-v_N^\star\|.
$$

若差为零结论显然；否则约去一项，得到

$$
\|v_N-v_N^\star\|\le\varepsilon_N/m.
$$

### 9.4 有限维最优解逼近无限维最优解

令 $e=v_N^\star-v^\star$。$\mathsf P_Nv^\star$ 在认证球内。利用有限维驻点条件、无限维驻点条件和梯度 Lipschitz 性，

$$
\begin{aligned}
m\|e\|^2
&\le\langle\nabla\mathcal J(v_N^\star)-\nabla\mathcal J(v^\star),e\rangle\\
&=\langle\nabla\mathcal J(v_N^\star)-\nabla\mathcal J(v^\star),\mathsf P_Nv^\star-v^\star\rangle\\
&\le L\|e\|\|(I-\mathsf P_N)v^\star\|.
\end{aligned}
$$

第二行因为 $v_N^\star-\mathsf P_Nv^\star\in \mathcal H_N$，其与有限维梯度的内积为零。于是

$$
\boxed{\|v_N-v^\star\|_2\le
\frac{\varepsilon_N}{m}+\frac Lm\|(I-\mathsf P_N)v^\star\|_2.}
$$

取任何 $N_j\to\infty$ 和 $\varepsilon_{N_j}\to0$，由于 $v^\star\in\ell^2$，右端趋零。因此整个检查点序列在 $\ell^2$ 中收敛，不只是子序列，不需另外假设有限分支收敛。

### 9.5 回到原系统的局部最优性

轨迹映射连续给出 $(X(v_N),U(v_N))\to(X(v^\star),U(v^\star))$ 于 $\ell^2\times\ell^2$，故任意固定前缀收敛。又由强凸性，

$$
\mathcal J(v)\ge\mathcal J(v^\star)+\frac m2\|v-v^\star\|^2.
$$

对附近任何满足原动力学的有限能量轨迹，令 $v_k=u_k-\mu_k(x_k)$。因为 $D\mu$ 有界，该逆变换在状态—控制 $\ell^2$ 拓扑下连续；足够接近的轨迹落入上述球内，并由同一递推唯一生成。因此所得结论是**原问题在相应可行状态—控制轨迹邻域中的严格局部最优性**。它不是原开环控制任意 $\ell^2$ 扰动下的光滑性，也不是全局最优性。

## 10. 何时停止：直接认证无限时域局部精度

有限前缀残差小，仍需考虑尚未开放的变量。对 $k\ge N$，因 $v_k=0$，相应全梯度分量为

$$
q_f(x_k)=\ell_u(x_k,\kappa_f(x_k))
+F_u(x_k,\kappa_f(x_k))^T\nabla T_f(F_f(x_k)).
$$

注意这里是对独立变量 $v_k$ 求导，故直接项是 $\ell_u$；反馈对状态的导数已经包含在 $T_f$ 中。第 7 节取 $H=0$ 给出 $\|\nabla T_f(z)\|\le C_T\|z\|$，其中 $C_T=C_rC_xC_d/(1-q_xq_d)$。若 $\|\kappa_f(z)\|\le K_f\|z\|$ 且 $\|F_u\|\le B_f$，可以取 $C_q=2\|R\|K_f+B_fC_TC_xq_x$，从而 $\|q_f(z)\|\le C_q\|z\|$。于是

$$
\left(\sum_{k=N}^\infty\|q_f(x_k)\|^2\right)^{1/2}
\le E_{\rm tail}(x_N):=\frac{C_qC_x}{\sqrt{1-q_x^2}}\|x_N\|.
$$

因此可以用有限计算得到完整无限梯度证书

$$
\boxed{\|\nabla\mathcal J(v_N)\|_2\le
\zeta_N:=\sqrt{(\|\widehat g_N\|_2+\delta_{g,N})^2+E_{\rm tail}(x_N)^2}.}
$$

尾界较粗时，可以显式计算更多反馈时刻的 $q_f$，再只对剩余后缀求几何界。这也可用于第 5.2 节的中心全梯度认证。

在已经证明的局部曲率条件下，

$$
\|v_N-v^\star\|_2\le\zeta_N/m,\qquad
0\le\mathcal J(v_N)-\mathcal J(v^\star)\le\zeta_N^2/(2m).
$$

第二式由 $\mathcal J(v^\star)\ge\mathcal J(v_N)+\langle g,v^\star-v_N\rangle+(m/2)\|v^\star-v_N\|^2$，对右侧线性二次式配方得到。它是该已认证局部最优点的性能差上界，不是与全局最优值的差。

当 $\zeta_N$ 达到用户指定精度即可返回策略，无需预设巨大的优化时域。第 9 节的 $\ell^2$ 轨迹收敛还给出 $x_{N_j}(v_{N_j})\to0$：其范数不超过 $\|X(v_{N_j})-X(v^\star)\|_2+\|x_{N_j}(v^\star)\|$。所以有限前缀门趋零时，完整证书也趋零；这不是从“末端始终在终端集”直接推出的。

## 11. 时域选择、总算法与工作量

### 11.1 总体状态机

算法的数据流只有一条：

$$
\text{终端准备}
\longrightarrow
\text{阶段一短块证书}
\longrightarrow
\text{阶段间局部入口}
\longrightarrow
\text{阶段二固定时域校正与几何延拓}
\longrightarrow
\text{完整无限梯度认证}.
$$

任一环节失败都返回该环节之前最后保存的合法对象。终端准备或短块搜索失败时尚未取得稳定尾策略；局部入口失败时已有可执行策略但没有局部最优证书；阶段二预算耗尽时返回最后接受策略及当前误差界。第 11.2—11.8 节按这个状态机展开。

### 11.2 算法主线：先取得可执行入口，再在同一目标上延拓

整个方法不是“不断增大一个带移动终端罚的有限时域问题”，而是下面两个职责不同的阶段。

1. **阶段一只解决可执行性。** 从给定初值出发，用真实非线性 rollout 逐块寻找 CLF 收缩，直到轨迹进入局部切换集。它交付的是“有限前缀 $+$ 稳定反馈尾部”的有限代价策略，不声称该策略已经驻定。
2. **阶段二解决局部优化与精度认证。** 先把阶段一策略写成参考反馈 $\mu_k$，再以稳定反馈扰动 $v_k$ 为优化变量。所有有限时域问题都是同一个真实无限时域目标 $\mathcal J$ 在嵌套子空间 $\mathcal H_N$ 上的限制；增加 $N$ 只是开放新的自由变量，不更换目标。

两阶段之间必须设置局部入口门。仅仅到达 $\mathcal X_f$ 只能说明稳定尾部可执行，不能推出第二阶段所需的光滑性、强凸性或局部最优性。

### 11.3 输入、离线准备与运行时状态

算法输入包括原始动力学 $F$、阶段代价 $\ell$、初值 $x_0$、最终无限梯度目标 $\varepsilon_{\rm goal}>0$、固定时域强迫序列

$$
\varepsilon_N=\frac{c_\varepsilon}{(N+1)^p},\qquad c_\varepsilon>0,\quad p>0,
$$

几何增长因子 $r>1$，以及计算预算。计算预算只决定软件何时返回当前结果，不进入理想算法的无穷延拓证明。

运行前先完成以下准备。

1. 在平衡点线性化并求 DARE，得到 $P$ 和局部反馈 $\kappa_f$；定义 $W(x)=x^TPx$、$V_f=\gamma W$。
2. 在 $\mathcal X_f=\{W\le\alpha_f\}$ 上验证真实非线性闭环不变性和 Lyapunov 下降，并建立尾部状态、导数及代价的几何界。由此才能有限计算 $T_f$、$\nabla T_f$ 及其误差 $E_0,E_1$。
3. 选择严格位于终端认证区内部的切换阈值 $0<\alpha_{\rm sw}<\alpha_f$，以及阶段一参数 $\bar\rho$、$h_{\max}$ 和中间增长系数。若使用有限精度 CLF 门，还要满足第 4.10 节的误差地板条件。

算法始终保存“最后一个已经通过全部相应门槛的对象”。阶段一保存包含最新合格块的整个真实前缀；阶段二保存最后一个经真实下降验证的无限策略。试探计算失败不能覆盖这些已保存对象。

### 11.4 阶段一：证书优先的短块扩展

本节只说明第一阶段在总状态机中的调用位置；完整求解、假设、定理与证明见第 4.11—4.15 节。

初始化 $N_0=0$、当前前缀为空、$z_0=x_0$ 和 $\beta_0=W(x_0)$。若 $W(x_0)\le\alpha_{\rm sw}$，可直接把空前缀与 $\kappa_f$ 尾部拼接，然后进入阶段间认证；否则对块编号 $b=0,1,\ldots$ 重复以下过程。

**第一步：搜索候选块。** 对整数长度 $h=1,\ldots,h_{\max}$，从当前终点 $z_b$ 出发调用第 4.2—4.7 节的两层求解器：先用最小违反量 SQP/SCP 找真实零违反块，再可选地在硬证书内改善性能。候选只有在重新用原动力学 rollout 后满足

$$
W(z_h)\le\bar\rho W(z_b),\qquad
W(z_t)\le cW(z_b),\quad 1\le t<h,
$$

才有资格提交。默认策略是按块长从短到长搜索并提交第一条通过的证书；若并行评估多个长度，也可以按预先规定的性能指标选取，但选择规则不能删除硬证书。遍历全部允许长度仍没有合格块时，返回 certificate_not_found。它只表示本次有限搜索没有找到证书，不表示原系统不可达。

**第二步：提交并保存见证。** 将合格块拼到当前前缀后，立即保存这条由真实动力学生成的见证轨迹，并令 $\beta_{b+1}=\bar\rho\beta_b$。此时已经得到

$$
W(x_{N_{b+1}})\le\beta_{b+1}.
$$

**第三步：可选地改善前缀性能。** 可以在当前固定前缀长度上做有限次 iLQR，但每个试探步必须同时满足：冻结的阶段一性能目标真实下降、真实 rollout 有定义、当前终点仍满足 $W(x_{N_{b+1}})\le\beta_{b+1}$。校正失败就恢复刚保存的块，不把性能改善失败解释成证书失败，也不要求改动后的历史轨迹继续满足过去每一个块端点的收缩式。

**第四步：判断是否切换。** 若当前终点尚未满足 $W(x_{N_{b+1}})\le\alpha_{\rm sw}$，就在新终点继续搜下一个块。若已经满足，则记最终前缀长度为 $N_0$，在其后接入 $\kappa_f$。至此得到第一条已认证的有限代价无限策略。有限命中依赖“每次都能在 $h_{\max}$ 内找到并验证合格块”的条件，而不是由局部求解器名称自动保证。

### 11.5 阶段间接口：把安全策略变成局部优化坐标

由阶段一保存的轨迹 $(\bar x_k,\bar u_k)_{k<N_0}$ 构造一次并固定参考反馈

$$
\mu_k(x)=
\begin{cases}
\bar u_k+K_k^{\rm ref}(x-\bar x_k),&k<N_0,\\
\kappa_f(x),&k\ge N_0,
\end{cases}
$$

并改用 $u_k=\mu_k(x_k)+v_k$。中心 $v=0$ 就是刚才保存的安全策略。随后必须在某个球 $\mathcal B_{r_{\mathrm B}}$ 上验证：真实轨迹映射适定且为 $C^2$、导数有界；同一无限目标满足

$$
0<mI\preceq D^2\mathcal J(v)\preceq LI;
$$

中心完整无限梯度具有可计算上界 $\eta$，且

$$
\eta<mr_{\mathrm B}/2.
$$

若这些条件未通过，返回阶段一的前缀加反馈尾部，并报告 feasible_policy_local_certificate_unavailable。此时已经有可执行的有限代价策略，但没有资格启动后续局部收敛和最优性定理。若入口通过，则下降子水平集严格位于球内，第二阶段不需要把该球当作原问题的硬约束。

### 11.6 阶段二外循环：开放时域并保持策略不变

初始化 $v=0$ 和 $N=\max\{1,N_0\}$。在任一当前时域，只允许前 $N$ 个扰动分量非零，即 $v\in\mathcal H_N$；$k\ge N$ 时执行固定反馈尾部。相应目标为

$$
\mathcal J_N(v)=\sum_{k=0}^{N-1}\ell(x_k,\mu_k(x_k)+v_k)+T_f(x_N).
$$

对当前 $N$，先运行第 11.7 节的固定时域校正器，直到有限前缀真实残差门通过。然后计算第 10 节的完整无限梯度上界 $\zeta_N$。

- 若 $\zeta_N\le\varepsilon_{\rm goal}$，停止并返回当前“优化前缀 $+$ 固定反馈尾部”、残差证书 $\zeta_N$，以及

  $$
  \|v-v^\star\|_2\le\zeta_N/m,
  \qquad
  0\le\mathcal J(v)-\mathcal J(v^\star)\le\zeta_N^2/(2m).
  $$

  这里的 $v^\star$ 是已认证邻域中的唯一严格局部最优点，不是全局最优解。
- 若 $\zeta_N>\varepsilon_{\rm goal}$，令

  $$
  N^+=\max\{N+1,\lceil rN\rceil\},
  $$

  并把 $v$ 零填充到 $\mathcal H_{N^+}$。由于

  $$
  \mathcal J_{N^+}(v,0)=\mathcal J_N(v),
  $$

  零填充前后的真实无限策略、状态轨迹和目标值完全相同；它只把时刻 $N,\ldots,N^+-1$ 的扰动从固定为零改成可优化变量。然后令 $N\leftarrow N^+$，降低阈值到 $\varepsilon_N$，再次调用固定时域校正器。

### 11.7 固定时域校正器：一次迭代的完整顺序

固定 $N$ 后，校正器重复以下步骤，直到残差门通过。

1. **真实 rollout。** 用 $u_k=\mu_k(x_k)+v_k$ 在原非线性系统上生成前缀，并从 $x_N$ 开始 rollout 固定尾部。所有 Jacobian 都在这条当前真实轨迹上重新计算。
2. **精化尾部和误差界。** 选择有限尾长 $H$，计算 $T_{f,H}$、$\nabla T_{f,H}$、伴随梯度估计 $\widehat g_N$、函数值区间及 $\delta_{g,N}$。若误差尚不足以判定残差或下降，就增加 $H$ 或算术精度后全部重算。
3. **先查停止门。** 若

   $$
   \|\widehat g_N\|_2+\delta_{g,N}\le\varepsilon_N,
   $$

   则当前固定时域校正完成。不能用小状态变化、小代价变化或大正则化后的微小步长替代这个门。
4. **构造方向。** 用第 8.2 节的 Gauss--Newton/Riccati 递推得到 iLQR 方向 $h$，并检查统一的定量下降和范数界。若模型不定、数值失败或方向检查不通过，则改用 $h=-\widehat g_N$；正则化只能改善模型，不能充当收敛证书。
5. **使方向判断可认证。** 令 $s=-\widehat g_N^Th$。精化计算直到

   $$
   \delta_g\le\varepsilon_N/4,
   \qquad
   \delta_g\|h\|\le s/4.
   $$

   精化后必须重算梯度、方向和 $s$。
6. **回溯并用真值区间验收。** 对试探点 $v^{\rm trial}=v+\alpha h$ 做原系统 rollout，检查它仍在认证球内且轨迹有定义，并缩小新旧目标值区间，直到可以验证

   $$
   \overline J(v+\alpha h)\le\underline J(v)-\alpha s/4.
   $$

   通过后才令 $v\leftarrow v+\alpha h$ 并保存策略；否则缩小 $\alpha$。在第 5.2 节条件下，总存在足够小的可接受步长；实际程序若因精度或预算无法找到，只能报告校正/数值证书失败并返回最后已保存策略。

固定时域过门后，$\|\mathsf P_N\nabla\mathcal J(v)\|\le\varepsilon_N$。这只控制已经开放的变量，所以随后还必须计算

$$
\zeta_N=\sqrt{(\|\widehat g_N\|+\delta_{g,N})^2+E_{\rm tail}(x_N)^2}
$$

来决定整个算法是否真正达到无限时域精度，不能只凭有限前缀残差停止。

### 11.8 完整伪代码与退出语义

    输入：F, ell, x0, 终端证书, epsilon_goal, epsilon_N, r, 计算预算

    离线准备：
      构造 P, kappa_f, X_f, T_f 的有限计算界和切换集。
      若终端不变性、Lyapunov 下降或尾部误差界无法验证：
          返回 terminal_certificate_unavailable。

    阶段一：安全前缀获取
      prefix <- empty; z <- x0; beta <- W(x0)
      while W(z) > alpha_sw:
          block_found <- false
          for h=1,...,h_max:
              用多组初值运行最小违反量 SQP/SCP。
              每个增量都通过真实 rollout 和违反量接受比筛选。
              若找到真实零违反块：
                  立即保存 certificate_block。
                  accepted_block <- certificate_block
                  可选运行硬约束性能 SQP；若成功则更新 accepted_block。
                  性能求解失败则保持 accepted_block=certificate_block。
                  block_found <- true; break。
          若 block_found=false：返回 certificate_not_found。
          prefix <- append(prefix, accepted_block)
          z <- 真实前缀的当前终点
          beta <- rho_bar * beta；保存整个合格 prefix。
          可选有限次前缀 iLQR；只接受真实下降且终点不超过 beta 的步。
          失败则恢复刚保存的合格 prefix；成功后更新 z。
      N0 <- length(prefix); safe_policy <- prefix + kappa_f tail

    阶段间认证：
      由 safe_policy 构造固定参考反馈 mu，令 v <- 0。
      验证可行球内轨迹适定性、C2 性、一致曲率 m,L 和中心全梯度界。
      若 eta >= m*r_B/2 或其他条件未认证：
          返回 safe_policy，状态 feasible_policy_local_certificate_unavailable。

    阶段二：同一无限目标上的时域延拓
      N <- max(1,N0)
      loop:
          repeat:                                      # 固定时域 corrector
              真实 rollout，计算有限尾、g_hat、delta_g 和目标值区间。
              若 norm(g_hat)+delta_g <= epsilon_N：break。
              计算并检查 Riccati 方向；失败则用负梯度方向。
              精化误差，回溯，并只接受经区间证明的真实下降。
              保存最后接受策略。
          计算完整无限梯度上界 zeta_N。
          若 zeta_N <= epsilon_goal：
              返回当前策略、zeta_N 和局部距离/性能差证书。
          若预算用尽：
              返回最后保存策略和当前误差界，状态 budget_reached。
          N <- max(N+1,ceil(r*N)); v <- zero_pad(v,N)

terminal_certificate_unavailable、certificate_not_found、feasible_policy_local_certificate_unavailable、校正/数值失败和 budget_reached 是不同状态。前两者尚未建立阶段一所需的稳定尾策略；局部入口失败时已有可执行策略但没有局部最优证书；预算用尽只表示计算停止。它们都不能写成 converged。

### 11.9 为什么这条流程能够延拓

这一版用固定几何增长即可保证理论延拓，无需把估计收敛率、候选兼容门、$N+1$ 回退和常数 corrector 次数一起加入主定理。其逻辑链是：阶段一的终端证书产生一个可执行中心；局部入口保证所有下降迭代留在统一可行球内；固定 $N$ 时，真实残差大于正阈值就产生可认证的真实下降，因此有限次通过门；零填充保持同一无限策略不变；几何增长使 $N_j\to\infty$。之后可以再研究残差驱动的增长优化。

### 11.10 可以证明到什么复杂度

精确或具有相对误差控制的计算下，充分下降与强凸性给出目标间隙一致线性收缩：存在 $q\in(0,1)$，

$$
\mathcal J_N(v^{i+1})-\mathcal J_N(v_N^\star)
\le q[\mathcal J_N(v^i)-\mathcal J_N(v_N^\star)].
$$

具体地，精确计算的充分下降常数为 $a$，而 $\|g_N\|^2\ge2m(\mathcal J_N-\mathcal J_N^\star)$，可取不超过一的收缩因子 $1-2ma$，必要时保守减小 $a$。有限精度方向相对误差和接受式只改变常数。

所有阶段起始间隙不超过 $\mathcal J(0)$。由梯度 Lipschitz 和强凸性，

$$
\|g_N(v)\|^2\le L^2\|v-v_N^\star\|^2
\le(2L^2/m)[\mathcal J_N(v)-\mathcal J_N(v_N^\star)].
$$

因此 $\varepsilon_N=c(N+1)^{-p}$ 对应 $O(\log N)$ 次校正。固定维数下一次 Riccati 及 rollout 需要 $O(N+H)$ 次阶段运算；尾误差几何衰减，在所需多项式精度下 $H=O(\log N)$。几何和给出至 $N_{\rm final}$ 的阶段运算上界

$$
\boxed{O(N_{\rm final}\log N_{\rm final}).}
$$

这是第二阶段、固定模型维数和统一局部常数下的实数运算量估计；不包含第一阶段非凸找块的最坏复杂度，也不包含高精度区间运算的位复杂度。不能称为整个任意非线性控制问题的普遍复杂度。

若另能认证每次零填充的真实残差与新阈值之比一致有界，则每阶段可有一致有界校正次数，进而得到 $O(N_{\rm final})$ 的加强结果。该条件并不由 Lyapunov 下降自动给出。

与每个整数时域都完整校正相比，几何增长避免累加全部 $1+2+\cdots+N$。与直接解最终大时域相比，两者可能同阶；本文不证明普遍更快。延拓的明确价值是无需预知所需时域、可提前按尾证书停止、始终保留可执行策略，并便于利用热启动。

## 12. 当前数值证据与实现状态

### 12.1 当前无输入约束代码的复现

运行 scripts/audit_continuation_claims.py，当前结果见 results/theory_audit_20260906_unconstrained.json。该脚本记录源码哈希、回归测试及无输入限制的终端门反例；额外的旧有界输入复现明确标记为 legacy，不混入当前模型。

当前 representative 的 control_limit=None，最大时域 900、bootstrap growth 1.25、每时域最多 80 次校正：

| 指标 | 无输入约束实测 |
|---|---:|
| 原算法退出状态 | maximum_horizon_reached |
| 首次终端命中时域 | 450 |
| 最终预条件残差 $\chi$ | 0.00628376 |
| 代码门阈值 $\tau$ | 0.03331483 |
| 原有限时域目标真实梯度最大分量 | 1.06986556 |
| 最坏梯度坐标 $k=0$ 的控制 | -41.50641498 |
| 无约束一步解析最优控制 | -34.77348589 |
| 该一步解的真实梯度 | 0（数值精度内） |

一步解和长时域控制均可超出旧幅值上限，确认当前模型不再裁剪输入。取消输入范围并没有自动修复旧代码的 $\chi$ 精度门；当前无约束驻点诊断直接使用真实伴随梯度。该次运行到达的是计算时域上限，不是无限时域最优性证书。

固定时域的大阻尼反例仍成立：初始阻尼 $10^{12}$ 可使求解器一步报告 converged，但真实梯度约为 $0.132893$。当前 16 项回归测试全部通过；其中新增检查覆盖默认无输入限制与无约束一步解析解。

2026-09-05 的旧审计和历史基准保持原样。那些结果包含非仿射输入或显式有界输入设置，只用于历史比较，不作为当前无输入约束问题的实验证据。

### 12.2 新路线的独立非线性原型

`scripts/stable_tail_prototype.py` 实现

$$
x^+=1.1x+0.02x^3+u,\qquad \ell=x^2+u^2,\qquad
\kappa_f(x)=-0.6x,\qquad x_0=0.04.
$$

稳定变量下 $x^+=0.5x+0.02x^3+v$，用 $\|v\|_2\le0.05$ 作为局部证明邻域，物理输入 $u\in\mathbb R$ 没有硬幅值约束。它使用固定反馈尾 rollout、真实伴随、Gauss–Newton Riccati 方向和严格下降区间验证；尾部变量不参与优化。该例从局部区域开始，**没有验证非局部 CLF 搜索阶段**。

在 $|x|\le1$ 上可用 $q_x=0.52,q_d=0.56$。记

$$
s=|x_0|+\frac{0.05}{\sqrt{1-q_x^2}}<1,\qquad
\Lambda=\frac{2(1.36s+0.6\cdot0.05)}{1-q_d}.
$$

由首次越界反证，卷积上界保证所有球内扰动均保持 $|x_k|\le s<1$；由轨迹界可推得 $|u_k|\le0.6s+0.05$，这只是该邻域内的派生估计，不是原问题输入上限，也不用于限幅。具体地，在任何假定的首次越界时刻之前，$|x_k|\le q_x^k|x_0|+\sum_{t<k}q_x^{k-1-t}|v_t|\le s$，与首次越界矛盾。相同卷积估计还给出 $\|X(v)\|_2\le|x_0|/\sqrt{1-q_x^2}+\|v\|_2/(1-q_x)$，故真实轨迹具有有限能量。

对方向 $h$，一阶变分为 $z_{k+1}=A_kz_k+h_k$、$z_0=0$，$|A_k|\le q_d$，故 $\|DXh\|_2\le\|h\|_2/(1-q_d)$。两个方向对应的二阶变分为 $y_{k+1}=A_ky_k+0.12x_kz_k^{(1)}z_k^{(2)}$、$y_0=0$。利用 $\|z^{(1)}z^{(2)}\|_2\le\|z^{(1)}\|_2\|z^{(2)}\|_2$，得到 $\|D^2X[h_1,h_2]\|_2\le0.12s\|h_1\|_2\|h_2\|_2/(1-q_d)^3$。系数对 $X$ 连续，且上述估计在略大于认证球的邻域仍成立，因此一、二阶变分是连续有界算子；多项式递推的 Taylor 余项与这些卷积界给出 $X$ 的 $C^2$ 性。$U=-0.6X+v$ 随之具有相同正则性。

真实 Hessian 中唯一非线性动力学余项为 $\sum_k p_{k+1}(0.12x_k)(\delta x_k)^2$，其中 $|p_k|\le\Lambda$。故可以取

$$
c_2=\frac{0.12s\Lambda}{(1-q_d)^2},\qquad
m=\frac2{1.36}-c_2>0,
$$

$$
L=2\left[\frac1{(1-q_d)^2}+\left(1+\frac{0.6}{1-q_d}\right)^2\right]+c_2.
$$

这给出了该原型的解析曲率依据，而非仅采样 Hessian。中心全梯度也由有限前缀与后缀界验证满足 $\eta<mr_{\mathrm B}/2$。原型的前缀使用 Gauss–Newton 曲率；尾部使用该标量例中正定的有限尾精确 Hessian，属于第 8 节允许的有界半正定终端模型。具体数值和每次接受步证书记录在 `results/stable_tail_prototype_20260906_unconstrained.json`。

2026-09-06 无输入约束原型运行得到 $m\approx1.42505583$，$\eta\approx0.02833457<mr_{\mathrm B}/2\approx0.03562640$。在 $N=1,2,4,8,16$ 上分别使用 $2,2,2,2,1$ 次校正；所有有限前缀残差达到 $10^{-12}$，$N=16$ 的完整无限梯度上界约为 $1.043\times10^{-8}$。九次接受更新均通过含无限尾余项的区间下降验证。观察到的少量校正不能推为全部时域的一致次数定理。

该原型证明新接口可以落为有限算法，并检验零填充恒等式、梯度和接受规则。它不代表当前二维摆算例已满足局部曲率或入口条件，也不支持全局收敛、普遍常数校正次数或普遍加速结论。

## 13. 后续迁移应按什么顺序

先把新路线作为独立算法实现和证明对象，不将旧 `V_f` 代码的实验标签换成新算法名称。

1. 原算法若继续维护，优先使用无约束真实梯度门并修正停滞状态；若保留终端硬约束，须单列其 KKT，并保留见证与最后可行策略。
2. 第一阶段实现短块最小松弛搜索、真实提交门和当前终点预算；独立测试失败路径。
3. 为当前二维模型建立稳定尾部的值、一阶、二阶敏感性界；确认新坐标下的真实轨迹适定且处于相应局部稳定域。
4. 实现真实尾价值接口、误差区间、稳定坐标 iLQR 和完整无限梯度证书。
5. 验证局部入口及曲率；若验证失败，明确报告所缺条件，不把数值收敛当作证明。
6. 在相同无输入约束模型、初值、实际尾目标和残差精度下，与大固定时域比较实际工作量。

本次交付已完成理论审计、修订算法及证明、旧算法反例复现和局部非线性原型。生产版 CLF—稳定尾部算法与二维摆模型的完整证书尚未实现，不能据本文提前宣称已完成。

## 14. 文献定位与证据边界

- [Noroozi 等：有限步 CLF](https://arxiv.org/abs/1908.09660) 提供有限步下降与块控制的背景；本文仍需单独验证每次在线找到的真实控制块。
- [Roulet 等，JMLR 2025：iLQR/DDP 收敛](https://jmlr.org/papers/v26/22-1271.html) 讨论相应条件下有限时域的驻点及局部收敛。不能直接转用为本文增长维数、终端筛选和无限时域结论。
- [Na 与 Anitescu：非线性动态规划敏感性衰减](https://arxiv.org/abs/1912.06734) 的一致条件为边界敏感性分析提供参考。本文的第 9 节采用同一目标的嵌套空间证明，没有把该文献当作未验证分支收敛的替代。

本笔记中的新定理依赖已逐项列出的条件，证明在文内给出；上述文献用于定位，而不是把相似题目的结论直接移植为证据。
