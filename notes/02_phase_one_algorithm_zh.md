# 第一阶段算法：一般非线性系统的见证保护型终端集获取

> 2026-09-06 研究范围更新：原问题的控制取值为 $u_k\in\mathbb R^m$，不施加输入范围约束。终端集、CLF 进度条件和局部校正邻域属于证书或算法条件，不是原问题的输入限制。本文保留历史路线的职责；当前修订及证明边界以 [06_review_and_revised_method_zh.md](06_review_and_revised_method_zh.md) 为准。

## 1. 研究问题

考虑一般离散时间非线性系统

$$
x_{k+1}=f(x_k,u_k),\qquad x_0=\bar x,\qquad u_k\in\mathbb R^m,
$$

以及固定初值的无折扣无限时域问题

$$
V_\infty(\bar x)
=
\inf_{\boldsymbol u}
\sum_{k=0}^{\infty}
\ell(x_k,u_k),
$$

其中

$$
\ell(x,u)=x^\top Qx+u^\top Ru,
\qquad Q\succ0,\quad R\succ0.
$$

系统只要求写成一般映射 $f(x,u)$，不要求关于控制仿射。阶段指标选择二次型，是为了提供正定性、局部 Riccati 结构以及清晰的有限代价推论；非凸性仍然来自非线性轨迹映射。

算法最终交付的不是一条实际存储到无穷远的控制序列，而是：

1. 一条真实动力学可行的有限控制前缀 $U_N$；
2. 对应状态轨迹 $X_N$；
3. 一个经过证明或验证的局部终端反馈 $\kappa_f$；
4. 一条可实现的无限时域策略

$$
\Pi_N(U_N)
=
(u_0,\ldots,u_{N-1})\oplus\kappa_f.
$$

## 2. 为什么仍然采用两阶段

两阶段不是为了形式复杂，而是因为两个任务的逻辑性质不同。

第一阶段只解决：

> 是否已经找到一条有限前缀，使其进入一个具有安全、可行、有限代价尾部的终端集合？

第二阶段才解决：

> 在不丢失这个安全尾部的前提下，能否逐渐增加时域并减小有限时域驻点误差？

因此第一阶段的结束条件只应是终端集合命中，不能再混入驻点、局部最优或全局最优结论。第二阶段的第一个任务是在首次命中时域上执行安全 corrector，并通过残差门。

定义两个不同的时域：

$$
N_{\mathrm{hit}}
=
\text{首次获得安全终端前缀的时域},
$$

$$
N_{\mathrm{cert}}
=
\text{首次同时满足终端成员关系和残差门的时域}.
$$

数值上二者可能相等，但逻辑上不能合并。

## 3. 基础系统假设

### 3.1 问题层假设

采用以下基础条件：

1. $f(0,0)=0$；
2. $f$ 在原点邻域二阶连续可微；
3. $Q\succ0$，$R\succ0$；
4. 原点线性化

$$
A_0=f_x(0,0),\qquad B_0=f_u(0,0)
$$

中的 $(A_0,B_0)$ 可稳定。

这些都是系统和性能指标本身的条件，不涉及“第几次 iLQR 迭代恰好满足某个性质”。

### 3.2 iLQR 所需的附加正则性

在性能优化所访问的轨迹邻域内，需要 $f$ 的一、二阶导数存在并具有适当的局部 Lipschitz 性。控制直接在 $\mathbb R^m$ 中优化。数值信赖域只约束本次增量以控制模型误差，不是原问题的控制取值范围。

这些条件服务于局部模型、ratio test 和 corrector 分析，不参与第一阶段的有限命中证明。

## 4. 局部终端三元组不是额外假设

### 4.1 DARE 构造

由 $(A_0,B_0)$ 可稳定、$Q\succ0$ 和 $R\succ0$，离散代数 Riccati 方程给出 $P\succ0$ 和 $K_f$。记

$$
A_c=A_0-B_0K_f,
$$

$$
S=Q+K_f^\top RK_f\succ0,
$$

则

$$
A_c^\top PA_c-P=-S.
$$

取任意 $\gamma>1$，定义

$$
\kappa_f(x)=-K_fx,
$$

$$
V_f(x)=\gamma x^\top Px,
$$

$$
\mathcal X_f
=
\{x:V_f(x)\le\alpha_f\}.
$$

### 4.2 非线性余项证明

由 $f$ 在原点二阶连续可微，

$$
f(x,-K_fx)=A_cx+r(x),
$$

并且在足够小的邻域内存在 $c_r>0$ 使

$$
\|r(x)\|\le c_r\|x\|^2.
$$

于是

$$
\begin{aligned}
&V_f(f(x,-K_fx))-V_f(x)+\ell(x,-K_fx)\\
&=-(\gamma-1)x^\top Sx
+2\gamma x^\top A_c^\top Pr(x)
+\gamma r(x)^\top Pr(x).
\end{aligned}
$$

第一项是严格负二次项，后两项分别为三阶和四阶小量。因此存在足够小但非零的 $\alpha_f$，使

$$
V_f(f(x,\kappa_f(x)))-V_f(x)
\le-\ell(x,\kappa_f(x))
$$

对所有 $x\in\mathcal X_f$ 成立。

该终端半径只用于保证非线性余项界及 Lyapunov 下降有效；无需再加入输入幅值条件。

由 $\ell\ge0$ 立即得到

$$
V_f(f(x,\kappa_f(x)))\le V_f(x),
$$

因此

$$
f(x,\kappa_f(x))\in\mathcal X_f.
$$

所以局部终端三元组的存在性可以从系统层条件推出。实际代码仍需使用解析余项界、区间界或其他可靠方法选择具体 $\alpha_f$；单纯在有限采样点上检查不能构成证明。

## 5. 第一阶段采用带裕量终端集

固定 $\delta\in(0,1)$，定义

$$
\mathcal X_f^\delta
=
\{x:V_f(x)\le(1-\delta)\alpha_f\}.
$$

显然

$$
\mathcal X_f^\delta\subset\operatorname{int}\mathcal X_f.
$$

对任意真实 rollout $X_N=(x_0,\ldots,x_N)$，定义

$$
T_\delta(X_N)
=
\min\{k\in\{0,\ldots,N\}:x_k\in\mathcal X_f^\delta\}.
$$

若集合为空，则记 $T_\delta(X_N)=+\infty$。

算法扫描整条轨迹，而不是只检查 $x_N$。一旦中间状态首先命中，就在该位置截断前缀并接入 $\kappa_f$。

## 6. 第一阶段最小的系统层条件

定义有限代价域

$$
\mathcal D_{\mathrm{fin}}
=
\left\{
x:\exists\boldsymbol u,\ 
J_\infty(x;\boldsymbol u)<\infty
\right\}.
$$

第一阶段所需的本质系统条件是

$$
\bar x\in\mathcal D_{\mathrm{fin}}.
$$

这比全局可控、全局稳定反馈或全局控制 Lyapunov 函数弱得多，只要求给定初值存在一条有限代价可行策略。

### 6.1 有限代价与有限步终端可达等价

在已经构造局部终端三元组的条件下，下列两件事等价：

1. $\bar x\in\mathcal D_{\mathrm{fin}}$；
2. 存在有限 $T$ 和可行控制前缀 $U_T^b$，使其真实 rollout 满足

$$
x_T^b\in\mathcal X_f^\delta.
$$

证明分为两个方向。

若存在有限代价无限轨迹，则

$$
\lambda_{\min}(Q)
\sum_{k=0}^{\infty}\|x_k\|^2
\le
J_\infty(\bar x;\boldsymbol u)
<\infty.
$$

因此

$$
x_k\to0.
$$

连续性给出

$$
V_f(x_k)\to0,
$$

所以轨迹必在有限步进入 $\mathcal X_f^\delta$。

反之，若有限前缀进入 $\mathcal X_f^\delta$，则接入 $\kappa_f$ 后

$$
\sum_{k=T}^{\infty}
\ell(x_k,\kappa_f(x_k))
\le V_f(x_T)<\infty.
$$

加上有限前缀代价，即得到有限无限时域代价。

该等价性只证明“存在一条轨迹”，不证明局部 iLQR 能找到它。

## 7. 为什么必须把可达见证作为独立输入

一般非线性有限时域控制序列优化是非凸的。只知道终端集合可达，不能推出从任意初始化出发的局部 iLQR 会进入正确吸引域。

如果希望算法本身具有有限命中保证，就必须采用以下至少一种接口：

1. 直接提供一条有限步可达见证 $U_{T_b}^b$；
2. 提供一个已知在 $\bar x$ 上可行并渐近稳定的基准控制器；
3. 提供一个对所研究状态域完备的可达性求解器。

当前算法采用第一种接口作为理论主形式，也允许第二种接口在线生成见证。

第二种接口的一种可验证系统层充分条件是：存在前向不变区域 $\mathcal D_b$、可行反馈 $\kappa_b$、正定函数 $W_b$ 和正定标量函数 $\alpha_b$，使 $\bar x\in\mathcal D_b$ 且

$$
W_b(f(x,\kappa_b(x)))-W_b(x)
\le-\alpha_b(\|x\|).
$$

在相关 $W_b$ 子水平集紧的条件下，求和得到

$$
\sum_{k=0}^{\infty}\alpha_b(\|x_k\|)<\infty,
$$

从而 $x_k\to0$，所以该反馈在有限步内生成终端可达见证。这仍是系统和基准控制器的性质，不是 iLQR 迭代性质。

给定候选控制序列后，代码重新使用原始非线性动力学 rollout，并检查：

1. 每个控制维数是否正确且数值有限；
2. 所有状态是否有限；
3. 动力学是否由同一 $f$ 递推；
4. 是否在有限步进入 $\mathcal X_f^\delta$。

因此见证不是“相信另一个优化器成功”，而是可以独立验证的有限数据。

若没有提供有效见证，算法只能报告 `reachability_witness_required`、`reachability_witness_invalid` 或 `reachability_witness_not_found`，不能报告系统不可行。

## 8. 第一阶段的双通道算法

### 8.1 证书通道

证书通道独立保存已经验证的

$$
(X_{T_b}^b,U_{T_b}^b).
$$

它只承担：

1. 有限步命中保证；
2. 原动力学一致性与轨迹有限性；
3. 失败时的安全回退策略。

任何性能优化步骤都不能覆盖该数据。

### 8.2 性能通道

性能通道从 $N=1$ 的真实非线性问题开始：

$$
\min_{u_0\in\mathbb R^m}
\ell(\bar x,u_0)+V_f(f(\bar x,u_0)).
$$

与旧版本不同，第一阶段后续不再用 $\kappa_f$ 生成性能延拓。给定已经得到的
$N$ 步真实性能轨迹，在当前终端状态 $x_N$ 上重新直接求解

$$
u_N^g
\in
\arg\min_{u\in\mathbb R^m}
\left\{
\ell(x_N,u)+V_f(f(x_N,u))
\right\},
$$

这里一般非线性 $f$ 的一步问题未必凸。若 $f(x,\cdot)$ 全局连续、$V_f\ge0$，则目标至少为 $\lambda_{\min}(R)\|u\|^2$，因而在 $\mathbb R^m$ 上连续且强制，最小值存在；这不保证局部求解器找到全局解。若只取得数值候选，下述延拓恒等式仍成立，但与任意基准控制的最优值比较仅在已验证目标比较时成立。

在控制仿射专门情形 $f(x,u)=a(x)+B(x)u$ 下，一步问题严格凸，直接解线性方程

$$
[R+\gamma B(x)^\top PB(x)]u^g=-\gamma B(x)^\top Pa(x)
$$

即可得到唯一全局解，无需截断或输入投影。

然后用真实动力学确定

$$
x_{N+1}^g=f(x_N,u_N^g).
$$

定义这一步的可计算延拓缺陷

$$
\Delta_N^g
:=
\ell(x_N,u_N^g)
+V_f(x_{N+1}^g)
-V_f(x_N).
$$

于是延拓前后的增广代价满足精确恒等式

$$
\widehat J_{N+1}(U_N\oplus u_N^g)
=
\widehat J_N(U_N)+\Delta_N^g.
$$

若一次从 $N$ 增长到 $M>N$，则在每一个新生成状态上重复上述一步求解，再进行
有限次 iLQR correction。每个 correction 仍使用原始非线性 rollout，并且只有在

$$
\widehat J_M(U_M^+)\le\widehat J_M(U_M)
$$

时才接受。因此 correction 只会改善下面的代价上界，不参与第一阶段有限命中的
必要条件。

在相同继承终端状态 $x_N$ 上，只要该一步 rollout 定义良好且采用全局一步最小值，一步直接最小化给出

$$
\Delta_N^g
\le
\ell(x_N,\kappa_f(x_N))
+V_f(f(x_N,\kappa_f(x_N)))
-V_f(x_N).
$$

所以它在“一步增广代价增长”这个指标上不劣于 $\kappa_f$。但在
$x_N\notin\mathcal X_f$ 时，两边都不一定非正，不能把这个比较偷换成全局
Lyapunov 下降。

随后可在一列有限的增长时域上局部改进

$$
\min_{U_N\in\mathbb R^{mN}}
\widehat J_N(U_N),
$$

其中

$$
\widehat J_N(U_N)
=
\sum_{k=0}^{N-1}\ell(x_k,u_k)+V_f(x_N).
$$

每个时域只允许有限次 ratio-globalized iLQR 尝试。性能通道的任务是尽量找到：

1. 更早的终端命中；
2. 更低的认证前缀代价；
3. 更适合后续 corrector 的初值。

它不承担无条件的可行性回退证明；该职责仍由独立见证承担。但直接一步缺陷允许
额外建立一个不依赖 iLQR 驻点收敛的条件性有限命中定理。

### 8.3 累计正缺陷与有限命中

定义截至性能时域 $N_j$ 的累计正延拓缺陷

$$
D_j^+
:=
\sum_{\text{已执行的一步延拓}}
[\Delta_k^g]_+,
\qquad
[a]_+:=\max\{a,0\}.
$$

由一步恒等式和 iLQR 接受准则，

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

由于 $Q\succ0$、$V_f\ge0$，

$$
\min_{0\le k<N_j}\|x_k\|^2
\le
\frac{\widehat J_{N_j}(U_{N_j})}
{N_j\lambda_{\min}(Q)}
\to0.
$$

而 $\mathcal X_f^\delta$ 包含原点邻域，所以某条有限的修正后 rollout 必有一个
状态进入 $\mathcal X_f^\delta$。算法扫描整条轨迹并在首次命中处截断，因而得到

$$
N_{\mathrm{hit}}<\infty.
$$

这个条件是充分条件，不是一步求解自动给出的事实。代码同时记录
$D_j^+$ 和 $\widehat J_{N_j}/N_j$；有限实验中二者没有表现出预期趋势时，只能
报告该条件未被数值支持，不能报告系统不可达。

### 8.4 第一阶段终止

第一阶段按以下逻辑运行：

```text
输入：一般非线性问题、终端三元组、裕量 delta、
      经验证的有限步可达见证、有限的性能工作预算

重新 rollout 并验证见证
独立保存见证及其首次命中前缀

求解真实非线性 N=1 性能问题
for 每个不超过见证时域的性能时域:
    扫描整条真实性能 rollout
    若命中 Xf_delta:
        截断并返回该性能前缀

    在每个新增状态上直接求解一步问题
    记录每个 Delta_g、累计正缺陷 D_plus 和 Jhat_N/N
    对直接延拓轨迹重新扫描首次命中

    尝试有限次 iLQR correction
    只接受真实增广代价不增加的 rollout
    每次接受后重新扫描整条 rollout
    correction 失败时保留上一个可行性能轨迹

若性能通道没有更早命中:
    返回独立保存的见证前缀
```

因为见证长度 $T_b$ 有限、性能时域个数有限、每个时域的尝试数有限，所以带见证
版本无条件有限结束。若不使用见证而让性能时域持续增长，则需要
$\widehat J_{N_j}=o(N_j)$ 或 $D_j^+=o(N_j)$ 才能推出有限命中。两个证明都不
需要“iLQR 至少接受一步”或“残差每次收缩”。

## 9. 第一阶段可以严格得到什么

第一阶段返回的轨迹满足

$$
x_{N_{\mathrm{hit}}}\in
\mathcal X_f^\delta
\subset\operatorname{int}\mathcal X_f.
$$

因此接入终端反馈后：

1. 有限前缀满足原始动力学；
2. 所有控制属于 $\mathbb R^m$，无需幅值检查；
3. 尾部状态始终留在 $\mathcal X_f$；
4. 尾部代价满足

$$
J_{\mathrm{tail}}
\le V_f(x_{N_{\mathrm{hit}}});
$$

5. 整条实现策略具有有限无限时域代价。

第一阶段不能推出：

1. 前缀满足有限时域 KKT 条件；
2. 前缀是局部极小值；
3. 前缀是全局极小值；
4. 前缀等于无限时域最优策略；
5. iLQR 对任意可达初值都能找到同一前缀。

直接一步性能通道还能条件性推出：

$$
D_j^+=o(N_j)
\Longrightarrow
\widehat J_{N_j}=o(N_j)
\Longrightarrow
\min_{0\le k<N_j}\|x_k\|\to0
\Longrightarrow
N_{\mathrm{hit}}<\infty.
$$

这条链只说明某条有限 rollout 中存在终端命中状态，不说明第一阶段终端状态序列
$x_{N_j}$ 本身收敛，也不说明命中前缀满足驻点条件。

## 10. 第一阶段假设为什么已经接近最少

不能完全删除可达性条件。若

$$
\bar x\notin\mathcal D_{\mathrm{fin}},
$$

则不存在任何有限代价无限时域策略，任何算法都不可能返回所要求的结果。

当前主条件

$$
\bar x\in\mathcal D_{\mathrm{fin}}
$$

等价于有限步进入 $\mathcal X_f^\delta$，因此既不是人为加强到全局可控，也不是对 iLQR 迭代行为的假设。

但存在性不能自动给出控制序列。为了让可执行算法也具有保证，必须把有限见证作为输入，或额外提供完备可达规划器。对一般连续非线性系统，若既不提供见证，又不提供完备规划器，却声称局部 iLQR 必然找到证书，就是逻辑跳步。

## 11. 第一阶段工作量

验证长度为 $T_b$ 的见证需要

$$
O(T_b)
$$

次动力学计算。

若性能通道访问

$$
1=N_0<N_1<\cdots<N_s\le T_b
$$

并且每个时域最多尝试 $K_b$ 次 iLQR，则额外 Riccati sweep 工作量为

$$
\sum_{j=0}^{s}O(K_bN_j).
$$

密集逐点时域最坏为

$$
O(K_bT_b^2),
$$

几何时域序列为

$$
O(K_bT_b).
$$

每个新增控制还需要求解一个低维真实一步问题。因此，从 $1$ 延拓到 $T_b$ 的
一步求解次数为 $O(T_b)$。控制仿射情形采用固定维数的正定线性求解时，该部分仍为线性工作量。一般非线性无界域上的网格与局部精修只能作为候选生成器，其成本和全局精度需另行说明。

性能工作量可以通过令 $K_b=0$ 关闭所有 iLQR correction，但直接一步延拓和独立
见证回退仍然保留。

## 12. 第一阶段代码字段和退出状态

结果中分别保存：

- `first_hit_horizon`：$N_{\mathrm{hit}}$；
- `certificate_source`：选中的第一阶段通道；
- `reachability_witness_horizon`：独立见证首次命中时域；
- `terminal_hit_margin`：$\delta$。
- `bootstrap_extension_mode`：第一阶段延拓方式，当前主算法为
  `one_step_greedy`；
- `bootstrap_extensions`：实际执行的一步延拓数；
- `bootstrap_cumulative_extension_defect`：$\sum_k\Delta_k^g$；
- `bootstrap_cumulative_positive_defect`：$D_j^+$；
- `bootstrap_average_positive_defect`：$D_j^+/$ 一步延拓数。

主要退出状态：

- `reachability_witness_required`：缺少见证，不表示不可行；
- `reachability_witness_invalid`：见证未通过真实 rollout 验证；
- `reachability_witness_not_found`：基准策略在有限搜索上限内未命中，不表示系统不可行。

第二阶段使用的 `certified_horizon`、残差字段和 corrector 退出状态见
`03_phase_two_algorithm_zh.md`。

## 13. 历史代表性实验中的第一阶段（不作为无输入范围版本结果）

以下数字来自旧输入受限实验，保留为历史记录，不能沿用为本次无输入范围问题的运行结果；无约束版本需要重新运行。理论不依赖该历史系统结构。

当前设置采用

$$
\delta=0.1.
$$

独立逆动力学基准控制器的真实 rollout 在

$$
T_b=601
$$

首次进入 $\mathcal X_f^\delta$。因此即使所有性能 iLQR 步骤失败，第一阶段仍具有有限回退前缀。

采用直接一步延拓、性能时域增长因子 $1.25$ 时，性能通道在

$$
N_{\mathrm{hit}}=489.
$$

该点同时通过第二阶段首个残差门，因此

$$
N_{\mathrm{cert}}=489.
$$

第一阶段接受 $16$ 次性能 correction；在相同增长因子和其他参数下，旧的
$\kappa_f$ 性能延拓接受 $23$ 次。两者的首次命中均为 $489$，所以该实验支持的
结论是“一步直接延拓减少了局部校正工作”，而不是“它必然更早命中”。

直接一步延拓记录的累计正缺陷为

$$
D^+=2.8295\times10^4,
$$

平均每个新增控制约为 $49.90$。在这个有限区间内，该量没有提供
$D_j^+=o(N_j)$ 的数值证据；有限命中事实来自实际 rollout 扫描，算法的无条件
回退仍来自 $T_b=601$ 的独立见证。这个负面诊断必须与工作量改进同时报告。

严格逐步的直接一步延拓，即每次只做 $N\mapsto N+1$ 并最多接受一次 iLQR，
得到

$$
N_{\mathrm{hit}}=N_{\mathrm{cert}}=491,
$$

共接受 $478$ 次性能 correction，累计 Riccati sweep 长度为 $120217$。它与
增长因子 $1.25$ 的结果说明：直接一步延拓改善种子并不消除密集时域的二次级
重复工作，第一阶段仍需要单独选择时域调度。

## 14. 第一阶段最终逻辑链

当前第一阶段的严格结论是

$$
\boxed{
\begin{aligned}
&f\text{ 局部二阶光滑}
+(A_0,B_0)\text{ 可稳定}
+Q,R\succ0\\
&\qquad\Longrightarrow
\text{存在非平凡局部终端三元组},\\
&\bar x\in\mathcal D_{\mathrm{fin}}
\Longleftrightarrow
\text{有限步可达 }\mathcal X_f^\delta,\\
&D_j^+=o(N_j)
\Longrightarrow
\text{直接一步性能通道有限命中},\\
&\text{提供并验证有限可达见证}
\Longrightarrow
N_{\mathrm{hit}}<\infty,\\
&N_{\mathrm{hit}}<\infty
+\text{终端下降}
\Longrightarrow
\text{可行且有限代价的无限时域策略}.
\end{aligned}
}
$$

第一条性能通道结论只使用实际增广代价上界，不使用 iLQR 残差收缩。第二条见证
结论则连 $D_j^+$ 的次线性条件也不需要。局部 iLQR 在第一阶段只负责性能改进，
在第二阶段才负责驻点校正；$\kappa_f$ 也只在命中后承担认证延拓和有限尾部。
