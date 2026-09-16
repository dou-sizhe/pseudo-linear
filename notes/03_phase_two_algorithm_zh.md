# 第二阶段算法：安全约束下的认证时域延拓 iLQR

> 2026-09-06 研究范围更新：原问题的控制取值为 $u_k\in\mathbb R^m$，不施加输入范围约束。终端集、CLF 进度条件和局部校正邻域属于证书或算法条件，不是原问题的输入限制。本文保留历史路线的职责；当前修订及证明边界以 [06_review_and_revised_method_zh.md](06_review_and_revised_method_zh.md) 为准。

本稿要求的真实梯度认证门属于算法规范；旧生产代码仍使用 $\chi_N$ 作为主门，尚未完成该认证逻辑的迁移。本次代码已取消默认输入限幅，二者应分开理解。

## 1. 本文档的作用

本文档只说明第二阶段：如何从第一阶段交付的安全终端前缀出发，在保持可实现
无限尾部的前提下，逐渐增加有限时域并减小驻点残差。

第一阶段的问题假设、终端三元组构造、有限可达见证和
$N_{\mathrm{hit}}<\infty$ 的证明见 `02_phase_one_algorithm_zh.md`。本文不重复
证明这些内容，只把第一阶段结论作为显式输入接口。

## 2. 第二阶段的输入与任务

系统和有限时域目标仍为

$$
x_{k+1}=f(x_k,u_k),
\qquad
\ell(x,u)=x^\top Qx+u^\top Ru,
$$

$$
\widehat J_N(U_N)
=
\sum_{k=0}^{N-1}\ell(x_k,u_k)+V_f(x_N).
$$

第一阶段交付：

1. 一条由原始非线性动力学 rollout 得到的有限前缀
   $(X_{N_{\mathrm{hit}}},U_{N_{\mathrm{hit}}})$；
2. 严格终端成员关系

$$
x_{N_{\mathrm{hit}}}\in
\mathcal X_f^\delta
\subset\operatorname{int}\mathcal X_f;
$$

3. 终端反馈 $\kappa_f$ 和终端 Lyapunov 函数 $V_f$，满足正不变性及

$$
V_f(f(x,\kappa_f(x)))-V_f(x)
\le-\ell(x,\kappa_f(x)),
\qquad x\in\mathcal X_f.
$$

第二阶段承担两个不同任务：

1. 在 $N_{\mathrm{hit}}$ 上获得第一个同时满足安全性和残差门的检查点
   $N_{\mathrm{cert}}$；
2. 从该检查点出发，沿同一局部解分支进行非递减时域延拓。

第二阶段不重新证明终端集合可达，也不声称把局部 iLQR 变成全局优化器。

## 3. 固定时域 iLQR corrector

给定名义控制序列 $\bar U_N$，先用原动力学生成

$$
\bar x_{k+1}=f(\bar x_k,\bar u_k).
$$

在名义轨迹上使用真实 Jacobian

$$
A_k=f_x(\bar x_k,\bar u_k),
\qquad
B_k=f_u(\bar x_k,\bar u_k).
$$

iLQR 保留动力学的一阶模型和价值函数的二次近似。反向 Riccati 递推给出前馈
方向 $d_k$ 与反馈矩阵 $L_k$，trial 控制为

$$
u_k^+
=
\bar u_k+\alpha d_k+L_k(x_k^+-\bar x_k).
$$

每条 trial 状态轨迹都由

$$
x_{k+1}^+=f(x_k^+,u_k^+)
$$

重新 rollout。控制直接取值于 $\mathbb R^m$，不作饱和或投影。算法不把线性化状态当作可行状态，不对状态轨迹作凸组合。

为处理局部模型失真，反向递推在控制 Hessian 上使用
Levenberg--Marquardt 正则化，并通过真实代价的 actual/predicted reduction
ratio 接受或拒绝整步。拒绝时保留上一条可行轨迹并增大正则化。

## 4. 从 $N_{\mathrm{hit}}$ 到 $N_{\mathrm{cert}}$

### 4.1 安全 corrector

从首次命中开始，所有 corrector trial 除了真实代价 ratio test 外，还必须满足

$$
x_N^+\in\mathcal X_f.
$$

第一阶段使用 $\mathcal X_f^\delta$ 为入口提供非零裕量，但该裕量不能保证后续校正始终远离终端边界，也不能保证安全区域含有原有限时域目标的驻点；反例见第 06 号修订稿。取消输入范围约束不消除这个终端边界问题。

若 trial 不满足终端成员关系，则整步拒绝。若所有允许的步长和正则化均失败，
第二阶段停止，但第一阶段保存的安全前缀不被覆盖。

### 4.2 两种残差必须分开

时域候选的工作量预测使用 Riccati 预条件残差

$$
\chi_N(U_N)
=
\max_{0\le k<N}\|d_k^0(U_N)\|,
$$

其中 $d_k^0$ 是在未正则化局部模型上计算的 iLQR 前馈方向；若对应控制
Hessian 非正定或奇异，则该诊断记为不可用。强迫容差取

$$
\tau_N=\frac{\tau_0}{(N+1)^p},
\qquad p>0.
$$

另外取趋零的真实梯度容差 $\varepsilon_N>0$。第一个同时满足

$$
x_N\in\mathcal X_f,
\qquad
\chi_N(U_N)\le\tau_N,
\qquad g_N^{\max}\le\varepsilon_N
$$

的时域定义为 $N_{\mathrm{cert}}$。其中 $g_N^{\max}$ 在下面定义；若仅通过 $\chi_N$ 门，记录为工作量检查点，不能称为原目标驻点认证。

$\chi_N$ 不是原控制坐标下的 KKT 梯度。代码还必须独立报告真实伴随梯度

$$
g_N^{\max}
=
\max_{0\le k<N}
\|2Ru_k+f_u(x_k,u_k)^\top p_{k+1}\|,
$$

其中

$$
p_N=2\gamma Px_N,
$$

$$
p_k
=
2Qx_k+f_x(x_k,u_k)^\top p_{k+1}.
$$

预条件残差只用于估计 corrector 工作量。无输入范围限制的原有限时域问题，其驻点判据为 $\nabla\widehat J_N=0$，实际精度门必须显式检查真实梯度范数。仅通过 $\chi_N$ 门不能称为原问题近似驻点。若将 $c_N(U_N)=V_f(x_N(U_N))-\alpha_f\le0$ 作为硬约束，则应另外报告 $\nabla\widehat J_N+\nu_N\nabla c_N=0$、$\nu_N\ge0$、$\nu_Nc_N=0$ 与可行性；该受终端约束问题的 KKT 点未必是原问题驻点。

### 4.3 入口失败时的结论

若安全 corrector 在 $N_{\mathrm{hit}}$ 上失败：

1. 返回 `trust_region_failure`；
2. 保留第一阶段得到的可行无限时域策略；
3. 令 `certified_horizon` 为空；
4. 不把优化失败改写成系统不可行；
5. 不声称已经获得有限时域驻点或局部最优解。

因此，在第一阶段的最小系统假设下可以保证已验证见证带来
$N_{\mathrm{hit}}<\infty$，但不能保证局部 iLQR 必在有限次后给出
$N_{\mathrm{cert}}<\infty$。

## 5. 终端反馈延拓

对 $M>N\ge N_{\mathrm{cert}}$，定义

$$
\mathcal E_N^M U_N
=
(u_0,\ldots,u_{N-1},
\kappa_f(x_N),\ldots,\kappa_f(x_{M-1})).
$$

所有新增状态均由原始非线性动力学 rollout。由于新增部分位于正不变终端集中，
扩展控制和状态可行，并且终端下降条件逐项给出

$$
\widehat J_M(\mathcal E_N^M U_N)
\le
\widehat J_N(U_N).
$$

这一不等式只证明扩展轨迹具有安全、非增的带终端罚有限时域代价。它不证明扩展
轨迹在新时域上已经驻点，因为新终端边界会产生延拓缺陷。

## 6. 候选时域和 corrector 预算

候选集合取为

$$
\mathcal H(N)
=
\{N+1\}
\cup
\{\lceil r_1N\rceil,\ldots,\lceil r_sN\rceil\},
\qquad r_i>1.
$$

对每个候选 $M$，先计算扩展序列
$\mathcal E_N^M U_N$ 的预条件残差 $\chi_M$，再用经验保守收缩率
$\widehat\theta\in(0,1)$ 估计达到 $\tau_M$ 所需的校正次数：

$$
\widehat K(M)
=
\left\lceil
\frac{\log(\chi_M/\tau_M)}
{\log(1/\widehat\theta)}
\right\rceil_+.
$$

选择满足

$$
\widehat K(M)\le K_{\mathrm{tar}}
$$

的最大候选；若没有几何候选通过预算，则保底尝试 $N+1$。

$\widehat K(M)$ 只是时域选择启发式。无论预测值多小，corrector 后都必须重新
检查真实终端成员关系、$\chi_M\le\tau_M$ 和真实梯度 $g_M^{\max}\le\varepsilon_M$。以下关于 $\chi$ 收缩的历史工作量估计只控制预条件残差门；若要同时保证任意给定的梯度容差在有限次内达到，还需证明两种残差在所访问邻域的定量关系，或使用第 06 号稿的真实梯度充分下降校正。

## 7. 第二阶段完整流程

```text
输入：第一阶段安全前缀、终端三元组、残差强迫序列、
      候选增长因子、每个时域的有限 corrector 预算

在 N_hit 上运行安全 iLQR corrector
若终端成员关系或残差门不能通过:
    返回第一阶段安全策略，并报告尚无 N_cert
否则:
    记录 N_cert

while N < N_max:
    构造候选时域集合 H(N)
    对每个候选 M:
        用终端反馈扩展并以原动力学 rollout
        计算扩展残差和预测 corrector 次数
    选择预算允许的最大 M
    在 M 上运行有限次安全 iLQR corrector
    每个 trial 都检查控制及状态数值有限、真实代价 ratio 和终端成员关系
    若 chi_M <= tau_M 且真实梯度 gmax_M <= epsilon_M:
        接受检查点并继续
    否则:
        保留最近的安全策略并报告对应失败状态
```

达到 $N_{\max}$ 只表示有限实现到达用户设定的计算上限，不能称为无限时域
收敛。

## 8. 有限校正次数需要哪些附加条件

终端 Lyapunov 下降本身不能保证延拓后残差很小，也不能保证 iLQR corrector
收缩。若要证明每个时域只需一致有界的校正次数，至少需要两个独立接口。

第一，延拓缺陷界。例如存在 $C_e>0$ 和 $\beta\in(0,1)$，使

$$
\chi_M(\mathcal E_N^M U_N)
\le
\chi_N(U_N)+C_e\beta^N.
$$

第二，时域一致的局部 corrector 收缩。即在相关局部分支的统一邻域内存在
$\theta\in(0,1)$，使每次被接受的安全 corrector 满足

$$
\chi_M(U_M^+)\le\theta\chi_M(U_M).
$$

还需要扩展初值始终落在该统一收缩邻域内。上述条件涉及解分支的强正则性、
远端边界影响衰减和局部吸引域，不能由 $(A_0,B_0)$ 可稳定或
$Q,R\succ0$ 单独推出。

在这些条件成立时，几何候选时域间 $\tau_N/\tau_M$ 有界，故每阶段所需校正次数
可以一致有界；若时域几何增长，总 Riccati sweep 工作量可为

$$
O(N_{\max}).
$$

若逐点采用 $N\leftarrow N+1$，即使每个时域只需常数次 sweep，累计工作量仍可能
达到

$$
O(N_{\max}^2).
$$

## 9. 渐近结论还需要哪些条件

由每个有限时域检查点通过残差门，只能得到一列有限维近似驻点。若要进一步得到
固定前缀收敛，需要：

1. 局部驻点分支在时域变化下具有强正则性；
2. 远端终端边界对任意固定前缀的影响随 $N$ 衰减；
3. 检查点保持在同一个局部吸引域内；
4. 数值误差和残差强迫 $\tau_N$ 足够小。

若还要把极限称为严格局部最优，需要相应的无限维二阶充分条件。一般非线性
iLQR continuation 不提供全局最优性；全局结论必须来自额外的凸性、精确动态
规划、梯度支配或可验证全局证书。

## 10. 正确性结论的层次

第二阶段必须按以下顺序报告：

$$
\text{第一阶段安全前缀}
\Longrightarrow
\text{可行且有限代价的回退策略},
$$

$$
\text{安全 corrector}
+\|\nabla\widehat J_N\|\le\varepsilon_N
\Longrightarrow
\text{有限时域近似驻点检查点},
$$

$$
\text{一致 corrector 收缩}
+\text{延拓缺陷衰减}
\Longrightarrow
\text{条件性的有界校正次数},
$$

$$
\text{强正则性}
+\text{远端影响衰减}
\Longrightarrow
\text{条件性的固定前缀收敛},
$$

$$
\text{无限维二阶充分条件}
\Longrightarrow
\text{严格局部最优}.
$$

任何一行都不能跳过。特别地，终端集合成员关系、小的预条件残差和到达有限
$N_{\max}$ 都不是全局无限时域最优性的证明。

## 11. 第二阶段代码字段和退出状态

结果中与第二阶段直接相关的字段包括：

- `certified_horizon`：$N_{\mathrm{cert}}$；
- `certificate_horizon`：为兼容旧结果保留，语义等于
  `certified_horizon`；
- `final_preconditioned_stationarity_residual`：最终 $\chi_N$；
- `final_forcing_tolerance`：最终 $\tau_N$；
- `final_adjoint_gradient_max`：最终原始伴随梯度诊断；
- `final_terminal_set_membership`：最终有限前缀是否仍满足终端成员关系。

主要退出状态包括：

- `trust_region_failure`：安全 corrector 未找到可接受步，但已有安全策略仍保留；
- `corrector_iteration_limit`：达到当前时域校正上限，不表示收敛或不可行；
- `maximum_horizon_reached`：达到有限计算上限，不表示无限时域收敛。

## 12. 历史代表性实验中的第二阶段（需重新运行无输入范围版本）

以下记录来自旧配置，不能作为本次无输入范围版本的数值结论。该历史设置采用 $\delta=0.1$。增长因子 $1.25$ 时，第一阶段在

$$
N_{\mathrm{hit}}=489
$$

获得安全前缀。该点的安全 corrector 同时通过残差门，因此

$$
N_{\mathrm{cert}}=489.
$$

随后认证 continuation 达到有限计算上限

$$
N_{\max}=900.
$$

最终诊断为

$$
\chi_{900}=1.402\times10^{-3},
$$

$$
\tau_{900}=3.331\times10^{-2}.
$$

这些数值只说明有限实现保持了终端证书并通过预条件残差门，不说明已经在数学
意义上达到无限时域极限，也不说明全局最优。

增长因子 $1.5$ 时，第一阶段在 $N_{\mathrm{hit}}=468$ 已获得安全尾部，但
安全 corrector 失败，所以 `certified_horizon` 为空。该失败不是系统不可行；
它说明当前局部 corrector 没有在给定预算和安全约束下建立
$N_{\mathrm{cert}}$。

## 13. 第二阶段最终逻辑链

$$
\boxed{
\begin{aligned}
&N_{\mathrm{hit}}<\infty
+\text{终端下降}
\Longrightarrow
\text{始终存在安全回退策略},\\
&x_N\in\mathcal X_f
+\|\nabla\widehat J_N\|\le\varepsilon_N
\Longrightarrow
\text{得到有限时域近似驻点检查点},\\
&\text{延拓缺陷衰减}
+\text{一致局部收缩}
\Longrightarrow
\text{有界 corrector 工作量},\\
&\text{强正则性}
+\text{远端影响衰减}
\Longrightarrow
\text{条件性的固定前缀收敛}.
\end{aligned}
}
$$

第一行是第一阶段证书在第二阶段中的不变安全底座；后面三行依次需要新的条件，
不能由第一行自动推出。
