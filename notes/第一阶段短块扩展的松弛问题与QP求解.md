# 第一阶段短块扩展的松弛问题与 QP 求解

本文可接在草稿的“假设2”之后，作为第一阶段短块控制的具体计算方法。符号沿用 $W_f(x)=\gamma x^\top Px$，以 $l$ 表示短块扩展次数，以 $j$ 表示固定块长下的内部求解迭代次数。

## 1. 短块控制的可行性子问题

在第 $l$ 次短块扩展中，固定当前轨迹终点
$$
z_0^l=\bar x_N,\qquad W_0^l:=W_f(z_0^l)>\alpha_f.
$$
对候选块长 $h\in\{1,\dots,h_{\max}\}$，记待求控制块为
$$
U_h^l=\operatorname{col}(w_0^l,\dots,w_{h-1}^l)\in\mathbb R^{mh}.
$$
给定 $U_h^l$ 后，由原始非线性动力学递推得到
$$
z_{t+1}^l=a(z_t^l)+B(z_t^l)w_t^l,\qquad t=0,\dots,h-1.
$$

为统一表示终点收缩条件和中间状态增长条件，定义
$$
b_t^l=
\begin{cases}
cW_0^l,&1\le t<h,\\
\bar\rho W_0^l,&t=h.
\end{cases}
$$
则一个合格控制块应满足
$$
W_f(z_t^l)-b_t^l\le0,\qquad t=1,\dots,h.
$$
当 $h=1$ 时，仅保留终点收缩条件，不存在中间状态增长条件。

引入公共松弛变量 $s\ge0$，构造可行性子问题
$$
\begin{aligned}
\min_{z_1^l,\dots,z_h^l,\;U_h^l,\;s}\quad&s\\
\mathrm{s.t.}\quad
&z_{t+1}^l=a(z_t^l)+B(z_t^l)w_t^l,
&&t=0,\dots,h-1,\\
&W_f(z_t^l)-b_t^l\le s,
&&t=1,\dots,h,\\
&s\ge0.
\end{aligned}
\tag{1}
$$
其中，$z_0^l$、$W_0^l$ 和 $b_t^l$ 在本次固定块长的求解过程中保持不变。

对任意控制块 $U_h^l$，定义其在真实动力学下的最大约束违反量
$$
\Phi_h^l(U_h^l)
:=
\max\left\{
0,\;
\max_{1\le t\le h}
\bigl[W_f(z_t^l(U_h^l))-b_t^l\bigr]
\right\}.
\tag{2}
$$
在固定控制块及其真实状态轨迹后，满足式（1）的最小松弛量恰为
$$
s=\Phi_h^l(U_h^l).
$$
因此，优化松弛量的实质是通过调整控制块，减小真实轨迹的最大约束违反量。第一阶段只需找到满足
$$
\Phi_h^l(U_h^l)=0
\tag{3}
$$
的控制块，而不要求先求得式（1）的全局最优解。

## 2. 动力学与约束的一阶近似

式（1）包含非线性动力学，通常不能直接作为二次规划求解。为此，在固定 $l$ 和 $h$ 后，采用逐次线性化方法构造局部子问题。

设第 $j$ 次内部迭代的名义控制块为
$$
\bar U_h^{l,j}
=
\operatorname{col}
(\bar w_0^{l,j},\dots,\bar w_{h-1}^{l,j}),
$$
并从固定起点 $z_0^l$ 出发，使用原始动力学计算名义状态
$$
\bar z_0^{l,j}=z_0^l,\qquad
\bar z_{t+1}^{l,j}
=
a(\bar z_t^{l,j})
+B(\bar z_t^{l,j})\bar w_t^{l,j}.
\tag{4}
$$

引入状态与控制增量
$$
\delta z_t\in\mathbb R^n,\qquad
\delta w_t\in\mathbb R^m.
$$
由于本次短块起点固定，有 $\delta z_0=0$。

记 $F(z,w):=a(z)+B(z)w$，并定义
$$
A_t^{l,j}
:=
a_x(\bar z_t^{l,j})
+
\sum_{q=1}^{m}
\bar w_{t,q}^{l,j}
D B_{\cdot q}(\bar z_t^{l,j}),
\qquad
B_t^{l,j}:=B(\bar z_t^{l,j}),
\tag{5}
$$
其中 $B_{\cdot q}(x)$ 表示 $B(x)$ 的第 $q$ 列，$D B_{\cdot q}(x)\in\mathbb R^{n\times n}$ 为该列向量关于状态的 Jacobian 矩阵。

在名义轨迹附近，一阶展开为
$$
\begin{aligned}
&F(\bar z_t^{l,j}+\delta z_t,\,
       \bar w_t^{l,j}+\delta w_t)\\
&\qquad\approx
F(\bar z_t^{l,j},\bar w_t^{l,j})
+A_t^{l,j}\delta z_t
+B_t^{l,j}\delta w_t.
\end{aligned}
$$
利用名义轨迹满足式（4），得到线性化动力学
$$
\boxed{
\delta z_0=0,\qquad
\delta z_{t+1}
=A_t^{l,j}\delta z_t+B_t^{l,j}\delta w_t,
\quad t=0,\dots,h-1.
}
\tag{6}
$$

另一方面，由 $W_f(z)=\gamma z^\top Pz$ 可得
$$
\nabla W_f(z)=2\gamma Pz,\qquad
\nabla^2W_f(z)=2\gamma P.
$$
因此
$$
\begin{aligned}
W_f(\bar z_t^{l,j}+\delta z_t)
={}&W_f(\bar z_t^{l,j})
+2\gamma(\bar z_t^{l,j})^\top P\delta z_t\\
&+\gamma\delta z_t^\top P\delta z_t.
\end{aligned}
\tag{7}
$$
舍去式（7）中的二阶增量项，定义一阶预测
$$
\widehat W_t^{l,j}(\delta z_t)
:=
W_f(\bar z_t^{l,j})
+2\gamma(\bar z_t^{l,j})^\top P\delta z_t.
\tag{8}
$$
于是式（1）中的不等式近似为
$$
\widehat W_t^{l,j}(\delta z_t)-b_t^l\le s,
\qquad t=1,\dots,h.
\tag{9}
$$

式（6）和式（9）均为关于状态增量、控制增量及松弛变量的线性约束。

## 3. 消去状态增量并构造二次规划

将控制增量堆叠为
$$
p:=\operatorname{col}(\delta w_0,\dots,\delta w_{h-1})
\in\mathbb R^{mh}.
$$
定义选择矩阵
$$
E_t\in\mathbb R^{m\times mh},\qquad E_tp=\delta w_t.
$$
进一步定义状态增量对控制增量的敏感度矩阵
$$
M_t^{l,j}\in\mathbb R^{n\times mh}.
$$
由式（6）递推可得
$$
M_0^{l,j}=0,\qquad
M_{t+1}^{l,j}
=A_t^{l,j}M_t^{l,j}+B_t^{l,j}E_t,
\tag{10}
$$
从而
$$
\delta z_t=M_t^{l,j}p.
\tag{11}
$$

将式（11）代入式（9），记
$$
v_t^{l,j}
:=W_f(\bar z_t^{l,j})-b_t^l,
\qquad
g_t^{l,j}
:=2\gamma(M_t^{l,j})^\top P\bar z_t^{l,j}
\in\mathbb R^{mh},
\tag{12}
$$
得到
$$
v_t^{l,j}+(g_t^{l,j})^\top p\le s,
\qquad t=1,\dots,h.
\tag{13}
$$

令
$$
v^{l,j}
=
\operatorname{col}(v_1^{l,j},\dots,v_h^{l,j}),
$$
$$
G^{l,j}
=
\begin{bmatrix}
(g_1^{l,j})^\top\\
\vdots\\
(g_h^{l,j})^\top
\end{bmatrix}
\in\mathbb R^{h\times mh}.
$$
则全部预测不等式可以写为
$$
G^{l,j}p+v^{l,j}\le s\mathbf1_h.
\tag{14}
$$

为限制一阶模型的使用范围，引入信赖域
$$
\|Dp\|_\infty\le\Delta_j,
\tag{15}
$$
其中 $D\in\mathbb R^{mh\times mh}$ 为正对角尺度矩阵，$\Delta_j>0$ 为当前信赖域半径。该条件约束本次控制增量，不改变原问题中 $w_t^l\in\mathbb R^m$ 的控制取值范围。

为保持可行性优先，采用字典序目标：先最小化预测松弛量，再在最小松弛水平内选择较小的控制增量。具体分为以下两步。

第一步求解
$$
\begin{aligned}
s_j^{\mathrm{pred}}
=\min_{p,s}\quad&s\\
\mathrm{s.t.}\quad
&G^{l,j}p+v^{l,j}\le s\mathbf1_h,\\
&s\ge0,\\
&-\Delta_j\mathbf1_{mh}
\le Dp\le
\Delta_j\mathbf1_{mh}.
\end{aligned}
\tag{16}
$$
这是线性规划，也可以视为 Hessian 为零的凸二次规划。由于 $p=0$ 配合足够大的 $s$ 总是可行，式（16）不会仅因当前轨迹违反收缩条件而不可行。

第二步固定最小预测松弛量，求解
$$
\boxed{
\begin{aligned}
\min_p\quad&\frac12p^\top H_jp\\
\mathrm{s.t.}\quad
&G^{l,j}p+v^{l,j}
\le s_j^{\mathrm{pred}}\mathbf1_h,\\
&-\Delta_j\mathbf1_{mh}
\le Dp\le
\Delta_j\mathbf1_{mh},
\end{aligned}
}
\tag{17}
$$
其中 $H_j\succ0$，例如可取 $H_j=I_{mh}$。式（17）是具有线性约束的严格凸二次规划。其作用是在保持最小预测违反量的前提下，选择较小的控制增量。

在浮点实现中，式（17）的松弛上界可改为
$$
s_j^{\mathrm{pred}}+\tau_s,\qquad \tau_s\ge0,
$$
以容纳两次求解之间的数值误差；后续预测违反量应根据实际所得增量重新计算。

上述两步分别确定“能够将预测违反量降低到什么程度”以及“采用哪个控制增量实现这一降低”。因此，$s$ 由优化器与控制增量共同确定，不需要另外指定松弛量的更新公式。

## 4. 真实动力学验证与内部迭代

记式（17）返回的控制增量为 $p_j$，对应分量为 $\delta w_t^j=E_tp_j$。构造试探控制
$$
w_t^{l,j,\mathrm{trial}}
=
\bar w_t^{l,j}+\delta w_t^j,
\tag{18}
$$
并从同一个起点重新进行真实动力学递推：
$$
z_0^{l,j,\mathrm{trial}}=z_0^l,
$$
$$
z_{t+1}^{l,j,\mathrm{trial}}
=
a(z_t^{l,j,\mathrm{trial}})
+B(z_t^{l,j,\mathrm{trial}})
w_t^{l,j,\mathrm{trial}}.
\tag{19}
$$

当前真实违反量为
$$
\phi_j:=\Phi_h^l(\bar U_h^{l,j}),
$$
试探控制的真实违反量为
$$
\phi_j^{\mathrm{trial}}
:=
\max\left\{
0,\;
\max_{1\le t\le h}
\left[
W_f(z_t^{l,j,\mathrm{trial}})-b_t^l
\right]
\right\}.
\tag{20}
$$
局部模型对试探控制的预测违反量为
$$
\widehat\phi_j(p_j)
:=
\max\left\{
0,\;
\max_{1\le t\le h}
\left[
v_t^{l,j}+(g_t^{l,j})^\top p_j
\right]
\right\}.
\tag{21}
$$
相应定义预测下降量和真实下降量
$$
\operatorname{pred}_j
:=\phi_j-\widehat\phi_j(p_j),
\qquad
\operatorname{ared}_j
:=\phi_j-\phi_j^{\mathrm{trial}}.
\tag{22}
$$
当 $\operatorname{pred}_j>0$ 时，定义下降比
$$
r_j:=
\frac{\operatorname{ared}_j}{\operatorname{pred}_j}.
\tag{23}
$$

若试探控制已经满足真实约束，即
$$
\phi_j^{\mathrm{trial}}=0,
\tag{24}
$$
则直接保存该控制块并结束本次可行性搜索。若初始名义控制块本身已经满足真实约束，也可直接保存，无须执行局部求解。

若真实约束尚未全部满足，则以违反量下降决定是否接受为内部迭代。给定 $0<\eta_{\mathrm{acc}}<1$，当真实递推有定义且
$$
\operatorname{pred}_j>0,\qquad
r_j\ge\eta_{\mathrm{acc}}
\tag{25}
$$
时，接受试探控制：
$$
\bar U_h^{l,j+1}=U_h^{l,j,\mathrm{trial}}.
$$
随后根据新的真实轨迹重新计算动力学 Jacobian、状态敏感度及约束线性化。

若式（25）不满足，则保留原名义控制块，并缩短信赖域，例如
$$
\Delta_{j+1}=\theta_{\mathrm{dec}}\Delta_j,
\qquad 0<\theta_{\mathrm{dec}}<1,
$$
重新构造并求解局部子问题。若预测与真实下降符合较好，则可以保持或适当增大信赖域。当前初值下没有取得有效进展时，可以更换初值或尝试下一个候选块长。

需要强调，式（16）得到
$$
s_j^{\mathrm{pred}}=0
$$
仅表示线性化模型预测可行，不能代替式（24）的真实验收。由式（7），即使暂不考虑动力学线性化误差，舍去的项
$$
\gamma\delta z_t^\top P\delta z_t\ge0
$$
也可能使真实二次函数值超过其一阶预测。因此，控制块的最终接受必须以原始非线性动力学递推得到的状态为依据。有限精度计算还应预留收缩裕量；若要求严格认证，则需要用可验证的误差界覆盖数值递推和函数值计算误差。

## 5. 控制块提交与第一阶段终止

对 $h=1,\dots,h_{\max}$ 依次执行上述搜索。一旦找到满足
$$
W_f(z_h^l)\le\bar\rho W_f(z_0^l),
\qquad
W_f(z_t^l)\le cW_f(z_0^l),\quad1\le t<h
\tag{26}
$$
的真实控制块，即将该控制块及其对应状态拼接到当前轨迹末尾，并更新
$$
N\leftarrow N+h,\qquad
z_0^{l+1}=z_h^l.
\tag{27}
$$
若
$$
W_f(z_0^{l+1})\le\alpha_f,
$$
则第一阶段结束；否则继续进行下一次短块扩展。

假设2保证满足式（26）的短块存在，但不自动保证局部逐次优化方法能够找到该块。因此，若要进一步证明上述数值算法有限终止，还需补充求解器能够在每个实际访问端点有限找到并可靠验收合格控制块的条件。若预定块长和初值均未找到合格块，只能判定本次搜索未成功，不能据此判定系统不可达。
