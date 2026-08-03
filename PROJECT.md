---
title: 项目级记忆
type: project-memory
status: active
updated: 2026-07-31
tags:
  - project/memory
  - ptc
aliases:
  - PROJECT
  - 项目记忆
related:
  - "[[01-项目管理/研究方案|研究方案]]"
  - "[[02-理论基础/理论推导教程|理论推导教程]]"
---

# PROJECT.md

> 项目名称：手征光子时间晶体中的无序诱导时间拓扑相  
> 英文工作题目：Disorder-Induced Temporal Topology in a Chiral Photonic Time Crystal  
> 当前阶段：单层传播矩阵已通过验收，干净双层 PTC 的 $k$-gap 已完成复现  
> 最近更新：2026-07-31

## 1. 项目级结论

本项目研究：保持有效手征结构的时间无序，能否把一个干净、平庸的光子时间晶体（PTC）驱动为非平庸时间拓扑相。

目标证据链为：

$$
\{n_m,\tau_m\}
\longrightarrow
\text{时间界面递推}
\longrightarrow
H_{\mathrm{eff}}
\longrightarrow
\{v_n,w_n\}
\longrightarrow
\gamma,\nu
\longrightarrow
\text{中隙时间界面态与因果场响应}.
$$

项目的核心不是“无序产生了局域峰”，而是同时证明：

1. 干净体系为平庸相；
2. 特定无序使体拓扑指标改变；
3. 同一参数窗口出现干净体系中不存在的中隙时间界面态；
4. 完整 PTC 时间传播计算重现有效模型结果；
5. 统计、有限尺寸和必要对照排除普通 Anderson 局域等非拓扑解释。

研究定位是 simulation-first 的数值拓扑光子学项目。已有理论用于定义模型和验收结果，主要贡献应来自物理无序方案、相图、界面态、时间场动力学、统计和对照，而不是发展新的随机拓扑数学理论。

## 2. 当前项目状态

### 已完成

- 建立了干净手征 PTC、SSH 模型和二者映射的理论框架。
- 推导了任意逐层、手征保持无序下的精确 PTC–SSH 映射。
- 推导了无序拓扑判据

  $$
  \gamma
  =
  \left\langle
  \ln\frac{n_{2,n}}{n_{1,n}}
  \right\rangle.
  $$

- 推导了准周期外部层调制的解析临界幅度和均匀随机外部层调制的隐式临界条件。
- 固定了单位胞、拓扑标签、物理假设、最小证据链、必要对照和六个月路线。
- 确定采用“理论推导—同步仿真”工作方式：每完成一层文献复现或本课题推导，就立即完成对应的最小数值实现、结果图和验收测试。
- 已用 Python 实现单层时间传播矩阵，并通过行列式、前后互逆、均匀介质合并和独立 Maxwell ODE 对照；当前示例的最大归一化状态相对误差约为 $1.6\times10^{-11}$。
- 已用 Yang Figure 1(c) 参数 $n_1=4,n_2=2,\bar t=1\,\mathrm{fs}$ 和周期矩阵 $M=P_2P_1$ 复现干净手征 PTC 色散及周期性 $k$-gap。
- 已将 `C:\Users\Lenovo\Desktop\PTC` 重构为唯一 Obsidian Vault，核心文档、文献、后续仿真和结果均在同一工作目录维护。
- 已初始化以 `main` 为默认分支的 Git 仓库；Markdown、Bases、项目配置和后续代码进入版本控制，PDF 文献由 Git LFS 管理。

### 尚未完成

- 除单层传播矩阵和干净 PTC $k$-gap 基准外，其余仿真模块尚未建立。
- 尚未用代码复现干净 SSH 和 PTC–SSH 色散等价关系。
- 尚未验证有限无序链的相图、中隙态、统计收敛和完整 TMM/Maxwell 响应。
- 尚未评估具体实验平台中的调制幅度、有限开关、色散和损耗。
- 最终论文命名尚未确定；在真正随机无序得到完整证据前，不使用强表述 “temporal topological Anderson phase”。

## 3. 核心科学问题与边界

主问题：

> 从干净平庸手征 PTC 出发，保持手征结构的时间无序能否使其体拓扑从 $\nu=0$ 变为 $\nu=1$，并产生新的时间界面态？

需要回答的子问题：

1. 候选拓扑窗口在哪里，解析预测与有限链结果是否一致？
2. 拓扑改变时是否出现 $\lambda\approx0$ 的中隙态？
3. 中隙态的局域长度、子晶格极化、物理电位移和放大响应是什么？
4. 有效 SSH 结论能否由原始 PTC 时间传播矩阵重现？
5. 结果在真正随机无序、不同尺寸和不同样本中是否稳定？

项目暂不以以下内容为主线：

- 随机拓扑的严格数学证明；
- 完整 bulk–edge correspondence 证明；
- K 理论和完整非 Hermitian 分类；
- 多种拓扑指标并行开发；
- 多种无序分布的大规模比较；
- 同时处理损耗、色散、连续开关和多种噪声；
- 具体空间器件的全波仿真。

## 4. 固定物理假设

精确无序 PTC–SSH 映射依赖以下条件：

1. 每个时间层满足

   $$
   \frac{\tau_m}{n_m}=\bar t.
   $$

2. 介质在目标频段内可视为无色散，$n_m$ 为实数。
3. 磁导率固定为 $\mu_0$。
4. 时间界面足够快，可使用分段常数模型。
5. 介质空间均匀，因此理想时间界面保持真实空间波矢 $k$。
6. 基础模型不含额外损耗、增益、导电项或材料记忆。
7. 常规三点递推要求 $\sin(kc\bar t)\ne0$；$\lambda=\pm1$ 的特殊点应回到原始传播矩阵处理。

若后续加入有限开关、色散、损耗或固定 $\tau$ 的折射率扰动，必须重新检查手征结构，不能默认继续使用精确映射和解析相界。

## 5. 固定符号、单位胞与拓扑约定

- 第 $n$ 个单位胞包含时间界面站点 $A_n,B_n$。
- 内部时间层 $(n_{1,n},\tau_{1,n})$ 位于 $A_n$ 与 $B_n$ 之间。
- 外部时间层 $(n_{2,n},\tau_{2,n})$ 位于 $B_n$ 与 $A_{n+1}$ 之间。
- $v_n$ 是 $A_n\leftrightarrow B_n$ 的胞内耦合。
- $w_n$ 是 $B_n\leftrightarrow A_{n+1}$ 的胞间耦合。
- SSH 站点是时间界面处的电位移变量，不是时间层。
- 干净体系始终采用

  $$
  v>w\Rightarrow\nu=0,\qquad
  v<w\Rightarrow\nu=1.
  $$

- 在当前 PTC 映射中

  $$
  n_2>n_1\Rightarrow\nu=0,\qquad
  n_1>n_2\Rightarrow\nu=1.
  $$

- 不得在推导、代码或作图中途交换单位胞边界或链的终止方式。
- 干净周期体系中 $q=\Omega T$；有效本征值与真实波矢满足

  $$
  \lambda=\cos(kc\bar t).
  $$

- $\Omega$、$q$、$k$ 和 $\lambda$ 是不同物理量，不得混用。

## 6. 采用的模型

### 6.1 干净 SSH 基准

$$
H(q)=
\begin{pmatrix}
0&v+we^{-iq}\\
v+we^{iq}&0
\end{pmatrix},
\qquad
\sigma_zH(q)\sigma_z=-H(q).
$$

干净 PTC 映射为

$$
v=\frac{n_2}{n_1+n_2},
\qquad
w=\frac{n_1}{n_1+n_2}.
$$

该公式只用于干净周期体系和代码基准，不得逐单位胞套用于无序系统。

### 6.2 任意手征保持无序的精确映射

令

$$
r_m=\frac{1}{n_m},
\qquad
s_m=r_{m-1}+r_m.
$$

时间界面递推形成广义本征值问题

$$
K\mathbf u=\lambda S\mathbf u,
$$

其中 $K$ 为实对称最近邻矩阵，$S$ 为正定对角矩阵。标准 Hermitian 有效矩阵必须写成

$$
H_{\mathrm{eff}}=S^{-1/2}KS^{-1/2},
\qquad
\boldsymbol\psi=S^{1/2}\mathbf u.
$$

不得直接把通常不对称的 $S^{-1}K$ 当作 SSH Hamiltonian。

定义

$$
r_{1,n}=\frac{1}{n_{1,n}},
\qquad
r_{2,n}=\frac{1}{n_{2,n}},
$$

$$
s_{A,n}=r_{2,n-1}+r_{1,n},
\qquad
s_{B,n}=r_{1,n}+r_{2,n}.
$$

精确有效耦合为

$$
v_n=
\frac{r_{1,n}}{\sqrt{s_{A,n}s_{B,n}}},
\qquad
w_n=
\frac{r_{2,n}}{\sqrt{s_{B,n}s_{A,n+1}}}.
$$

有效波函数还原为物理电位移时使用

$$
D(t_m)=u_m=\frac{\psi_m}{\sqrt{s_m}}.
$$

### 6.3 无序拓扑判据

零模递推为

$$
a_{n+1}=-\frac{v_n}{w_n}a_n.
$$

有限链典型衰减指数为

$$
\gamma_N=
\frac{1}{N}\sum_{n=1}^{N}
\ln\left|\frac{v_n}{w_n}\right|.
$$

当前终止方式下

$$
\gamma<0\Rightarrow\nu=1,\qquad
\gamma>0\Rightarrow\nu=0,\qquad
\gamma=0\Rightarrow\text{相变}.
$$

对物理 PTC，局部归一化因子在长链中望远镜消去：

$$
\gamma=
\left\langle\ln n_{2,n}\right\rangle
-
\left\langle\ln n_{1,n}\right\rangle.
$$

因此相界由折射率的几何平均决定，而不是算术平均：

$$
n_{2,\mathrm{typ}}=n_{1,\mathrm{typ}}.
$$

局域长度的典型估计为

$$
\xi_{\mathrm{edge}}\simeq\frac{1}{|\gamma|}.
$$

临界附近 $\xi_{\mathrm{edge}}$ 发散，必须做更长链和有限尺寸分析。

## 7. 当前采用的无序方案

### 主方案 A：准周期外部层调制

固定内部层：

$$
n_{1,n}=n_{1,0},
\qquad
\tau_{1,n}=\bar t n_{1,0}.
$$

调制外部层：

$$
n_{2,n}
=n_{2,0}+\Delta\cos(2\pi\alpha n+\phi),
\qquad
\tau_{2,n}=\bar t n_{2,n},
$$

并要求 $0\le\Delta<n_{2,0}$。干净平庸条件为 $n_{2,0}>n_{1,0}$。

解析临界幅度为

$$
\Delta_c=
2\sqrt{n_{1,0}(n_{2,0}-n_{1,0})}.
$$

正折射率范围内存在相变要求

$$
\frac{n_{2,0}}{2}<n_{1,0}<n_{2,0}.
$$

这一路线应准确称为 “quasiperiodicity-induced topology”。

### 主方案 B：真正随机外部层调制

$$
n_{2,n}=n_{2,0}(1+\rho\xi_n),
\qquad
\xi_n\sim U[-1,1],
\qquad
0\le\rho<1,
$$

并同步设置 $\tau_{2,n}=\bar t n_{2,n}$。相界由

$$
\ln\frac{n_{2,0}}{n_{1,0}}+L(\rho_c)=0
$$

隐式给出，其中

$$
L(\rho)=
\frac{
(1+\rho)\ln(1+\rho)
-
(1-\rho)\ln(1-\rho)
-
2\rho
}{2\rho}.
$$

正折射率均匀无序能诱导相变的必要条件为

$$
\frac{n_{2,0}}{n_{1,0}}<\frac{e}{2}.
$$

若相变只在 $\rho\rightarrow1$ 时出现，应把干净参数移近相界，而不是使用接近零折射率的可疑层。

### 理论验证参数

首个无量纲基准采用

$$
n_{1,0}=2.0,\qquad n_{2,0}=2.1.
$$

对应

$$
v_0\approx0.51220,\qquad
w_0\approx0.48780,
$$

为靠近相界的干净平庸体系。解析预测：

- 准周期调制 $\Delta_c\approx0.89443$；
- 均匀随机调制 $\rho_c\approx0.518$；
- $\rho=0.70$ 时 $\gamma\approx-0.0489$，$\xi_{\mathrm{edge}}\approx20.5$ 个单位胞；
- 首轮随机链建议从 $N=100$ 开始，并比较 $N=50,200,400$。

这些数值仅用于理论和代码验证，不等于已选定实验材料参数。

### 不作为物理主线的方案

- 抽象 SSH 中直接设 $w_n=w_0+V\cos(2\pi\alpha n+\phi)$：从 $v_0>w_0$ 出发通常需 $V>2v_0$，并伴随近零键和符号变化，不作为首选物理 PTC 方案。
- 两类时间层使用相同相对无序：其对数平均修正相消，主要作为已有拓扑态鲁棒性对照，通常不能诱导体相改变。
- 只随机化 $n_m$ 而保持 $\tau_m$ 不变：会破坏 $\tau_m/n_m=\bar t$，仅作为手征破坏对照。

## 8. 数值工作流

### 阶段 1：干净 Maxwell/PTC 基线

- 已从 $(D,B)^\mathsf T$ 一阶系统构造单层传播矩阵。
- 已验证单层矩阵行列式为 1、前后传播互逆、均匀介质合并，并与独立 Maxwell ODE 积分一致。
- 已构造双层周期矩阵 $M=P_2P_1$ 并得到 Floquet 色散。
- 已用 Yang Figure 1(c) 参数复现干净 $k$-gap。

### 阶段 2：干净 PTC–SSH 复现

- 复现干净 SSH 能带、winding number 和开边界零模。
- 验证 $v>w$ 为平庸、$v<w$ 为拓扑。
- 验证 SSH 色散与 PTC/TMM Floquet 色散一致。
- 固化单位胞、边界终止、相位和归一化约定。

### 阶段 3：任意无序精确映射

- 从折射率序列构造 $K,S$ 和 $H_{\mathrm{eff}}$。
- 验证 $H_{\mathrm{eff}}$ 实对称、零对角、最近邻和手征对称。
- 用短随机链比较广义本征值与 $H_{\mathrm{eff}}$ 本征值。
- 验证 $\psi_m=\sqrt{s_m}D_m$ 的场还原。

### 阶段 4：相图和有限链拓扑

- 先扫描 $\gamma_N$ 或 $\nu_N$，再计算有限链谱。
- 准周期扫描包含 $\Delta,\phi,n_{1,0}/n_{2,0}$。
- 随机扫描包含 $\rho$、样本数和系统尺寸。
- 用扭曲边界 winding number 交叉检查 $\gamma_N$。
- 选择干净平庸、无序拓扑和强无序三个代表点。

### 阶段 5：物理 PTC 与完整证据闭环

- 构造平庸区与无序诱导拓扑区的时间拼接。
- 寻找 $\lambda\approx0$ 的谱界面态。
- 用 $D_m=\psi_m/\sqrt{s_m}$ 与原始 TMM/Maxwell 场形比较。
- 单独计算因果输入场的未来演化、局域、增长、衰减和放大。
- 完成随机统计、有限尺寸、必要对照和一种现实因素测试。

原则：理论和仿真逐层同步推进；每层先得到最小可视结果并通过数值验收，再进入下一层。正式结论阶段先确认体拓扑和中隙态，再优化局域峰或放大图。

## 9. 计划代码结构

当前采用 Python。`propagation_matrix.py` 是可导入的单层传播与验证模块，`propagation_matrix` 是对应的命令行入口；`clean_ptc_k_gap.py` 负责干净双层 PTC 色散、$k$-gap 和 $\cos(\Omega T)=\mathrm{Tr}(M)/2$ 判据图。后续模块仍按下面的职责划分逐步建立。

```text
04-仿真/
  propagation_matrix
  propagation_matrix.py
  propagation_matrix_validation.png
  clean_ptc_k_gap.py
  clean_ptc_k_gap.png
  effective_ssh/
    clean_ssh
    disordered_mapping
    topological_invariant
    quasiperiodic_scan
    random_disorder_scan
    edge_spectrum
  ptc_transfer_matrix/
    layer_propagator
    clean_chiral_ptc
    ptc_ssh_mapping
    disordered_ptc
    temporal_field_dynamics
  controls/
    clean_trivial_interface
    equal_relative_disorder
    internal_layer_disorder
    ordinary_anderson
    clean_topological_robustness
    chiral_breaking
  data/
  tests/
05-结果/
  figures/
06-实验记录/
  每日记录/
```

推荐从 Python 的可复用模块加少量 notebook 开始：模块负责矩阵构造、扫描、统计和测试，notebook 只负责探索与作图。不要把核心算法只保存在 notebook 单元格中。

## 10. 预期运行流程

当前可执行入口：

```powershell
python "04-仿真\propagation_matrix"
python "04-仿真\clean_ptc_k_gap.py"
```

第一个入口运行单层传播矩阵基准并保留数值验收；第二个入口默认使用 Yang Figure 1(c) 参数，仅复现并保存干净 PTC 色散、$k$-gap 与 $\cos(\Omega T)=\mathrm{Tr}(M)/2$ 判据图。后续代码应形成固定顺序：

1. 运行干净 SSH 和干净 PTC 基准测试；
2. 运行短链精确映射一致性测试；
3. 运行准周期粗扫并定位候选窗口；
4. 在代表点计算谱、态、IPR 和子晶格极化；
5. 运行真正随机粗扫；
6. 在临界区增加样本数和尺寸；
7. 运行完整 PTC/TMM 对照；
8. 运行因果场演化和必要对照；
9. 从保存的数据统一生成论文图。

后续仍需补充依赖版本、统一配置格式和各阶段正式入口。

## 11. 数据与可复现性规范

每次扫描必须保存：

- 完整物理参数与无量纲参数；
- 单位胞和边界终止约定；
- 系统尺寸；
- 准周期 $\alpha,\phi$ 或随机种子；
- 实际使用的折射率、持续时间和有效耦合序列；
- $\gamma_N$、单样本整数 $\nu_N$ 和拓扑样本比例；
- 本征值、中隙态位置、IPR 和子晶格极化；
- 物理场局域长度、峰值放大和其他可观测量；
- 样本数、均值、标准误或置信区间；
- 代码版本、数据格式版本和生成时间。

随机扫描建议：

- 粗扫至少 100 个样本；
- 临界附近和主图参数使用 500 至 1000 个样本；
- 比较 $N=50,100,200,400$；
- 不对整个二维参数空间盲目使用最高样本数。

长链传输矩阵出现溢出或病态时，改用 QR、SVD、散射矩阵递推或逐步归一化；不要继续直接相乘大量矩阵。

## 12. 最低测试要求

### 代数与结构

- 单层传播矩阵满足 $\det P_j=1$。
- $H_{\mathrm{eff}}$ 为实对称、零对角、最近邻矩阵。
- $\|\Gamma H_{\mathrm{eff}}\Gamma+H_{\mathrm{eff}}\|$ 接近数值零。
- 谱满足 $\lambda\leftrightarrow-\lambda$ 成对。
- 短链广义本征值与标准化后的 $H_{\mathrm{eff}}$ 本征值一致。

### 干净极限

- 干净 SSH 数值能谱与解析解一致。
- $v>w$ 得到 $\nu=0$，$v<w$ 得到 $\nu=1$。
- 开边界零模的位置与终止方式一致。
- 无序关闭时精确回到干净结果。
- PTC 干净色散与手征 PTC 文献结果一致。
- SSH 色散与 TMM Floquet 色散一致。

### 无序与物理场

- $\gamma_N$ 与折射率对数平均一致。
- $\gamma_N$ 与扭曲边界 winding number 的标签一致。
- 周期环中的望远镜边界项严格消失。
- $D_m=\psi_m/\sqrt{s_m}$ 的场还原与原始方程一致。
- 时间界面处 $D,B$ 连续。
- 短链传输矩阵与直接场演化一致。
- 不同数值精度和时间步长下结果收敛。

## 13. 最小但充分的验收标准

只有同时满足以下条件，才支持“无序诱导时间拓扑”：

1. 干净基线为 $\nu=0$，无目标中隙态。
2. 有限无序窗口内单样本 $\nu_N=1$，拓扑样本比例随尺寸收敛。
3. $\gamma_N$ 与扭曲边界 winding number 一致。
4. 同一窗口出现新的 $\lambda\approx0$ 中隙态。
5. 中隙态具有强子晶格极化，并在界面两侧局域。
6. 中隙态的局域长度趋势与 $1/|\gamma|$ 一致。
7. 完整 PTC/TMM 在对应参数窗口重现本征值、场形和时间界面响应。
8. 随机样本、有限尺寸和误差统计支持结论。
9. 必要对照不能产生相同的“拓扑指标改变 + 中隙态 + PTC 响应”组合。

IPR 增大、局域峰、放大、接近中隙的孤立态或弱无序下已有边界态继续存在，均不能单独证明无序诱导拓扑。

## 14. 必要对照

1. **干净平庸–平庸界面**：应无目标中隙态。
2. **等相对无序**：令两类层使用相同相对分布，理论预测体相标签不变。
3. **内部层无序**：只随机化 $n_1$，理论预测体系更平庸。
4. **普通时间 Anderson 局域**：可以有较大 IPR 或局域峰，但 $\nu$ 不改变。
5. **手征破坏**：固定 $\tau$ 而随机化 $n$，检查谱钉扎和拓扑指标退化。
6. **干净拓扑态鲁棒性**：从 $n_1>n_2$ 出发加等相对无序，区分“保护已有态”和“创造新态”。
7. **准周期与真正随机**：决定最终是否可以使用 “topological Anderson” 命名。

## 15. 风险与止损条件

出现以下任一情况时应调整方案或降低主张：

- 准周期模型有拓扑窗口，但真正随机无序没有：改称 quasiperiodicity-induced topology。
- 有效模型有结果，但精确 PTC/TMM 无法重现：检查手征条件、归一化和边界；仍失败则降低为有效时间晶格结果或更换平台。
- 只有局域而拓扑指标不改变：按普通时间 Anderson 局域或无序放大处理。
- 结果只存在于单个尺寸、单个相位或少数样本。
- 拓扑指标改变与中隙态区间明显不一致。
- 只有破坏手征结构的极端参数才能产生目标态。
- 所需折射率或调制幅度超出合理物理范围。
- 理论推导连续超过两周仍未产生可执行测试或结果图：暂停扩展推导，返回数值复现。

## 16. 文档优先级与长期维护

当前 Vault 核心文档：

1. `PROJECT.md`：项目级长期记忆和当前采用方案；后续任务优先查阅。
2. `02-理论基础/理论推导教程.md`：2026-07-30 理论推导终稿，是公式、物理约定和验收细节的权威来源。
3. `01-项目管理/研究方案.md`：2026-07-24 simulation-first 项目方案，是研究定位、里程碑、证据组织和风险管理的主要来源。
4. `03-文献/PDF/`：本地参考文献原文。
5. `03-文献/文献笔记/`：结构化文献笔记，与 PDF 一一关联。
6. `00-项目主页.md`：Obsidian 中的统一工作入口。

`C:\Users\Lenovo\Desktop\PTC` 本身就是唯一工作 Vault。文档、后续代码、数据索引和结果均以该目录为权威来源，不再维护嵌套 Vault 副本。

当项目方案与理论教程冲突时，以较新的理论教程为物理模型基准。特别是：

- 不再把干净耦合公式逐胞套用于无序 PTC；
- 不再把准周期胞间键调制作为物理 PTC 首选；
- 主线使用外部时间层折射率调制，并同步调节持续时间以保持 $\tau_m/n_m=\bar t$；
- 相同相对无序被归类为鲁棒性对照，而不是诱导机制。

后续每次出现以下长期变化时更新本文件：

- 模型、单位胞、物理假设或主无序方案改变；
- 确定编程语言、目录、配置格式或运行入口；
- 解析预测被数值结果确认或否定；
- 验收标准、主图结构或论文主张改变；
- 新方案取代旧方案；
- 确定实验平台和现实参数。

### Git 版本控制约定

- 默认分支为 `main`。
- PDF 文献使用 Git LFS，避免普通 Git 历史被大二进制文件占满。
- 跟踪项目文档、Bases、模板、可复现代码和必要的 Obsidian 配置。
- 不跟踪 `.obsidian/workspace.json`、下载的主题和插件、缓存、虚拟环境、临时文件及大规模生成数据。
- 稳定的小型结果图可以进入 Git；原始扫描数据应放在已忽略的数据目录，并通过实验记录保存其生成参数和位置。
- 提交应对应可解释的项目状态；不要把一次大范围扫描产生的临时文件与代码修改混在同一提交中。

不要在本文件记录临时调试过程、一次性扫描结果或已废弃的细节。此类内容应放在实验记录、issue、notebook 或数据元数据中。

## 17. 当前开放问题

- 具体实验平台及可实现的 $n,\tau,\bar t$ 范围是什么？
- 有限上升时间、色散和损耗下，手征结构能保持到什么程度？
- 随机外部层方案在完整 PTC/TMM 中是否产生与解析相界一致的中隙态？
- 强无序下是否存在重入平庸区、近零时间层或额外非拓扑局域？
- 因果输入场能否以足够信噪比激发候选时间界面态？
- 最终应使用 “disorder-induced temporal topology”、 “quasiperiodicity-induced temporal topology” 还是 “temporal topological Anderson phase”？
- 项目代码采用 Python、MATLAB 还是 Mathematica；当前推荐 Python，但尚未形成正式决策。
- 项目方案引用的 Gao 等准周期光子 SSH 原文尚未加入本地 `03-文献/PDF/`，正式写作前必须补齐并复核。
