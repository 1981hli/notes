> 2022.10.03

# Lagrangian mechanics

- stationary-action 原理 / least action 原理

- Lagrange 力学描述的是 $(M,L)$ 系统 / (构象空间,Lagrange 函数)系统; $L=T-V$

- 构象空间与相空间有何异同?

- stationary-action 原理要求作用量泛函在 stationary 点出不变. 这个约束使得可以由 Lagrange 方程导出运动方程.

## 1 Introduction

- 通过选择适当的广义坐标, 可以避免考虑约束力

- 多粒子系统的广义坐标 $\mathbf{r_1}=(x_1,y_1,z_1)$, $\mathbf{r_2}=(x_2,y_2,z_2)$, ...; 广义速度 $\mathbf{v_1}=\frac{d\mathbf{r_1}}{dt}$, $\mathbf{v_2}=\frac{d\mathbf{r_2}}{dt}$, ...; 满足 Newton 定律
  $$\sum\mathbf{F}=m\frac{d^2 \mathbf{r}}{dt^2}$$

### 1.1 Lagrangian

- Lagrange 力学使用能量而不是力.

- Lagrange 力学的中心是 Lagrange 量; 它给出系统的全部动力学; Lagrange 量具有能量的量纲; 只要给出正确的运动方程, Lagrange 量可以是任何函数

- 粒子系统的非相对论 Lagrangian:
  $$L=T-V=\frac{1}{2}\sum_{k=1}^{N}m_k v_k^2$$
  即粒子系统的总能量

- 系统的动能只是粒子速度的函数, 与粒子的位置, 时间无关: $T=T(\mathbf{v_1},\mathbf{v_2},\cdots)$

- 系统的势能反映粒子之间的相互作用的能量.
  保守力(如 Newton 引力)下, 势能只与粒子的位矢有关: $v=v(\mathbf{r_1},\mathbf{r_2}, \cdots)$
  非保守力时, 势能由势能函数给出(如电磁势能), 此时与速度有关 $v=v(\mathbf{r _1},\mathbf{r_2} ,\cdots,\mathbf{v_1} ,\mathbf{v_2} ,\cdots)$;
  外场与时间有关时 $v=v(\mathbf{r _1},\mathbf{r_2} ,\cdots,\mathbf{v_1} ,\mathbf{v_2} ,\cdots,t)$

- 以上形式的 Lagrange 不适用于相对论的情形; 耗散系统的情况下还要引入额外的项.

- 完整约束(holonomic constraint) 只与粒子的位置有关: $f_1(\mathbf{r},t)=0,f_2(\mathbf{r},t)=0,\cdots,f_C(\mathbf{r},t)=0$

- 完整约束的方程给出粒子的可能的路径, 但未给出某时刻粒子的位置和速度.

- 非完整约束取决于粒子的位矢的各级时间导数 $\dot{\mathbf{r}},\ddot{\mathbf{r}},\cdots$

- Lagrange 力学只适用于完整约束的系统.

- 3 类非完整约束: 约束方程不可积; 约束有不等式; 有复杂的非保守力(如摩擦力)

- 含有非完整约束的系统有特别的处理方式, 比如回到 Newton 力学.

- 受到含时约束或者含时外力的影响, Langrange 显含时间;
  如果动能, 势能都不显含时间, Lagrange不显含时间;
  所有情况下, Lagrange 都隐含时间.

- 第一型Lagrange方程:
  $$\frac{\partial L}{\partial \mathbf{r}_k}-\frac{d}{dt}\frac{\partial L}{\partial \mathbf{\dot{r}}_k}+\sum_{i=1}^C\lambda_i\frac{\partial f_i}{\partial r_k}=0;$$
  $$\tag{第一型Lagrange方程}\;$$
  每个约束都贡献一个 Lagrange 乘子 $\lambda_i$;
  其中
  $$\frac{\partial}{\partial\mathbf{r}_k}=\left(\frac{\partial}{\partial x_k},\frac{\partial}{\partial y_k},\frac{\partial}{\partial z_k}\right),\quad\frac{\partial}{\partial\dot{\mathbf{r}_k}}=\left(\frac{\partial}{\partial \dot{x}_k},\frac{\partial}{\partial \dot{y}_k},\frac{\partial}{\partial \dot{z}_k}\right)$$

- Lagrange方程含有 $3N+C$ 个方程, 比Newton理论多出来$C$个约束方程;
  Lagrangian中的坐标和速度是独立变量;
  每一个约束都抵扣一个坐标, 故独立的坐标数为 $n=3N-C$.

- 以广义坐标表达位置坐标;
  广义坐标 n-tuple $\mathbf{q}=(q_1,q_2,\cdots,q_n)$;
  $\mathbf{r_k}=\mathbf{r_k}(\mathbf{q},t)=\left( x_k(\mathbf{q},t), y_k(\mathbf{q},t), z_k(\mathbf{q},t) \right)$;
  $\mathbf{q}$ 是构象空间中的一点.

- 广义速度 $\dot{q}_j=\frac{dq_j}{dt}$;
  以广义速度表达空间速度
  $$\mathbf{v_k}=\sum_{j=1}^n\frac{\partial\mathbf{r_k}}{\partial q_j}\dot{q}_j+\frac{\partial \mathbf{r_k}}{\partial t},$$
  其中 $k$ 是粒子编号, $j$ 是广义坐标编号.
  **其实我不太清楚这个公式怎么得到的?**

- 系统的能量决定于系统的广义坐标, 广义速度, 时间: $T=T(\mathbf{q},\dot{\mathbf{q}},t)$

- 第二型 Euler-Lagrange 方程 $$\frac{d}{dt}\left( \frac{\partial L}{\partial \dot{q}_j} \right)=\frac{\partial L}{\partial q_j}$$
  $$\tag{第二型Lagrange方程}\;$$

- 两型 Lagrange 方程的区别是, 第一型以空间坐标表达, 第二型以广义坐标表达.

## 2 从 Newton 力学到 Lagrange 力学

### 2.1 Newton 定律

- 单粒子的 Newton 第二定律 $\mathbf{F}=m\mathbf{a}$; 在3维空间中, 它是3个耦合的2阶常微分方程; Newton 定律在笛卡尔坐标中形式较简单.

- 曲线坐标系 $\xi=\left(\xi^1,\xi^2,\xi^3\right)$ 下, 以张量形式表示的Newton定律
  $$\tag{2.1-1}
  \begin{aligned}
  F^a&=m\left ( \frac{d^2\xi ^a}{dt^2} +\Gamma^a_{}{}_{bc}\frac{d\xi^b}{dt} \frac{d\xi^c}{dt}  \right )\\&=g^{ak}\left ( \frac{d}{dt}\frac{\partial T}{\partial\dot\xi^k}-\frac{\partial T}{\partial\xi^k}\right )
  \end{aligned}$$
  $$\tag{\color{red}这个公式是怎么来的呢?}\;$$

- 其中动能为 $$ T=\frac{1}{2}mg_{bc}\frac{\mathrm{d}\xi^b}{\mathrm{d}x}\frac{\mathrm{d}\xi^c}{\mathrm{d}x} $$

- 曲线坐标与广义坐标不是一回事.

- 合外力为 0 时, 上面方程的解是测地线.

- 合外力=合约束力+合非约束力 $$\mathbf{F}=\mathbf{C}+\mathbf{N}$$

- 如果有约束, 曲线坐标系不是独立的, 它们受到约束方程的限制.

### 2.2 D'Alembert 原理

- D'Alembert 原理: 虚位移的虚功为零.  $$ \sum_{k=1}^{N}\left(\mathbf{N}_k+\mathbf{C}_k-m_k\mathbf{a}_k\right)\cdot\delta\mathbf{r}_k=0\tag{2.2-1} $$
  上式可以解释为 $$非约束力的虚功+约束力的虚功+??=0$$

- 其中. 虚位移与约束力的方向垂直, 它的虚功为零
  $$\sum_{k=1}^{N}\mathbf{C}_k\cdot\delta\mathbf{r}_k=0$$
  考虑到如上这点, D'Alembert 原理简化为 $$ \sum_{k=1}^{N}\left(\mathbf{N}_k-m_k\mathbf{a}_k\right)\cdot\delta\mathbf{r}_k=0 \tag{2.2-3} $$

- D'Alembert 原理使得我们研究运动方程的时候只需要考虑非约束力. 

### 2.3 由 D'Alembert 原理给出运动方程

- 思路就是分别把D'Alembert原理式(2.2-3)左边的各项表示成广义坐标的函数.

- 虚位移 $\delta \mathbf{r}_k=\left(\delta x_k,\delta y_k,\delta z_k\right)$ 可以表达成广义坐标的函数, 于是
  $$ \delta \mathbf{r}_k=\sum_{j=1}^{n}\frac{\partial \mathbf{r}_k}{\partial q_j}\delta q_j $$
  此式形式上类似取全微分 $\mathrm{d}\mathbf{r}_k$. 
  根据虚位移的意义, 它是某一时刻的可能位移, 它不会是时间的函数, 因此 $\delta \mathbf{r}_k$ 不含有对时间的变分.

- (2.2-3)中非约束力的虚功可以写成
  $$\begin{aligned} \sum_{k=1}^{N}\mathbf{N}_k\cdot\delta\mathbf{r}_k&=\sum_{k=1}^{N}\mathbf{N}_k\cdot \left(\sum_{j=1}^{n}\frac{\partial \mathbf{r}_k}{\partial q_j}\delta q_j \right)\\ &=\sum_{j=1}^{n}\left(\sum_{k=1}^{N}\mathbf{N}_k\cdot\frac{\partial \mathbf{r}_k}{\partial q_j}\right)\delta q_j\\ &=\sum_{j=1}^{n}Q_j\delta q_j \end{aligned} \tag{2.3-2}$$
  其中广义力 $$Q_j=\sum_{k=1}^{N}\mathbf{N}_k\cdot \frac{\partial \mathbf{r}_k}{\partial q_j}$$

- D'Alembert 原理的最后一项我不知道怎么往下推导的 **(而这一步是得到后面的广义运动方程的关键! 不搞明白是不行的!)** 
  $$\begin{aligned} \sum_{k=1}^N{m_k}\mathbf{a}_k\cdot \delta \mathbf{r}_k&=\sum_{k=1}^N{m_k}\mathbf{a}_k\cdot \left( \sum_{j=1}^n{\frac{\partial \mathbf{r}_k}{\partial q_j}}\delta q_j \right)\\ &=\sum_{j=1}^n{\left( \sum_{k=1}^N{m_k}\mathbf{a}_k\cdot \frac{\partial \mathbf{r}_k}{\partial q_j} \right)}\delta q_j\\ &\xlongequal{?}\sum_{j=1}^n{\left( \frac{\mathrm{d}}{\mathrm{d}t}\frac{\partial T}{\partial \dot{q}^j}-\frac{\partial T}{\partial q_j} \right)}\delta q_j\\ \end{aligned} \tag{2.3-3}$$

- 把 (2.3-2), (2.3-3) 代入 (2.2-3) 中, 得
  $$
  \sum_{j=1}^n{\left( Q_j-\left( \frac{\mathrm{d}}{\mathrm{d}t}\frac{\partial T}{\partial \dot{q}^j}-\frac{\partial T}{\partial q_j} \right) \right)}\delta q_j=0
  \tag{2.3-4}
  $$

  各虚位移 $\delta q_j$ 独立且非零, 故上式要求
  $$
  \tag{2.3-5}
  Q_j=\frac{\mathrm{d}}{\mathrm{d}t}\frac{\partial T}{\partial \dot{q}^j}-\frac{\partial T}{\partial q_j}
  $$
  $$\tag{Lagrange eq.}\;$$
  此式即广义运动方程或 Lagrange 方程.

- (2.3-5) 与非约束力的 Newton 定律等价. 
  D'Alembert 原理中没有约束力.
  (2.3-5) 中的广义力有非约束力导出.
  广义力可以是非保守力.

### 2.4 Euler-Lagrange 方程和 Hamilton 原理

- 对于依赖于速度的非保守力, 可以找到含位置和速度的势能函数$V(\mathbf{r},\mathbf{v})$

- 如果广义力可以由势得到, 即
  $$
  \tag{2.4-1}
  Q_j=\frac{\mathrm{d}}{\mathrm{d}t}\frac{\partial V}{\partial \dot{q}^j}-\frac{\partial V}{\partial q_j}
  $$
  将其代入(2.3-5)中, 即可得到第二型Lagrange方程/Euler-Lagrange方程 
  $$
  \tag{2.4-2}
  \frac{\partial L}{\partial q_j}-\frac{\mathrm{d}}{\mathrm{d}t}\frac{\partial L}{\partial \dot{q}_j}=0
  $$
  $$\tag{Euler-Lagrange eq.}\;$$
  其中的 Lagrangian 为 $L=T-V$

- 对于非保守力的情形, 只有找到势的条件下, Euler-Lagrange 方程才适用. 而 Lagrange 方程则没有这个限制, 因为它含由广义力. Lagrange 方程的适用范围比 Euler-Lagrange 方程的适用范围更广.

- 也可以从变分法得到 Euler-Lagrange 方程.
  $$ \tag{2.4-3} \delta L=\sum_{j=1}^n{\left( \frac{\partial L}{\partial q_j}\delta q_j+\frac{\partial L}{\partial \dot{q}_j}\delta \dot{q}_j \right)} $$
  其中 $$ \delta \dot{q}_j=\delta \frac{\mathrm{d}q_j}{\mathrm{d}t}=\frac{\mathrm{d}\left( \delta q_j \right)}{\mathrm{d}t} $$

- $\delta L$ 对时间的积分
  $$
  \tag{2.4-5}
  \begin{aligned}
    \int_{t_1}^{t_2}{\delta L\mathrm{d}t}&=\int_{t_1}^{t_2}{\sum_{j=1}^n{\left( \frac{\partial L}{\partial q_j}\delta q_j+\frac{\partial L}{\partial \dot{q}_j}\delta \dot{q}_j \right)}\mathrm{d}t}\\
    &=\int_{t_1}^{t_2}{\sum_{j=1}^n{\left( \frac{\partial L}{\partial q_j}\delta q_j+\frac{\partial L}{\partial \dot{q}_j}\frac{\mathrm{d}\delta q_j}{\mathrm{d}t} \right)}\mathrm{d}t}\\
    &=\int_{t_1}^{t_2}{\sum_{j=1}^n{\left( \frac{\partial L}{\partial q_j}\delta q_j+\frac{\mathrm{d}}{\mathrm{d}t}\left( \frac{\partial L}{\partial \dot{q}_j}\delta q_j \right) -\left( \frac{\mathrm{d}}{\mathrm{d}t}\frac{\partial L}{\partial \dot{q}_j} \right) \delta q_j \right)}\mathrm{d}t}\\
    &=\sum_{j=1}^n{\left( \frac{\partial L}{\partial \dot{q}_j}\delta q_j \right) _{t_1}^{t_2}}+\int_{t_1}^{t_2}{\sum_{j=1}^n{\left( \frac{\partial L}{\partial q_j}-\left( \frac{\mathrm{d}}{\mathrm{d}t}\frac{\partial L}{\partial \dot{q}_j} \right) \right)}\delta q_j\mathrm{d}t}\\
    &=0+\int_{t_1}^{t_2}{\sum_{j=1}^n{\left( \frac{\partial L}{\partial q_j}-\left( \frac{\mathrm{d}}{\mathrm{d}t}\frac{\partial L}{\partial \dot{q}_j} \right) \right)}\delta q_j\mathrm{d}t}\\
  \end{aligned}
  $$

- 上式中. 初态和末态的变分为零, 即 $\delta q_j(t_1)=\delta q_j(t_1)=0$. 各个自由度 $\delta q_j$ 相互独立, 要使 $\int_{t_1}^{t_2}{\sum_{j=1}^n{\left( \frac{\partial L}{\partial q_j}-\left( \frac{\mathrm{d}}{\mathrm{d}t}\frac{\partial L}{\partial \dot{q}_j} \right) \right)}\delta q_j\mathrm{d}t}=0$, 只有当系数都等于 0. 运动方程为
  $$
  \frac{\partial L}{\partial q_j}-\left( \frac{\mathrm{d}}{\mathrm{d}t}\frac{\partial L}{\partial \dot{q}_j} \right)=0
  \tag{2.4-6}
  $$
  这就是 Hamilton 原理 $$\int_{t_1}^{t_2}\delta L\;\mathrm{d}t=0$$

- Lagrangian 的时间积分称为作用量 $$S=\int_{t_1}^{t_2}L\;\mathrm{d}t$$ 因此 Hamilton 原理就是 $\delta S=0$.

- Hamilton 原理也称为最小作用量原理. 但是作用量积分只需是 stationary, 未必是极大或者极小.

### 2.5 Lagrange 乘子和约束 

- 多粒子的 Lagrangian 对笛卡尔坐标的变分原理
  $$
  \int_{t_1}^{t_2}{\sum_{k=1}^n{\left( \frac{\partial L}{\partial \boldsymbol{r}_k}-\left( \frac{\mathrm{d}}{\mathrm{d}t}\frac{\partial L}{\partial \dot{\boldsymbol{r}}_k} \right) \right) \cdot}\delta \boldsymbol{r}_k\mathrm{d}t}=0
  $$

- 对于各坐标非独立的情况, Hamilton 原理依然有效, 但是不能直接得出运动方程了. 此时可以引入 Lagrange 乘子, 构造新的 Lagrangian
  $$
  L^{\prime}=L\left( \boldsymbol{r}_1,\boldsymbol{r}_2,\cdots ;\dot{\boldsymbol{r}}_1,\dot{\boldsymbol{r}}_2,\cdots ;t \right) +\sum_{i=1}^C{\lambda _i\left( t \right) f_i\left( \boldsymbol{r}_k,t \right)}
  $$

- 对其作变分得
  $$
  \int_{t_1}^{t_2}{\delta L^{\prime}}\mathrm{d}t=\int_{t_1}^{t_2}{\sum_{k=1}^N{\left( \frac{\partial L}{\partial \boldsymbol{r}_k}-\frac{\mathrm{d}}{\mathrm{d}t}\frac{\partial L}{\partial \dot{\boldsymbol{r}}_k}+\sum_{i=1}^C{\lambda _i\frac{\partial f_i}{\partial \boldsymbol{r}_k}} \right) \cdot}\delta \boldsymbol{r}_k\mathrm{d}t}=0
  $$
  引入了 Lagrange 乘子以后, $\delta \boldsymbol{r}_k$ 又可以视为独立的自由度了. 于是有
  $$
  \tag{2.5-4}
  \begin{aligned}
    &\frac{\partial L^{\prime}}{\partial \mathbf{r}_k}-\frac{\mathrm{d}}{\mathrm{d}t}\frac{\partial L^{\prime}}{\partial \dot{\mathbf{r}}_k}=0\\
    \Rightarrow \quad &\frac{\partial L}{\partial \mathbf{r}_k}-\frac{\mathrm{d}}{\mathrm{d}t}\frac{\partial L}{\partial \dot{\mathbf{r}}_k}+\sum_{i=1}^C{\lambda _i}\frac{\partial f_i}{\partial \mathbf{r}_k}=0\\
  \end{aligned}
  $$

- 新 Lagrangian 对乘子的变分给出约束方程 
  $$
  \frac{\partial L^{\prime}}{\partial \mathbf{r}_k}-\frac{\mathrm{d}}{\mathrm{d}t}\frac{\partial L^{\prime}}{\partial \dot{\mathbf{r}}_k}=0 \quad\Rightarrow\quad f_i(\boldsymbol{r}_{k},t)=0
  $$

- 如果系统的保守力由势能 $V(\mathbf{r}_k)$ 的导数给出, 将 $L=T-V$ 代入运动方程得
  $$
  \tag{2.5-6}
  \begin{aligned}
  &\frac{\partial L}{\partial \mathbf{r}_k}-\frac{\mathrm{d}}{\mathrm{d}t}\frac{\partial L}{\partial \dot{\mathbf{r}}_k}+\sum_{i=1}^C{\lambda _i}\frac{\partial f_i}{\partial \mathbf{r}_k}=0\\
  \Rightarrow\quad & \underset{-\mathbf{F}_{k}}{\underline{\frac{\partial T}{\partial \mathbf{r}_{k}}-\frac{\mathrm{d}}{\mathrm{d} t} \frac{\partial T}{\partial \dot{\mathbf{r}}_{k}}}}+\underset{\mathbf{N}_{k}}{\underline{-\frac{\partial V}{\partial \mathbf{r}_{k}}}}+\underset{\mathbf{C}_k}{\underline{\sum_{i=1}^{C} \lambda_{i} \frac{\partial f_{i}}{\partial \mathbf{r}_{k}}}}=0\\
  \Rightarrow\quad & -合外力+非约束力+约束力=0
  \end{aligned}
  $$

## 3 Lagrangian 的性质

### 3.1 非唯一性

- 系统的 Lagrangian 不唯一. $L$ 与 $L^\prime=aL+b$ 描述相同的运动.

- 若限定变分的各路径两端固定, 则 Lagrangian 的不唯一性可以相差某个函数对时间的全微分
  $$
  L'(\mathbf{q},\dot{\mathbf{q}},t)=L(\mathbf{q},\dot{\mathbf{q}},t)+\frac{\mathrm{d}f(\mathbf{q},t)}{\mathrm{d}t}
  $$
  $L'(\mathbf{q},\dot{\mathbf{q}},t)$ 与 $L(\mathbf{q},\dot{\mathbf{q}},t)$ 给出相同的运动方程, 因为它们的作用量满足关系
  $$
  \begin{aligned}
  S'[\mathbf{q}]&=\int_{t_0}^{t_1}L'(\mathbf{q}(t),\dot{\mathbf{q}}(t),t)\mathrm{d}t\\
  &=\int_{t_0}^{t_1}L(\mathbf{q}(t),\dot{\mathbf{q}}(t),t)\mathrm{d}t+\int_{t_0}^{t_1}\frac{\mathrm{d}f(\mathbf{q}(t),t)}{\mathrm{d}t}\mathrm{d}t\\
  &=S[\mathbf{q}]+f(\mathbf{q}(t_1),t_1)-f(\mathbf{q}(t_0),t_0)
  \end{aligned}
  $$
  求变分时, $\mathbf{q}(t)$ 的端点 $\mathbf{q}(t_0)$ 和 $\mathbf{q}(t_1)$ 都是固定不动的. 就是说最后两项对 $\mathbf{q}(t)$ 的变分是 0,
  $$
  \frac{\delta}{\delta\mathbf{q}(t)}\left(f(\mathbf{q}(t_1),t_1)-f(\mathbf{q}(t_0),t_0)\right)=0
  $$
  因而有 $\frac{\delta}{\delta \mathbf{q}(t)}S'[\mathbf{q}]=\frac{\delta}{\delta\mathbf{q}(t)} S[\mathbf{q}]$.

### 3.2 点变换下的不变性

- 广义坐标变换 $\mathbf{q}\to \mathbf{s}, \mathbf{q}=\mathbf{q}(\mathbf{s},t)$ 下, Lagrangian 的变换为 $L\left(\mathbf{q}(\mathbf{s},t),\dot{\mathbf{q}}(\mathbf{s},\dot{\mathbf{s}},t),t\right)=L'(\mathbf{s},\dot{\mathbf{s}},t)$, Lagrange 方程的形式不变
  $$\frac{\mathrm{d}}{\mathrm{d}t}\frac{\partial L'}{\partial \dot{s}_{i}}=\frac{\partial L'}{\partial s_i}$$

- 可以用坐标变换的方法化简运动方程.

### 3.3 循环坐标和守恒动量

- 广义动量的定义 $p_i=\frac{\partial L}{\partial \dot{q}_i}$ .

- Lagrangian 中不含某个广义坐标, 根据运动方程, 与此广义坐标对应的广义动量不随时间变化
  $$\dot{p}_{i}=\frac{\mathrm{d}}{\mathrm{d}t}\frac{\partial L}{\partial \dot{q}_{i}}=\frac{\partial L}{\partial q_i}=0$$
  该广义动量对时间的积分是一个常数, 是系统的守恒量.
  这是 Noether 定理的一个特例. 
  此坐标称为循环坐标.

### 3.4 能量

- 能量的定义 $$E=\sum_{i=1}^{n}\dot{q}_i\frac{\partial L}{\partial \dot{q}_i}-L\tag{3.4-1}$$

- 构象空间做坐标变换时, 能量不变 $E(\mathbf{q},\dot{\mathbf{q}},t)=E(\mathbf{Q},\dot{\mathbf{Q}},t)$, $\mathbf{q}$ and $\mathbf{Q}$ 是构象空间的一个点的两个坐标.
  可以理解为能量是定义在构象流形上的函数.

- 系统是 closed $\Leftrightarrow$ Lagrangian 不依赖时间. 此时能量是运动积分.

- 求 Lagrangian 对时间 $t$ 的全导数得出
  $$
  \tag{3.4-2}
  \begin{aligned}
    \frac{\mathrm{d}}{\mathrm{d}t}E&=\frac{\mathrm{d}}{\mathrm{d}t}\left( \sum_{i=1}^n{\dot{q}_i\frac{\partial L}{\partial \dot{q}_i}}-L \right)\\
    &=\sum_{i=1}^n{\left( \frac{\mathrm{d}}{\mathrm{d}t}\dot{q}_i \right) \frac{\partial L}{\partial \dot{q}_i}}+\underset{[1]}{\underline{\sum_{i=1}^n{\dot{q}_i\frac{\mathrm{d}}{\mathrm{d}t}\frac{\partial L}{\partial \dot{q}_i}}}}-\frac{\mathrm{d}}{\mathrm{d}t}L\\
    &=\sum_{i=1}^n{\frac{\partial \dot{q}_i}{\partial t}\frac{\partial L}{\partial \dot{q}_i}}+\sum_{i=1}^n{\dot{q}_i\frac{\partial L}{\partial q_i}}-\frac{\mathrm{d}}{\mathrm{d}t}L\\
    &=\underset{3}{\underline{\sum_{i=1}^n{\frac{\partial \dot{q}_i}{\partial t}\frac{\partial L}{\partial \dot{q}_i}}}}+\underset{2}{\underline{\sum_{i=1}^n{\dot{q}_i\frac{\partial L}{\partial q_i}}}}-\left( \frac{\partial L}{\partial t}+\underset{2}{\underline{\sum_{i=1}^n{\frac{\partial q_i}{\partial t}\frac{\partial L}{\partial q_i}}}}+\underset{3}{\underline{\sum_{i=1}^n{\frac{\partial \dot{q}_i}{\partial t}\frac{\partial L}{\partial \dot{q}_i}}}} \right)\\
    &=-\frac{\partial L}{\partial t}\\
  \end{aligned}
  $$
  对 [1] 使用了 Euler-Lagrange 方程.

- 由上式可见, 如果 Lagrangian 不含时间, 即 $\frac{\partial L}{\partial t}=0$, 则系统能量 $E(\mathbf{q}(t),\dot{\mathbf{q}}(t),t)$ 守恒.

- k-homogeneous 函数: $f(s\mathbf{v})=s^{k}f(\mathbf{v})$.
  Euler's homogeneous 函数定理:
  $$kf(x_1,\cdots,x_n)=\sum_{i=1}^{n}x_{i}\frac{\partial f}{\partial x_i}(x_1,\cdots,x_n)$$
  动能是广义坐标的 2-homogeneous 函数
  $$2T=\sum_{i=1}^{n}\dot{q}_i\frac{\partial T}{\partial \dot{q}_i}$$

- 如果 $V$ 只是 $\mathbf{r}$ 的函数, 不是速度的函数, 则
  $$\tag{3.4-5}\sum_{i=1}^{n}\dot{q}_{i}\frac{\partial L}{\partial \dot{q}_{i}}=\sum_{i=1}^{n}\dot{q}_{i}\frac{\partial (T-V)}{\partial \dot{q}_{i}}=\sum_{i=1}^{n}\dot{q}_{i}\frac{\partial T}{\partial \dot{q}_{i}}=2T$$

- 由上式得
  $$\begin{aligned}
  &\sum_{i=1}^{n}\dot{q}_{i}\frac{\partial L}{\partial \dot{q}_{i}}-L=2T-L\\
  \Rightarrow & E=2T-(T-V)\\
  \Rightarrow & E=T+V
  \end{aligned}$$
  即系统得能量等于动能与势能的和.

### 3.5 力学相似性

### 3.6 相互作用粒子

- 一个系统可分成两个无相互作用的子系统时, 系统的 Lagrangian 可以写成两个子系统的 Lagrangian 的和;
  如果有相互作用, 则加入两个子系统的相互作用项.

## 4 例

## 5 扩展到包含非守恒力

## 6 其它内容

<!--stackedit_data:
eyJoaXN0b3J5IjpbNjAzMDYzNjgwLC0xMTYzNDgwNzk1XX0=
-->