# 感知RIS论文结构

硬件名字：Sensing RIS

论文名字：Sensing RIS: Enabling Dimension-Independent CSI Acquisition for Beamforming

1. Introduction
   - RIS是未来通信技术，有...好处；RIS是由...构成，通过调相的方式提供波束赋形增益；实现可靠波束赋形需要准确的CSI，因此信道估计是RIS系统中的重要问题。
   - 信道估计已经在传统MIMO系统中被广泛研究，然而RIS系统面临两个新挑战。一是RIS不具备信号处理能力。二是待估计参数量增加，导频开销大大增加。
   - A. Prior Contributions
     - 现有系统的信道估计和波束赋形是割裂的。现有信道估计方法可以被大致分为以下几类...现有波束赋形方法可以被大致分为以下几类...
     - 总的来说，现有方法(1)都是对信道的数学特征进行提取，(2)两个步骤割裂，导致导频开销只能被降低至O(N)。然而对于大规模的RIS阵列（举例），导频开销仍然很高，是否存在一种方案使得导频开销与RIS单元数无关？
   - B. Our Contributions
     - 为了解决RIS系统导频开销过大的问题，受启发于电磁干涉，我们提出了一种基于IRF的信道估计方案。
       - 我们提出了一种新的信号发送方法，并揭示了RIS系统中的电磁干涉现象。区别于以往分别发信号，在同时发信号的系统中，我们推导了RIS单元上信号幅度的时域变化情况，并将此命名为IRF。
       - To exploit this IRF, we employ sensing RIS...并进一步提出了相应的几种算法来实现RIS系统信道估计，把两个步骤联合在一起，使得导频开销可以降至O(1)。
       - 我们给出了此系统中的CRLB的表达式，并给出了其渐近结果。我们通过数值仿真比较了几种参数估计方法与CRLB的关系。
       - 仿真结果...
2. System Model
   - RIS-Aided SISO System
     - Beamforming design
     - Channel estimation
3. 信号发送/Signal Model for IRF（目前文章中使用的标题：干涉随机场，Interference Random Field）
   - IRF图解
   - IRF如何辅助信道估计？
4. Sensing RIS-based Channel Estimation
   - Hardware Architecture of Sensing RIS
   - 算法1：LS算法/算法2：精确度太低，MMSE算法/算法3：复杂度太高，von Mises-EM算法
   - MISO系统到MIMO系统的扩展方法（各种扩展的讨论）
5. Performance Analysis
   - CRLB
6. Simulation Results
7. Conclusion
