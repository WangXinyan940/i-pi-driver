# i-pi-driver

作为非常实用的PIMD engine，开源项目i-PI一直在疯狂更新，保持着极高的活跃度。I-PI不提供量子化学和分子力学计算，只负责进行Path Integral动力学积分运算，势能的计算则交给了其他通用软件。目前，i-PI仅支持与CP2K、Quantum-Espresso和Lammps，而本项目则为i-PI提供了一个通用的Python Socket Client，用以读取任意QM或MM软件的输出文件并传递给i-PI，以承担i-PI与QM或MM软件间的通信工作。

项目中的BaseDriver类实现了与i-PI server的基本交互，包括接收盒子及坐标信息以及发送坐标对应的力、能量与维里。用户使用时应继承driver.py中的BaseDriver类，并重新定义其中的grad函数实现能量及梯度的计算。其形式为：
```
class NewDriver(BaseDriver):
    def grad(self, crd):
        #do the calculation
        ...
        return energy, grad
```

其中crd、energy与grad均为SI单位。driver.py中提供了常数用以实现常用单位到SI单位的转换，示例如下：

```
# constants from driver.py
BOHR = 5.291772108e-11  # Bohr -> m
ANGSTROM = 1e-10  # angstrom -> m
AMU = 1.660539040e-27  # amu -> kg
FEMTO = 1e-15
PICO = 1e-12
EH = 4.35974417e-18  # Hartrees -> J
EV = 1.6021766209e-19  # eV -> J
H = 6.626069934e-34  # Planck const
KB = 1.38064852e-23  # Boltzmann const
MOLE = 6.02214129e23
KJ = 1000.0 # kJ -> J
KCAL = 4184.0 # kcal -> J

energy_si = 12.5 * EH # convert Hartree to J
force_si = np.array([2.0,1.0,1.0]) * (EH / BOHR) # convert Hartree / Bohr to Newton
energy_amu = energy_si / EH # convert J to Hartree
force_amu = force_si / (EH / BOHR) # convert Newton to Hartree / Bohr
```

代码中所有物理量均为SI单位，建议从输出文件中读取数值之后立即转化为SI单位以便于计算。项目同时提供了与Gaussian的接口代码（class GaussDriver），并在 example/h2/ 中提供了使用Gaussian对氢气分子做Path Integral Langevin Dynamics的基本示例。
作为科学计算程序，本项目尽量遵循了PEP8规则，未遵循的部分多半为人力所不能及，或者是单纯的懒。不过反正八成只有我一个人用，就凑合着看吧。
