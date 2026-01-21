import math

class CementProductionModel:
    def __init__(self):
        # 模型初始化
        self.constants = {}
        self.parameters = {}
        self.control_variables = {}
        self.variables = {}
        self.section_configs = {} # 三段工艺配置
        self.cells = [] # 有限体积单元列表
        self.converged = False # 收敛标志
        self.max_iterations = 1000 # 最大迭代次数
        self.convergence_threshold = 1e-4
        self.dt = 0.08  # 时间步长（平衡反应时间）
        self.dz = 1.0  # 空间步长
        
    def define_constants(self):
        # 定义固定常数
        # 物理常数
        self.constants['R'] = 8.314  # 气体常数，J/(mol·K)
        self.constants['sigma'] = 5.67e-8  # 斯特藩-玻尔兹曼常数，W/(m²·K⁴)
        
        # 热力学标准条件
        self.constants['T0'] = 298.15  # 标准热力学温度，K
        self.constants['P0'] = 101325  # 标准大气压，Pa
        self.constants['Tref'] = 1200  # 参考温度，K
        
        # 反应焓变（单位：J/mol）
        self.constants['reaction_enthalpies'] = {
            'r1': 179170,  # CaCO3 → CaO + CO2 (吸热)
            'r2': -75220,  # 2CaO + SiO2 → C2S (放热)
            'r3': 9920,    # CaO + C2S → C3S (吸热)
            'r4': -1540,   # 3CaO + Al2O3 → C3A (放热)
            'r5': 170880,  # 4CaO + Al2O3 + Fe2O3 → C4AF (吸热)
            'r6': -565976, # 2CO + O2 → 2CO2 (强放热)
            'r7': -41168,  # CO + H2O → CO2 + H2 (放热)
            'r8': -483640, # 2H2 + O2 → 2H2O (强放热)
            'r9': -221044, # 2C + O2 → 2CO (放热)
            'r10': 131298, # C + H2O → CO + H2 (吸热)
            'r11': 172466  # C + CO2 → 2CO (吸热)
        }
        
        # 固体物质属性
        self.constants['solid_properties'] = {
            'CaCO3': {
                'thermal_conductivity': 2.248,  # W/(K·m)
                'density': 2.71,  # g/cm³
                'molar_mass': 100.09  # g/mol
            },
            'CaO': {
                'thermal_conductivity': 30.1,  # W/(K·m)
                'density': 3.34,  # g/cm³
                'molar_mass': 56.08  # g/mol
            },
            'SiO2': {
                'thermal_conductivity': 3.4,  # W/(K·m)
                'density': 2.65,  # g/cm³
                'molar_mass': 60.09  # g/mol
            },
            'Al2O3': {
                'thermal_conductivity': [12.8, 14.5],  # W/(K·m)
                'density': 3.99,  # g/cm³
                'molar_mass': 101.96  # g/mol
            },
            'Fe2O3': {
                'thermal_conductivity': [0.3, 0.37],  # W/(K·m)
                'density': 5.25,  # g/cm³
                'molar_mass': 159.69  # g/mol
            },
            'C2S': {
                'thermal_conductivity': 3.45,  # W/(K·m)
                'density': 3.35,  # g/cm³
                'molar_mass': 172.24  # g/mol
            },
            'C3S': {
                'thermal_conductivity': 3.35,  # W/(K·m)
                'density': 3.13,  # g/cm³
                'molar_mass': 228.32  # g/mol
            },
            'C3A': {
                'thermal_conductivity': 3.74,  # W/(K·m)
                'density': 3.04,  # g/cm³
                'molar_mass': 270.19  # g/mol
            },
            'C4AF': {
                'thermal_conductivity': 3.17,  # W/(K·m)
                'density': [3.7, 3.9],  # g/cm³
                'molar_mass': 485.97  # g/mol
            }
        }
        
        # 气体物质属性
        self.constants['gas_properties'] = {
            'CO2': {
                'thermal_conductivity': {300: 16.77, 1000: 70.78},  # 10⁻³ W/(K·m)
                'molar_mass': 44.01,  # g/mol
                'viscosity': {300: 15.0, 1000: 41.18},  # μPa·s
                'diffusion_volume': 16.3  # cm³
            },
            'N2': {
                'thermal_conductivity': {300: 25.97, 1000: 65.36},  # 10⁻³ W/(K·m)
                'molar_mass': 28.014,  # g/mol
                'viscosity': {300: 17.89, 1000: 41.54},  # μPa·s
                'diffusion_volume': 18.5  # cm³
            },
            'O2': {
                'thermal_conductivity': {300: 26.49, 1000: 71.55},  # 10⁻³ W/(K·m)
                'molar_mass': 31.998,  # g/mol
                'viscosity': {300: 20.65, 1000: 49.12},  # μPa·s
                'diffusion_volume': 16.3  # cm³
            },
            'Ar': {
                'thermal_conductivity': {300: 17.84, 1000: 43.58},  # 10⁻³ W/(K·m)
                'molar_mass': 39.948,  # g/mol
                'viscosity': {300: 22.74, 1000: 55.69},  # μPa·s
                'diffusion_volume': 16.2  # cm³
            },
            'CO': {
                'thermal_conductivity': {300: 25, 600: 43.2},  # 10⁻³ W/(K·m)
                'molar_mass': 28.010,  # g/mol
                'viscosity': {300: 17.8, 1000: 29.1},  # μPa·s
                'diffusion_volume': 18  # cm³
            },
            'C_sus': {
                'thermal_conductivity': 0.1,  # W/(m·K)
                'molar_mass': 12.011,  # g/mol
                'viscosity': 1e-5,     # Pa·s
                'diffusion_volume': 15.9  # cm³
            },
            'H2O': {
                'thermal_conductivity': {300: 609.50, 1000: 95.877},  # 10⁻³ W/(K·m)
                'molar_mass': 18.015,  # g/mol
                'viscosity': {300: 853.74, 1000: 37.615},  # μPa·s
                'diffusion_volume': 13.1  # cm³
            },
            'H2': {
                'thermal_conductivity': {300: 193.1, 1000: 459.7},  # 10⁻³ W/(K·m)
                'molar_mass': 2.016,  # g/mol
                'viscosity': {300: 8.938, 1000: 20.73},  # μPa·s
                'diffusion_volume': 6.12  # cm³
            }
        }
        
        # 摩尔热容
        self.constants['molar_heat_capacity'] = {
            'CaO': {'C0': 71.69, 'C1': -3.08, 'C2': 0.22, 'T_range': (200, 1800)},
            'SiO2': {'C0': 58.91, 'C1': 5.02, 'C2': 0, 'T_range': (844, 1800)},
            'Al2O3': {'C0': 233.004, 'C1': -19.59, 'C2': 0.94, 'T_range': (200, 1800)},
            'Fe2O3': {'C0': 103.9, 'C1': 0, 'C2': 0, 'T_range': (1650, 1800)},
            'C2S': {'C0': 199.6, 'C1': -2.33, 'C2': 0, 'T_range': (200, 1800)},
            'C3S': {'C0': 333.92, 'C1': 0, 'C2': 0, 'T_range': (298, 1800)},
            'C3A': {'C0': 260.58, 'C1': 9.58/2, 'C2': 0, 'T_range': (298, 1800)},
            'C4AF': {'C0': 374.43, 'C1': 36.4, 'C2': 0, 'T_range': (298, 1863)},
            'CO2': {'C0': 25.98, 'C1': 43.61, 'C2': -1.49, 'T_range': (298, 1500)},
            'N2': {'C0': 27.31, 'C1': 5.19, 'C2': -1.55e-4, 'T_range': (298, 1500)},
            'O2': {'C0': 25.82, 'C1': 12.63, 'C2': -0.36, 'T_range': (298, 1100)},
            'Ar': {'C0': 20.79, 'C1': 0, 'C2': 0, 'T_range': (298, 1500)},
            'CO': {'C0': 26.87, 'C1': 6.94, 'C2': -0.08, 'T_range': (298, 1500)},
            'C_sus': {'C0': -0.45, 'C1': 35.53, 'C2': -1.31, 'T_range': (298, 1500)},
            'H2S': {'C0': 30.89, 'C1': 7.86, 'C2': 0.25, 'T_range': (298, 1500)},
            'H2O': {'C0': 28.95, 'C1': -0.58, 'C2': 0.19, 'T_range': (298, 1500)}
        }
        
        # 标准生成焓
        self.constants['standard_enthalpy'] = {
            # 固体
            'CaCO3': -1207.6 * 1000 ,  # kJ/mol
            'CaO': -634.92 * 1000 ,  # kJ/mol
            'SiO2': -910.94 * 1000 ,  # kJ/mol
            'Al2O3': -1675.7 * 1000 ,  # kJ/mol
            'Fe2O3': -825.50 * 1000 ,  # kJ/mol
            'C2S': -2256 * 1000 ,  # kJ/mol
            'C3S': -2881 * 1000 ,  # kJ/mol
            'C3A': -3582 * 1000 ,  # kJ/mol
            'C4AF': -4870 * 1000 ,  # kJ/mol
            # 气体
            'CO2': -393.51 * 1000 ,  # kJ/mol
            'N2': 0,  # kJ/mol
            'O2': 0,  # kJ/mol
            'Ar': 0,  # kJ/mol
            'CO': -110.522 * 1000 ,  # kJ/mol
            'C_sus': 0,  # kJ/mol 
            'H2O': -241.82 * 1000 ,  # kJ/mol
            'H2': 0  # kJ/mol
        }
        
        # WSGG模型系数
        self.constants['wsgg_coefficients'] = {
            'K1': [0.055, 0.88, 10, 135],
            'K2': [0.012, -0.021, -1.6, -35],
            'C1': [
                [0.358, 0.392, 0.142, 0.0798],  # C1_j,1
                [0.0731, -0.212, -0.0831, -0.0370],  # C1_j,2
                [-0.0466, 0.0191, 0.0148, 0.0023]  # C1_j,3
            ],
            'C2': [
                [0.165, 0.291, 0.348, 0.0866],  # C2_j,1
                [-0.0554, 0.644, -0.294, -0.106],  # C2_j,2
                [0.0930, -0.209, 0.0662, 0.0305]  # C2_j,3
            ],
            'C3': [
                [0.0598, 0.0784, -0.122, -0.0127],  # C3_j,1
                [0.0028, -0.197, 0.118, 0.0169],  # C3_j,2
                [-0.0256, 0.0662, -0.0295, -0.0051]  # C3_j,3
            ]
        }
        
        # 反应速率系数
        self.constants['reaction_rate_coefficients'] = {
            'r1': {'kr': 1e8 * 1000 , 'n': 0, 'EA': 175.7 * 1000 , 'alpha1': 1, 'alpha2': 0, 'alpha3': 0, 'beta2': 0, 'unit': 'kg/(m³·s)'},
            'r2': {'kr': 1e9 * 1000 , 'n': 0, 'EA': 240 * 1000 , 'alpha1': 2, 'alpha2': 1, 'alpha3': 0, 'beta2': 0, 'unit': 'kg/(m³·s)'},  # 从1e7提升到1e9
            'r3': {'kr': 1e10 * 1000 , 'n': 0, 'EA': 420 * 1000 , 'alpha1': 1, 'alpha2': 1, 'alpha3': 0, 'beta2': 0, 'unit': 'kg/(m³·s)'},  # 从1e9提升到1e10
            'r4': {'kr': 1e9 * 1000 , 'n': 0, 'EA': 310 * 1000 , 'alpha1': 3, 'alpha2': 1, 'alpha3': 0, 'beta2': 0, 'unit': 'kg/(m³·s)'},   # 从1e8提升到1e9
            'r5': {'kr': 1e9 * 1000 , 'n': 0, 'EA': 330 * 1000 , 'alpha1': 4, 'alpha2': 1, 'alpha3': 1, 'beta2': 0, 'unit': 'kg/(m³·s)'},   # 从1e8提升到1e9
            'r6': {'kr': 7.0e4, 'n': 0, 'EA': 66.5 * 1000 , 'alpha1': 1, 'alpha2': 1, 'alpha3': 0, 'beta2': 0, 'unit': 'mol/(m³·s)'},
            'r7': {'kr': 2.8e6, 'n': 0, 'EA': 83.7 * 1000 , 'alpha1': 1, 'alpha2': 1, 'alpha3': 0, 'beta2': 0, 'unit': 'mol/(m³·s)'},
            'r8': {'kr': 1.4e6, 'n': 0.5, 'EA': 295.5 * 1000 , 'alpha1': 1, 'alpha2': 1, 'alpha3': 0, 'beta2': 0, 'unit': 'mol/(m³·s)'},
            'r9': {'kr': 8.8e11, 'n': 0, 'EA': 239 * 1000 , 'alpha1': 0.5, 'alpha2': 0.5, 'alpha3': 0, 'beta2': 0, 'unit': 'mol/(m³·s)'},
            'r10': {'kr': 2.6e8, 'n': 0, 'EA': 237 * 1000 , 'alpha1': 0, 'alpha2': 0, 'alpha3': 0, 'beta2': 0.6, 'unit': 's⁻¹'},
            'r11': {'kr': 3.1e6, 'n': 0, 'EA': 215 * 1000 , 'alpha1': 0, 'alpha2': 0, 'alpha3': 0, 'beta2': 0.4, 'unit': 's⁻¹'}
        }
        
        # 化学反应式
        self.constants['reactions'] = {
            'r1': {'name': 'CaCO3分解', 'equation': 'CaCO3 → CaO + CO2'},
            'r2': {'name': 'C2S形成', 'equation': '2CaO + SiO2 → C2S'},
            'r3': {'name': 'C3S形成', 'equation': 'CaO + C2S → C3S'},
            'r4': {'name': 'C3A形成', 'equation': '3CaO + Al2O3 → C3A'},
            'r5': {'name': 'C4AF形成', 'equation': '4CaO + Al2O3 + Fe2O3 → C4AF'},
            'r6': {'name': 'CO燃烧', 'equation': '2CO + O2 → 2CO2'},
            'r7': {'name': '水煤气变换反应', 'equation': 'CO + H2O → CO2 + H2'},
            'r8': {'name': 'H2燃烧', 'equation': '2H2 + O2 → 2H2O'},
            'r9': {'name': 'C氧化生成CO', 'equation': '2C + O2 → 2CO'},
            'r10': {'name': 'C与水蒸气反应', 'equation': 'C + H2O → CO + H2'},
            'r11': {'name': 'Boudouard反应', 'equation': 'C + CO2 → 2CO'}
        }
        
        # 化学计量矩阵 (组分数 × 反应数)
        # 反应物系数为负，产物为正，常数
        # 组分顺序: CaCO3, CaO, SiO2, Al2O3, Fe2O3, C2S, C3S, C3A, C4AF, C_sus, CO2, CO, O2, H2O, H2
        self.constants['stoichiometric_matrix'] = {
            'components': ['CaCO3', 'CaO', 'SiO2', 'Al2O3', 'Fe2O3', 'C2S', 'C3S', 'C3A', 'C4AF', 'C_sus', 'CO2', 'CO', 'O2', 'H2O', 'H2'],
            'matrix': [
                # r1  r2  r3  r4  r5  r6  r7  r8  r9  r10 r11
                [-1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],  # CaCO3
                [+1, -2, -1, -3, -4,  0,  0,  0,  0,  0,  0],  # CaO
                [0,  -1,  0,  0,  0,  0,  0,  0,  0,  0,  0],  # SiO2
                [0,   0,  0, -1, -1,  0,  0,  0,  0,  0,  0],  # Al2O3
                [0,   0,  0,  0, -1,  0,  0,  0,  0,  0,  0],  # Fe2O3
                [0,  +1, -1,  0,  0,  0,  0,  0,  0,  0,  0],  # C2S
                [0,   0, +1,  0,  0,  0,  0,  0,  0,  0,  0],  # C3S
                [0,   0,  0, +1,  0,  0,  0,  0,  0,  0,  0],  # C3A
                [0,   0,  0,  0, +1,  0,  0,  0,  0,  0,  0],  # C4AF
                [0,   0,  0,  0,  0,  0,  0,  0, -2, -1, -1],  # C_sus
                [+1,  0,  0,  0,  0, +2, +1,  0,  0,  0, -1],  # CO2
                [0,   0,  0,  0,  0, -2, -1,  0, +2, +1, +2],  # CO
                [0,   0,  0,  0,  0, -1,  0, -1, -1,  0,  0],  # O2
                [0,   0,  0,  0,  0,  0, -1, +2,  0, -1,  0],  # H2O
                [0,   0,  0,  0,  0,  0, +1,  0,  0, +1,  0]   # H2
            ]
        }
    
    def define_parameters(self):
        # 定义系统参数
        
        # 物理参数
        self.parameters['a_omega'] = 0.05  # 休止角比例系数
        self.parameters['b_omega'] = 0.6  # 休止角截距常数，rad
        self.parameters['mu_g'] = 4.5e-5  # 气体平均粘度，Pa·s
        self.parameters['k_s'] = 3.5  # 固体热导率，W/(K·m)
        self.parameters['d_p'] = 3e-5  # 颗粒直径，m
        self.parameters['epsilon_s'] = 0.9  # 固体发射率
        
        # 系统几何参数
        self.parameters['equipment'] = {
            'preheater': {
                'length': 4,  # 长度，m
                'segments': 4,  # 分段数量
                'radius': 1  # 半径，m
            },
            'calciner': {
                'length': 5,  # 长度，m
                'segments': 5,  # 分段数量
                'radius': 1  # 半径，m
            },
            'kiln': {
                'length': 7,  # 长度，m
                'segments': 7,  # 分段数量
                'radius': 1,  # 半径，m
                'angle': math.radians(2)  # 倾角，rad
            }
        }
        
        # 计算各设备的空间步长
        self.parameters['equipment']['preheater']['dz'] = self.parameters['equipment']['preheater']['length'] / self.parameters['equipment']['preheater']['segments']
        self.parameters['equipment']['calciner']['dz'] = self.parameters['equipment']['calciner']['length'] / self.parameters['equipment']['calciner']['segments']
        self.parameters['equipment']['kiln']['dz'] = self.parameters['equipment']['kiln']['length'] / self.parameters['equipment']['kiln']['segments']
    
    def define_control_variables(self):
        # 定义控制量：连续进料导入预热器入口
        
        # 1. 固体进料（生料）
        self.control_variables['solid_feed'] = {
            'temperature': 600,  # K
            'total_rate': 666.7 * 1000  / 3600,  # 总进料速率，g/s
            'composition': {
                'CaCO3': 0.76,  # 质量分数
                'SiO2': 0.11,
                'Al2O3': 0.03,
                'Fe2O3': 0.02,
                'C2S': 0.02
            }
        }
        
        # 2. 气体进料（风）
        self.control_variables['gas_feed'] = {
            'total_rate': 15,  # 总风量，mol/s
            'initial_velocity': 0.5,  # 初始风速，m/s
            'temperature': 1300,  # K （提升至 1300K）
            'composition': {
                'O2': 0.21,  # 摩尔分数
                'N2': 0.70,
                'H2O': 0.09
            }
        }
        
        # 3. 燃料量（煤粉）
        self.control_variables['fuel'] = {
            'rate': 30 * 1000 / 3600,  # 燃料速率，g/s
            'composition': {
                'C_sus': 0.70,  # 质量占比
                'CO': 0.18,
                'CO2': 0.10,
                'H2O': 0.02
            }
        }
        
        # 4. 窑转速
        self.control_variables['omega'] = 0.4189  # 窑转速，rad/s
        
        # 5. 操作参数
        self.control_variables['excess_air_ratio'] = 1.1  # 空气过剩系数

    def configure_sections(self):
        # 三段工艺差异化参数化配置
        # 仅写「规则标识」，工艺差异与计算逻辑完全解耦
        
        self.section_configs['preheater'] = {
            'energy_eq': None,  # 无反应和公式
            'reactions': None,  # 无反应
            'solid_velocity': 0,  # 固相流速为0
            'gas_velocity': 'given_parameter',  # 气相流速为给定参数
            'heat_transfer': None,  # 无需对流换热计算
            'internal_energy_constraint': None,  # 无需内能密度约束
            'solid_concentration_constant': True,  # 固相浓度和摩尔数恒定不变
            'gas_concentration_constant': True,  # 气相浓度和摩尔数不变
            # 移除微分方程和代数方程
        }
        
        self.section_configs['calciner'] = {
            'energy_eq': 'formula_42',  # 能量方程单方程调用公式42
            'reactions': ['r1', 'r6', 'r7', 'r8', 'r9', 'r10', 'r11'],  # 仅r1分解反应
            'solid_velocity': 'formula_71',  # 固相流速调用公式71
            'heat_transfer': 'formula_35',  # 对流换热调用公式35
            'internal_energy_constraint': 'formula_65',  # 内能密度约束调用公式65
            'solid_concentration_constant': False,  # 固相浓度变化
            'gas_concentration_constant': False,  # 气相浓度变化
            # 移除微分方程和代数方程
        }
        
        self.section_configs['kiln'] = {
            'energy_eq': ['formula_44', 'formula_45'],  # 能量方程气固双方程调用公式44、45
            'reactions': ['r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'r7', 'r8', 'r9', 'r10', 'r11'],  # r2~r11熟料反应
            'solid_velocity': 'formula_13',  # 固相流速调用公式13
            'heat_transfer': 'formula_39',  # 对流换热调用公式39
            'internal_energy_constraint': ['formula_2', 'formula_3'],  # 内能密度约束调用公式2、3
            'solid_concentration_constant': False,  # 固相浓度变化
            'gas_concentration_constant': False,  # 气相浓度变化
            # 移除微分方程和代数方程
        }
    
    def spatial_discretization(self):
        # 有限体积空间离散
        # 将预热器→分解炉→回转窑沿物料流动方向做一维连续有限体积离散
        # 建立全局数组，统一存储所有单元的核心变量
        
        # 从参数中获取各设备的分段数量
        preheater_segments = self.parameters['equipment']['preheater']['segments']  # 预热器分段数量
        calciner_segments = self.parameters['equipment']['calciner']['segments']    # 分解炉分段数量
        kiln_segments = self.parameters['equipment']['kiln']['segments']       # 回转窑分段数量
        
        # 计算总分段数量
        total_segments = preheater_segments + calciner_segments + kiln_segments
        
        # 初始化全局变量数组
        self.variables = { #状态变量子字典，存储需要动态求解的变量
            # 状态变量 (向量形式：x=[C, U^]^T)
            'state_variables': {
                'C': [],  # 各相组分浓度向量，每个元素包含该段的所有组分浓度
                'U_hat': []  # 单位体积内能向量，每个元素包含该段的内能密度
            },
            
            # 代数变量 (向量形式：y=[Tg, Ts, P]^T)
            'algebraic_variables': { #代数变量子字典，存储实时计算的辅助变量
                'Tg': [],  # 气体温度向量，每个元素包含该段的气体温度
                'Ts': [],  # 固体温度向量，每个元素包含该段的固体温度
                'P': []  # 系统压力向量，每个元素包含该段的压力
            },
            
            # 设备类型标识，用于确定每个段所属的设备
            'device_type': []  # 每个元素为 'preheater', 'calciner', 或 'kiln'
        }
        
        # 获取所有组分名称
        #从化学计量矩阵中获取所有组分名称，用于后续初始化各单元的浓度字典
        components = self.constants['stoichiometric_matrix']['components']
        
        # 初始化每个分段的变量
        for i in range(total_segments):
            # 确定设备类型
            if i < preheater_segments:
                device = 'preheater'
            elif i < preheater_segments + calciner_segments:
                device = 'calciner'
            else:
                device = 'kiln'
            
            # 初始化状态变量
            # C: 各相组分浓度向量，初始化为0
            segment_C = {component: 0.0 for component in components}
            self.variables['state_variables']['C'].append(segment_C)
            
            # U_hat: 单位体积内能向量，初始化为0
            #将当前单元的初始内能密度添加到全局状态变量的U_hat列表中
            self.variables['state_variables']['U_hat'].append(0.0)
            
            # 初始化代数变量
            # Tg: 气体温度向量，初始化为300K
            self.variables['algebraic_variables']['Tg'].append(600.0)
            # Ts: 固体温度向量，初始化为300K
            self.variables['algebraic_variables']['Ts'].append(1300.0)
            
            # P: 系统压力向量，初始化为标准大气压
            self.variables['algebraic_variables']['P'].append(self.constants['P0'])
            
            # 记录设备类型
            self.variables['device_type'].append(device)
        
        # 初始化单元列表，每个单元对应一个分段
        self.cells = []
        for i in range(total_segments): #再次遍历所有离散单元，构建单元字典并添加到cells列表
            # 计算每个单元的delta_z
            if i < preheater_segments:
                # 预热器：长4m，分为4段
                delta_z = 4.0 / preheater_segments
            elif i < preheater_segments + calciner_segments:
                # 分解炉：长5m，分为5段
                delta_z = 5.0 / calciner_segments
            else:
                # 回转窑：长7m，分为7段
                delta_z = 7.0 / kiln_segments
            
            cell = {
                'index': i, #存储单元索引（从 0 开始）
                'device_type': self.variables['device_type'][i], #存储单元的设备类型
                'delta_z': delta_z,
                'state_variables': { #单元的状态变量子字典
                    'C': self.variables['state_variables']['C'][i], #引用全局C列表中当前单元的浓度字典
                    'U_hat': self.variables['state_variables']['U_hat'][i]
                },
                'algebraic_variables': { #单元的代数变量子字典
                    'Tg': self.variables['algebraic_variables']['Tg'][i],
                    'Ts': self.variables['algebraic_variables']['Ts'][i],
                    'P': self.variables['algebraic_variables']['P'][i]
                }
            }
            self.cells.append(cell) #将当前单元字典添加到cells列表中
        
        print(f"空间离散完成：{total_segments} 个单元，其中预热器 {preheater_segments} 个，分解炉 {calciner_segments} 个，回转窑 {kiln_segments} 个")
    
    # 通用计算子函数库 - 仅定义公式签名，不实现具体内容
    # 基础定义与代数约束方程
    def formula_1(self, V_hat_g, V_hat_s):
        # 1气体与固体体积分数守恒
        # 输入：V^g（气体相的单位体积分数，公式67输出）、V^s（固体相的单位体积体积分数，公式66输出）
        # 输出：体积分数守恒为1
        return V_hat_g + V_hat_s
    
    def formula_2(self, P, H_g, V_g):
        # 2回转窑气体内能密度定义
        # 输入：P（压力，变量）、H_g（气体单位体积焓，公式64输出）、V_g（气体相的单位体积分数，公式67输出）
        # 输出：U_g（气体内能密度）
        return H_g - P * V_g
    
    def formula_3(self, H_s):
        # 3回转窑固体内能密度定义
        # 输入：H_s（固体单位体积焓，公式63输出）
        # 输出：U_s（固体内能密度）
        return H_s
    
    def formula_4(self, Hf_i, n_i, integral_enthalpy):
        # 4气/固相总焓计算
        # 输入：Hf,i（标准生成焓，常数）、ni（各组分摩尔数列表，变量，公式57输出）、integral_enthalpy（各组分积分焓列表，公式5输出）
        # 输出：H（气体或固体的总焓）
        # 公式：H(T,P,n) = Σ n_i (ΔHf,i + ∫T0到T cp,i(τ)dτ)
        total_enthalpy = 0.0
        for hf, n, integral in zip(Hf_i, n_i, integral_enthalpy):
            total_enthalpy += n * (hf + integral)
        return total_enthalpy
    
    def formula_5(self, C0, C1, C2, T0, T):
        # 5组分i积分焓
        # 输入：C0,C1,C2（单个组分i的热容系数，常数）、T0（常数）、T（温度，变量）
        # 输出：单个组分i的摩尔热容积分值
        # 公式：∫T0到T cp,i(τ)dτ = C0(T-T0) + 1/2 C1(T²-T0²) + 1/3 C2(T³-T0³)
        return C0 * (T - T0) + 0.5 * C1 * (T**2 - T0**2) + (1/3) * C2 * (T**3 - T0**3)
    
    def formula_6(self, C0, C1, C2, T):
        # 6摩尔热容
        # 输入：C0,C1,C2（组分i的热容系数，常数）、T（温度，变量）
        # 输出：cp,i（组分i的摩尔热容）
        return C0 + C1 * T + C2 * T**2
    
    def formula_7(self, R, T, P, n_i):
        # 7气体体积计算
        # 输入：R（气体常数，常数）、T（温度，变量）、P（压力，变量）、ni（各组分摩尔数列表/字典，变量，公式57输出）
        # 输出：Vg（气体体积）
        # 公式：V_g = (RT/P) * Σn_i
        total_moles = sum(n_i.values()) if isinstance(n_i, dict) else sum(n_i)
        return (R * T / P) * total_moles
    
    def formula_8(self, M_i, rho_i, n_i):
        # 8固体体积计算
        # 输入：Mi（各组分摩尔质量列表/字典，常数）、ρi（各组分密度列表/字典，常数）、ni（各组分摩尔数列表/字典，变量，公式57输出）
        # 输出：Vs（固体体积）
        # 公式：V_s = Σ(M_i/ρ_i)n_i
        total_volume = 0.0
        
        # 辅助函数：获取数值，如果是列表则取第一个元素
        def get_value(val):
            if isinstance(val, list) and len(val) > 0:
                return val[0]
            return val
        
        if isinstance(M_i, dict) and isinstance(rho_i, dict) and isinstance(n_i, dict):
            # 所有输入都是字典，遍历共同的组分
            for component in M_i:
                if component in rho_i and component in n_i:
                    m = get_value(M_i[component])
                    rho = get_value(rho_i[component])
                    n = get_value(n_i[component])
                    total_volume += (m / rho) * n
        elif hasattr(M_i, '__iter__') and hasattr(rho_i, '__iter__') and hasattr(n_i, '__iter__'):
            # 所有输入都是可迭代对象，按顺序遍历
            for M, rho, n in zip(M_i, rho_i, n_i):
                m = get_value(M)
                rho_val = get_value(rho)
                n_val = get_value(n)
                total_volume += (m / rho_val) * n_val
        else:
            # 单个组分情况
            m = get_value(M_i)
            rho_val = get_value(rho_i)
            n_val = get_value(n_i)
            total_volume += (m / rho_val) * n_val
        return total_volume
    
    # 流动与传输方程
    def formula_9(self, V_delta, A_g):
        # 9水力直径计算
        # 输入：VΔ（单段总体积，参数，公式58输出）、Ag（气体横截面积，公式78输出，变量）
        # 输出：DH（水力直径）
        # 公式：D_H = 4VΔ / A_g
        return 4 * V_delta / A_g
    
    def formula_10(self, DH, mu_g, rho_g, dP_dz):
        # 10气体轴向速度
        # 输入：DH（公式9输出）、μg（气体粘度，常数）、ρg（气体密度，公式11输出）、∂zP（压力轴向梯度，公式63输出）
        # 输出：vg（气体速度）
        # 公式：vg = (2 / 0.316) * (DH^5 * |dP/dz| / (mu_g * rho_g^3))^(1/7)
        term = (DH**5 * abs(dP_dz)) / (mu_g * rho_g**3)
        return (2 / 0.316) * (term ** (1/7))
    
    def formula_11(self, M_i, C_i_g):
        # 11气体密度计算
        # 输入：Mi（组分i摩尔质量，常数）、Ci,g（组分i气体组分浓度，变量）
        # 输出：ρg（气体密度）
        # 公式：ρi = Σ M_i * C_i
        if isinstance(M_i, dict) and isinstance(C_i_g, dict):
            total_density = 0.0
            for component in C_i_g:
                if component in M_i:
                    # 处理M_i是嵌套字典的情况
                    if isinstance(M_i[component], dict) and 'molar_mass' in M_i[component]:
                        molar_mass = M_i[component]['molar_mass']
                    else:
                        molar_mass = M_i[component]
                    total_density += molar_mass * C_i_g[component]
            # 将密度从g/m³转换为kg/m³，以匹配后续公式中的单位使用
            return total_density / 1000
        else:
            total_density = sum(M_i * C_i for M_i, C_i in zip(M_i, C_i_g))
            # 将密度从g/m³转换为kg/m³，以匹配后续公式中的单位使用
            return total_density / 1000
    
    def formula_12(self, a_omega, b_omega, omega):
        # 12休止角计算
        # 输入：aω（比例系数）、bω（截距常数）、ω（窑转速，控制量）
        # 输出：ξ（休止角）
        return a_omega * omega + b_omega
    
    def formula_13(self, omega, psi, xi, r_c, Lc, phi_z):
        # 13回转窑固体速度
        # 输入：ω（窑转速，控制量）、ψ（窑倾角）、ξ（公式10输出）、rc（窑内半径，参数）、Lc（弦长，公式14输出）、ϕ(z)（床层坡度角，公式73输出）
        # 输出：vs（固体速度）
        # 公式：vs = ω * (ψ + phi_z * cos(xi)) / sin(xi) * (2*rc / sin(Lc/(2*rc)))
        return omega * (psi + phi_z * math.cos(xi)) / math.sin(xi) * (2 * r_c / math.sin(Lc / (2 * r_c)))
    
    def formula_14(self, r_c, theta_z):
        # 14弦长
        # 输入：rc（参数）、θ(z)（填充角，公式18输出）
        # 输出：Lc（弦长）
        # 公式：Lc = 2 * r_c * math.sin(theta_z / 2)
        return 2 * r_c * math.sin(theta_z / 2)
    
    def formula_15(self, V_delta, C_i_g_t):
        # 15迭代过程气体组分摩尔数
        # 输入：VΔ（窑炉单段总体积）、Ci,g(t)（t时刻气体组分i的摩尔浓度）
        # 输出：ni,g(t)（t时刻气体组分i的摩尔数）
        return C_i_g_t * V_delta
    
    def formula_16(self, r_c, theta_z):
        # 16固体横截面积
        # 输入：rc（参数）、θ(z)（填充角，公式18输出）
        # 输出：As（固体横截面积）
        # 公式：A_s = r_c**2 / 2 * (theta_z - math.sin(theta_z))
        return r_c**2 / 2 * (theta_z - math.sin(theta_z))
    
    def formula_17(self, r_c):
        # 17窑总横截面积，参数
        # 输入：rc（参数）
        # 输出：At（窑总横截面积）
        # 公式：A_t = math.pi * r_c**2
        return math.pi * r_c**2
    
    def formula_18(self, theta):
        # 18填充率计算公式
        # 输入：θ（填充角，变量）
        # 输出：η（填充因子）
        # 公式：η = (θ - sinθ)/(2π)
        return (theta - math.sin(theta)) / (2 * math.pi)
    
    def _calculate_fill_angle(self, A_s, A_t, eta):
        # 辅助函数：使用牛顿迭代法求解填充角
        # 输入：As（固体横截面积）、At（窑总横截面积）、η（填充因子）
        # 输出：θ(z)（填充角）
        # 公式：η = A_s / A_t = (θ - sinθ)/(2π)
        
        # 1. 确定使用的填充率
        if eta is None or eta < 0:
            # 计算实际填充率
            if A_t == 0:
                return 0.0
            actual_eta = A_s / A_t
        else:
            # 使用传入的填充率
            actual_eta = eta
        
        # 2. 使用牛顿迭代法求解填充角 θ
        # 定义函数 f(θ) = (θ - sinθ)/(2π) - η = 0
        # 导数 f'(θ) = (1 - cosθ)/(2π)
        
        # 初始猜测值
        theta = math.pi  # 初始猜测为180度
        tolerance = 1e-6  # 收敛容差
        max_iter = 100  # 最大迭代次数
        
        for i in range(max_iter):
            # 计算f(theta)
            f_theta = (theta - math.sin(theta)) / (2 * math.pi) - actual_eta
            
            # 计算f'(theta)
            df_dtheta = (1 - math.cos(theta)) / (2 * math.pi)
            
            # 检查导数是否为0
            if abs(df_dtheta) < 1e-10:
                break
            
            # 牛顿迭代更新
            delta_theta = -f_theta / df_dtheta
            theta += delta_theta
            
            # 检查收敛
            if abs(delta_theta) < tolerance:
                break
        
        # 3. 返回计算得到的填充角
        return theta
    
    def formula_19(self, delta_z, V_s_t_plus_dt):
        # 19迭代后的固体横截面积公式
        # 输入：Δz（窑炉单段长度）、Vs(t+Δt)（迭代后的固体总体积，公式15输出）
        # 输出：As(t+Δt)（迭代后的固体横截面积）
        # 公式：A_s_t_plus_dt = V_s_t_plus_dt / delta_z
        return V_s_t_plus_dt / delta_z
    
    def formula_20(self, T):
        # 20气体热导率
        # 输入：T（气体温度）
        # 输出：kg（气体热导率）
        # 公式：kg = 0.026 * (T / 300)**0.8
        return 0.026 * (T / 300)**0.8
    
    def formula_21(self, vg, C_i_g, Tg, P, components, xj, cg, dC_i_g_dz, i):
        # 21气体组分i物质通量
        # 输入：vg（公式10输出）、Ci,g（组分i气体组分浓度，变量）、Tg（气体温度）、P（压力）、components（组分列表）、xj（各组分摩尔分数）、cg（气体总摩尔浓度）、∂zCi,g（组分i浓度梯度，公式23输入）、i（当前计算的组分索引）
        # 输出：Ni,g（气体组分i物质通量）
        # 公式：N_i,g = v_g * C_i,g - D_i,g * ∂_z C_i,g
        # 计算有效扩散系数
        D_i_g = self.formula_22(Tg, P, components, xj, cg, i)
        return vg * C_i_g - D_i_g * dC_i_g_dz
    
    def formula_22(self, T, P, components, xj, cg, i):
        # 22气体有效扩散系数
        # 输入：T（温度，变量）、P（压力，变量）、components（组分列表）、xj（各组分摩尔分数列表，公式61输出）、cg（气体总摩尔浓度，变量，公式77输出）、i（当前计算的组分索引）
        # 输出：Di,g（气体有效扩散系数）
        # 公式：D_i,g = (Σ_{j≠i} x_j / (c_g D_ij))^(-1)
        sum_term = 0.0
        
        # 计算二元扩散系数矩阵
        D_ij = []
        for comp_i in components:
            row = []
            for comp_j in components:
                # 获取组分i和j的摩尔质量和扩散体积
                M_i = self.constants['gas_properties'][comp_i]['molar_mass']
                M_j = self.constants['gas_properties'][comp_j]['molar_mass']
                V_i = self.constants['gas_properties'][comp_i]['diffusion_volume']
                V_j = self.constants['gas_properties'][comp_j]['diffusion_volume']
                # 计算二元扩散系数
                D = self.formula_24(T, P, M_i, M_j, V_i, V_j)
                row.append(D)
            D_ij.append(row)
        
        # 计算有效扩散系数，使用指定的组分索引i
        for j in range(len(xj)):
            if j != i and cg != 0 and D_ij[i][j] != 0:
                sum_term += xj[j] / (cg * D_ij[i][j])
        
        if sum_term == 0:
            return 0.0
        
        return 1.0 / sum_term
    
    def formula_23(self, C_i_g_current, C_i_g_prev, delta_z):
        # 23气体组分i浓度轴向梯度
        # 输入：Ci,g(z_k)（当前段气体组分i浓度，变量）、Ci,g(z_k-1)（前一段气体组分i浓度，变量）、Δz（段长度）
        # 输出：∂zCi,g（气体组分i浓度轴向梯度）
        # 公式：∂_z C_i,g(z_k) = (C_i,g(z_k) - C_i,g(z_k-1)) / Δz
        return (C_i_g_current - C_i_g_prev) / delta_z
    
    def formula_24(self, T, P, M_i, M_j, V_i, V_j):
        # 24二元扩散系数
        # 输入：T（温度，变量）、P（压力，变量）、Mi,Mj（组分i和j的摩尔质量，常数）、Vi,Vj（组分i和j的扩散体积，常数）
        # 输出：Dij（二元扩散系数）
        # 公式：D_ij = 0.00143 T^1.75 / (P M_ij^(1/2) [(V_i)^(1/3) + (V_j)^(1/3)]^2)
        # 计算平均摩尔质量
        M_ij = 2 / (1/M_i + 1/M_j)
        
        # 计算扩散系数
        numerator = 0.00143 * (T ** 1.75)
        denominator = P * (M_ij ** 0.5) * ((V_i ** (1/3) + V_j ** (1/3)) ** 2)
        
        if denominator == 0:
            return 0.0
        
        return numerator / denominator
    
    def formula_25(self, vs, C_i_s):
        # 25固体组分i物质通量
        # 输入：vs（公式13输出）、Ci,s（固体组分i浓度，变量）
        # 输出：Ni,s（固体组分i物质通量）
        # 公式：N_i,s = v_s * C_i,s
        return vs * C_i_s
    
    # 反应动力学与质量守恒方程
    def formula_26(self):
        # 26反应速率公式
        # 输出：反应速率公式字符串
        # 公式：r_j = k_r * T^n * e^(-E_A/RT) * ∏ P_i^β_i * ∏ C_i^α_i
        return "r_j = k_r * T^n * e^(-E_A/RT) * ∏ P_i^β_i * ∏ C_i^α_i"
    
    def _calculate_reaction_rate(self, reaction_id, T, P_i_dict, C_i_dict):
        # 辅助函数：计算反应速率
        # 输入：reaction_id（反应标识）、T（温度）、P_i_dict（组分分压字典）、C_i_dict（组分浓度字典）
        # 输出：rj（反应速率）
        
        # 获取反应速率系数
        reaction_coeff = self.constants['reaction_rate_coefficients'][reaction_id]
        kr = reaction_coeff['kr']
        n = reaction_coeff['n']
        EA = reaction_coeff['EA']
        R = self.constants['R']
        
        # 获取反应级数
        alpha1 = reaction_coeff['alpha1']
        alpha2 = reaction_coeff['alpha2']
        alpha3 = reaction_coeff['alpha3']
        beta2 = reaction_coeff['beta2']
        
        # 根据反应ID确定反应物顺序，按反应物从左到右的顺序匹配对应的指数
        reaction_reactants = {
            'r1': ['CaCO3'],  # CaCO3 → CaO + CO2
            'r2': ['CaO', 'SiO2'],  # 2CaO + SiO2 → C2S
            'r3': ['CaO', 'C2S'],  # CaO + C2S → C3S
            'r4': ['CaO', 'Al2O3'],  # 3CaO + Al2O3 → C3A
            'r5': ['CaO', 'Al2O3', 'Fe2O3'],  # 4CaO + Al2O3 + Fe2O3 → C4AF
            'r6': ['CO', 'O2'],  # 2CO + O2 → 2CO2
            'r7': ['CO', 'H2O'],  # CO + H2O → CO2 + H2
            'r8': ['H2', 'O2'],  # 2H2 + O2 → 2H2O
            'r9': ['C_sus', 'O2'],  # 2C + O2 → 2CO
            'r10': ['C_sus', 'H2O'],  # C + H2O → CO + H2
            'r11': ['C_sus', 'CO2']  # C + CO2 → 2CO
        }
        
        # 获取当前反应的反应物列表
        reactants = reaction_reactants.get(reaction_id, [])
        
        # 计算分压乘积项 ∏ P_i^β_i
        product_P = 1.0
        # beta2对应第二个反应物的分压指数（如果第二个反应物是气体）
        if len(reactants) > 1:
            second_reactant = reactants[1]
            if second_reactant in P_i_dict:
                product_P *= P_i_dict[second_reactant] ** beta2
        
        # 计算浓度乘积项 ∏ C_i^α_i
        product_C = 1.0
        # 按反应物从左到右的顺序匹配对应的alpha指数
        for i, reactant in enumerate(reactants):
            if reactant in C_i_dict:
                if i == 0:
                    product_C *= C_i_dict[reactant] ** alpha1
                elif i == 1:
                    product_C *= C_i_dict[reactant] ** alpha2
                elif i == 2:
                    product_C *= C_i_dict[reactant] ** alpha3
        
        # 计算反应速率，将活化能EA从kJ/mol转换为J/mol以匹配气体常数R的单位
        rj = kr * (T ** n) * math.exp(-(EA) / (R * T)) * product_P * product_C
        return rj
    
    def formula_27(self):
        # 27反应源项计算公式
        # 输出：反应源项公式字符串
        # 公式：[R_s; R_g] = v * r
        return "[R_s; R_g] = v * r"
    
    def _calculate_reaction_source(self, T, P_i_dict, C_i_dict):
        # 辅助函数：计算反应源项
        # 输入：T（温度）、P_i_dict（组分分压字典）、C_i_dict（组分浓度字典）
        # 输出：Rs,Rg（固/气相生成速率）
        
        # 获取化学计量矩阵
        stoichiometric_matrix = self.constants['stoichiometric_matrix']['matrix']
        
        # 获取所有反应标识
        reactions = list(self.constants['reaction_rate_coefficients'].keys())
        
        # 计算所有反应的速率
        r = []
        for reaction_id in reactions:
            rj = self._calculate_reaction_rate(reaction_id, T, P_i_dict, C_i_dict)
            r.append(rj)
        
        # 确定固体和气体组分的分界点（根据化学计量矩阵定义，前9个为固体，后6个为气体）
        solid_gas_split = 9
        
        # 计算固体相生成速率
        Rs = []
        for i in range(solid_gas_split):
            R = 0.0
            # 确保反应速率r的长度与化学计量矩阵的列数匹配
            for j in range(min(len(r), len(stoichiometric_matrix[i]))):
                R += stoichiometric_matrix[i][j] * r[j]
            Rs.append(R)
        
        # 计算气体相生成速率
        Rg = []
        for i in range(solid_gas_split, len(stoichiometric_matrix)):
            R = 0.0
            # 确保反应速率r的长度与化学计量矩阵的列数匹配
            for j in range(min(len(r), len(stoichiometric_matrix[i]))):
                R += stoichiometric_matrix[i][j] * r[j]
            Rg.append(R)
        
        return Rs, Rg
    
    def formula_28(self, dN_i_s_dz, R_s_i):
        # 28固体质量守恒
        # 输入：∂zNi,s（固体组分i通量梯度，公式29输出）、Rs,i（公式27输出）
        # 输出：∂tCi,s（固体组分i浓度时间变化率）
        # 公式：∂_t C_i,s = -∂_z N_i,s + R_s,i
        return -dN_i_s_dz + R_s_i
    
    def formula_29(self, N_i_s_current, N_i_s_prev, delta_z):
        # 29固体物质通量轴向梯度
        # 输入：Ni,s（公式25输出，空间离散点数据）通过相邻段Ni,s差值计算
        # 输出：∂zNi,s（固体组分i物质通量轴向梯度）
        # 公式：∂_z N_i,s(z_k) = (N_i,s(z_k) - N_i,s(z_k-1)) / Δz
        return (N_i_s_current - N_i_s_prev) / delta_z
    
    def formula_30(self, dN_i_g_dz, R_g_i):
        # 30气体质量守恒
        # 输入：∂zNi,g（气体组分i通量梯度，公式31输出）、Rg,i（公式27输出）
        # 输出：∂tCi,g（气体组分i浓度时间变化率）
        # 公式：∂_t C_i,g = -∂_z N_i,g + R_g,i
        return -dN_i_g_dz + R_g_i
    
    def formula_31(self, N_i_g_current, N_i_g_prev, delta_z):
        # 31气体物质通量轴向梯度
        # 输入：Ni,g（公式21输出，空间离散点数据），通过相邻段Ni,g差值计算
        # 输出：∂zNi,g（气体组分i物质通量轴向梯度）
        # 公式：∂_z N_i,g(z_k) = (N_i,g(z_k) - N_i,g(z_k-1)) / Δz
        return (N_i_g_current - N_i_g_prev) / delta_z
    
    # 能量守恒与热传递方程
    def formula_32(self, cp_g, mu_g, k_g):
        # 32分解炉普朗特数计算
        # 输入：cp,g（气体比热容，公式81输出）、μg（常数）、kg（气体热导率，公式20输出）
        # 输出：Pr（普朗特数）
        # 公式：Pr = (cp,g * mu_g) / k_g
        return (cp_g * mu_g) / k_g
    
    def formula_33(self, rho_g, vg, De, mu_g):
        # 33轴向雷诺数
        # 输入：ρg（公式11输出）、vg（公式10输出）、De（有效直径，公式34输出）、μg（常数）
        # 输出：ReD（轴向雷诺数）
        # 公式：Re_D = (rho_g * vg * De) / mu_g
        return (rho_g * vg * De) / mu_g
    
    def formula_34(self, r_c, theta):
        # 34有效直径
        # 输入：rc（参数）、θ（填充角，公式18输出）
        # 输出：De（有效直径）
        # 公式：De = 2 * r_c * (math.pi - theta/2 + math.sin(theta)/2) / (math.pi - theta/2 + math.sin(theta/2))
        numerator = math.pi - theta/2 + math.sin(theta)/2
        denominator = math.pi - theta/2 + math.sin(theta/2)
        return 2 * r_c * numerator / denominator
    
    def formula_35(self, k_g, dp, ReD, Pr):
        # 35分解炉对流换热量
        # 输入：kg（气体热导率，公式20输出）、dp（颗粒直径，常数）、ReD（轴向雷诺数，公式33输出）、Pr（公式32输出）
        # 输出：Qgscv（分解炉对流换热量）
        # 公式：Qgs^cv = (k_g / dp) * 0.3 * ReD^0.6 * Pr^0.33
        return (k_g / dp) * 0.3 * (ReD ** 0.6) * (Pr ** 0.33)
    
    def formula_36(self, rho_g, De, mu_g, omega):
        # 36旋转雷诺数
        # 输入：ρg（公式11输出）、De（有效直径，公式34输出）、μg（粘度，常数）、ω（控制量）
        # 输出：Reω（旋转雷诺数）
        # 公式：Re_omega = (rho_g * omega * De**2) / mu_g
        return (rho_g * omega * De**2) / mu_g
    
    def formula_37(self, ReD, Re_omega, eta):
        # 37努塞尔数
        # 输入：ReD（公式33输出）、Reω（公式36输出）、η（公式18输出）
        # 输出：Nu（努塞尔数）
        # 公式：Nu = 0.46 * ReD^0.535 * Re_omega^0.104 * eta^(-0.341)
        return 0.46 * (ReD ** 0.535) * (Re_omega ** 0.104) * (eta ** (-0.341))
    
    def formula_38(self, k_g, De, Nu):
        # 38回转窑对流换热系数
        # 输入：kg（气体热导率，公式20输出）、De（公式34输出）、Nu（公式37输出）
        # 输出：β（对流换热系数）
        # 公式：β = (k_g / De) * Nu
        return (k_g / De) * Nu
    
    def formula_39(self, Ags, beta, Tg, Ts):
        # 39回转窑对流换热量
        # 输入：Ags（气固换热面积，公式40输出）、β（回转窑对流换热系数，公式38输出）、Tg,Ts（气\固温度，变量）
        # 输出：Qgscv（回转窑对流换热量）
        # 公式：Qgs^cv = Ags * beta * (Tg - Ts)
        return Ags * beta * (Tg - Ts)
    
    def formula_40(self, Lc, delta_z):
        # 40气固换热面积
        # 输入：Lc（公式14输出）、Δz（段长度，参数）
        # 输出：Ags（气固换热面积）
        # 公式：Ags(k) ≈ Lc(z_k)·Δz = [2r_c sin(θ_k/2)]Δz
        # 由于Lc已经是2r_c sin(θ_k/2)，所以直接相乘
        return Lc * delta_z
    
    def formula_41(self, sigma, Ags, epsilon_gs, Tg, Ts):
        # 41辐射换热量
        # 输入：σ（斯特藩-玻尔兹曼常数，5.67×10−8）、Ags（公式40输出）、ϵgs（气固发射率，变量，公式80输出）、Tg,Ts（气温度，变量）
        # 输出：Qgsrad（辐射换热量）
        # 公式：Qgs^rad = sigma * Ags * epsilon_gs * (Tg^4 - Ts^4)
        return sigma * Ags * epsilon_gs * (Tg**4 - Ts**4)
    
    def formula_42(self, H_tilde_s_zk, H_tilde_g_zk, Q_tilde_g_zk, H_tilde_s_zk_prev, H_tilde_g_zk_prev, Q_tilde_g_zk_prev, delta_z, Qgscv, Qgsrad, V_delta, Jsg=0.0):
        # 42分解炉能量平衡
        # 输入：H~s(zk), H~g(zk), Q~g(zk)（当前段焓通量密度和热传导）
        #       H~s(zk-1), H~g(zk-1), Q~g(zk-1)（前一段焓通量密度和热传导）
        #       Δz（段长度）、Qgscv,Qgsrad（对流辐射换热）、VΔ（单段总体积）、Jsg（反应焓变）
        # 输出：∂tU^c（分解炉内能密度变化率）
        # 公式：∂_t Û_c = -∂_z (H~_g + H~_s + Q~_g) - (Qgs^rad + Qgs^cv) / VΔ + Jsg / VΔ
        
        # 计算当前段的总焓通量和热传导
        total_current = H_tilde_g_zk + H_tilde_s_zk + Q_tilde_g_zk
        # 计算前一段的总焓通量和热传导
        total_prev = H_tilde_g_zk_prev + H_tilde_s_zk_prev + Q_tilde_g_zk_prev
        # 计算梯度 ∂_z (H~_g + H~_s + Q~_g)
        gradient = (total_current - total_prev) / delta_z
        
        return -gradient - (Qgsrad + Qgscv) / V_delta + Jsg 
    
    def formula_43(self, T_i):
        # 43温度轴向梯度
        # 输入：Ti（气体或固体温度，空间离散点数据，包含当前段和前一段的温度）
        # 输出：∂Ti/∂z（梯度）
        # 公式：∂Ti/∂z (z_k) = (Ti(z_k) - Ti(z_k-1)) / Δz
        # 这里假设T_i是一个包含当前段和前一段温度的列表
        if len(T_i) < 2:
            return 0.0
        return (T_i[-1] - T_i[-2]) / self.dz
    
    def formula_44(self, H_tilde_s_zk, Q_tilde_s_zk, H_tilde_s_zk_prev, Q_tilde_s_zk_prev, delta_z, Qgsrad, Qgscv, V_delta, Jsg):
        # 44回转窑固体能量平衡
        # 输入：H~s(zk), Q~s(zk)（当前段固体焓通量密度和热传导）
        #       H~s(zk-1), Q~s(zk-1)（前一段固体焓通量密度和热传导）
        #       Δz（段长度）、Qgsrad,Qgscv（对流辐射换热）、VΔ（单段总体积）、Jsg（反应焓变）
        # 输出：∂tU^s（固体内能密度变化率）
        # 公式：∂_t Û_s = -∂_z (H~_s + Q~_s) + (Qgs^rad + Qgs^cv) / VΔ - J_sg
        
        # 计算当前段的总焓通量和热传导
        total_current = H_tilde_s_zk + Q_tilde_s_zk
        # 计算前一段的总焓通量和热传导
        total_prev = H_tilde_s_zk_prev + Q_tilde_s_zk_prev
        # 计算梯度 ∂_z (H~_s + Q~_s)
        gradient = (total_current - total_prev) / delta_z
        
        return -gradient + (Qgsrad + Qgscv) / V_delta + Jsg
    
    def formula_45(self, H_tilde_g_zk, Q_tilde_g_zk, H_tilde_g_zk_prev, Q_tilde_g_zk_prev, delta_z, Qgsrad, Qgscv, V_delta, Jsg):
        # 45回转窑气体能量平衡
        # 输入：H~g(zk), Q~g(zk)（当前段气体焓通量密度和热传导）
        #       H~g(zk-1), Q~g(zk-1)（前一段气体焓通量密度和热传导）
        #       Δz（段长度）、Qgsrad,Qgscv（对流辐射换热）、VΔ（单段总体积）、Jsg（反应焓变）
        # 输出：∂tU^g（气体内能密度变化率）
        # 公式：∂_t Û_g = -∂_z (H~_g + Q~_g) - (Qgs^rad + Qgs^cv) / VΔ + J_sg
        
        # 计算当前段的总焓通量和热传导
        total_current = H_tilde_g_zk + Q_tilde_g_zk
        # 计算前一段的总焓通量和热传导
        total_prev = H_tilde_g_zk_prev + Q_tilde_g_zk_prev
        # 计算梯度 ∂_z (H~_g + Q~_g)
        gradient = (total_current - total_prev) / delta_z
        
        return -gradient + (Qgsrad + Qgscv) / V_delta + Jsg
    
    def formula_46(self, cj_i, T):
        # 46权重系数
        # 输入：cj,i（温度相关系数，变量）、T（气体温度，变量）
        # 输出：aj（权重系数，j=0∼4）
        # 公式：a0 = 1 - Σa_j(j=1到4), aj = Σcj,i(T/Tref)^(i-1)(i=1到3)
        # 从constants获取Tref
        Tref = self.constants['Tref']
        
        a = []
        # cj_i的结构是 [3 rows][4 columns]，对应 i=1-3, j=1-4
        # 所以我们需要转置，j从0到3（对应j=1到4）
        for j in range(4):
            aj = 0.0
            for i in range(3):
                # 确保索引不越界
                if i < len(cj_i) and j < len(cj_i[i]):
                    aj += cj_i[i][j] * (T / Tref) ** (i)
            a.append(aj)
        
        # 计算a0
        a0 = 1 - sum(a)
        a.insert(0, a0)
        
        return a
    
    def formula_47(self, xH2O, xCO2):
        # 47吸收系数
        # 输入：xH2O,xCO2（摩尔分数，公式61输出）
        # 输出：kj（吸收系数，j=1∼4）
        # 公式：kj = K1,j + K2,j * (xH2O / xCO2)
        # 从constants获取wsgg_coefficients
        wsgg_coeffs = self.constants['wsgg_coefficients']
        K1_j = wsgg_coeffs['K1']
        K2_j = wsgg_coeffs['K2']
        
        k = []
        for K1, K2 in zip(K1_j, K2_j):
            if xCO2 != 0:
                kj = K1 + K2 * (xH2O / xCO2)
            else:
                kj = K1
            k.append(kj)
        return k
    
    def formula_48(self, xH2O, xCO2):
        # 48温度相关系数cj,i
        # 输入：xH2O,xCO2（摩尔分数，公式61输出）
        # 输出：cj,i（温度相关系数）
        # 公式：cj,i = C1,j,i + C2,j,i*(xH2O/xCO2) + C3,j,i*(xH2O/xCO2)^2
        # 从constants获取wsgg_coefficients
        wsgg_coeffs = self.constants['wsgg_coefficients']
        
        # 计算xH2O/xCO2
        if xCO2 != 0:
            xH2O_xCO2 = xH2O / xCO2
        else:
            xH2O_xCO2 = 0.0
        
        # 获取C1, C2, C3系数
        C1_j_i = wsgg_coeffs.get('C1', [[0.0]*3 for _ in range(4)])
        C2_j_i = wsgg_coeffs.get('C2', [[0.0]*3 for _ in range(4)])
        C3_j_i = wsgg_coeffs.get('C3', [[0.0]*3 for _ in range(4)])
        
        cj_i = []
        for C1, C2, C3 in zip(C1_j_i, C2_j_i, C3_j_i):
            cj = []
            for c1, c2, c3 in zip(C1, C2, C3):
                cj_val = c1 + c2 * xH2O_xCO2 + c3 * (xH2O_xCO2 ** 2)
                cj.append(cj_val)
            cj_i.append(cj)
        return cj_i
    
    def formula_49(self, r_c):
        # 49平均射线长度
        # 输入：rc（参数）
        # 输出：Sm（平均射线长度）
        # 公式：Sm = 0.95 × 2rc
        return 0.95 * 2 * r_c
    
    def formula_50(self, Hf_i, integral_enthalpy):
        # 50摩尔焓
        # 输入：Hf,i（组分i标准生成焓，常数）、积分焓（组分i的积分焓，公式5输出）
        # 输出：hi（组分i摩尔焓）
        # 公式：hi = ΔHf,i + ∫(T0到T) cp,i(τ)dτ
        return Hf_i + integral_enthalpy
    
    def formula_51(self, Q_tilde_i_zk, Q_tilde_i_zk_minus_1, delta_z):
        # 51热传导轴向的变化率
        # 输入：Q~i(zk)（当前热传导）、Q~i(zk−1)（前一段热传导）、Δz（段长度）
        # 输出：∂zQ（热传导的梯度）
        # 公式：∂_z Q~i(z_k) = (Q~i(z_k) - Q~i(z_k-1)) / Δz
        return (Q_tilde_i_zk - Q_tilde_i_zk_minus_1) / delta_z
    
    def formula_52(self, reaction_rates, reaction_enthalpies):
        # 52反应焓变计算（支持多个反应）
        # 输入：reaction_rates（反应速率字典，键为反应ID，值为反应速率）、reaction_enthalpies（反应焓变字典，键为反应ID，值为反应焓变）
        # 输出：Jsg（总反应焓变）
        # 公式：J_sg = Σ (r_j × ΔH_rj(T))
        total_enthalpy = 0.0
        for reaction_id, rate in reaction_rates.items():
            if reaction_id in reaction_enthalpies:
                total_enthalpy += rate * reaction_enthalpies[reaction_id]
        return total_enthalpy
    
    def formula_53(self, T, P, xH2O, xCO2, r_c):
        # 53气体发射率
        # 输入：T（温度）、P（压力）、xH2O,xCO2（摩尔分数）、r_c（参数）
        # 输出：ϵg（气体发射率）
        # 公式：ϵg = Σ(a_j(1 - e^(-kj*Sm*P*(xH2O+xCO2))))
        
        # 1. 计算平均射线长度（公式49）
        S_m = self.formula_49(r_c)
        
        # 2. 计算温度相关系数cj,i（公式48）
        cj_i = self.formula_48(xH2O, xCO2)
        
        # 3. 计算权重系数aj（公式46）
        a_j = self.formula_46(cj_i, T)
        
        # 4. 计算吸收系数kj（公式47）
        k_j = self.formula_47(xH2O, xCO2)
        
        # 5. 计算气体发射率（累加计算，j从0到4）
        epsilon_g = 0.0
        for j in range(5):  # j=0到4
            # 获取当前j对应的a_j和k_j
            a = a_j[j] if j < len(a_j) else 0.0
            k = k_j[j-1] if j > 0 and j-1 < len(k_j) else (k_j[0] if j == 0 and len(k_j) > 0 else 0.0)
            
            # 计算指数项
            exponent = -k * S_m * P * (xH2O + xCO2)
            # 限制指数值范围，避免溢出
            exponent = max(exponent, -700)  # math.exp(-700) 已经非常接近0
            exponent = min(exponent, 700)   # math.exp(700) 会溢出，但这里应该是负数
        
            
            # 累加计算气体发射率
            epsilon_g += a * (1 - math.exp(exponent))
        
        return epsilon_g
    
    def formula_54(self, k_g, dT_i_dz):
        # 54气体热传导
        # 输入：kg（气体热导率，公式20输出）、∂Ti/∂z（温度梯度，公式43输出）
        # 输出：Q~g（气体热传导）
        # 公式：Q~_g = -k_g ∂T_g/∂z
        return -k_g * dT_i_dz
    
    def formula_55(self, dT_i_dz):
        # 55固体热传导
        # 输入：∂Ti/∂z（温度梯度，公式43输出）
        # 输出：Q~s（固体热传导）
        # 公式：Q~_s = -k_s ∂T_s/∂z
        k_s = self.parameters['k_s']  # 固体热导率，W/(K·m)
        return -k_s * dT_i_dz
    
    def formula_56(self, m_dot_i_s_0, M_i):
        # 56入口边界初始化固体组分摩尔数
        # 输入：m˙i,s,0（入口固体组分i的质量流量）、Mi（固体组分i的摩尔质量）
        # 输出：ni,s,0（初始时刻（t=0）固体组分i的摩尔数）
        return m_dot_i_s_0  / M_i
    
    def formula_57(self, V_delta, C_i_s_t):
        # 57迭代过程固体组分摩尔数
        # 输入：VΔ（窑炉单段总体积，公式58输出）、Ci,s(t)（t时刻固体组分i的摩尔浓度）
        # 输出：ni,s(t)（t时刻固体组分i的摩尔数）
        # 公式：n_i,s(t) = C_i,s(t)·VΔ
        return C_i_s_t * V_delta
    
    def formula_57_gas(self, V_delta, C_i_g_t):
        # 57-气体迭代过程气体组分摩尔数
        # 输入：VΔ（窑炉单段总体积，公式58输出）、Ci,g(t)（t时刻气体组分i的摩尔浓度）
        # 输出：ni,g(t)（t时刻气体组分i的摩尔数）
        # 公式：n_i,g(t) = C_i,g(t)·VΔ
        return C_i_g_t * V_delta
    
    def formula_58(self, A_t, delta_z):
        # 58单段总体积
        # 输入：At（窑总横截面积，公式17输出）、Δz（窑炉单段长度）
        # 输出：VΔ（窑炉单段总体积）
        # 公式：VΔ = A_t·Δz
        return A_t * delta_z
    
    def formula_59(self, V_delta, n_i_s_0):
        # 59初始化固体组分i浓度
        # 输入：VΔ（窑炉单段总体积，公式58输出）、ni,s,0（初始时刻固体组分i的摩尔数，公式56输出）
        # 输出：Ci,s,0（初始固体组分i的摩尔浓度）
        # 公式：C_i,s,0 = n_i,s,0 / VΔ
        return n_i_s_0 / V_delta
    
    def formula_60(self, x_i_g_0, P_in, T, R):
        # 60初始化气体组分i浓度
        # 输入：xi,g,0（入口气体组分i的摩尔分数）、Pin（入口边界压力）、T（入口风温）、R（气体常数）
        # 输出：Ci,g,0（初始气体组分i的摩尔浓度）
        # 公式：C_i,g,0 = x_i,g,0 · P_in / (R·T)
        return x_i_g_0 * P_in / (R * T)
    
    def formula_61(self, C_i_g_t, cg_t):
        # 61迭代过程中更新摩尔分数
        # 输入：Ci,g(t)（迭代时刻气体组分i的浓度）、cg(t)（迭代时刻气体总浓度）
        # 输出：xi,g(t)（迭代时刻气体组分i的摩尔分数）
        # 公式：xi,g(t) = ni / ng = Ci,g(t) / cg(t)
        if cg_t == 0:
            return 0.0
        return C_i_g_t / cg_t
    
    def formula_62(self, delta_z, P_zk, P_zk_minus_1):
        # 62压力梯度
        # 输入：Δz（窑炉单段长度）、P(zk)（当前段k的压力）、P(zk−1)（前一段k−1的压力）
        # 输出：∂zP(zk)（当前段k的轴向压力梯度）
        # 公式：dP/dz = (P_zk - P_zk_minus_1) / delta_z
        return (P_zk - P_zk_minus_1) / delta_z
    
    def formula_63(self, H_s, V_delta):
        # 63固体单位体积焓
        # 输入：Hs（固体总焓，公式4输出）、VΔ（单段总体积，参数）
        # 输出：H^s（固体单位体积焓）
        return H_s / V_delta
    
    def formula_64(self, H_g, V_delta):
        # 64气体单位体积焓
        # 输入：Hg（气体总焓，公式4输出）、VΔ（单段总体积，参数）
        # 输出：H^g（气体单位体积焓）
        return H_g / V_delta
    
    def formula_65(self, P, H_hat_g, H_hat_s, V_g):
        # 65分解炉混合内能密度定义
        # 输入：P（压力，变量）、H^g（气体单位体积焓，公式64输出）、H^s（固体单位体积焓，公式63输出）、V^g（公式67输出）
        # 输出：U^c（内能密度）
        return H_hat_s + H_hat_g * V_g - P * V_g
    
    def formula_66(self, V_s, V_delta):
        # 66固体体积分数
        # 输入：Vs（固体体积，公式8输出）、VΔ（单段总体积，参数）
        # 输出：V^s（固体相体积分数）
        return V_s / V_delta
    
    def formula_67(self, V_g, V_delta):
        # 67气体体积分数
        # 输入：Vg（气体体积，公式7输出）、VΔ（单段总体积，参数）
        # 输出：V^g（气体相体积分数）
        return V_g / V_delta
    
    def formula_68(self, h_i_s, N_i_s):
        # 68固体焓通量密度
        # 输入：hi,s（固体组分i摩尔焓，公式50输出）、Ni,s（固体组分i通量，公式25输出）
        # 输出：H~s（固体焓通量密度）
        # 公式：H~_s = Σ (h_i,s · N_i,s)
        H_tilde_s = 0.0
        for h, N in zip(h_i_s, N_i_s):
            H_tilde_s += h * N
        return H_tilde_s
    
    def formula_69(self, h_g, N_g):
        # 69气体焓通量密度
        # 输入：hg（气体组分i摩尔焓，公式50输出）、Ng（气体组分i通量，公式21输出）
        # 输出：H~g（气体焓通量密度）
        # 公式：H~_g = Σ (h_i,g · N_i,g)
        H_tilde_g = 0.0
        for h, N in zip(h_g, N_g):
            H_tilde_g += h * N
        return H_tilde_g
    
    def formula_70(self, H, delta_z):
        # 70焓通量密度梯度
        # 输入：H（气体或固体焓通量密度，公式68、69输出）、Δz（窑炉单段长度）
        # 输出：∂zH（气体或固体焓通量密度梯度）
        # 公式：∂_z H(z_k) = (H(z_k) - H(z_k-1)) / Δz
        # 这里假设H是一个包含当前段和前一段焓通量密度的列表
        if len(H) < 2:
            return 0.0
        return (H[-1] - H[-2]) / delta_z
    
    def formula_71(self, vg):
        # 71分解炉固体轴向速度
        # 输入：vg气体轴向速度（公式10输出）
        # 输出：vs（分解炉固体轴向速度，等于气体速度）
        # 公式：vs = vg
        return vg
    
    def formula_72(self, r_c, theta_z):
        # 72床层高度
        # 输入：rc（参数）、θ(z)（填充角，公式18输出）
        # 输出：h(z)（床层高度）
        # 公式：h(z) = r_c * (1 - math.cos(theta_z / 2))
        return r_c * (1 - math.cos(theta_z / 2))
    
    def formula_73(self, dh_z_dz):
        # 73床层坡度角
        # 输入：∂h(z)/∂z（床层高度轴向梯度，公式74输出）
        # 输出：ϕ(z)（床层坡度角）
        # 公式：phi(z) = math.atan(dh_z_dz)
        return math.atan(dh_z_dz)
    
    def formula_74(self, h_z):
        # 74床层高度轴向梯度
        # 输入：h(z)（公式72输出，空间离散点数据，包含当前段和前一段的床层高度）
        # 输出：∂h(z)/∂z（床层高度轴向梯度）
        # 公式：∂h(z)/∂z = (h(z_k) - h(z_{k-1})) / Δz
        if len(h_z) < 2:
            return 0.0
        return (h_z[-1] - h_z[-2]) / self.dz
    
    def formula_75(self, x_i_g_0, P_in, T, R, V_delta):
        # 75入口边界初始化气体组分摩尔数
        # 输入：xi,g,0（入口气体组分i的摩尔分数）、Pin（入口边界压力）、T（入口风温）、R（气体常数）、VΔ（窑炉单段总体积）
        # 输出：ni,g(0)（初始时刻气体组分i的摩尔数）
        return (x_i_g_0 * P_in * V_delta) / (R * T)
    
    def formula_76(self, M_i, M_j):
        # 76平均摩尔质量
        # 输入：Mi（组分i摩尔质量，常数）、Mj（组分j摩尔质量，常数）
        # 输出：Mij（平均摩尔质量）
        # 公式：M_ij = 2 / (1/M_i + 1/M_j)
        return 2 / (1/M_i + 1/M_j)
    
    def formula_77(self, C_i_g):
        # 77气体总摩尔浓度
        # 输入：Ci,g（气体i组分摩尔浓度）
        # 输出；cg（气体总摩尔浓度）
        # 公式：c_g = Σ C_j,g
        if isinstance(C_i_g, dict):
            return sum(C_i_g.values())
        else:
            return sum(C_i_g)
    
    def formula_78(self, A_t, A_s):
        # 78气体横截面积 #删去
        # 输入：At（窑总横截面积，公式17输出）、As（固体横截面积，公式16输出）
        # 输出：Ag（气体横截面积）
        # 公式：A_g = A_t - A_s
        return A_t - A_s
    
    def formula_79(self, A_s, A_t):
        # 79气体截面积
        # 输入：As（固体横截面积）、At（窑总横截面积）
        # 输出：Ag（气体横截面积）
        # 公式：A_g = A_t - A_s
        return A_t - A_s
    
    def formula_80(self, epsilon_g, epsilon_s):
        # 80气固发射率
        # 输入：ϵg（气体发射率，公式53输出）、ϵs（固体发射率，参数）
        # 输出：ϵgs（气固发射率）
        # 公式：ϵgs = ϵg + ϵs - ϵgϵs
        return epsilon_g + epsilon_s - epsilon_g * epsilon_s
    
    def formula_81(self, cp_i_list, M_i_list, n_i_list):
        # 81气体比热容
        # 输入：cp,i_list（各气体组分摩尔热容列表，公式 6 输出）、M_i_list（各气体组分摩尔质量列表，常数）、n_i_list（各气体组分摩尔数列表，公式57_gas输出）
        # 输出： cp,g（气体比热容）
        # 公式：cp,g = Σ(ni*cp,i) / Σ(ni*Mi)
        numerator = 0.0
        denominator = 0.0
        for cp_i, M_i, n_i in zip(cp_i_list, M_i_list, n_i_list):
            numerator += n_i * cp_i
            denominator += n_i * M_i
        
        if denominator == 0:
            return 0.0
        
        return numerator / denominator
    
    def solve(self):
        # 微分 - 代数同步求解核心循环
        # 外层时间迭代主循环，循环条件：未收敛且未达到最大迭代步数
        iteration = 0 #初始化迭代计数器，记录当前迭代次数，从 0 开始
        while not self.converged and iteration < self.max_iterations:
            # 存储前一时刻核心变量用于收敛判断
            prev_core_vars = self._get_core_variables()
            
            # 【关键修复】每个时间步开始时，先进行物料传递，确保浓度能够沿管道流动
            # 从预热器→分解炉→回转窑，逐个单元传递浓度和热量
            for i in range(len(self.cells) - 1):
                # 所有相邻单元都进行物料转移（这样确保浓度沿着系统流动）
                self._transfer_between_cells(self.cells[i], self.cells[i+1])
            
            # 内层单元遍历循环：沿物料流动方向依次遍历所有有限体积单元
            for i, cell in enumerate(self.cells):
                # ① 识别单元所属设备，根据单元索引获取单元所属的工艺区段
                section = self._get_cell_section(i)
                
                # 预热器单元不需要计算
                if section == 'preheater':
                    continue
                
                # 获取设备配置
                config = self.section_configs[section].copy()
                
                # 传入设备类型，以便获取设备特定的执行规则
                #存储当前单元的工艺类型
                config['device_type'] = section
                
                # ② 读取匹配：对应设备的「差异化配置」+ 全局「通用公共规则」
                #整合通用规则和当前设备的差异化规则，形成完整执行清单
                execution_rules = self._get_execution_rules(config)
                
                # ③ 代数计算：按配置调用对应公式，计算所有代数变量
                algebraic_vars = self._perform_algebraic_calculations(cell, execution_rules)
                
                # ④ 微分求解：按配置调用对应公式，代入代数变量，求解浓度变化率和内能变化率
                dC_dt, dU_dt = self._solve_differential_equations(cell, execution_rules, algebraic_vars)
                
                # ⑤ 同步更新：新值 = 原值 + 变化率 × Δt
                # 同一个时间步 Δt 内，同一个单元里，代数计算和微分计算同时做，无先后顺序，结果互相代入
                self._update_variables(cell, dC_dt, dU_dt)

            
            # 单元边界衔接：确保三段工艺无缝衔接
            #将前一段工艺最后一个单元的变量传递给后一段工艺的第一个单元
            self._handle_section_transitions()
            
            # 收敛判定：计算全局核心微分变量的迭代差值
            current_core_vars = self._get_core_variables() #获取当前迭代后的核心变量
            # 至少进行10次迭代后才开始检查收敛，确保反应充分进行
            if iteration >= 10 and self.check_convergence(prev_core_vars, current_core_vars): #对比迭代前后的核心变量，判断差值是否均小于收敛阈值（1e-6）
                self.converged = True
                print(f"模型在 {iteration} 次迭代后达到稳态")
                break
            
            iteration += 1
        
        # 结果输出
        if self.converged:
            self.output_results()
        else:
            print("达到最大迭代步数，模型未收敛")
    
    def _get_cell_section(self, cell_index):
        # 根据单元索引判断所属工艺区段
        # 前一区段最后一个单元无缝衔接后一区段第一个单元
        return self.variables['device_type'][cell_index] #根据索引获取该单元的工艺类型
    
    def _get_execution_rules(self, config):
        # 整合通用公共规则和当前设备的差异化配置，生成单元的完整执行规则清单
        # 通用公式统一计算，差异化公式根据设备分别计算
        # 差异化部分直接使用传入的config（来自section_configs），不再使用硬编码的device_specific_rules
        # 确保使用差异化配置中的能量方程、反应速率类型、对流换热、固相速度、内能密度约束求解
        
        #定义通用公共规则字典，键为规则名称，值为对应的公式标识
        # 通用公共规则（所有设备都需要计算的公式）
        common_rules = {
            'mass_balance': 'formula_28',  # 质量守恒公式
            'reaction_source': 'formula_27',  # 反应源项公式
            'species_flux': 'formula_21',  # 物种通量公式
            'enthalpy_flux': 'formula_68',  # 焓通量公式
            'gas_temperature': 'formula_53',  # 气体温度相关公式
            'solid_temperature': 'formula_50',  # 固体温度相关公式
            'pressure_calculation': 'formula_62',  # 压力计算相关公式
            'velocity_calculation': 'formula_10',  # 速度计算相关公式
        }
        
        # 直接合并通用规则和差异化配置，差异化配置优先级更高
        # 这样可以确保使用section_configs中的能量方程、反应速率类型、对流换热、固相速度、内能密度约束求解
        execution_rules = common_rules.copy() #复制通用规则，作为基础规则
        execution_rules.update(config) #用当前设备的差异化配置 覆盖通用规则
        
        return execution_rules #返回整合后的完整执行规则清单，供后续计算调用
    
    def _perform_algebraic_calculations(self, cell, rules):
        # 执行代数计算
        # 通用公式统一计算，差异化公式根据设备分别计算
        algebraic_vars = {} #存储当前单元计算出的所有代数变量
        
        index = cell['index'] #获取当前单元的索引
        device_type = cell['device_type']
      
#一、基础几何与填充参数计算  
        # 1. 计算填充率和填充角（只计算一次，供后续所有计算使用）
        # 使用实际设备半径
        r_c = self.parameters['equipment'][device_type]['radius']
        
        # 计算窑总横截面积（公式17）
        A_t = self.formula_17(r_c)
        
        # 初始填充率设置为0.13
        eta = 0.13
        
        # 使用迭代方法计算填充角
        # 1. 初始计算固体横截面积（使用初始填充率）
        # 计算填充角（使用辅助函数）
        theta = self._calculate_fill_angle(0, A_t, eta)
        
        # 2. 使用计算得到的填充角计算固体横截面积（公式16）
        A_s = self.formula_16(r_c, theta)
        
        # 3. 迭代更新填充率和填充角，直到收敛
        max_iter = 10
        tolerance = 1e-6
        for _ in range(max_iter):
            # 更新填充率（公式：η = A_s / A_t）
            new_eta = A_s / A_t
            
            # 检查收敛
            if abs(new_eta - eta) < tolerance:
                break
            
            # 更新填充率
            eta = new_eta
            
            # 重新计算填充角
            theta = self._calculate_fill_angle(A_s, A_t, eta)
            
            # 重新计算固体横截面积
            A_s = self.formula_16(r_c, theta)
        
        # 4. 使用填充率计算公式验证（公式18）
        formula_eta = self.formula_18(theta)
        # 确保填充率一致性
        if abs(formula_eta - eta) > tolerance:
            # 如果不一致，以公式计算结果为准，重新计算填充角
            theta = self._calculate_fill_angle(0, A_t, formula_eta)
            A_s = self.formula_16(r_c, theta)
            eta = formula_eta
        
        # 将填充角和填充率存储在代数变量中，供后续计算使用
        algebraic_vars['theta'] = theta
        algebraic_vars['eta'] = eta
        algebraic_vars['A_t'] = A_t
        algebraic_vars['A_s'] = A_s
        algebraic_vars['r_c'] = r_c
     
#二、初始温度与压力设置   
        # 获取当前时刻的状态变量
        state_vars = cell['state_variables'] #获取单元的状态变量字典
        C = state_vars['C']  # 组分浓度向量
        U_hat = state_vars['U_hat']  # 单位体积内能向量
        
        # 初始温度和压力（如果是第一次迭代，使用默认值）
        # 判断单元的代数变量中是否已存在气固温度;若存在，直接读取上一次的温度值
        if 'Tg' in cell['algebraic_variables'] and 'Ts' in cell['algebraic_variables']:
            Tg = cell['algebraic_variables']['Tg']
            Ts = cell['algebraic_variables']['Ts']
        else:
            # 第一次迭代根据设备类型设置合理的初始温度
            if device_type == 'preheater':
                Tg = 390.0  # 预热器初始气体温度
                Ts = 390.0  # 预热器初始固体温度，对应固体进料温度
            elif device_type == 'calciner':
                Tg = 800.0  # 分解炉初始气体温度
                Ts = 700.0  # 分解炉初始固体温度
            else:  # kiln
                Tg = 1000.0  # 回转窑初始气体温度
                Ts = 900.0   # 回转窑初始固体温度
        
        # 限制温度范围，避免数值溢出
        Tg = max(Tg, 300.0)  # 最低温度300K
        Tg = min(Tg, 2000.0)  # 最高温度2000K
        Ts = max(Ts, 300.0)  # 最低温度300K
        Ts = min(Ts, 2000.0)  # 最高温度2000K
        
        if 'P' in cell['algebraic_variables']:
            P = cell['algebraic_variables']['P']
        else:
            P = self.constants['P0']  # 默认压力
  
#三、组分摩尔焓计算      
        # 获取所有组分列表
        stoichiometric_matrix = self.constants['stoichiometric_matrix'] #从化学计量矩阵中获取组分信息
        all_components = stoichiometric_matrix['components'] #提取所有组分名称
        solid_components = all_components[:9]
        gas_components = all_components[9:]
        
        # 1. 计算各组分的摩尔焓
        h_i_s = []  # 固体组分i摩尔焓列表
        h_i_g = []  # 气体组分i摩尔焓列表
        
        T0 = self.constants['T0']
        
        for component in solid_components: #遍历所有固体组分，逐个计算摩尔焓
            # 获取组分的热容系数
            if component in self.constants['molar_heat_capacity']: #判断该组分是否有热容系数数据
                coeffs = self.constants['molar_heat_capacity'][component]
            else:
                coeffs = {'C0': 0.0, 'C1': 0.0, 'C2': 0.0}
            
            C0 = coeffs['C0']
            C1 = coeffs['C1']
            C2 = coeffs['C2']
            
            # 计算积分焓（公式5），固体使用固体温度Ts
            integral_enthalpy = self.formula_5(C0, C1, C2, T0, Ts)
            
            # 获取标准生成焓
            Hf_i = self.constants['standard_enthalpy'].get(component, 0.0)
            
            # 计算摩尔焓（公式50）
            h_i = self.formula_50(Hf_i, integral_enthalpy)
            h_i_s.append(h_i) #将固体组分的摩尔焓添加到h_i_s列表
        
        for component in gas_components:
            # 获取组分的热容系数
            if component in self.constants['molar_heat_capacity']:
                coeffs = self.constants['molar_heat_capacity'][component]
            else:
                coeffs = {'C0': 0.0, 'C1': 0.0, 'C2': 0.0}
            
            C0 = coeffs['C0']
            C1 = coeffs['C1']
            C2 = coeffs['C2']
            
            # 计算积分焓（公式5），气体使用气体温度Tg
            integral_enthalpy = self.formula_5(C0, C1, C2, T0, Tg)
            
            # 获取标准生成焓
            Hf_i = self.constants['standard_enthalpy'].get(component, 0.0)
            
            # 计算摩尔焓（公式50）
            h_i = self.formula_50(Hf_i, integral_enthalpy)
            h_i_g.append(h_i)  # 将气体组分的摩尔焓添加到h_i_g列表

#四、物质通量计算        
        # 4. 计算各组分的通量
        # 获取速度，使用差异化配置中的固相速度公式
        
        # 计算压力梯度（公式62）
        delta_z = cell['delta_z']
        if index > 0:
            prev_cell = self.cells[index-1]
            prev_P = prev_cell['algebraic_variables']['P']
        else:
            prev_P = P
        dP_dz = self.formula_62(delta_z, P, prev_P)
        
        # 计算水力直径（公式9）
        # 计算气体横截面积（公式78）
        A_g = self.formula_78(A_t, A_s)
        # 计算单段总体积（公式58）
        V_delta = self.formula_58(A_t, delta_z)
        DH = self.formula_9(V_delta, A_g)
        
        # 计算气体密度（公式11）
        gas_C = {component: C[component] for component in gas_components}
        rho_g = self.formula_11(self.constants['gas_properties'], gas_C)
        
        # 获取气体粘度
        mu_g = self.parameters['mu_g']
        
        # 计算气体轴向速度（公式10）
        # 避免除以零错误，如果rho_g为零或非常小，使用默认值
        if rho_g <= 1e-10 or mu_g <= 1e-20:
            vg = 0.084  # 使用默认气体速度
        else:
            vg = self.formula_10(DH, mu_g, rho_g, dP_dz)
        
        # 获取固相速度，使用差异化配置中的公式
        solid_velocity_rule = rules.get('solid_velocity', 'formula_71')
        if isinstance(solid_velocity_rule, str) and solid_velocity_rule == 'formula_71':
            # 分解炉固相流速：vs = vg（公式71）
            vs = self.formula_71(vg)
        elif isinstance(solid_velocity_rule, str) and solid_velocity_rule == 'formula_13':
            # 回转窑固相流速（公式13）
            # 完整使用公式13计算
            omega = self.control_variables.get('omega', 0.4189)  # 窑转速
            psi = self.parameters['equipment']['kiln']['angle']  # 窑倾角
            
            # 计算休止角（公式12）
            a_omega = self.parameters.get('a_omega', 0.05)
            b_omega = self.parameters.get('b_omega', 0.6)
            xi = self.formula_12(a_omega, b_omega, omega)
            
            # 从预计算的代数变量中获取窑内半径、填充角
            r_c = algebraic_vars['r_c']
            theta = algebraic_vars['theta']
            
            # 计算弦长（公式14），使用预计算的填充角
            Lc = self.formula_14(r_c, theta)
            
            # 计算床层高度轴向梯度和床层坡度角
            # 1. 计算当前单元的床层高度（公式72）
            h_current = self.formula_72(r_c, theta)
            
            # 2. 获取前一个单元的床层高度
            h_prev = h_current  # 默认情况：如果是第一个单元，前一个单元的床层高度等于当前单元
            if index > 0: #若不是第一个单元
                prev_cell = self.cells[index-1]
                prev_theta = prev_cell['algebraic_variables'].get('theta', theta) #读取前一个单元的填充角
                h_prev = self.formula_72(r_c, prev_theta) #计算前一个单元的床层高度
            
            # 3. 构建床层高度列表，包含当前段和前一段的床层高度
            h_z = [h_prev, h_current]
            
            # 4. 使用公式74计算床层高度轴向梯度
            dh_z_dz = self.formula_74(h_z)
            
            # 5. 使用公式73计算床层坡度角
            phi_z = self.formula_73(dh_z_dz)
            
            # 计算固体速度（公式13）
            vs = self.formula_13(omega, psi, xi, r_c, Lc, phi_z)
        else:
            # 其他情况，预热器,使用默认值
            vs = algebraic_vars.get('vs', 0.0)
        
        # 计算固体组分通量（公式25）
        N_i_s = [] #初始化固体组分通量列表
        for component in solid_components: #遍历所有固体组分
            C_i_s = C[component] #读取该组分的浓度
            N_i = self.formula_25(vs, C_i_s)
            N_i_s.append(N_i) #将通量添加到N_i_s列表
        
        # 计算气体组分通量（公式21）
        N_i_g = []
        # 准备气体通量计算所需的额外参数
        components = gas_components #气体组分列表
        
        # 计算气体总摩尔浓度（公式77）
        cg = self.formula_77(C)
        
        # 计算各组分摩尔分数（公式61）
        xj = [] #初始化摩尔分数列表
        for component in gas_components: #初始化摩尔分数列表
            C_i_g_t = C[component] #读取该组分的浓度
            x_i = self.formula_61(C_i_g_t, cg)
            xj.append(x_i) #将摩尔分数添加到xj列表
        
        for i, component in enumerate(gas_components): #遍历所有气体组分，获取索引
            C_i_g = C[component]
            
            # 计算组分i浓度梯度（公式23）
            # 获取前一个单元的浓度
            if index > 0: #若不是第一个单元
                prev_cell = self.cells[index-1] #获取前一个单元
                prev_C = prev_cell['state_variables']['C']  # 获取前一个单元
                prev_C_i_g = prev_C[component]  # 读取前一个单元该组分的浓度
                delta_z = cell['delta_z'] #读取当前单元的空间步长
                dC_i_g_dz = self.formula_23(C_i_g, prev_C_i_g, delta_z)
            else:
                dC_i_g_dz = 0.0 #初始化浓度梯度为 0
            
            N_i = self.formula_21(vg, C_i_g, Tg, P, components, xj, cg, dC_i_g_dz, i)
            N_i_g.append(N_i)
   
#五、焓通量密度计算     
        # 3. 计算焓通量密度
        # 固体焓通量密度（公式68）
        H_tilde_s = self.formula_68(h_i_s, N_i_s)
        algebraic_vars['H_tilde_s'] = H_tilde_s #存储固体焓通量密度
        
        # 气体焓通量密度（公式69）
        H_tilde_g = self.formula_69(h_i_g, N_i_g)
        algebraic_vars['H_tilde_g'] = H_tilde_g
  
#六、气体比热容与对流换热计算      
        # 4. 计算气体组分的摩尔数
        V_delta = self.formula_58(algebraic_vars.get('A_t', 3.1416), cell['delta_z'])
        
        # 5. 计算各气体组分的摩尔热容（公式6）
        cp_i_list = [] #初始化气体组分摩尔热容列表
        M_i_list = [] #初始化气体组分摩尔质量列表
        n_i_list = [] #初始化气体组分摩尔数列表
        
        for component in gas_components: #遍历所有气体组分
            # 获取组分的热容系数
            if component in self.constants['molar_heat_capacity']: #判断是否有热容系数数据
                coeffs = self.constants['molar_heat_capacity'][component] #读取热容系数
            else:
                coeffs = {'C0': 0.0, 'C1': 0.0, 'C2': 0.0}
            
            C0 = coeffs['C0']
            C1 = coeffs['C1']
            C2 = coeffs['C2']
            
            # 计算摩尔热容（公式6），气体使用气体温度Tg
            cp_i = self.formula_6(C0, C1, C2, Tg)
            cp_i_list.append(cp_i)
            
            # 获取摩尔质量
            if component in self.constants['gas_properties']: #判断该组分是否有摩尔质量数据
                M_i = self.constants['gas_properties'][component]['molar_mass'] #读取摩尔质量
            else:
                M_i = 0.0
            M_i_list.append(M_i) #添加到摩尔质量列表
            
            # 计算气体组分摩尔数（公式57_gas）
            C_i_g_t = C[component] #读取该组分的浓度
            n_i = self.formula_57_gas(V_delta, C_i_g_t)
            n_i_list.append(n_i) #添加到摩尔数列表
        
        # 6. 计算气体比热容（公式81）
        cp_g = self.formula_81(cp_i_list, M_i_list, n_i_list)
        algebraic_vars['cp_g'] = cp_g #存储气体比热容
        
        # 7. 计算普朗特数和对流换热（使用差异化配置中的对流换热公式）
        # 使用公式81计算出的cp_g
        mu_g = self.parameters['mu_g']  # 气体平均粘度，Pa·s
        # 计算气体热导率（公式20），使用气体温度Tg
        k_g = self.formula_20(Tg)
        # 计算普朗特数
        Pr = self.formula_32(cp_g, mu_g, k_g)
        algebraic_vars['Pr'] = Pr #存储普朗特数
        
        # 使用差异化配置中的对流换热公式
        heat_transfer_rule = rules.get('heat_transfer', 'formula_35') #从执行规则中获取对流换热公式
        if isinstance(heat_transfer_rule, str) and heat_transfer_rule == 'formula_35':
            # 分解炉对流换热（公式35）
            # 完整使用公式35计算
            dp = self.constants.get('d_p', 3e-5)  # 颗粒直径
            
            # 计算气体密度（公式11）
            # 只传递气体组分的浓度
            gas_C = {component: C[component] for component in gas_components}
            rho_g = self.formula_11(self.constants['gas_properties'], gas_C)
            
            # 计算有效直径（公式34）
            # 从预计算的代数变量中获取窑内半径、填充角
            r_c = algebraic_vars['r_c']
            theta = algebraic_vars['theta']
            De = self.formula_34(r_c, theta)
            
            # 计算轴向雷诺数（公式33）
            ReD = self.formula_33(rho_g, vg, De, mu_g)
            
            # 计算对流换热量（公式35）
            Qgscv = self.formula_35(k_g, dp, ReD, Pr)
            algebraic_vars['Qgscv'] = Qgscv #存储分解炉对流换热量
            algebraic_vars['ReD'] = ReD #存储轴向雷诺数
        elif isinstance(heat_transfer_rule, str) and heat_transfer_rule == 'formula_39':
            # 回转窑对流换热（公式39）
            # 完整使用公式39计算
            # 计算弦长（公式14）
            # 从预计算的代数变量中获取窑内半径、填充角
            r_c = algebraic_vars['r_c']
            theta = algebraic_vars['theta']
            Lc = self.formula_14(r_c, theta)
            
            # 计算气固换热面积（公式40）
            Ags = self.formula_40(Lc, cell['delta_z'])
            
            # 计算有效直径（公式34）
            De = self.formula_34(r_c, theta)
            
            # 计算气体密度（公式11）
            rho_g = self.formula_11(self.constants['gas_properties'], C)
            
            # 计算轴向雷诺数（公式33）
            ReD = self.formula_33(rho_g, vg, De, mu_g)
            
            # 计算旋转雷诺数（公式36）
            Re_omega = self.formula_36(rho_g, De, mu_g, self.control_variables.get('omega', 0.4189))
            
            # 计算努塞尔数（公式37）
            Nu = self.formula_37(ReD, Re_omega, eta)
            
            # 计算对流换热系数（公式38）
            beta = self.formula_38(k_g, De, Nu)
            
            # 计算对流换热量（公式39）
            Qgscv = self.formula_39(Ags, beta, Tg, Ts)
            algebraic_vars['Qgscv'] = Qgscv #存储回转窑对流换热量
            algebraic_vars['Ags'] = Ags
            algebraic_vars['beta'] = beta #存储对流换热系数
            algebraic_vars['ReD'] = ReD
        
#七、内能密度约束计算
        # 9. 计算内能密度约束，使用差异化配置中的公式
        internal_energy_constraint = rules.get('internal_energy_constraint', 'formula_65')
        
        # 计算固体和气体的总焓（公式4）
        # 准备总焓计算所需参数
        all_components = self.constants['stoichiometric_matrix']['components'] #提取所有组分名称
        solid_components = all_components[:9]
        gas_components = all_components[9:]
        
        # 计算各组分的积分焓和摩尔数
        Hf_i_s = [] #固体组分标准生成焓列表
        n_i_s = [] #固体组分摩尔数列表
        integral_enthalpy_s = [] #固体组分积分焓列表
        
        Hf_i_g = []
        n_i_g = []
        integral_enthalpy_g = []
        
        T0 = self.constants['T0']
        
        # 计算单段总体积（公式58），用于计算摩尔数
        V_delta = self.formula_58(algebraic_vars['A_t'], cell['delta_z'])
        
        for component in solid_components: #遍历固体组分，准备总焓计算参数
            # 固体组分
            if component in self.constants['molar_heat_capacity']: #判断是否有热容系数数据
                coeffs = self.constants['molar_heat_capacity'][component] #读取热容系数
            else: 
                coeffs = {'C0': 0.0, 'C1': 0.0, 'C2': 0.0}
            
            C0 = coeffs['C0']
            C1 = coeffs['C1']
            C2 = coeffs['C2']
            
            integral_h = self.formula_5(C0, C1, C2, T0, Ts)
            integral_enthalpy_s.append(integral_h) #添加到固体积分焓列表
            
            Hf_i_s.append(self.constants['standard_enthalpy'].get(component, 0.0)) #添加固体标准生成焓
            
            # 使用公式57计算固体组分摩尔数
            C_i_s_t = C[component]  # 读取固体组分浓度
            n_i = self.formula_57(V_delta, C_i_s_t)  # 计算固体组分摩尔数
            n_i_s.append(n_i) #添加到固体摩尔数列表
        
        for component in gas_components:
            # 气体组分
            if component in self.constants['molar_heat_capacity']:
                coeffs = self.constants['molar_heat_capacity'][component]
            else:
                coeffs = {'C0': 0.0, 'C1': 0.0, 'C2': 0.0}
            
            C0 = coeffs['C0']
            C1 = coeffs['C1']
            C2 = coeffs['C2']
            
            integral_h = self.formula_5(C0, C1, C2, T0, Tg)
            integral_enthalpy_g.append(integral_h)
            
            Hf_i_g.append(self.constants['standard_enthalpy'].get(component, 0.0))
            
            # 使用公式57_gas计算气体组分摩尔数
            C_i_g_t = C[component]
            n_i = self.formula_57_gas(V_delta, C_i_g_t)
            n_i_g.append(n_i)
        
        # 计算固体总焓（公式4）
        H_s = self.formula_4(Hf_i_s, n_i_s, integral_enthalpy_s)
        
        # 计算气体总焓（公式4）
        H_g = self.formula_4(Hf_i_g, n_i_g, integral_enthalpy_g)
        
        # 计算单段总体积
        r_c = self.parameters['equipment'][device_type]['radius']
        A_t = self.formula_17(r_c)
        V_delta = self.formula_58(A_t, cell['delta_z'])
        
        # 计算固体单位体积焓（公式63）
        H_hat_s = self.formula_63(H_s, V_delta)
        
        # 计算气体单位体积焓（公式64）
        H_hat_g = self.formula_64(H_g, V_delta)
        
        # 计算气体和固体体积分数
        # 计算气体体积（公式7）
        #构建气体组分摩尔质量字典
        gas_molar_mass = {comp: self.constants['gas_properties'][comp]['molar_mass'] for comp in gas_components if comp in self.constants['gas_properties']}
        gas_moles = {comp: C[comp] for comp in gas_components}#构建气体组分摩尔数字典
        V_g = self.formula_7(self.constants['R'], Tg, P, gas_moles)
        
        # 计算固体体积（公式8）
        solid_molar_mass = {comp: self.constants['solid_properties'][comp]['molar_mass'] for comp in solid_components if comp in self.constants['solid_properties']}
        solid_density = {comp: self.constants['solid_properties'][comp]['density'] for comp in solid_components if comp in self.constants['solid_properties']}
        solid_moles = {comp: C[comp] for comp in solid_components}
        V_s = self.formula_8(solid_molar_mass, solid_density, solid_moles)
        
        # 计算体积分数
        V_s_fraction = self.formula_66(V_s, V_delta)
        V_g_fraction = self.formula_67(V_g, V_delta)
        
        # 验证体积分数守恒（公式1）
        volume_fraction_sum = self.formula_1(V_g_fraction, V_s_fraction)
        
        # 如果体积分数和不为1，进行调整以满足守恒条件
        if abs(volume_fraction_sum - 1.0) > 1e-10:
            # 以固体体积分数为基准，调整气体体积分数
            V_g_fraction = 1.0 - V_s_fraction
            # 或者以气体体积分数为基准，调整固体体积分数
            # V_s_fraction = 1.0 - V_g_fraction
        
        if isinstance(internal_energy_constraint, str): #若内能约束为单公式
            if internal_energy_constraint == 'formula_65':
                # 分解炉内能密度约束（公式65）
                U_hat = self.formula_65(P, H_hat_g, H_hat_s, V_g_fraction)
                algebraic_vars['U_hat'] = U_hat #存储分解炉内能密度
        elif isinstance(internal_energy_constraint, list): #若内能约束为多公式
            # 回转窑内能密度约束（气体和固体）
            if len(internal_energy_constraint) >= 2: #若约束列表包含至少两个公式
                # 气体内能密度约束（公式2）
                U_g = self.formula_2(P, H_hat_g, V_g_fraction)
                algebraic_vars['U_g'] = U_g #存储气体内能密度
                
                # 固体内能密度约束（公式3）
                U_s = self.formula_3(H_hat_s)
                algebraic_vars['U_s'] = U_s
                
                # 总内能密度
                U_hat = U_g + U_s
                algebraic_vars['U_hat'] = U_hat #存储回转窑总内能密度
        
        # 10. 设置其他必要的代数变量
        algebraic_vars['Tg'] = Tg #存储气体温度 Tg
        algebraic_vars['Ts'] = Ts #存储固体温度 Ts
        algebraic_vars['P'] = P 
        algebraic_vars['vg'] = vg
        algebraic_vars['vs'] = vs
        algebraic_vars['A_t'] = A_t
        algebraic_vars['A_s'] = A_s
        algebraic_vars['eta'] = eta
        algebraic_vars['theta'] = theta #存储填充角 θ
        
        # 11. 更新单元的代数变量
        cell['algebraic_variables'].update(algebraic_vars)
        #返回当前单元的所有代数变量
        return algebraic_vars
    
    def _solve_differential_equations(self, cell, rules, algebraic_vars):
        # 求解微分方程：当前有限体积单元、单元执行规则、前序代数计算得到的变量
        # 计算浓度变化率 dC/dt 和内能变化率 dU/dt
        dC_dt = {} #存储各组分的浓度时间变化率
        dU_dt = 0.0  # 初始化为浮点数，因为U_hat是标量值
        
        index = cell['index'] #获取当前单元的索引
        
        # 获取当前单元的状态变量
        state_vars = cell['state_variables'] #获取单元的状态变量字典
        C = state_vars['C']  # 组分浓度向量
        U_hat = state_vars['U_hat']  # 单位体积内能向量
        
        # 获取单元所属设备类型
        section = cell['device_type']
        
        # 获取从代数计算得到的约束变量
        Tg = algebraic_vars.get('Tg', 1000.0)  # 气体温度，从代数计算得到
        Ts = algebraic_vars.get('Ts', 1000.0)  # 固体温度，从代数计算得到
        P = algebraic_vars.get('P', self.constants['P0'])  # 压力，从代数计算得到
        vg = algebraic_vars.get('vg', 0.084)  # 气体速度，从代数计算得到
        vs = algebraic_vars.get('vs', 0.0)  # 固体速度，从代数计算得到
        
        # 获取所有组分列表
        stoichiometric_matrix = self.constants['stoichiometric_matrix'] #从化学计量矩阵中获取组分信息
        all_components = stoichiometric_matrix['components'] #提取所有 15 个组分的名称列表
        solid_components = all_components[:9]
        gas_components = all_components[9:]

#一、反应源项计算        
        # 1. 计算反应源项（公式27）
        # 基于质量守恒定律，结合反应速率和生成速率计算
        
        # 根据差异化配置获取当前设备的反应列表
        reactions = rules.get('reactions', [])
        
        # 如果没有指定反应列表，根据设备类型设置默认反应
        if not reactions: #判断反应列表是否为空
            if section == 'calciner':
                reactions = ['r1', 'r6', 'r7', 'r8', 'r9', 'r10', 'r11']  # 分解炉仅r1分解反应
            elif section == 'kiln':
                reactions = ['r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'r7', 'r8', 'r9', 'r10', 'r11']  # 回转窑r2~r5熟料反应
            else:
                reactions = []  # 其他设备无反应
        
        # 复制原始反应速率系数字典
        original_reactions = self.constants['reaction_rate_coefficients'].copy()
        
        try: #异常捕获起始，确保无论计算是否出错，都能恢复原始反应列表
            # 根据设备类型过滤反应速率系数
            if reactions: #若反应列表非空
                #筛选出当前设备的反应速率系数
                filtered_reactions = {r: original_reactions[r] for r in reactions if r in original_reactions}
                #更新反应速率系数字典，仅保留当前设备的反应
                self.constants['reaction_rate_coefficients'] = filtered_reactions
            
            P_i_dict = {}  # 组分分压字典，从代数计算得到
            C_i_dict = C  # 组分浓度字典，当前状态变量
            Rs, Rg = self._calculate_reaction_source(Tg, P_i_dict, C_i_dict)
        finally: #异常捕获结束，确保必须执行的恢复操作
            # 恢复原有反应列表
            self.constants['reaction_rate_coefficients'] = original_reactions
   
#二、浓度变化率计算     
        # 2. 计算通量和通量梯度
        
        # 固体组分处理
        for i, component in enumerate(solid_components): #遍历所有固体组分，逐个计算浓度变化率
            # 计算固体通量（公式25）
            N_i_s = self.formula_25(vs, C[component])
            
            # 计算通量梯度（公式29）
            # 获取前一个单元的通量
            prev_N_i_s = N_i_s  # 默认边界条件：入口处梯度为0
            if index > 0: #若不是第一个单元
                prev_cell = self.cells[index-1] #获取前一个单元
                prev_C = prev_cell['state_variables']['C']  #读取前一个单元的浓度字典
                prev_vs = prev_cell['algebraic_variables'].get('vs', vs)  #读取前一个单元的固体速度
                prev_N_i_s = self.formula_25(prev_vs, prev_C[component])
            
            # 获取轴向步长
            delta_z = cell['delta_z']
            
            # 计算通量梯度
            dN_i_s_dz = self.formula_29(N_i_s, prev_N_i_s, delta_z)
            
            # 计算浓度变化率（公式28）
            dC_dt[component] = self.formula_28(dN_i_s_dz, Rs[i])
        
        # 气体组分处理
        # 气体组分列表
        components = gas_components
        
        # 计算气体总摩尔浓度（公式77）
        cg = self.formula_77(C)
        
        # 计算各组分摩尔分数（公式61）
        xj = []
        for component in gas_components: #遍历所有气体组分
            C_i_g_t = C[component] #读取该组分的浓度
            x_i = self.formula_61(C_i_g_t, cg)
            xj.append(x_i) #将摩尔分数添加到列表
        
        # 初始浓度梯度为0
        dC_i_g_dz = 0.0
        
        for i, component in enumerate(gas_components): #遍历所有气体组分，逐个计算浓度变化率
            # 计算气体通量（公式21）
            N_i_g = self.formula_21(vg, C[component], Tg, P, components, xj, cg, dC_i_g_dz, i)
            
            # 计算通量梯度（公式31）
            # 获取前一个单元的通量
            prev_N_i_g = N_i_g  # 默认边界条件：入口处梯度为0
            if index > 0:
                prev_cell = self.cells[index-1] #获取前一个单元
                prev_C = prev_cell['state_variables']['C'] #读取前一个单元的浓度字典
                prev_vg = prev_cell['algebraic_variables'].get('vg', vg) #读取前一个单元的气体速度
                prev_N_i_g = self.formula_21(prev_vg, prev_C[component], Tg, P, components, xj, cg, dC_i_g_dz, i)
            
            # 获取轴向步长
            delta_z = cell['delta_z']
            
            # 计算通量梯度
            dN_i_g_dz = self.formula_31(N_i_g, prev_N_i_g, delta_z)
            
            # 计算浓度变化率（公式30）
            dC_dt[component] = self.formula_30(dN_i_g_dz, Rg[i])

#三、内能变化率计算            
        # 3. 计算内能变化率 dŨ/dt
        # 基于能量守恒定律，整合各种能量项
        # 使用执行规则中指定的能量方程，根据差异化配置计算
        H_tilde_s = algebraic_vars.get('H_tilde_s', 0.0)  # 从代数变量中读取固体焓通量密度
        H_tilde_g = algebraic_vars.get('H_tilde_g', 0.0)  # 气体焓通量密度
        
        # 从执行规则中获取能量方程、执行规则中的差异化配置
        energy_eq = rules.get('energy_eq', 'formula_42')
        
        # 提前计算通用参数，避免重复计算
        sigma = self.constants.get('sigma', 5.67e-8)
        epsilon_s = self.constants.get('epsilon_s', 0.9)
            
        # 根据能量方程类型计算内能变化率
        if isinstance(energy_eq, str):
            # 单方程情况（如分解炉）
            if energy_eq == 'formula_42':
                # 分解炉能量平衡：∂_t Û_c = -∂_z (H~_g + H~_s + Q~_g) - (Qgs^rad + Qgs^cv) / VΔ
                
                # 计算气体总摩尔浓度（公式77）
                cg = self.formula_77(C)
                
                # 计算气体发射率（公式53）
                # 获取气体组分摩尔分数
                xH2O = C.get('H2O', 0.0) / cg if cg > 0 else 0.0
                xCO2 = C.get('CO2', 0.0) / cg if cg > 0 else 0.0
                r_c = self.parameters['equipment']['calciner']['radius']
                epsilon_g = self.formula_53(Tg, P, xH2O, xCO2, r_c)
                
                # 计算气固发射率（公式80）
                epsilon_gs = self.formula_80(epsilon_g, epsilon_s)
                
                # 计算气固换热面积（公式40）
                # 使用预计算的填充角和相关参数
                theta = algebraic_vars.get('theta', math.pi)  # 使用代数计算中预计算的填充角
                r_c = algebraic_vars.get('r_c', r_c)  # 使用代数计算中预计算的半径
                Lc = self.formula_14(r_c, theta)
                Ags = self.formula_40(Lc, cell['delta_z'])
                
                # 计算辐射换热量
                Qgsrad = self.formula_41(sigma, Ags, epsilon_gs, Tg, Ts)
                
                # 计算气体温度梯度（公式43）
                current_Tg = Tg
                prev_Tg = current_Tg
                if index > 0: 
                    prev_cell = self.cells[index-1] 
                    prev_Tg = prev_cell['algebraic_variables'].get('Tg', current_Tg) 
                Tg_list = [prev_Tg, current_Tg]
                dTg_dz = self.formula_43(Tg_list)
                
                # 计算气体热传导（公式54），使用气体温度梯度
                k_g = self.formula_20(Tg)
                Q_tilde_g = self.formula_54(k_g, dTg_dz)
                
                # 计算固体温度梯度（公式43）
                current_Ts = Ts
                prev_Ts = current_Ts
                if index > 0: 
                    prev_cell = self.cells[index-1] 
                    prev_Ts = prev_cell['algebraic_variables'].get('Ts', current_Ts) 
                Ts_list = [prev_Ts, current_Ts]
                dTs_dz = self.formula_43(Ts_list)
                
                # 从代数变量中获取对流换热量
                Qgscv = algebraic_vars.get('Qgscv', 0.0)
                
                # 计算单段总体积（公式58）
                A_t = algebraic_vars.get('A_t', 0.0)  # 从代数变量中获取总横截面积
                V_delta = self.formula_58(A_t, cell['delta_z'])
                
                # 计算反应焓变（公式52）
                # 分解炉包含的反应：r1、r6、r7、r8、r9、r10、r11
                reaction_rates = {}
                reaction_enthalpies = self.constants['reaction_enthalpies']
                
                # 计算所有相关反应的速率
                for reaction_id in reactions:
                    P_i_dict = {}  # 组分分压字典
                    C_i_dict = C  # 组分浓度字典
                    # 对于不同反应，使用不同的温度
                    if reaction_id in ['r1']:  # 固体反应使用固体温度
                        T_react = Ts
                    else:  # 气体反应使用气体温度
                        T_react = Tg
                    # 计算反应速率
                    reaction_rates[reaction_id] = self._calculate_reaction_rate(reaction_id, T_react, P_i_dict, C_i_dict)
                
                # 计算总反应焓变
                Jsg = self.formula_52(reaction_rates, reaction_enthalpies)
                
                # 获取前一段的焓通量密度和热传导
                H_tilde_s_zk_prev = H_tilde_s
                H_tilde_g_zk_prev = H_tilde_g
                Q_tilde_g_zk_prev = Q_tilde_g
                if index > 0:
                    prev_cell = self.cells[index-1]
                    prev_algebraic_vars = prev_cell['algebraic_variables']
                    H_tilde_s_zk_prev = prev_algebraic_vars.get('H_tilde_s', H_tilde_s)
                    H_tilde_g_zk_prev = prev_algebraic_vars.get('H_tilde_g', H_tilde_g)
                    # 前一段的热传导需要重新计算
                    prev_Tg = prev_cell['algebraic_variables'].get('Tg', Tg)
                    prev_dTg_dz = self.formula_43([prev_cell['algebraic_variables'].get('Tg', Tg), Tg])
                    prev_k_g = self.formula_20(prev_Tg)
                    Q_tilde_g_zk_prev = self.formula_54(prev_k_g, prev_dTg_dz)
                
                # 计算分解炉内能变化率，添加反应焓变Jsg
                dU_dt = self.formula_42(
                    H_tilde_s, H_tilde_g, Q_tilde_g, 
                    H_tilde_s_zk_prev, H_tilde_g_zk_prev, Q_tilde_g_zk_prev, 
                    cell['delta_z'], Qgscv, Qgsrad, V_delta, Jsg
                )
        elif isinstance(energy_eq, list):
            # 多方程情况（如回转窑气固双方程）
            if len(energy_eq) >= 2 and energy_eq[0] == 'formula_44' and energy_eq[1] == 'formula_45':
                # 固体能量平衡和气体能量平衡
                
                # 计算气体总摩尔浓度（公式77）
                cg = self.formula_77(C)
                
                # 计算气体发射率（公式53）
                xH2O = C.get('H2O', 0.0) / cg if cg > 0 else 0.0
                xCO2 = C.get('CO2', 0.0) / cg if cg > 0 else 0.0
                r_c = self.parameters['equipment']['kiln']['radius']
                epsilon_g = self.formula_53(Tg, P, xH2O, xCO2, r_c)
                epsilon_gs = self.formula_80(epsilon_g, epsilon_s)
                
                # 计算气固换热面积（公式40）
                # 使用预计算的填充角和相关参数
                theta = algebraic_vars.get('theta', math.pi)  # 使用代数计算中预计算的填充角
                r_c = algebraic_vars.get('r_c', r_c)  # 使用代数计算中预计算的半径
                Lc = self.formula_14(r_c, theta)
                Ags = self.formula_40(Lc, cell['delta_z'])
                Qgsrad = self.formula_41(sigma, Ags, epsilon_gs, Tg, Ts)  # 辐射换热量
                
                # 从代数变量中获取对流换热量
                Qgscv = algebraic_vars.get('Qgscv', 0.0)
                
                # 计算单段总体积（公式58）
                A_t = algebraic_vars.get('A_t', 0.0)  # 从代数变量中获取总横截面积
                V_delta = self.formula_58(A_t, cell['delta_z'])
                
                # 计算气体温度梯度（公式43）
                current_Tg = Tg
                prev_Tg = current_Tg
                if index > 0:
                    prev_cell = self.cells[index-1] 
                    prev_Tg = prev_cell['algebraic_variables'].get('Tg', current_Tg) 
                Tg_list = [prev_Tg, current_Tg]
                dTg_dz = self.formula_43(Tg_list)
                
                # 计算固体温度梯度（公式43）
                current_Ts = Ts
                prev_Ts = current_Ts
                if index > 0:
                    prev_cell = self.cells[index-1] 
                    prev_Ts = prev_cell['algebraic_variables'].get('Ts', current_Ts) 
                Ts_list = [prev_Ts, current_Ts]
                dTs_dz = self.formula_43(Ts_list)
                
                # 计算热传导
                k_g = self.formula_20(Tg)
                Q_tilde_g = self.formula_54(k_g, dTg_dz)  # 使用气体温度梯度
                Q_tilde_s = self.formula_55(dTs_dz)  # 使用固体温度梯度
                
                # 计算反应焓变（公式52）
                # 根据设备类型和能量方程类型确定需要计算的反应
                reaction_rates = {}
                reaction_enthalpies = self.constants['reaction_enthalpies']
                
                # 计算所有相关反应的速率
                for reaction_id in reactions:
                    P_i_dict = {}  # 组分分压字典
                    C_i_dict = C  # 组分浓度字典
                    # 对于不同反应，使用不同的温度
                    if reaction_id in ['r1', 'r2', 'r3', 'r4', 'r5']:  # 固体反应使用固体温度
                        T_react = Ts
                    else:  # 气体反应使用气体温度
                        T_react = Tg
                    # 计算反应速率
                    reaction_rates[reaction_id] = self._calculate_reaction_rate(reaction_id, T_react, P_i_dict, C_i_dict)
                
                # 计算总反应焓变
                Jsg = self.formula_52(reaction_rates, reaction_enthalpies)
                
                # 获取前一段的焓通量密度和热传导
                H_tilde_s_zk_prev = H_tilde_s
                H_tilde_g_zk_prev = H_tilde_g
                Q_tilde_s_zk_prev = Q_tilde_s
                Q_tilde_g_zk_prev = Q_tilde_g
                if index > 0:
                    prev_cell = self.cells[index-1]
                    prev_algebraic_vars = prev_cell['algebraic_variables']
                    H_tilde_s_zk_prev = prev_algebraic_vars.get('H_tilde_s', H_tilde_s)
                    H_tilde_g_zk_prev = prev_algebraic_vars.get('H_tilde_g', H_tilde_g)
                    # 前一段的热传导需要重新计算
                    prev_Ts = prev_cell['algebraic_variables'].get('Ts', Ts)
                    prev_Tg = prev_cell['algebraic_variables'].get('Tg', Tg)
                    prev_dTs_dz = self.formula_43([prev_cell['algebraic_variables'].get('Ts', Ts), Ts])
                    prev_dTg_dz = self.formula_43([prev_cell['algebraic_variables'].get('Tg', Tg), Tg])
                    prev_k_g = self.formula_20(prev_Tg)
                    Q_tilde_g_zk_prev = self.formula_54(prev_k_g, prev_dTg_dz)
                    Q_tilde_s_zk_prev = self.formula_55(prev_dTs_dz)
                
                # 计算固体能量变化率
                dU_dt_solid = self.formula_44(
                    H_tilde_s, Q_tilde_s, 
                    H_tilde_s_zk_prev, Q_tilde_s_zk_prev, 
                    cell['delta_z'], Qgsrad, Qgscv, V_delta, Jsg
                )
                
                # 计算气体能量变化率
                dU_dt_gas = self.formula_45(
                    H_tilde_g, Q_tilde_g, 
                    H_tilde_g_zk_prev, Q_tilde_g_zk_prev, 
                    cell['delta_z'], Qgsrad, Qgscv, V_delta, Jsg
                )
                
                # 综合固体和气体的内能变化率
                dU_dt = (dU_dt_solid + dU_dt_gas) / 2
        else:
            # 默认使用分解炉能量方程
            # 使用当前段的值作为前一段的值
            dU_dt = self.formula_42(
                H_tilde_s, H_tilde_g, 0.0, 
                H_tilde_s, H_tilde_g, 0.0, 
                cell['delta_z'], 0.0, 0.0, cell['delta_z']
            )
        
        return dC_dt, dU_dt #返回浓度变化率字典dC_dt和内能变化率dU_dt 
 
#一、变量更新方法   
    def _update_variables(self, cell, dC_dt, dU_dt):
        # 同步更新单元变量：当前单元、浓度变化率、内能变化率
        # 采用「代数 - 微分同步计算、结果互代」模式
        # 同一时间步 Δt 内，同一单元里，代数计算和微分计算结果互相代入
        index = cell['index'] #获取当前单元的索引，用于同步更新全局变量数组
          
        # 获取当前单元的状态变量和代数变量
        state_vars = cell['state_variables']
        algebraic_vars = cell['algebraic_variables']
        
        # 1. 更新状态变量：浓度 C 和单位体积内能 Ũ
        # 按公式「新值 = 原值 + 变化率 × Δt」
        
        # 更新浓度变量 C
        for component in dC_dt: #遍历所有有浓度变化率的组分
            # 同步更新状态变量
            state_vars['C'][component] += dC_dt[component] * self.dt #按公式更新浓度（原值 + 变化率 ×dt）
            # 确保浓度非负
            if state_vars['C'][component] < 0: #判断浓度是否小于 0
                state_vars['C'][component] = 0.0
            # 同步更新全局数组
            self.variables['state_variables']['C'][index][component] = state_vars['C'][component]
        
        # 更新内能变量 Ũ
        state_vars['U_hat'] += dU_dt * self.dt #按公式更新单位体积内能（原值 + 变化率 ×dt）
        # 同步更新全局数组
        self.variables['state_variables']['U_hat'][index] = state_vars['U_hat']
        
        # 2. 根据更新后的状态变量，重新计算并更新代数变量
        # 确保代数变量 T、P 随状态变量更新同步修正，保持变量耦合一致性
        
        # 重新获取所有组分列表
        stoichiometric_matrix = self.constants['stoichiometric_matrix'] #从化学计量矩阵获取组分信息
        all_components = stoichiometric_matrix['components'] #所有组分列表
        solid_components = all_components[:9] 
        gas_components = all_components[9:]
        
        # 重新计算温度（基于能量守恒和比热容）
        # 温度变化 = 内能变化 / (质量 × 比热容)
        # 分别计算固体和气体的质量、比热容和温度变化
        
        # 获取当前单元的设备类型
        device_type = cell['device_type']
        
        # 1. 计算固体相关参数
        # 计算固体质量
        solid_mass = 0.0
        for component in solid_components: #遍历所有固体组分
            concentration = state_vars['C'][component] #获取浓度值
            if concentration > 0:
                molar_mass = self.constants['solid_properties'][component]['molar_mass']
                # 浓度单位是mol/m³，质量 = 浓度 × 摩尔质量 × 体积
                solid_mass += concentration * molar_mass
        
        # 计算固体比热容（使用平均比热容）
        solid_cp = 800.0  # J/(kg·K)，固体平均比热容
        
        # 2. 计算气体相关参数
        # 计算气体质量
        gas_mass = 0.0
        for component in gas_components:
            concentration = state_vars['C'][component]
            if concentration > 0:
                molar_mass = self.constants['gas_properties'][component]['molar_mass']
                # 浓度单位是mol/m³，质量 = 浓度 × 摩尔质量 × 体积
                gas_mass += concentration * molar_mass
        
        # 计算气体比热容（使用平均比热容）
        gas_cp = 1000.0  # J/(kg·K)，气体平均比热容
        
        # 3. 计算温度变化
        # 温度变化 = 内能变化 / (质量 × 比热容)
        # 假设内能变化dU_dt已经考虑了体积因素，单位是J/m³·s
        dt = self.dt  # 时间步长
        
        # 确保dU_dt是实数，处理复数情况
        if isinstance(dU_dt, complex):
            dU_dt = dU_dt.real
        dU_dt = float(dU_dt)
        
        # 更新气体温度
        if gas_mass > 0 and gas_cp > 0:
            # 计算气体内能变化
            gas_energy_change = dU_dt * dt  # J/m³
            # 温度变化 = 能量变化 / (质量 × 比热容)
            Tg_change = gas_energy_change / (gas_mass * gas_cp)  # K
            # 确保Tg_change是实数，处理复数情况
            if isinstance(Tg_change, complex):
                Tg_change = Tg_change.real
            Tg_change = float(Tg_change)
            new_Tg = algebraic_vars['Tg'] + Tg_change
            # 确保new_Tg是实数，处理复数情况
            if isinstance(new_Tg, complex):
                new_Tg = new_Tg.real
            new_Tg = float(new_Tg)
            # 确保温度在合理范围内
            if new_Tg < 300.0:
                new_Tg = 300.0
            if new_Tg > 2000.0:
                new_Tg = 2000.0
            algebraic_vars['Tg'] = new_Tg
        
        # 更新固体温度
        if solid_mass > 0 and solid_cp > 0:
            # 计算固体内能变化
            solid_energy_change = dU_dt * dt  # J/m³
            # 温度变化 = 能量变化 / (质量 × 比热容)
            Ts_change = solid_energy_change / (solid_mass * solid_cp)  # K
            # 确保Ts_change是实数，处理复数情况
            if isinstance(Ts_change, complex):
                Ts_change = Ts_change.real
            Ts_change = float(Ts_change)
            new_Ts = algebraic_vars['Ts'] + Ts_change
            # 确保new_Ts是实数，处理复数情况
            if isinstance(new_Ts, complex):
                new_Ts = new_Ts.real
            new_Ts = float(new_Ts)
            # 确保温度在合理范围内
            if new_Ts < 300.0:
                new_Ts = 300.0
            if new_Ts > 2000.0:
                new_Ts = 2000.0
            algebraic_vars['Ts'] = new_Ts
        
        # 重新计算压力 P（基于理想气体状态方程）
        # 压力与气体浓度和温度成正比
        gas_total_concentration = self.formula_77(state_vars['C']) #计算气体总摩尔浓度
        
        # 避免除零错误的压力更新逻辑
        denominator = gas_total_concentration - sum(dC_dt.values()) * self.dt #计算分母（更新前的气体总浓度）
        if denominator <= 0: #判断分母是否为零或负数
            # 分母为零或负数时，使用更安全的压力更新方式
            new_P = algebraic_vars['P'] * (1 + 0.01 * self.dt)  # 缓慢增加压力
        else:
            new_P = algebraic_vars['P'] * (gas_total_concentration / denominator)  # 按浓度比值更新压力
        
        # 限制压力范围：压力限制在标准大气压的 0.5~1.5 倍
        algebraic_vars['P'] = max(min(new_P, self.constants['P0'] * 1.5), self.constants['P0'] * 0.5)  # 压力范围限制
        
        # 3. 同步更新全局数组中的代数变量
        self.variables['algebraic_variables']['Tg'][index] = algebraic_vars['Tg']
        self.variables['algebraic_variables']['Ts'][index] = algebraic_vars['Ts']
        self.variables['algebraic_variables']['P'][index] = algebraic_vars['P']
        
        # 4. 确保变量耦合一致性
        # 将更新后的代数变量同步到单元字典
        cell['algebraic_variables'].update(algebraic_vars)
        
        # 5. 同步更新全局设备类型标识
        self.variables['device_type'][index] = cell['device_type']
    
#二、相邻单元传递方法
    def _transfer_between_cells(self, cell1, cell2):
        # 相邻单元连续传递：前一个单元、后一个单元
        # 前一区段最后一个单元 = 后一区段第一个单元
        # 复制单元1的状态变量到单元2
        cell2['state_variables']['C'] = cell1['state_variables']['C'].copy() #复制浓度字典
        cell2['state_variables']['U_hat'] = cell1['state_variables']['U_hat'] #复制内能字典
        
        # 复制单元1的代数变量到单元2
        # 复制气固温度
        if 'Tg' in cell1['algebraic_variables']:
            cell2['algebraic_variables']['Tg'] = cell1['algebraic_variables']['Tg']
        if 'Ts' in cell1['algebraic_variables']:
            cell2['algebraic_variables']['Ts'] = cell1['algebraic_variables']['Ts']
        if 'P' in cell1['algebraic_variables']:
            cell2['algebraic_variables']['P'] = cell1['algebraic_variables']['P']

#三、预热器进料补给方法
    def _feed_preheater_inlet(self):
        # 预热器入口按初始进料值「补料」
        # 判断是否存在单元，且第一个单元为预热器（进料入口）
        if self.cells and self.cells[0]['device_type'] == 'preheater':
            # 获取预热器第一个单元
            inlet_cell = self.cells[0]
            
            # 设置初始进料温度，区分气固温度
            inlet_cell['algebraic_variables']['Tg'] = self.control_variables['gas_feed']['temperature']
            inlet_cell['algebraic_variables']['Ts'] = self.control_variables['solid_feed']['temperature']
            
            # 设置初始进料浓度
            # 计算单段总体积（公式58）
            A_t = self.formula_17(self.parameters['equipment']['preheater']['radius'])
            delta_z = inlet_cell['delta_z']
            V_delta = self.formula_58(A_t, delta_z)
            
            # 固体进料
            total_solid_rate = self.control_variables['solid_feed']['total_rate'] #获取固体总进料速率
            for component, mass_fraction in self.control_variables['solid_feed']['composition'].items(): #遍历固体进料的各组分
                # 计算该组分的质量流量（总速率 × 质量分数）
                mass_rate = total_solid_rate * mass_fraction
                # 获取摩尔质量
                molar_mass = self.constants['solid_properties'][component]['molar_mass']
                # 入口边界初始化固体组分摩尔数（公式56）
                n_i_s_0 = self.formula_56(mass_rate, molar_mass)
                # 初始化固体组分i浓度（公式59）
                C_i_s_0 = self.formula_59(V_delta, n_i_s_0)
                # 设置浓度
                inlet_cell['state_variables']['C'][component] = C_i_s_0
            
            # 气体进料
            gas_temp = self.control_variables['gas_feed']['temperature'] #获取气体进料温度
            gas_pressure = self.constants['P0']
            R = self.constants['R']
            for component, mole_fraction in self.control_variables['gas_feed']['composition'].items(): #遍历气体进料的各组分
                # 初始化气体组分i浓度（公式60）
                C_i_g_0 = self.formula_60(mole_fraction, gas_pressure, gas_temp, R)
                # 设置浓度
                inlet_cell['state_variables']['C'][component] = C_i_g_0
            # 燃料进料处理
            fuel = self.control_variables['fuel']
            fuel_rate = fuel['rate']  # g/s
            fuel_composition = fuel['composition']
            
            # 计算燃料各成分的浓度
            for component, mass_fraction in fuel_composition.items():
                # 计算该组分的质量流量（总速率 × 质量分数）
                mass_rate = fuel_rate * mass_fraction  # g/s
                
                # 获取摩尔质量
                if component in self.constants['solid_properties']:
                    molar_mass = self.constants['solid_properties'][component]['molar_mass']
                elif component in self.constants['gas_properties']:
                    molar_mass = self.constants['gas_properties'][component]['molar_mass']
                else:
                    # 使用平均摩尔质量
                    molar_mass = 28.0  # g/mol，近似值
                
                # 计算摩尔流量
                mole_rate = mass_rate / molar_mass  # mol/s
                
                # 计算浓度 (mol/m³) = 摩尔流量 / 体积流量
                # 体积流量 = 气体进料速率 / 气体密度
                # 使用单段总体积计算浓度
                fuel_concentration = mole_rate / V_delta  # mol/m³
                
                # 添加到浓度字典中
                if component in inlet_cell['state_variables']['C']:
                    inlet_cell['state_variables']['C'][component] += fuel_concentration
                else:
                    inlet_cell['state_variables']['C'][component] = fuel_concentration
 
 
#四、区段衔接处理方法 :处理预热器→分解炉、分解炉→回转窑的区段边界传递      
    def _handle_section_transitions(self):
        # 处理不同区段之间的边界衔接
        preheater_cells = [cell for cell in self.cells if cell['device_type'] == 'preheater'] #筛选所有预热器单元
        calciner_cells = [cell for cell in self.cells if cell['device_type'] == 'calciner']
        kiln_cells = [cell for cell in self.cells if cell['device_type'] == 'kiln']
        
        # 1. 预热器最后一个单元 → 分解炉第一个单元
        if preheater_cells and calciner_cells: #判断预热器和分解炉均有单元
            preheater_last = preheater_cells[-1] #获取预热器最后一个单元
            calciner_first = calciner_cells[0]#获取分解炉第一个单元
            # 调用传递方法:传递「升温后温度 + 原始浓度」
            self._transfer_between_cells(preheater_last, calciner_first)
        
        # 2. 分解炉最后一个单元 → 回转窑第一个单元
        if calciner_cells and kiln_cells: #判断分解炉和回转窑均有单元
            calciner_last = calciner_cells[-1] #获取分解炉最后一个单元
            kiln_first = kiln_cells[0] #获取回转窑第一个单元
            # 传递「反应后温度 + 浓度」
            self._transfer_between_cells(calciner_last, kiln_first)
        
        # 3. 回转窑最后一个单元 → 作为出料值
        if kiln_cells: #判断回转窑有单元
            kiln_last = kiln_cells[-1] #获取回转窑最后一个单元
            # 这里可以添加出料值的处理逻辑
            # 例如：存储到结果变量中或输出到文件
            pass
    
#五、核心变量获取方法
    def _get_core_variables(self):
        # 获取全局核心微分变量用于收敛判断
        # 温度、主物料浓度、内能
        core_vars = { #构建核心变量字典
            'gas_temperature': self.variables['algebraic_variables']['Tg'],
            'solid_temperature': self.variables['algebraic_variables']['Ts'],  # 区分气固温度
            'main_concentration': [cell['state_variables']['C']['CaCO3'] for cell in self.cells],  # 主物料浓度:以CaCO3为例
            'internal_energy': self.variables['state_variables']['U_hat'] #单位体积内能
        }
        return core_vars
  
#六、收敛判定方法  
    def check_convergence(self, prev_vars, current_vars):
        # 收敛判定:上一迭代核心变量\当前迭代核心变量
        # 所有核心变量迭代差值均小于收敛阈值时，判定系统达到稳态
        for var_name in prev_vars: #遍历所有核心变量类型
            for prev_val, curr_val in zip(prev_vars[var_name], current_vars[var_name]): #遍历各单元的变量值
                if abs(curr_val - prev_val) > self.convergence_threshold: #判断当前值与上一值的差值是否大于阈值
                    return False
        return True
    
#七、结果输出方法
    def output_results(self):
        # 结果整理与输出
        # 输出稳态下各单元的温度、浓度、反应转化率分布
        # 输出三段关键工艺指标（出口温度、产量、熟料成分占比）
        
        # 1. 找到回转窑的最后一个单元作为出料口
        kiln_cells = [cell for cell in self.cells if cell['device_type'] == 'kiln']
        if not kiln_cells:
            print("未找到回转窑单元")
            return
        
        # 出料口是回转窑的最后一个单元
        outlet_cell = kiln_cells[-1]
        outlet_index = outlet_cell['index']
        
        # 2. 获取出口温度（固体温度）
        outlet_temperature = outlet_cell['algebraic_variables']['Ts']
        
        # 3. 计算熟料组分占比
        # 获取固体组分浓度
        solid_components = self.constants['stoichiometric_matrix']['components'][:9]
        outlet_concentrations = outlet_cell['state_variables']['C']
        
        # 计算各熟料组分的摩尔数
        clinker_components = ['CaO', 'C2S', 'C3S', 'C3A', 'C4AF'] #熟料矿物组分列表
        clinker_moles = {} 
        for component in clinker_components: #循环遍历熟料组分
            clinker_moles[component] = outlet_concentrations.get(component, 0.0)#获取组分浓度并存入字典
        
        # 计算所有固体组分的总摩尔浓度（相对于所有固体组分）
        total_solid_moles = 0.0
        for component in solid_components:
            if component in outlet_concentrations:
                total_solid_moles += outlet_concentrations[component]
        
        # 如果total_solid_moles为0，则使用clinker_moles的总和作为备选
        if total_solid_moles == 0:
            total_solid_moles = sum(clinker_moles.values())
        
        # 计算各组分占比
        clinker_ratio = {}
        for component in clinker_components: #循环遍历熟料组分
            if total_solid_moles > 0:
                clinker_ratio[component] = (clinker_moles[component] / total_solid_moles) * 100
            else:
                clinker_ratio[component] = 0.0
        
        # 4. 打印结果
        print("=" * 60)
        print("出料口工艺指标")
        print("=" * 60)
        print(f"出口温度: {outlet_temperature:.2f} K")
        
        print("\n熟料及原料组分占比 (%):")
        for component in solid_components:
            ratio = clinker_ratio.get(component, 0)
            if component in ['CaO', 'C2S', 'C3S', 'C3A', 'C4AF']:
                print(f"{component} (熟料): {ratio:.2f}%")
            else:
                if component in outlet_concentrations:
                    conc = outlet_concentrations[component]
                    if total_solid_moles > 0:
                        ratio = (conc / total_solid_moles) * 100
                    else:
                        ratio = 0.0
                    print(f"{component} (原料): {ratio:.2f}%")
        print(f"\n总固体摩尔浓度: {total_solid_moles:.6f} mol/m³")
        print("=" * 60)

# 八、主运行模型
def main():
    # 创建模型实例
    model = CementProductionModel()
    
    # 全局参数初始化
    model.define_constants() #常数
    model.define_parameters() #参数
    model.define_control_variables() #控制变量

    
    # 三段工艺差异化参数化配置
    model.configure_sections()
    
    # 有限体积空间离散
    model.spatial_discretization()
    
    # 初始进料设置，确保第一个单元有正确的初始值
    model._feed_preheater_inlet()
    
    # 【调试】打印进料后第一个单元的浓度
    if model.cells:
        first_cell = model.cells[0]
        solid_components = model.constants['stoichiometric_matrix']['components'][:9]
    
    # 启动求解
    model.solve()

if __name__ == "__main__":
    main()