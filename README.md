# 无约束优化算法——坐标轮换法、鲍威尔法(直接求解法)

“坐标轮换法（Univariate_Search）”

坐标轮换法是一种沿坐标方向轮流进行搜索的寻优方法，即将多个变量的寻优问题轮流转化为单变量的优化问题。每轮寻优仅改变一个坐标维度上的坐标值（变量），依次循环直至达到优化目标或终止算法，因此这种方法又称作变量轮换法。在该方法进行搜索的过程中，由于设计变量的更新方向是确定的，因此无需计算目标函数的导数，只需要根据每个设计变量的取值计算目标函数的数值信息。

“鲍威尔法（Powell）”

作为直接求解法的一种，鲍威尔法直接利用函数值来构造共轭方向进行最优搜索，因此鲍威尔法又称为共轭方向法。通过构造共轭向量系来达到最快接近全局最优点的目的，与坐标轮换法相似，该方法不需要对目标函数的导数进行计算，只需要不断更新坐标值并计算目标函数值信息。

本项目的目标函数给定基于Sympy符号计算库，对于多为自变量的给定规则，约定为x1、x2、x3......xn，目标函数计算法则依据python计算符号规定。

'Univariate_Search'、'Powell'分别给出了坐标轮换法类模板和鲍威尔法类模板，以及主函数中包含了示例的算例，该所搜算法类的最优步长确定依据黄金分割法进行一维最优搜索。基于这些模板可以对不同的目标函数在不同的维度问题上进行求解，并且包含了多种终止准则（点距准则、目标值准则、迭代次数准则）。
