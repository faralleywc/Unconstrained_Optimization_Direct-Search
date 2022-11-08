from math import *
from sympy import *

class univariate_search():
    def __init__(self, function, start_point, threshold, times, step_range, gold_precision):
        self.function = sympify(function)
        self.start_point = start_point
        self.threshold = threshold
        self.times = times
        self.step_range = step_range
        self.gold_precision = gold_precision
        self.alpha = Symbol('alpha')
        self.dim = len(start_point)
    
    #  将符号化表达的x1、x2......xn赋以同名的变量，用于计算过程的取用
    def symbol_variable(self):
        X = []
        for i in range(self.dim):
            vars()['x%s'% (i+1)] = Symbol('x%s'% (i+1))
            X.append(symbols('x%s'% (i+1)))       
        return X

    # 将字符串x1、x2......xn存放在符号列表中
    def symbol_save(self):
        # 以符号表达方式存储自变量
        dim = len(self.start_point)
        symbol_x = []
        for i in range(dim):
            xi = 'x' + str(i+1)
            symbol_x.append(Symbol(xi))
        return symbol_x
    
    # 基于点距离准则和迭代次数准则执行算法迭代
    def algorithm_run(self):
        # 调用建立好的符号列表（[x1,x2,x3,...,xn]）
        symbol_x = self.symbol_save()
        # 设置点距初值和迭代次数初值
        time = 0
        # 迭代过程中的坐标，以起始点为初值
        x = self.start_point
        distance = float('inf')
        while (distance > self.threshold) and (time < self.times):
            # 保存上一代坐标值用于做终止准则判断
            x_last = x[:]
            # 选出第i个坐标进行轮换
            for i in range(self.dim):
                target = self.function
                # 只改变第i个坐标值，其余坐标保持不变
                for j in range(self.dim):
                    if j == i:
                        target = target.subs(symbol_x[j], x[j] + self.alpha)
                    elif j != i:
                        target = target.subs(symbol_x[j], x[j])
                
                #  根据黄金分割法确定最优步长并更新坐标值
                best_step = self.golden_ratio(target)
                x[i] = x[i] + best_step                

            # 完成一轮坐标更新后，判断是否达到终止准则
            res  = [i - j for i, j in zip(x_last, x)]
            error = 0
            for item in res:
                error += item ** 2
            distance = sqrt(error)
            time += 1
        
        return x

    #  计算在最优坐标下的最优目标函数值
    def target_calculate(self, x):
        X = self.symbol_variable()
        y = self.function.evalf(subs = dict(zip(X, x)))

        return y

    # 黄金分割法一维搜索
    def golden_ratio(self, target):
        a = -self.step_range
        b = self.step_range
        e = self.gold_precision
        ratio = float(GoldenRatio) - 1
        G1 = a + (1 - ratio) * (b - a)
        G2 = a + (b - a) * ratio
        while True:
            # 如果搜索区间长度在精度范围内，直接得到最优步长
            if (b - a) <= e:
                midpoint = (a + b) / 2
                best_step = midpoint
                break
            # 如果搜索区间长度在精度范围外，黄金分割缩小搜索区间
            elif (b - a) > e:
                # G2分割点更优，保留b点并使G2替代G1，构造新的搜索区间
                # print(target.subs(self.alpha, G1))
                if target.subs(self.alpha, G1) >= target.subs(self.alpha, G2):
                    a = G1
                    G1 = G2
                    G2 = a + ratio * (b - a)
                # G1分割点更优，保留a点并使G1替代G2，构造新的搜索区间
                elif target.subs(self.alpha, G1) <= target.subs(self.alpha, G2):
                    b = G2
                    G2 = G1
                    G1 = a + (1 - ratio) * (b - a)

        return best_step
        
if __name__ == "__main__":
    # 给定原始条件
    function = 'x1 ** 2 + x2 ** 2 - x1 * x2 - 10 * x1 - 4 * x2 + 60'
    start_point = [0, 0]
    threshold = 1e-1
    times = 100
    step_range = 1
    gold_precision = 1e-5
    # 得到结果
    myalgorithm = univariate_search(function, start_point, threshold, times, step_range, gold_precision)
    best_x = myalgorithm.algorithm_run()
    best_y = myalgorithm.target_calulate(best_x)
