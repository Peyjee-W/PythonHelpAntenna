from sympy import symbols, diff, sin, cos, exp, I, pi

# 定义符号变量
r, theta, k, I_L, eta = symbols('r theta k I_L eta', real=True)
j = I  # j 是虚数单位

# E_r 分量
E_r = eta * I_L / (2 * pi * r**2) * cos(theta) * (1 + 1 / (j * k * r)) * exp(-j * k * r)

# E_theta 分量
E_theta = j * eta * k * I_L / (4 * pi * r) * sin(theta) * (1 + 1 / (j * k * r) - 1 / (k**2 * r**2)) * exp(-j * k * r)

# 计算散度的每一项
# 1. 计算径向分量部分的散度项
div_E_r = (1 / r**2) * diff(r**2 * E_r, r)

# 2. 计算纬向分量部分的散度项
div_E_theta = (1 / (r * sin(theta))) * diff(E_theta * sin(theta), theta)

# 散度 ∇⋅E 的总和
div_E = div_E_r + div_E_theta

# 展开并简化
div_E_simplified = div_E.simplify()

# 在终端输出散度的结果
print("电场的散度为:", div_E_simplified)
