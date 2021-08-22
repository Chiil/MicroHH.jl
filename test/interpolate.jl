using Interpolations
using Tullio

## 
x = 0:0.1:1 |> collect
f = exp.(x)

interp_lin = LinearInterpolation(x, f)

x_new = 0:(1//3):1 |> collect
f_new = interp_lin(x_new)

## 
y = copy(x)
g = zeros(size(x, 1), size(y, 1))
@tullio g[i, j] = exp(x[i])

interp_lin_2d = LinearInterpolation((x, y), g)

y_new = copy(x_new)
g_new = interp_lin_2d(x_new, y_new)

