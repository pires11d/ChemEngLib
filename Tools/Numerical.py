"""Module that contains all Numerical Methods 
i.e. derivation, integration and root finding algorithms."""


def dfdx(f, xo, dx):
    df = (f(xo+dx)-f(xo-dx))/dx
    return df

def dfdx2(f, xo, dx):
    df2 = (dfdx(f,xo+dx,dx)-dfdx(f,xo-dx,dx))/dx
    return df2

def Newton(f, xo, dx, error=1e-3):
    error1 = 1000.0
    x0 = xo
    while error1 > error:
        x1 = x0 - f(x0)/dfdx(f,x0,dx)
        error1 = abs(x1 - x0)
        x0 = x1
    return x1

def qNewton(f, xo, dx, error=1e-3):
    error1 = 1000.0
    x0 = xo
    while error1 > error:
        x1 = x0 - dfdx(f, x0, dx) / dfdx2(f, x0, dx)
        error1 = abs(x1 - x0)
        x0 = x1
    return x1

def Trapezoidal(f, x_list):
    dx = (x_list[-1] - x_list[0]) / (len(x_list) - 1)
    acc = 0
    for i, x in enumerate(x_list):
        if i < len(x_list)-1:
            acc += (f(x)+f(x_list[i+1])) * dx / 2
    return acc

def Trapezoidal_list(f, x_list):
    dx = (x_list[-1] - x_list[0]) / (len(x_list) - 1)
    acc = 0
    acc_list = []
    for i, x in enumerate(x_list):
        if i < len(x_list)-1:
            acc += (f(x)+f(x_list[i+1])) * dx / 2
            acc_list.append(acc)
    return acc_list

def Simpson(f, x_list):
    dx = (x_list[-1] - x_list[0]) / (len(x_list) - 1)
    acc = 0
    for i, x in enumerate(x_list):
        if i < len(x_list) - 1:
            acc += (f(x)+4*f((x+x_list[i+1])/2)+f(x_list[i+1])) * dx / 6
    return acc

def Simpson_list(f, x_list):
    dx = (x_list[-1] - x_list[0]) / (len(x_list) - 1)
    acc = 0
    acc_list = []
    for i, x in enumerate(x_list):
        if i < len(x_list) - 1:
            acc += (f(x)+4*f((x+x_list[i+1])/2)+f(x_list[i+1])) * dx / 6
            acc_list.append(acc)
    return acc_list

def Euler(dydt, yo, t_list):
    F = dydt
    dt = (t_list[-1] - t_list[0]) / (len(t_list) - 1)
    y_list = []
    yi = yo
    for i, t in enumerate(t_list):
        y_list.append(yi)
        if i < len(t_list)-1:
            yi = yi + dt*(F(t, yi))
    return y_list

def RungeKutta(dydt, yo, t_list):
    F = dydt
    dt = (t_list[-1] - t_list[0]) / (len(t_list) - 1)
    y_list = []
    yi = yo
    for i, t in enumerate(t_list):
        y_list.append(yi)
        if i < len(t_list)-1:
            k1 = F(t, yi)
            k2 = F(t+dt/2,yi+k1*dt/2)
            k3 = F(t+dt/2,yi+k2*dt/2)
            k4 = F(t+dt/2,yi+k3*dt)
            yi = yi + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
    return y_list