from numericaltools import integrate

x = 1

def func(x):
    y = 2*x + 1
print(integrate(func, 0, 6, 7))