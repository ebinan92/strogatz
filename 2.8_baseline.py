import math
import matplotlib.pyplot as plt
import numpy as np

DELTA_T = 0.5
MAX_T = 5.0

t0 = 0.0
x0 = 1.0 #x(0)

def f(x):
    return -x

def getX1byEuler(t,x):
    check = 1.0/DELTA_T
    if not check.is_integer():
        return print("select DELTA_T so that 1.0 / DELTA_T is an integer")
    while t < 1.0:
        x = f(x)*DELTA_T + x
        t += DELTA_T
        #print(x,t)
    return print("Euler:x(1)={}".format(x))

def getX1byRunge(t,x):
    check = 1.0/DELTA_T
    if not check.is_integer():
        return print("select DELTA_T so that 1.0 / DELTA_T is an integer")
    while t < 1.0:
        k1 = f(x)*DELTA_T
        k2 = f(x+1/2*k1)*DELTA_T
        k3 = f(x+1/2*k2)*DELTA_T
        k4 = f(x+k3)*DELTA_T
        x = x + 1/6*(k1+2*k2+2*k3+k4)
        t += DELTA_T
        #print(x,t)
    return print("Runge_Kutta:x(1)={}".format(x))

#Euler
def Euler(t,x):
    x_hist = []
    t_hist = []
    while t < MAX_T:
        x = f(x)*DELTA_T + x
        t += DELTA_T
        x_hist.append(x)
        t_hist.append(t)
        #print(x,t)
    return plt.plot(t_hist, x_hist,color="skyblue",label="euler")

#FixEuler
def FixEuler(t,x):
    x_hist = []
    t_hist = []
    while t < MAX_T:
        x_tilda = f(x)*DELTA_T + x
        x = x + 1/2*(f(x)+f(x_tilda))*DELTA_T
        t += DELTA_T
        x_hist.append(x)
        t_hist.append(t)
        #print(x,t)
    return plt.plot(t_hist, x_hist,color="green",label="fixeuler")

#E = |x_tilda(t)-x(t)|
def NumError(t,x):
    t_hist = []
    E_hist = []
    while t < MAX_T:
        E = abs(f(x)*DELTA_T +x - pow(math.e,-t))
        x = f(x)*DELTA_T + x
        t += DELTA_T
        E_hist.append(E)
        t_hist.append(t)
        #print(x,t)

    '''log Error'''
    #E_hist = [math.log(e) for e in E_hist]
    #t_hist = [math.log(e) for e in t_hist]
    return plt.plot(t_hist, E_hist,color="red",label="error")

# Runge_kutta
def Runge_kutta(t,x):
    x_hist = []
    t_hist = []
    while t < MAX_T:
        k1 = f(x)*DELTA_T
        k2 = f(x+1/2*k1)*DELTA_T
        k3 = f(x+1/2*k2)*DELTA_T
        k4 = f(x+k3)*DELTA_T
        x = x + 1/6*(k1+2*k2+2*k3+k4)
        t += DELTA_T
        x_hist.append(x)
        t_hist.append(t)
        #print(x,t)
    return plt.plot(t_hist, x_hist,color="black",label="runge_kutta")


#Exact solution
def ExactSolution():
    t = np.linspace(0,MAX_T, num=100)
    x = pow(math.e,-t)
    return plt.plot(t, x,color="orange",label="exact")

'''execute'''
getX1byEuler(t0,x0)
getX1byRunge(t0,x0)

plt.subplot(2, 1, 2)
NumError(t0,x0)
plt.xlabel("t")
plt.ylabel("E(t)")
plt.legend()

plt.subplot(2, 1, 1)
Euler(t0,x0)
FixEuler(t0,x0)
Runge_kutta(t0,x0)
ExactSolution()
plt.xlabel("t")
plt.ylabel("x(t)")
plt.legend()
plt.show()