# -*- coding: utf-8 -*-
"""
Daisyworld Model

@author: Sivanessha Kalithasan
"""
from pylab import *
import numpy as np
import matplotlib.pylab as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import scipy.optimize as optimize


#Set up of parameters
T_opt = 22.5
gamma = 0.3
sigma = 5.6704*(10**(-8))
S = 917
q = 20
A_b = 0.25
A_x = 0.5
A_w = 0.75
L = 1.2

#Intial Condition 
y0 = np.array([0.1,0.1])
tmax = 101
dt = 1
tspan = np.arange(0,tmax,dt)

#planetary albedo
def A_p(alpha_b, alpha_w):
    return (alpha_w * A_w + alpha_b * A_b + (1-alpha_w-alpha_b)*A_x) 

#planetary temperature
def T_e(alpha_b, alpha_w): 
    z1 = S*L*(1-A_p(alpha_b, alpha_w))
    z2 = (z1/sigma)**(1/4) - 273.15
    return z2

#local temperature of white daisy
def T_w(alpha_b, alpha_w):
    z1 = q*(A_p(alpha_b, alpha_w)-A_w) + T_e(alpha_b, alpha_w)
    return z1
    

#local temperature of black daisy
def T_b(alpha_b, alpha_w):
    z1 = q*(A_p(alpha_b, alpha_w)-A_b) + T_e(alpha_b, alpha_w)
    return z1     

#growth rates of white daisy
def Beta_w(alpha_b, alpha_w):
    z1 = np.max([0, 1-(((T_opt - T_w(alpha_b, alpha_w))/17.5)**2)])  
    return z1

#growth rates of black daisy
def Beta_b(alpha_b, alpha_w):
    z1 = np.max([0, 1-(((T_opt - T_b(alpha_b, alpha_w))/17.5)**2)])  
    return z1

#growth rate of area of white daisy
def dy1dt(alpha_b, alpha_w):
    z1 = alpha_w*( (1-alpha_b-alpha_w)*Beta_w(alpha_b, alpha_w)  - gamma  )
    return z1

#growth rate of area of black daisy
def dy2dt(alpha_b, alpha_w):
    z1 = alpha_b*( (1-alpha_b-alpha_w)*Beta_b(alpha_b, alpha_w)  - gamma  )
    return z1

#function to calculate the growth rate of area of black and white daisy
def func(t,y): 
    dydt = np.array([dy1dt(y[1], y[0]), dy2dt(y[1], y[0])  ])
    return dydt

#solving for the differential equation(growth rate) using solva_ivp function
sol = solve_ivp(func, [tspan[0], tspan[-1]], y0, t_eval=tspan, method = 'RK45')

figure()
plot(sol.t, sol.y[0],label="white")
plot(sol.t, sol.y[1],label="black",color="black")
plot(sol.t, 1 - sol.y[0]-sol.y[1], label="bare ground")
xlabel("time")
ylabel("Area covered")
title("Area covered by organisms against time at Luminosity = 1.2")
legend()
savefig("figdata1.png")
show()


figure()
plot(sol.t, T_e(sol.y[1], sol.y[0]))
xlabel("time")
ylabel("Temperature")
title("Global Temperature against time at Luminosity = 1.2")
savefig("figdata2.png")
show()








"""Part 2 varied Luminosity"""
l = np.arange(0.60 , 1.42, 0.02)#varying the luminosity
solution_w = np.array([])
solution_b = np.array([])

#Calculating solution for varying luminosity
for i in range(len(l)):
    L = l[i]
    solution = solve_ivp(func, [tspan[0], tspan[-1]], y0, method = 'RK45')
    solution_w = r_[solution_w,solution.y[0][-1]]
    solution_b = r_[solution_b,solution.y[1][-1]]
   
figure()
plot(l, solution_w, label="white")
plot(l, solution_b, label="black", color="black")
plot(l, 1 - solution_w-solution_b, label="empty")
xlabel("Luminosity")
ylabel("Area covered")
legend()
title("Area covered by organisms at varied Luminosity")
savefig("figdata3.png")
show()

#planetary temperature with L as a parameter
def T_e_L(alpha_b, alpha_w, L): 
    z1 = S*L*(1-A_p(alpha_b, alpha_w))
    z2 = (z1/sigma)**(1/4) - 273.15
    return z2

figure()
plot(l, T_e_L(solution_b, solution_w, l), label="Global Temperature")
plot(l, T_e_L(np.zeros_like(l), np.zeros_like(l), l), ".b", label="Global Tmeperature without any organisms")
xlabel("Luminosity")
ylabel("Temperature")
legend()
title("Temperature at varied Luminosity")
savefig("figdata4.png")
show()








"""part 3"""
"""Periodic luminosity"""
l = 0.8*(np.sin(1/10*tspan))**2 + 0.6 #periodic luminosity function
figure()
plot(tspan,l)
xlabel("Time")
ylabel("Luminosity")
title("Luminosity aganist time")
savefig("figdata5.png")
show()

solution_w = np.array([])
solution_b = np.array([])

for i in range(len(l)):
    L = l[i]
    solution = solve_ivp(func, [tspan[0], tspan[-1]], y0, method = 'RK45')
    solution_w = r_[solution_w,solution.y[0][-1]]
    solution_b = r_[solution_b,solution.y[1][-1]]

figure()
plot(tspan, solution_w, label="white")
plot(tspan, solution_b, label="black", color="black")
plot(tspan, 1 - solution_w-solution_b, label="empty")
xlabel("time")
ylabel("Area covered")
legend()
title("Area covered by organisms at periodic Luminosity")
savefig("figdata6.png")
show()

figure()
plot(tspan, T_e_L(solution_b, solution_w, l), label="Global Temperature")
plot(tspan, T_e_L(np.zeros_like(l), np.zeros_like(l), l), ".b", label="Global Tmeperature without any organisms")
xlabel("time")
ylabel("Temperature")
legend()
title("Temperature at periodic Luminosity")
savefig("figdata7.png")
show()






"""More organisms(4 organims)"""

#Intial Condition
A_1 = 0.35
A_2 = 0.65
y0 = np.array([0.1,0.1,0.1,0.1])
tmax = 501
dt = 1
tspan = np.arange(0,tmax,dt)

#planetary albedo
def A_p(alpha_b, alpha_w, alpha_1, alpha_2):
    return (alpha_w * A_w + alpha_b * A_b + alpha_1*A_1 + alpha_2*A_2 +(1-alpha_w-alpha_b-alpha_1-alpha_2)*A_x) 

#planetary temperature
def T_e(alpha_b, alpha_w, alpha_1, alpha_2): 
    z1 = S*L*(1-A_p(alpha_b, alpha_w, alpha_1, alpha_2))
    z2 = (z1/sigma)**(1/4) - 273.15
    return z2

#local temperature of white daisy
def T_w(alpha_b, alpha_w, alpha_1, alpha_2):
    z1 = q*(A_p(alpha_b, alpha_w, alpha_1, alpha_2)-A_w) + T_e(alpha_b, alpha_w, alpha_1, alpha_2)
    return z1
    

#local temperature of black daisy
def T_b(alpha_b, alpha_w,alpha_1, alpha_2):
    z1 = q*(A_p(alpha_b, alpha_w,alpha_1, alpha_2)-A_b) + T_e(alpha_b, alpha_w,alpha_1, alpha_2)
    return z1     

#local temperature of daisy_1
def T_1(alpha_b, alpha_w, alpha_1, alpha_2):
    z1 = q*(A_p(alpha_b, alpha_w, alpha_1, alpha_2)-A_1) + T_e(alpha_b, alpha_w, alpha_1, alpha_2)
    return z1

#local temperature of white daisy_2
def T_2(alpha_b, alpha_w, alpha_1, alpha_2):
    z1 = q*(A_p(alpha_b, alpha_w, alpha_1, alpha_2)-A_2) + T_e(alpha_b, alpha_w, alpha_1, alpha_2)
    return z1

#growth rates of white daisy
def Beta_w(alpha_b, alpha_w,alpha_1, alpha_2):
    z1 = np.max([0, 1-(((T_opt - T_w(alpha_b, alpha_w,alpha_1, alpha_2))/17.5)**2)])  
    return z1

#growth rates of black daisy
def Beta_b(alpha_b, alpha_w,alpha_1, alpha_2):
    z1 = np.max([0, 1-(((T_opt - T_b(alpha_b, alpha_w,alpha_1, alpha_2))/17.5)**2)])  
    return z1

#growth rates of daisy_1
def Beta_1(alpha_b, alpha_w,alpha_1, alpha_2):
    z1 = np.max([0, 1-(((T_opt - T_1(alpha_b, alpha_w,alpha_1, alpha_2))/17.5)**2)])  
    return z1

#growth rates of daisy_2
def Beta_2(alpha_b, alpha_w,alpha_1, alpha_2):
    z1 = np.max([0, 1-(((T_opt - T_2(alpha_b, alpha_w,alpha_1, alpha_2))/17.5)**2)])  
    return z1

#growth rate of area of white daisy
def dy1dt(alpha_b, alpha_w,alpha_1, alpha_2):
    z1 = alpha_1*( (1-alpha_b-alpha_w-alpha_1-alpha_2)*Beta_1(alpha_b, alpha_w,alpha_1, alpha_2)  - gamma  )
    return z1

#growth rate of area of black daisy
def dy2dt(alpha_b, alpha_w,alpha_1, alpha_2):
    z1 = alpha_2*( (1-alpha_b-alpha_w-alpha_1-alpha_2)*Beta_2(alpha_b, alpha_w,alpha_1, alpha_2)  - gamma  )
    return z1

#growth rate of area of white daisy
def dywdt(alpha_b, alpha_w,alpha_1, alpha_2):
    z1 = alpha_w*( (1-alpha_b-alpha_w-alpha_1-alpha_2)*Beta_w(alpha_b, alpha_w,alpha_1, alpha_2)  - gamma  )
    return z1

#growth rate of area of black daisy
def dybdt(alpha_b, alpha_w,alpha_1, alpha_2):
    z1 = alpha_b*( (1-alpha_b-alpha_w-alpha_1-alpha_2)*Beta_b(alpha_b, alpha_w,alpha_1, alpha_2)  - gamma  )
    return z1

#function to calculate the growth rate of area of black and white daisy
def func(t,y): 
    dydt = np.array([dywdt(y[1], y[0], y[2], y[3]), dybdt(y[1], y[0], y[2], y[3]), dy1dt(y[1], y[0], y[2], y[3]), dy2dt(y[1], y[0], y[2], y[3]) ])
    return dydt

#solving for the differential equation(growth rate) using solva_ivp function
sol = solve_ivp(func, [tspan[0], tspan[-1]], y0, t_eval=tspan, method = 'RK45')

figure()
plot(sol.t, sol.y[0],label="white")
plot(sol.t, sol.y[1],label="black",color="black")
plot(sol.t, sol.y[2],label="daisy 1")
plot(sol.t, sol.y[3],label="daisy 2")
plot(sol.t, 1 - sol.y[0]-sol.y[1]-sol.y[2]-sol.y[3], label="bare ground")
xlabel("time")
ylabel("Area covered")
title("Area covered by organisms against time at Luminosity = 1.2")
legend()
savefig("figdata8.png")
show()


figure()
plot(sol.t, T_e(sol.y[1], sol.y[0],sol.y[2], sol.y[3]) )
xlabel("time")
ylabel("Temperature")
title("Global Temperature against time at Luminosity = 1.2")
savefig("figdata9.png")
show()




"""Part 3 varied Luminosity with more organisms"""
l = np.arange(0.60 , 1.42, 0.02)#varying the luminosity
solution_w = np.array([])
solution_b = np.array([])
solution_1 = np.array([])
solution_2 = np.array([])

#Calculating solution for varying luminosity
for i in range(len(l)):
    L = l[i]
    solution = solve_ivp(func, [tspan[0], tspan[-1]], y0, method = 'RK45')
    solution_w = r_[solution_w,solution.y[0][-1]]
    solution_b = r_[solution_b,solution.y[1][-1]]
    solution_1 = r_[solution_1,solution.y[2][-1]]
    solution_2 = r_[solution_2,solution.y[3][-1]]
   
figure()
plot(l, solution_w, label="white")
plot(l, solution_b, label="black", color="black")
plot(l, solution_1, label="daisy 1")
plot(l, solution_2, label="daisy 2")
plot(l, 1 - solution_w-solution_b-solution_1-solution_2, label="empty")
xlabel("Luminosity")
ylabel("Area covered")
legend()
title("Area covered by organisms at varied Luminosity")
savefig("figdata10.png")
show()

#planetary temperature with L as a parameter
def T_e_L(alpha_b, alpha_w, alpha_1, alpha_2, L): 
    z1 = S*L*(1-A_p(alpha_b, alpha_w, alpha_1, alpha_2))
    z2 = (z1/sigma)**(1/4) - 273.15
    return z2

figure()
plot(l, T_e_L(solution_b, solution_w,solution_1,solution_2, l), label="Global Temperature")
plot(l, T_e_L(np.zeros_like(l), np.zeros_like(l),np.zeros_like(l),np.zeros_like(l), l), ".b", label="Global Tmeperature without any organisms")
xlabel("Luminosity")
ylabel("Temperature")
legend()
title("Temperature at varied Luminosity")
savefig("figdata11.png")
show()




