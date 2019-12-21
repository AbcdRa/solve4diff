# -*- coding: utf-8 -*-
'''
Created on 9 dec. 2019 y.

@author: Michail
'''
import numpy as np
import scipy.optimize
import scipy.constants
import math
import sys

def equations(p): 
#E, k0 = p
    E = p[0]
    k0 = p[1]
    eq = [0.,0.,0.,0.,0.];
    for i in range(0,5):
        eq[i] = np.log(betta[i]/(T_max[i]-T0)) + (E/(scipy.constants.R*T_max[i])) - np.log(k0)
    return [eq[0], eq[1], eq[2], eq[3], eq[4]]

def equations1(p):
    if(type(p) == np.float64):
        n = p
    else:
        n = p[0]
    eq = [0.,0.,0.,0.,0.];
    for i in range(0,5):
        eq[i] = v[i] - (1 - x[i])*n*(k_avg*tetta[i])**(n-1)*(k_avg + tetta[i]*k_dot_avg)
    return [eq[0], eq[1], eq[2], eq[3], eq[4]]
    

num = 1
while(True):
    print('Запуск', num)
    #Начальные данные
    betta = [0.0,0.,0.,0.,0.];
    T_max = [0.,0.,0.,0.,0.];
    v = [0.,0.,0.,0.,0.];
    tetta = [0.,0.,0.,0.,0.];
    x = [0.,0.,0.,0.,0.];
    T0 = 296.0 #К
    
    
    for i in range(0,5):
        betta[i] = input('Введите betta для ' + str(i+1) + ' скорости\n')
        betta[i] = float(betta[i])

    for i in range(0,5):
        T_max[i] = input('Введите T_max для ' + str(i+1) +' скорости\n')
        T_max[i] = float(T_max[i])

    for i in range(0,5):   
        v[i] = input("Введите " +str(i+1) + " скорость x'\n")
        v[i] = float(v[i])

    for i in range(0,5):
        tetta[i]  = input ('Введите ' + str(i+1) + ' время tetta\n')
        tetta[i]  = float(tetta[i] )

    for i in range(0,5):   
        x[i]  = input ('Введите ' + str(i+1) + ' степень превращения x\n')
        x[i]  = float(x[i])
    
    print('\n' )
    print('-------------------------------------------') 
    
    result = scipy.optimize.root(equations, (1e5, 1e18, 0, 0, 0), method = "lm")
    E = result.get('x')[0]
    k0 = result.get('x')[1]
    
#    E = float(input("Enter E = "))
#    k0 = float(input("Enter k0 = "))

    print('E =', '{:e}'.format(E))
    print('k0 = ', '{:e}'.format(k0))
    print('Проверка:', equations((E, k0)), '\n') 
    error = np.linalg.norm( np.array( equations((E,k0)) ) )
    print("Средняя ошибка: " + str(error) +'\n' )
    
    k = [0.0,0.,0.,0.,0.]
    
    for i in range(0,5):
        k[i] = k0*np.exp( -E/(scipy.constants.R*T_max[i]) )
        
    k = np.array(k)
    k_avg = np.mean(k)
    print("k = " + str(k))
    print("k_avg = " + str(k_avg))
    error = np.linalg.norm(k-k_avg)
    print("Среднеквадратичная ошибка : " + str( error ) )
    print("DEBUG Относительная ошибка: " + str(error/k_avg*20) + "%\n")
        
    k_dot = [0.0,0.,0.,0.,0.]
    for i in range(0,5):
        k_dot[i] = k_avg*betta[i]*E/(scipy.constants.R*(T_max[i]**2)) 
    
    k_dot = np.array(k_dot)
    k_dot_avg = np.mean(k_dot)
    print("k' = " + str(k_dot))
    print("k'_avg= " + str(k_dot_avg))
    error = np.linalg.norm(k_dot-k_dot_avg) 
    print("Среднеквадратичная ошибка : " + str(error) )
    print("DEBUG Относительная ошибка: " + str(error/k_dot_avg*20) + "%\n")
    
    result = scipy.optimize.root(equations1, (0, 0, 0, 0, 0), method = "lm")
    n = result.get('x')[0]
    print("n = " + str(n) +'\n')
    print('Проверка:', equations1(n)) 
    error = np.linalg.norm( np.array( equations1(n) ))
    print("Средняя ошибка: ", error )
    print("DEBUG Относительная ошибка: " + str(error/n*20) + "%\n")
    num = num + 1
    print('-------------------------------------------')
    print('\n', '\n' )
