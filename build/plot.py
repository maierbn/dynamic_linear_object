#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

filename = "Fall3Reibung.csv"
if len(sys.argv) >= 2:
  filename = sys.argv[1]
  
# load data
with open(filename,"r") as f:
  data = f.readlines()
  
def points(coefficients, parameters, versch, theta0):
  def N(i, h, s):
    if (i-1)*h < s < i*h:
      return (s - (i - 1) * h) / h
    
    elif i*h <= s < (i+1)*h:
      return ((i + 1) * h - s) / h
      
    return 0
    
  def Theta(h,s):
    i1 = (int)(s / h)
    i2 = i1+1
    
    c1 = 0
    c2 = 0
    
    if 0 <= i1 < len(coefficients):
      c1 = coefficients[i1]
    if 0 <= i2 < len(coefficients):
      c2 = coefficients[i2]
      
    N1 = N(i1, h, s)
    N2 = N(i2, h, s)
    
    return c1*N1 + c2*N2
    
  [n, Rflex, L, rho, mu, simpson_iter, ver, reib, winkelkontrolle, t00, t01, t02, vx1, vy1, vx2, vy2] = map((float),parameters)
  h = L/n
  
  n_steps = 1000
  ds = L/n_steps
  x = versch[0]
  y = versch[1]
  points = []
  points.append([x,y])
  
  for i in range(n_steps+1):
    s = i*ds
    theta = Theta(h,s)
    x += np.cos(theta)*ds
    y += np.sin(theta)*ds
    if i % 4 == 0:
      points.append([x,y])
  
  return points
  
parameters = []
t_list = []
versch_list = []
theta0_list = []
coefficients_list = []
points_list = []

for i,line in enumerate(data):
  values = line.split(",")
  if "\n" in values[-1]:
    values[-1] = values[-1][0:values[-1].index("\n")]
    
  if i == 0:
    parameters = values
    [n, Rflex, L, rho, mu, simpson_iter, ver, reib, winkelkontrolle, t00, t01, t02, vx1, vy1, vx2, vy2] = map(float,parameters)
    description = "parameters: \n n: {}, Rflex: {}, L: {}, rho: {}, mu: {}, \n simpson_iter: {}, verschiebung: {}, reib: {}, winkelkontrolle: {}, \n winkel:       t02 * t^2 + t01 * t + t00  mit  t00: {}, t01: {}, t02: {},\n verschiebung: v2 * t^2 + v1 * t          mit vx1: {}, vy1: {}, vx2: {}, vy2: {}".format(n, Rflex, L, rho, mu, simpson_iter, ver, reib, winkelkontrolle, t00, t01, t02, vx1, vy1, vx2, vy2)
    continue
      
  elif len(values) == 1:
    print("duration of simulation: {}s".format(values[0]))
    break  
  
  t = (float)(values[0])
  versch = map(float,values[1:3])
  theta0 = (float)(values[3])
  coefficients = map(float,values[4:])
  
  t_list.append(t)
  versch_list.append(versch)
  theta0_list.append(theta0)
  coefficients_list.append(coefficients)
  
  #print(coefficients)
  
  points_list.append(points(coefficients, parameters, versch, theta0))
  
show_plots = True
#if len(sys.argv) > 1:
#  show_plots = False;

fig = plt.figure()
ax = plt.axes()
#ax.grid(which='both')

ax.plot([x for [x,y] in versch_list], [y for [x,y] in versch_list], "o-") 

for i,t in enumerate(t_list):
  if i % 10 == 0:
    ax.plot([x for [x,y] in points_list[i]], [y for [x,y] in points_list[i]], color=(0.5,0.5,0.5), label=t)

plt.title(description)
plt.axes().set_aspect('equal', 'datalim')
ax.set_xlabel("x")
ax.set_ylabel("y")
#plt.show()

# animate
def update_line(i, line1, t, n):
  xlist = [x for [x,y] in points_list[i]]
  ylist = [y for [x,y] in points_list[i]]
  line1.set_data(xlist, ylist)
  
  t.set_text('t='+str(t_list[i]))
  
  return [line1,t]

fig1, ax1 = plt.subplots()

line1, = ax1.plot([], [], color = "r")
ax1.set_xlim(ax.get_xlim())
ax1.set_ylim(ax.get_ylim())
ax1.set_aspect('equal')
title = ax1.text(0.1,ax1.get_ylim()[1]*0.8,"t")
line_ani = animation.FuncAnimation(fig1, update_line, len(points_list), fargs=[line1,title,len(points_list)],
                                   interval=2e3/len(points_list), blit=True)
plt.show()
