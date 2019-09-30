#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Visualize simulation results.
# Some settings such as figure size and stride of timesteps are hardcoded and can be adjusted depending on the scenario
#
# usage:
# ./plot.py <file.csv>

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# parse filename from command line argument
filename = "experiment3.csv"
if len(sys.argv) >= 2:
  filename = sys.argv[1]
  
# load data
with open(filename,"r") as f:
  data = f.readlines()
  
# get all the points of the state of the one-dimensional object for one time point, the points make up the line of the object
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
    
  [n, Rflex, L, rho, mu, ver, reib, winkelkontrolle, t00, t01, t02, vx1, vy1, vx2, vy2] = map((float),parameters)
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

# loop over data (each line is one simulation result)
for i,line in enumerate(data):
  values = line.split(",")
  
  # parse values in line
  if "\n" in values[-1]:
    values[-1] = values[-1][0:values[-1].index("\n")]
  
  # in first line, parse paramreters  
  if i == 0:
    parameters = values
    [n, Rflex, L, rho, mu, ver, reib, winkelkontrolle, t00, t01, t02, vx1, vy1, vx2, vy2] = map(float,parameters)
    description = "parameters: \n n: {}, Rflex: {}, L: {}, rho: {}, mu: {}, \n verschiebung: {}, reib: {}, winkelkontrolle: {}, \n winkel:       t02 * t^2 + t01 * t + t00  mit  t00: {}, t01: {}, t02: {},\n verschiebung: v2 * t^2 + v1 * t          mit vx1: {}, vy1: {}, vx2: {}, vy2: {}".format(n, Rflex, L, rho, mu, ver, reib, winkelkontrolle, t00, t01, t02, vx1, vy1, vx2, vy2)
    continue
  
  # if line contains only one entry, it is the last line
  elif len(values) == 1:
    break  
  
  # parse values and append to lists
  t = (float)(values[0])
  versch = map(float,values[1:3])
  theta0 = (float)(values[3])
  coefficients = map(float,values[4:])
  
  t_list.append(t)
  versch_list.append(versch)
  theta0_list.append(theta0)
  coefficients_list.append(coefficients)
  
  # append points of object line to list
  points_list.append(points(coefficients, parameters, versch, theta0))
  
show_plots = True
#if len(sys.argv) > 1:
#  show_plots = False;

# -----------------------------------------------
# create static plot that shows states of object
plt.rcParams.update({'font.size': 14})
plt.rcParams['lines.linewidth'] = 2

#fig = plt.figure(figsize=(10,3))  # exp2
#fig = plt.figure(figsize=(2,4))   # exp1
fig = plt.figure(figsize=(4,4))    # white2, tube
plt.tight_layout()
ax = plt.axes()
#ax.grid(which='both')

t_stride = 100   # white2
t_stride = 50   # tube
#t_stride = 10   # exp2
#t_stride = 20  # exp1

#ax.plot([-1.2,1.2],[0,0],"o",color="w")  # exp1
#ax.plot([0,0],[-1.2,1.2],"o",color="w")  # exp2

# start point
ax.plot([x for [x,y] in versch_list[::t_stride]], [y for [x,y] in versch_list[::t_stride]], "o-")

for i,t in enumerate(t_list):
  if i % t_stride == 0:
    #ax.plot([x for [x,y] in points_list[i]], [y for [x,y] in points_list[i]], color=(0.5,0.5,0.5), label=t)
    ax.plot([x for [x,y] in points_list[i]], [y for [x,y] in points_list[i]], label=t)

#ax.plot([x for [x,y] in points_list[0]], [y for [x,y] in points_list[0]], color='r', label=t)


# end points
ax.plot([points_list[i][-1][0] for i in range(len(points_list))], [points_list[i][-1][1] for i in range(len(points_list))], "--", color='grey')

#plt.title(description)
plt.axes().set_aspect('equal', 'datalim')
#ax.set_xlabel("x")
#ax.set_ylabel("y")
plt.savefig("a.png")
#plt.show()

# -----------------------------------------------
# create animation
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
                                   interval=t_list[-1]*1e3/len(points_list), blit=True)
plt.show()
