# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 00:07:05 2019

@author: Karnika
"""

from gurobipy import*
import os
import xlrd
from scipy import spatial
from sklearn.metrics.pairwise import euclidean_distances
import math


book = xlrd.open_workbook(os.path.join("data (2).xlsx"))

                ##################  PARAMETERS  ##################
                                
N=[]                                                                           ## Set of all Nodes
cij={}                                                                         ## Travelling Distance
tij={}                                                                         ## Travel time
K=[]                                                                           ## Set of Vehicles
Vehicle_Renting_Cost = {'V1':5000,'V2':5000,'V3':15000,'V4':20000}
T=[1,2,3]                                                                      ## set of time period
Available_time_per_period={1:15000,2:15000,3:15000,4:28800,5:0}
Q={'V1':3000000,'V2':3000000,'V3':15000000,'V4':20000000}                                                  ## Fleet Size
P=[1,2,3]                                                                      ## Set of Different Product
h={1:0.00017,2:0.00017,3:0.00017}                                              ## Holding Cost of products
U={1:300,2:300,3:400}                                                          ## Holding Capacity of Product
D={}                                                                           ## Demand of RATM
Iit_1={}                                                                       ## initial inventory level
co={}                                                                          ## Co-ordinates of DEPOT and RATM
sk=5                                                                           ## Average speed of Vehicle
Notes={1:2000,2:500,3:100}                                                     ## Various Currency notes

              ##################  DATA IMPORT  ##################
sh = book.sheet_by_name("Sheet1")
i = 1
while True:
    try:
        sp = sh.cell_value(i,0)
        N.append(sp)
        sp1=sh.cell_value(i,1)
        sp2=sh.cell_value(i,2)
        sp3=(sp1,sp2)
        co[sp]=sp3
        D[sp]=sh.cell_value(i,7)
        sp4=sh.cell_value(i,8)
        K.append(sp4)
        i = i + 1
    except IndexError:
        break

for i in range(len(K)):
    if K[i]=='':
        K=K[:i]
        break
        
i = 1
for k in N:
    j = 3
    for l in P:
        Iit_1[k,l] = sh.cell_value(i,j)
        j += 1
    i += 1
N_Dash=N[1:]
def calculate_dist(x1, x2):
    eudistance = spatial.distance.euclidean(x1, x2)    
    return(eudistance)    
for i in N:
    for j in N:
        cij[i,j]=int(round(calculate_dist(co[i],co[j])*10))
        tij[i,j]=cij[i,j]//sk
        
             ##################     MODEL   ##################       
             
             
m = Model("Basic model Of Inventory Routing Problem")

m.modelSense=GRB.MINIMIZE

            ##################   VARIABLES   ##################
        
yijkt  = m.addVars(N,N,K,T,vtype=GRB.BINARY ,name='Y_ijkt')                    ## Route variable Binary
Iitp   = m.addVars(N,T,P,  vtype=GRB.INTEGER,name='I_itp' )                    ## Inventory level of each node at the end of period
qiktp  = m.addVars(N,K,T,P,vtype=GRB.INTEGER,name='q_iktp')                    ## amount of quantity delivered or pickup in period t by vehicle k
Uikt   = m.addVars(N,K,T,  vtype=GRB.INTEGER,name='U_ikt' )                    ## Subtour elimination variable
Vkt    = m.addVars(K,T,    vtype=GRB.BINARY ,name='V_kt'  )                    ## Vehicle k used in period t
xikt   = m.addVars(N,K,T,  vtype=GRB.BINARY ,name='X_ikt' )                    ## binary variable visiting variable to node
zikt   = m.addVars(N,K,T,  vtype=GRB.INTEGER,name='Z_ikt' )                    ## Integer variable visiting Variable to node
f      = m.addVars(N,T,P,  vtype=GRB.INTEGER,name='f_itp' )                    ## No of Various Currency meeting demand of node i in period t

            #############   OBJECTIVE FUNCTION   #############
            
            
m.setObjective(
#               sum(h[p]*Iit_1[i,p] for i in N for p in P)+
               sum(h[p]*Iitp[i,t,p]*Notes[p] for i in N for t in T for p in P ) +                      
               sum(Vehicle_Renting_Cost[k]*Vkt[k,t] for k in K for t in T )
               )                                                               ## inimize total inventory cost and hiring cost of Vehicle

            ################## CONSTRAINTS ##################
###Constraint 1:
for i in N_Dash:
    for t in T:
        for p in P:
            if t == 1:
                m.addConstr(Iit_1[i,p]  + sum(qiktp[i,k,t,p]  for k in K )- f[i,t,p] == Iitp[i,t,p]  )
            else:
                m.addConstr(Iitp[i,t-1,p] +sum(qiktp[i,k,t,p]  for k in K) - f[i,t,p] == Iitp[i,t,p] )
###Constraint 2:            
for i in N_Dash:
    for t in T:            
        m.addConstr(sum(f[i,t,p]*Notes[p] for p in P)==D[i])                                                                  ## Inventory balance constraint at each node at end of period t
####Constraint 3:-                
for i in N_Dash:
    for t in T:
        for p in P:
            m.addConstr(Iitp[i,t,p] <= U[p])
                                                                               ## capacity constraint 
####Constraint 4:-
for i in N_Dash:
    for k in K:
        for t in T:
            for p in P:
                m.addConstr(qiktp[i,k,t,p]<=U[p]*zikt[i,k,t])                  ## delivered or pickup amount less than capacity RATM
                if t == 1:
                    m.addConstr(qiktp[i,k,t,p] + Iit_1[i,p] <= U[p])
                else:
                    m.addConstr(qiktp[i,k,t,p] + Iitp[i,t-1,p] <= U[p])


####Constraint 5:-                    
for k in K:
     for t in T:
         m.addConstr(sum(qiktp[i,k,t,p]*Notes[p] for i in N_Dash for p in P) <= Q[k]*zikt[1,k,t])
                                                                               ## Sum of all pickup and deliverd items is less than vehicle in any route
##Constraint 6:-
for i in N_Dash:
    for t in T:
        m.addConstr(sum(zikt[i,k,t] for k in K)<=1)
                                                                               ## Node must be visited by only one vehicle 
#####Constraint 7:-                
for i in N:
    for k in K:
        for t in T:
            m.addConstr(sum(yijkt[i,j,k,t] for j in N  if j!=i )==zikt[i,k,t]) 
            m.addConstr(sum(yijkt[j,i,k,t] for j in N  if i!=j )==zikt[i,k,t])                                  
                                                                               ## Route Constraint
####Constraint 8:-
                                                                               
            #############  Subtour Elimination Constraint #############
            
for i in N_Dash:
    for t in T:
        for k in K:
            for j in N:
                if i!=j:
                    m.addConstr((Uikt[i,k,t]-Uikt[j,k,t] + Q[k]*yijkt[i,j,k,t])<= Q[k] - sum(qiktp[j,k,t,p]*Notes[p] for p in P))  
                                                                               ## Subtour elimination constraint
###Constraint 9:-
                    
for i in N_Dash:
    for t in T:
        for k in K:
            m.addConstr(Uikt[i,k,t]<=Q[k]) and m.addConstr(Uikt[i,k,t]>=sum(qiktp[i,k,t,p]*Notes[p] for p in P))
                                                                               ## Subtour elimination constraint
### ARK Constraints:-
###Constraint 10:-
for i in N:
    for j in N:
        for k in K:
            for t in T:
                if i!=j:
                    m.addConstr(yijkt[i,j,k,t] <= xikt[i,k,t])
                                                                               ## If node is visited by vehicle k then xikt =1 binary
###Constraint 11:-
for i in N_Dash:
    for k in K:
        for t in T:
            m.addConstr(xikt[i,k,t]<=Vkt[k,t])
                                                                               ## Vehicle k used in period t
###Constraint 12:-
for k in K:
    for t in T:
        m.addConstr(sum(yijkt[i,j,k,t]*tij[i,j] for i in N for j in N if i!=j)+1800*sum(zikt[i,k,t] for i in N)<=Available_time_per_period[t])
                                                                               ## Time limiting constraint on vehicle                                                                                                                            
            
time=m.runtime
m.write('MTZ.lp')
m.setParam('TimeLimit',3600)
m.optimize()

                #############   RESULTS     #############
                
for v in m.getVars():
    if v.x > 0.01:
        print(v.varName, v.x)
print('Objective:',round(m.objVal,2))
print('time',m.runtime)

                ############# PROGRAME END  #############
