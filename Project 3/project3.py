# -*- coding: utf-8 -*-
"""
Created on Sun Jan  1 01:58:12 2017

@author: LLH
"""
import re
import numpy as np
import matplotlib.pyplot as plot
import random
from scipy.stats import bernoulli as ber

### Prject_3
class proj_3:
    ##[Testing Faith, with K-means Cluster](refined version with standardization)
    def q1():
        f=open('C:/Users/LLH/Desktop/EE 511/assignment/3/old-faithful.txt')
        predata=f.readlines()
        f.close
        dura=[];time=[]
        for i in predata:
            k=re.split(r'[\s]*',i)
            try:
                int(k[0])
            except ValueError:
                continue
            else:
                dura.append(float(k[1]))
                time.append(float(k[2]))
        ## K-means cluster
        # initial centers
        c=input('Please initial centers:format:x1 y1 x2 y2=')
        c=np.array([float(i) for i in c.split(' ')]);
        c1=c[0:2];c2=c[2:4]
        cr=float(input('Please input convergence criterion='))
        #standardization
        md=np.mean(dura);sd=np.std(dura)
        mt=np.mean(time);st=np.std(time)
        dura=(dura-md)/sd
        time=(time-mt)/st
        c1=(c1-np.array([md,mt]))/np.array([sd,st])
        c2=(c2-np.array([md,mt]))/np.array([sd,st])
        k=0
        # calculate distances,classify points and find new centers
        while 1:
            k+=1
            sx1=0;sy1=0;sx2=0;sy2=0;n1=0;n2=0
            c1_d=[];c1_t=[];c2_d=[];c2_t=[]
            for i in range(len(dura)):
                if ((dura[i]-c1[0])**2+(time[i]-c1[1])**2)\
                <((dura[i]-c2[0])**2+(time[i]-c2[1])**2):
                    sx1+=dura[i];sy1+=time[i];n1+=1
                    c1_d.append(dura[i]);
                    c1_t.append(time[i]);
                else:
                    sx2+=dura[i];sy2+=time[i];n2+=1
                    c2_d.append(dura[i]);
                    c2_t.append(time[i]);
            nc1_d=sx1/n1;nc1_t=sy1/n1
            nc2_d=sx2/n2;nc2_t=sy2/n2
            dist1=(nc1_d-c1[0])**2+(nc1_t-c1[1])**2
            dist2=(nc2_d-c2[0])**2+(nc2_t-c2[1])**2
            if dist1+dist2<=cr:
                break
            elif k>=1000:
                print('Warning: calculating times over 1k,\
                please change convergence criterion.')
                break
            else:
                c1=[nc1_d,nc1_t];c2=[nc2_d,nc2_t]
        #drawing scatter       
        c1_d=np.array(c1_d)*sd+md;c2_d=np.array(c2_d)*sd+md
        c1_t=np.array(c1_t)*st+mt;c2_t=np.array(c2_t)*st+mt
        c1=np.array(c1)*np.array([sd,st])+np.array([md,mt])
        c2=np.array(c2)*np.array([sd,st])+np.array([md,mt])
        nc1_d=nc1_d*sd+md;nc2_d=nc2_d*sd+md
        nc1_t=nc1_t*st+mt;nc2_t=nc2_t*st+mt
        fz=15
        plot.scatter(c1_d,c1_t,s=50,c='red',linewidth=0)
        plot.xticks(np.arange(1,5.5,0.5))
        plot.yticks(np.arange(40,100,10))
        plot.scatter(c2_d,c2_t,s=50,c='blue',linewidth=0)
        plot.xticks(np.arange(1,5.5,0.5))
        plot.yticks(np.arange(40,100,10))
        plot.scatter([nc1_d,nc2_d],[nc1_t,nc2_t],s=150,c='black')
        plot.xticks(np.arange(1,5.5,0.5))
        plot.yticks(np.arange(40,100,10))
        plot.xlabel('duration',fontsize=fz)
        plot.ylabel('waiting time',fontsize=fz)
        plot.grid()
        plot.show()       
        print('previous centers=',c1,c2)
        print('final centers=',[nc1_d,nc1_t],[nc2_d,nc2_t])
        print('calculating times=',k)       

    ##[Generating Mixed Samples]
    def q2():
        k=ber.rvs(0.4,size=1000);dataset=[]
        for i in k:
            if i==1:
                dataset.append(random.gauss(-1,1))
            else:
                dataset.append(random.gauss(1,1))
        c=input('Please initial center(1D-2 centers)\
        :format:x1 x2=')
        c=c.split(' ')
        c1=float(c[0]);c2=float(c[1])
        cr=float(input('Please input convergence criterion='))
        k=0;data1=[];data2=[]
        while 1:
            k+=1
            sx1=0;sx2=0;n1=0;n2=0
            for i in dataset:
                if (i-c1)**2<(i-c2)**2:
                    sx1+=i;n1+=1
                    data1.append(i)
                else:
                    sx2+=i;n2+=1
                    data2.append(i)
            nc1=sx1/n1;nc2=sx2/n2
            dist=(nc1-c1)**2+(nc2-c2)**2
            if dist<cr:
                break
            elif k>=1000:
                print('Warning: calculating times over 1k,\
                please change convergence criterion.')
                break
            else:
                c1=nc1;c2=nc2
        #drawing histogram    
        fz=15
        d_g=np.arange(-5,5,0.01)
        g1=(1/(2*np.pi))*np.exp(-(d_g-(-1))**2/(2*1))
        g2=(1/(2*np.pi))*np.exp(-(d_g-1)**2/(2*1))
        plot.plot(d_g,0.4*g1,c='r')
        plot.plot(d_g,0.6*g2,c='g')
        plot.plot(d_g,0.4*g1+0.6*g2,c='black')
        plot.figure()
        plot.hist(data1,bins=np.arange(-5,5,0.1),color='red')
        plot.hist(data2,bins=np.arange(-5,5,0.1),color='blue')
        plot.xlabel('Value',fontsize=fz)
        plot.ylabel('Frequency',fontsize=fz)
        plot.show()
        print('previous centers=',c1,c2)
        print('final centers=',nc1,nc2)
        print('calculating times=',k)