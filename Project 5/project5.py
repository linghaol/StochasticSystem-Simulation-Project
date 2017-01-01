# -*- coding: utf-8 -*-
"""
Created on Sun Jan  1 02:06:15 2017

@author: LLH
"""
import numpy as np
import matplotlib.pyplot as plot

### Project 5
class proj_5:
    ##[Monte Carlo]
    def q1_1():
        # K=50 runs of Pi-estimations with n=100 samples
        pi_est=[]
        for k in range(50):
            s=np.random.rand(100,2)
            count=0
            for i in s:
                if i[0]**2+i[1]**2<=1:
                    count+=1
            pi_est.append(count/100)
        plot.hist(pi_est)
        plot.xlabel("Pi")
    def q1_2():
        # Using different n, plot variance of estimations
        n=100
        v=[]
        while 1:
            pi_est2=[]       
            for k in range(50):
                s=np.random.rand(n,2)
                count=0
                for i in s:
                    if i[0]**2+i[1]**2<=1:
                        count+=1
                pi_est2.append(count/n)
            v.append(np.var(pi_est2))
            if n==1000:
                plot.plot(np.arange(100,1100,100),v)
                plot.xlabel("number of samples",fontsize=13)
                plot.ylabel("Variance",fontsize=13)
                break
            else:
                n+=100
    def q1_3():
        # introduction use
        from mpl_toolkits.mplot3d import Axes3D
        x=np.arange(0,1.005,0.005)
        y=np.arange(0,1.005,0.005)
        X,Y=np.meshgrid(x,y)
        G=np.abs(4*X-2)*np.abs(4*Y-2)
        n=100
        s=np.random.rand(n,3)
        s[:,2]=4*s[:,2]
        fig=plot.figure()
        ax=fig.add_subplot(111, projection='3d')
        ax.scatter3D(s[:,0],s[:,1],s[:,2],c='blue',s=50)
        ax.plot_surface(X,Y,G,color='yellow')
        # the 1st way
        inte=[];stder=[]
        n=100
        while 1:
            sample=[]
            for j in range(50):
                s=np.random.rand(n,3)
                s[:,2]=4*s[:,2]
                count=0
                for i in s:
                    if i[2]<=np.abs(4*i[0]-2)*np.abs(4*i[1]-2):
                        count+=1
                sample.append(count/n*4)
            inte.append(np.mean(sample))
            stder.append(np.std(sample))
            if n==10000:
                print('the 1st way:')
                print('integral=',inte)
                print('standard error=',stder,'\n')
                break
            else:
                n=n*10
        # the 2nd way
        n=100
        inte2=[];stder2=[]
        while 1:
            ss=np.random.rand(n,2)
            GG=np.abs(4*ss[:,0]-2)*np.abs(4*ss[:,1]-2)
            inte2.append(np.mean(GG))
            stder2.append(np.sqrt((np.sum((GG-np.mean(GG))**2)/(n-1)/n)))
            if n==10000:
                print('the 2nd way:')
                print('integral=',inte2)
                print('standard error=',stder2)
                break
            else:
                n=n*10
    ##[Variance Reduction Methods for Monte Carlo]
    def q2():
        #draw the functions
        from mpl_toolkits.mplot3d import Axes3D
        x=np.arange(0,1.005,0.005)
        y=np.arange(0,1.005,0.005)
        X,Y=np.meshgrid(x,y)
        G1=np.exp(5*np.abs(X-0.5)+5*np.abs(Y-0.5))       
        G3=np.abs(4*X-2)*np.abs(4*Y-2)
        fig=plot.figure()
        ax=fig.add_subplot(111, projection='3d')
        ax.plot_surface(X,Y,G1,color='yellow')
        fig=plot.figure()
        ax=fig.add_subplot(111, projection='3d')
        ax.plot_surface(X,Y,G3,color='green')
        xx=np.arange(-1,1.005,0.005)
        yy=np.arange(-1,1.005,0.005)
        XX,YY=np.meshgrid(xx,yy)
        G2=np.cos(np.pi+5*XX+5*YY) 
        fig=plot.figure()
        ax=fig.add_subplot(111, projection='3d')
        ax.plot_surface(XX,YY,G2,color='blue')
        #stratification: 4 subintervals
        #a,c
        def iv_a(n_acin,n_acout):
            sub_in=np.exp(5*np.abs(n_acin[:,0]-0.5)+5*np.abs(n_acin[:,1]-0.5))
            sub_out=np.exp(5*np.abs(n_acout[:,0]-0.5)+5*np.abs(n_acout[:,1]-0.5))
            itg_a=0.25*np.sum(sub_in)/np.shape(sub_in)[0]+0.75*np.sum(sub_out)/np.shape(sub_out)[0]
            var_a=1/16*np.var(sub_in)/(np.shape(sub_in)[0]-1)+\
            9/16*np.var(sub_out)/(np.shape(sub_out)[0]-1)
            return itg_a,var_a
        def iv_c(n_acin,n_acout):
            sub_in=np.abs(4*n_acin[:,0]-2)*np.abs(4*n_acin[:,1]-2)
            sub_out=np.abs(4*n_acout[:,0]-2)*np.abs(4*n_acout[:,1]-2)
            itg_c=0.25*np.sum(sub_in)/np.shape(sub_in)[0]+0.75*np.sum(sub_out)/np.shape(sub_out)[0]
            var_c=1/16*np.var(sub_in)/(np.shape(sub_in)[0]-1)+\
            9/16*np.var(sub_out)/(np.shape(sub_out)[0]-1)
            return itg_c,var_c
        sample_ac=np.random.rand(1000,2)
        suba=np.exp(5*np.abs(sample_ac[:,0]-0.5)+5*np.abs(sample_ac[:,1]-0.5))
        itg_sa=np.mean(suba)
        var_sa=np.var(suba)/(np.shape(suba)[0]-1)
        subc=np.abs(4*sample_ac[:,0]-2)*np.abs(4*sample_ac[:,1]-2)
        itg_sc=np.mean(subc)
        var_sc=np.var(subc)/(np.shape(subc)[0]-1)
        n_acin=[]
        n_acout=[]
        for i in sample_ac:
            if ((i[0]>0.25)&(i[0]<0.75))&((i[1]>0.25)&(i[1]<0.75)):
                n_acin.append(i)
            else:
                n_acout.append(i)
        n_acin=np.array(n_acin)
        n_acout=np.array(n_acout)
        ra=iv_a(n_acin,n_acout)
        rc=iv_c(n_acin,n_acout)
        print('result_a=',ra)
        print('result_c=',rc)
        #b
        def iv_b(n_bin,n_bout):
            sub_in=np.cos(np.pi+5*n_bin[:,0]+5*n_bin[:,1])
            sub_out=np.cos(np.pi+5*n_bout[:,0]+5*n_bout[:,1])
            itg_c=1*np.sum(sub_in)/np.shape(sub_in)[0]+3*np.sum(sub_out)/np.shape(sub_out)[0]
            var_c=1/16*np.var(sub_in)/(np.shape(sub_in)[0]-1)+\
            9/16*np.var(sub_out)/(np.shape(sub_out)[0]-1)
            return itg_c,var_c           
        sample_b=2*np.random.rand(1000,2)-1
        subb=np.cos(np.pi+5*sample_b[:,0]+5*sample_b[:,1])
        itg_sb=np.mean(subb)*4
        var_sb=np.var(subb)/(np.shape(subb)[0]-1)
        n_bin=[]
        n_bout=[]
        for i in sample_b:
            if ((i[0]>-0.5)&(i[0]<0.5))&((i[1]>-0.5)&(i[1]<0.5)):
                n_bin.append(i)
            else:
                n_bout.append(i)
        n_bin=np.array(n_bin)
        n_bout=np.array(n_bout)       
        rb=iv_b(n_bin,n_bout)
        print('result_b=',rb)

        #importance sampling
               
    ##[MCMC for Optimization]
    def q3_1():
        #2 dimension case
        x=np.arange(-500,503,2)
        y=np.arange(-500,503,2)
        X,Y=np.meshgrid(x,y)
        f=418.9829*2-(X*np.sin(np.sqrt(np.abs(X)))+
        Y*np.sin(np.sqrt(np.abs(Y))))
        plot.contour(X,Y,f)
        plot.colorbar()
        plot.xticks(np.arange(-500,600,100))
        plot.yticks(np.arange(-500,600,100))
    def q3_2(s,tstep):
        # simulated annealing
        from random import random as rd
        from random import choice
        from scipy.constants import k
        def cost(s):
            c=418.9829*2-(s[0]*np.sin(np.sqrt(np.abs(s[0])))+
            s[1]*np.sin(np.sqrt(np.abs(s[1]))))
            return c 
        def neighbor(s,c):
            ax=choice([-5,0,5])
            ay=choice([-5,0,5])
            nx=s[0]+ax
            ny=s[1]+ay
            if np.abs(nx)>500:
                nx=nx/np.abs(nx)*500
            if np.abs(ny)>500:
                ny=ny/np.abs(ny)*500
            if cost([nx,ny])-c<5:
                return np.round(1000*np.random.rand(2)-500)
            else:
                return np.array([nx,ny])
        def accept_prob(c,n_c,T,k):
            a=min([1,np.exp((n_c-c)/(k*T))])
            return a
        def cooling(T,t):
            return T*0.9*np.exp(-t)  #exponential
#            return 1/2*t**(-2)+1/2*t**(-1)  #polynomial  
#            return T/(1+5*np.log(1+t))  #logarithmic
        T=1.0
        c=cost(s)
        past=np.array([s])
        T=cooling(T,tstep)
        i=1
        while i<=100:
            n_s=neighbor(s,c)
            n_c=cost(n_s)
            ap=accept_prob(c,n_c,T,k)
            if ap>rd():
                past=np.append(past,[s],axis=0)
                s=n_s
                c=n_c
            i+=1
#        return c # used for q3_3
        return past # used for q3_4
    def q3_3():
        r20=[]
        for i in range(20):
            r20.append(proj_5.q3_2([0,0],i+1))
        plot.subplot(2,2,1)
        plot.hist(r20)
        plot.title('t=20')
        r50=[]
        for i in range(50):
            r50.append(proj_5.q3_2([0,0],i+1))
        plot.subplot(2,2,2)
        plot.hist(r50)
        plot.title('t=50')
        r100=[]
        for i in range(100):
            r100.append(proj_5.q3_2([0,0],i+1))
        plot.subplot(2,2,3)
        plot.hist(r100)
        plot.title('t=100')            
        r1000=[]
        for i in range(1000):
            r1000.append(proj_5.q3_2([0,0],i+1))
        plot.subplot(2,2,4)
        plot.hist(r1000)
        plot.title('t=1000')
    def q3_4():
        fp=proj_5.q3_2([0,0],1)
        proj_5.q3_1()
        plot.scatter(fp[:,0],fp[:,1])