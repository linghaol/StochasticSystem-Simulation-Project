# -*- coding: utf-8 -*-
"""
Created on Sun Jan  1 02:06:13 2017

@author: LLH
"""
import numpy as np
import re
import matplotlib.pyplot as plot

### Project_4
class proj_4:
    ##EM-estimation algorithm           
    def em(sph,switch):
    #EM algorithm
        from numpy import pi
        import numpy.linalg as lin
        para=input("initial 2-D parameters(format:alpha1 alpha2\
         mu1x1 mu1x2 mu2x1 mu2x2)\n=")
        para=para.split(" ")
        a1=float(para[0]);a2=float(para[1])
        mu1=np.array([float(para[2]),float(para[3])])      
        mu2=np.array([float(para[4]),float(para[5])])
        sigmax1=input("initial dist1 sigma(format:sigma11 sigma12\
         sigma21 sigma22)\n=")
        sigmax1=sigmax1.split(" ")
        sigma1=np.array([[float(sigmax1[0]),float(sigmax1[1])],
                           [float(sigmax1[2]),float(sigmax1[3])]])
        sigmax2=input("initial dist2 sigma(format:sigma11 sigma12\
         sigma21 sigma22)\n=")
        sigmax2=sigmax2.split(" ")
        sigma2=np.array([[float(sigmax2[0]),float(sigmax2[1])],
                           [float(sigmax2[2]),float(sigmax2[3])]])
        k=1;
        ori=sph
        while 1:
            if switch==1:
                #set noise standard deviation
                stdN=0.6/k**2
                noise=np.array(stdN*np.random.randn(sph.shape[0],2))               
                sph=noise+sph                
            fyH1=[];fyH2=[];fyH=[];sig1i=[];sig2i=[]
            for i in range(len(sph)):
                fyH1.append(1/(np.sqrt((2*pi)**2*lin.det(sigma1)))*np.exp(
                -0.5*np.dot(np.dot((sph[i,:]-mu1),lin.inv(sigma1)),(sph[i,:]-mu1).T)))
                fyH2.append(1/(np.sqrt((2*pi)**2*lin.det(sigma2)))*np.exp(
                -0.5*np.dot(np.dot((sph[i,:]-mu2),lin.inv(sigma2)),(sph[i,:]-mu2).T)))                
            fyH1=np.array(fyH1);fyH2=np.array(fyH2)
            fyH=a1*fyH1+a2*fyH2
            pz1yh=a1*fyH1/fyH;pz1yh=pz1yh.reshape(len(sph),1)
            pz2yh=a2*fyH2/fyH;pz2yh=pz2yh.reshape(len(sph),1)
            a1_n=sum(pz1yh)/len(sph)
            a2_n=sum(pz2yh)/len(sph)
            mu1_n=sum(pz1yh*sph)/sum(pz1yh)
            mu2_n=sum(pz2yh*sph)/sum(pz2yh)
            for i in range(len(sph)):
                sig1i.append(pz1yh[i]*(sph[i,:]-mu1).reshape(2,1)*(sph[i,:]-mu1))
                sig2i.append(pz2yh[i]*(sph[i,:]-mu2).reshape(2,1)*(sph[i,:]-mu2))
            sigma1_n=sum(sig1i)/sum(pz1yh)
            sigma2_n=sum(sig2i)/sum(pz2yh)
            k+=1
            jud1=(a1_n-a1<0.0001)&(sum(mu1_n-mu1)<0.0001/2)&(sum(sum(sigma1_n-sigma1))<0.001/4)
            jud2=(a2_n-a2<0.0001)&(sum(mu2_n-mu2)<0.0001/2)&(sum(sum(sigma2_n-sigma2))<0.001/4)
            if jud1&jud2:
                result=[a1_n,a2_n,mu1_n,mu2_n,sigma1_n,sigma2_n]
                print("converged parameters:(alpha1,alph2,mu1,mu2,sigma1,sigma2)\n="
                ,result)
                print("calculation times=",k)
                break
            else:
                a1=a1_n;a2=a2_n
                mu1=mu1_n;mu2=mu2_n
                sigma1=sigma1_n;sigma2=sigma2_n
        #contour and scatter    
        x=np.arange(np.floor(np.min(ori[:,0]))-5,np.floor(np.max(ori[:,0]))+5,1)
        y=np.arange(np.floor(np.min(ori[:,1]))-5,np.floor(np.max(ori[:,1]))+5,1)
        X,Y=np.meshgrid(x,y)
        Z=[]
        for i in range(len(x)):
            for j in range(len(y)):
                xy=np.array([x[i],y[j]])
                Z.append((a1_n*1/(np.sqrt((2*pi)**2*lin.det(sigma1_n)))*np.exp(
                -0.5*np.dot(np.dot((xy-mu1_n),lin.inv(sigma1_n)),(xy-mu1_n).T)))+\
                (a2_n*1/(np.sqrt((2*pi)**2*lin.det(sigma2_n)))*np.exp(
                -0.5*np.dot(np.dot((xy-mu2_n),lin.inv(sigma2_n)),(xy-mu2_n).T))))                 
        Z=np.array(Z).reshape(len(x),len(y)).T
        plot.contour(X, Y, Z)
        plot.scatter(ori[:,0],ori[:,1],c="red",linewidth=0)
        plot.show()
        #EM clustering
        set1=[];set2=[]
        for l in ori:
            p1=(a1_n*1/(np.sqrt((2*pi)**2*lin.det(sigma1_n)))*np.exp(
                -0.5*np.dot(np.dot((l-mu1_n),lin.inv(sigma1_n)),(l-mu1_n).T)))
            p2=(a2_n*1/(np.sqrt((2*pi)**2*lin.det(sigma2_n)))*np.exp(
                -0.5*np.dot(np.dot((l-mu2_n),lin.inv(sigma2_n)),(l-mu2_n).T)))
            if p1>=p2:
                set1.append(l)
            else:
                set2.append(l)
        set1=np.array(set1);set2=np.array(set2);
        plot.figure()
        plot.scatter(set1[:,0],set1[:,1],s=50,c="blue",linewidth=0)
        plot.scatter(set2[:,0],set2[:,1],s=50,c="red",linewidth=0)
        line=[];lx=np.arange(1,4.5,0.5)
        for n in lx:
            line.append(-60/3.5*n+(60/3.5+100))
        plot.plot(lx,line)   
    ##[EM]
    def q1(switch):
        #2D-Random Number Generator for GMM
        n_sam=150 #number of samples
        #spherical gaussian, mean[m1,m2],standard deviation[s1,s2]
        m1=[-2,-2];m2=[2,2]
        s1=np.array([1,1]);s2=np.array([1,1])
        gau1=np.sqrt(s1)*np.random.randn(n_sam,2)+np.array([m1[0],m1[1]])
        gau2=np.sqrt(s2)*np.random.randn(n_sam,2)+np.array([m2[0],m2[1]])
        sph=np.vstack([gau1,gau2])
        #run EM
        proj_4.em(sph,switch)         
    ##[Testing Faith Again]
    def q2(switch):
        #import dataset "old-faithful"
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
        dataset=np.array([dura,time]).T
        #run EM
        proj_4.em(dataset,switch)
    ##[Noise in GMM-EM]
    def q3(switch):
        #import dataset "old-faithful"
        proj_4.q2(switch)