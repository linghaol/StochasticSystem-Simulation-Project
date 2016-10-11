# -*- coding: utf-8 -*-

### EE511 Project
import numpy as np
import matplotlib.pyplot as plot
import random
from scipy import stats
global fs;fs=11

### Project_1
class proj_1():
    def q1():
        ## 1.[Adding Coins]
        # 100 simulated fair Bernoulli trials and hitogram
        a1=[int(np.round(random.uniform(0,1))) for i in range(100)]
        plot.style.use('ggplot');
        plot.subplot(1,3,1)
        plot.hist(a1,np.arange(-1.5,3),color='blue');plot.grid(1);
        plot.title('Q1.1',fontsize=fs+2)
        plot.xticks(np.arange(2),('0','1'),fontsize=fs);plot.yticks(fontsize=fs)
        plot.xlabel('Value',fontsize=fs);plot.ylabel('Number',fontsize=fs)
        plot.show()
    
        #success in 5 fair Bernoulli trials, 100 samples
        a2=[sum([int(np.round(np.random.uniform(0,1))) for i in range(5)]) 
        for j in range(100)]
        plot.subplot(1,3,2);
        plot.hist(a2,bins=np.arange(min(a2)-0.5,max(a2)+1.5),color='red');
        plot.title('Q1.2',fontsize=fs+2)
        plot.xticks(fontsize=fs);plot.yticks(fontsize=fs)
        plot.xlabel('Value',fontsize=fs);plot.ylabel('Number',fontsize=fs)
        
        #number of trials before the first success, 100 samples
        trial_number=20
        a3=[[np.round(random.uniform(0,1)) for i in range(trial_number)] 
        for j in range(100)]
        for k in range(100):
            a3[k]=a3[k].index(1)
        plot.subplot(1,3,3);
        plot.hist(a3,bins=np.arange(min(a3)-0.5,max(a3)+1.5),color='black');
        plot.title('Q1.3',fontsize=fs+2)
        plot.xticks(fontsize=fs);plot.yticks(fontsize=fs)
        plot.xlabel('Value',fontsize=fs);plot.ylabel('Number',fontsize=fs)
        plot.show()
        return
    
    def q2():
        ## 2.[Coin Limits]
        bt=['k=2','k=5','k=10','k=30','k=50']
        b=bk=[2,5,10,30,50]
        for i in range(5):
            b[i]=[sum([int(np.round(random.uniform(0,1))) for j in range(bk[i])]) 
            for k in range(300)]
            plot.subplot(2,3,i+1)
            plot.hist(b[i],bins=np.arange(min(b[i])-0.5,max(b[i])+1.5),color='green')
            plot.xticks(fontsize=fs-2);plot.yticks(fontsize=fs-2)
            plot.xlabel('Value('+bt[i]+')',fontsize=fs-1);
            plot.ylabel('Number',fontsize=fs-1)
        plot.show()
        return
    
    def q3():
        ## 3.[Bootstraps]
        import pandas as pd
        n_x=6;k_t=1000
        c=pd.read_table('C:/Users/LLH/Desktop/EE 511/assignment/1/data/NJGAS.dat',
                        header=-1,index_col=0)
        boot=[np.mean([random.choice(c.index) for i in range(n_x)])
        for k in range(k_t)]
        plot.hist(boot[0:1000],bins=70,color='purple')
        plot.title('sampling of mean (1000 times)',fontsize=fs+2)
        plot.xticks(fontsize=fs);plot.yticks(fontsize=fs)
        plot.xlabel('Mean',fontsize=fs);plot.ylabel('Number',fontsize=fs)
        plot.show()
        m=np.mean(boot);s=np.std(boot)
        ci_min=m-1.96*s/np.sqrt(k_t);ci_max=m+1.96*s/np.sqrt(k_t)
        print('95% bootstrap confidence interval: (',ci_min,',',ci_max,')')
        return
       

### Project_2
class proj_2():
    def q1():
        ## 1.[Waiting]
        print('Please input the number of samples:')
        n=int(input())
        s_u=[random.uniform(0,1) for i in range(n)]
        s_e=[];chiv=[];sam_sum=[];sam_num=[]
        lamda=5
        for i in range(n):
            s_e.append(np.log(1-s_u[i])/(-lamda))  
        avg=np.average(s_e)
        var=np.var(s_e)
        print('average of sample set after inverse:',avg,
              '\nvariance of sample set after inverse:',var)
        print('please input the number of bins:')
        nbins=int(input())
        inter=list(np.arange(0,np.max(s_e)/nbins*(nbins+1),np.max(s_e)/nbins))
        inter.remove(0)
        ec=stats.expon.cdf(inter,scale=1/lamda)
        es=[n*ec[0]]
        for i in range(1,len(ec)):
            es.append(n*(ec[i]-ec[i-1]))
        for i in range(len(inter)):
            seq=[]
            for j in range(len(s_e)):
                if s_e[j]<inter[i]:
                    seq.append(j)
            sam_sum.append(len(seq))
        sam_num=[sam_sum[0]]
        for i in range(1,len(sam_sum)):
            sam_num.append(sam_sum[i]-sam_sum[i-1])
        for i in range(nbins):
            chiv.append((sam_num[i]-es[i])**2/es[i])
        print('Current degrees of freedom :',nbins-1,
        '\nChi-square test statistic :',sum(chiv))
        #tips: 0.05 critical value of doff=9 : 16.92
        
    def q2_s():
        ## 2.[Counting] (Slow Version)
        lamda=5;s_exp=[];s_cum=[];interv=[];interv_num=[];sumoft=0
        while 1:
            temp=np.log(1-random.uniform(0,1))/(-lamda)
            sumoft=sumoft+temp
            if sumoft<1000:
               s_exp.append(temp)
            else:
               break
        cumu=0
        for i in range(len(s_exp)):
            cumu=cumu+s_exp[i]
            s_cum.append(cumu)
        for i in range(1,1001):
            for j in range(len(s_cum)):
                if s_cum[j]>i:
                    interv.append(j)
                    break
        for i in range(1,len(interv)):
            interv_num.append(interv[i]-interv[i-1])
        plot.hist(interv_num)
        plot.show()        
            
    def q2():
        ## 2.[Counting] (High-speed Version) (no 'for' circulation)
        lamda=5;interv_num=[];sumoft=0;count=0;lmt=1
        while 1:
            temp=np.log(1-random.uniform(0,1))/(-lamda)
            if sumoft+temp<1000:
               sumoft=sumoft+temp
               while sumoft>lmt:
                   interv_num.append(count)
                   count=0
                   lmt=lmt+1
               count=count+1               
            else:
               break
        plot.hist(interv_num,
                  bins=np.arange(np.min(interv_num)-1.5,
                  np.max(interv_num)+2.5),)
        plot.show()  
        
    def q3():
        ## 3.[Networking]
        import networkx as nx
        G1=nx.Graph();G2=nx.Graph()
        edgs1=[];edgs2=[]
        for i in range(1,51):
            for j in range(i+1,51):
                p=random.uniform(0,1)
                q=random.uniform(0,1)
                if p<=0.02:
                    edgs1.append((i,j))
                if q<=0.09:
                    edgs2.append((i,j))
        G1.add_nodes_from(range(1,51))
        G2.add_nodes_from(range(1,51))
        G1.add_edges_from(edgs1)
        G2.add_edges_from(edgs2)
        deg1=[G1.degree(i) for i in range(1,51)]
        deg2=[G2.degree(j) for j in range(1,51)]
        plot.subplot(2,2,1)
        nx.draw(G1,node_size=300,with_labels=True,
                node_color='r')
        plot.subplot(2,2,2)
        nx.draw(G2,node_size=300,with_labels=True,
                node_color='g')
        plot.subplot(2,2,3)       
        plot.hist(deg1,
                  bins=np.arange(np.min(deg1)-1.5,
                                 np.max(deg1)+2.5),
                                 color='r');        
        plot.subplot(2,2,4)
        plot.hist(deg2,
                  bins=np.arange(np.min(deg2)-1.5,
                                 np.max(deg2)+2.5),
                                 color='g');         
    ##tips:
           # p*=log(n)/n threshod of how tight of loose of the graph
           # We need to discuss this parameter in the project Q3.4
           # 

    
    def q4():
        ## 4.[Network at larger scale]
        import networkx as nx    
        G3=nx.Graph()
        edgs3=[]
        for i in range(1,251):
            for j in range(i+1,251):
                r=random.uniform(0,1)
                if r<=0.08:
                    edgs3.append((i,j)) 
        G3.add_nodes_from(range(1,251))
        G3.add_edges_from(edgs3)
        deg3=[G3.degree(k) for k in range(1,251)]            
        plot.subplot(1,2,1)
        nx.draw(G3,node_size=300,node_color='b')
        plot.subplot(1,2,2)
        plot.hist(deg3,
                  bins=np.arange(np.min(deg3)-1.5,
                                 np.max(deg3)+2.5),
                                 color='b');                 

 
                    
            
        

        







     
   

