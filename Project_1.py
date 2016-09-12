### EE511 Project_1
import numpy as np
import matplotlib.pyplot as plot
import random
fs=11

## 1.[Adding Coins]
# 100 simulated fair Bernoulli trials and hitogram
a1=[int(np.round(np.random.uniform(0,1))) for i in range(100)]
plot.style.use('ggplot');
plot.subplot(1,3,1)
plot.hist(a1,np.arange(-1.5,3),color='blue');plot.grid(1);
plot.title('Q1.1',fontsize=fs+2)
plot.xticks(np.arange(2),('0','1'),fontsize=fs);plot.yticks(fontsize=fs)
plot.xlabel('Value',fontsize=fs);plot.ylabel('Number',fontsize=fs)

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
a3=[[np.round(np.random.uniform(0,1)) for i in range(trial_number)] 
for j in range(100)]
for k in range(100):
    a3[k]=a3[k].index(1)
plot.subplot(1,3,3);
plot.hist(a3,bins=np.arange(min(a3)-0.5,max(a3)+1.5),color='black');
plot.title('Q1.3',fontsize=fs+2)
plot.xticks(fontsize=fs);plot.yticks(fontsize=fs)
plot.xlabel('Value',fontsize=fs);plot.ylabel('Number',fontsize=fs)
plot.show()

## 2.[Coin Limits]
bt=['k=2','k=5','k=10','k=30','k=50']
b=bk=[2,5,10,30,50]
for i in range(5):
    b[i]=[sum([int(np.round(np.random.uniform(0,1))) for j in range(bk[i])]) 
    for k in range(300)]
    plot.subplot(2,3,i+1)
    plot.hist(b[i],bins=np.arange(min(b[i])-0.5,max(b[i])+1.5),color='green')
    plot.xticks(fontsize=fs-2);plot.yticks(fontsize=fs-2)
    plot.xlabel('Value('+bt[i]+')',fontsize=fs-1);
    plot.ylabel('Number',fontsize=fs-1)
plot.show()

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
