import numpy as np
import numpy.linalg as la
import csv
with open('example.dat', 'r') as f:
    data = csv.reader(f, delimiter = ' ')
    A = np.array(list(data)).astype('float')
print('Original Data:')
print()
print(A)
A = A.T

def PCA_algo(A):
    X = A
    Y = 1/np.sqrt(X.shape[1] - 1) * X.T
    U, sigma, Vt = la.svd(Y, full_matrices=False)
    new_X = Vt.dot(X)
    print('Data in new PCA coordinates:')
    print()
    print(new_X.T.round(2))
    print('===========================')
    eigenv = sigma**2
    totalvar=eigenv.sum()
    for i in range(len(sigma)):
        print('Variance captured in PC'+str(i+1)+': ', str(np.round((eigenv[i]/totalvar*100),2))+'%',', Accumulative: ',str(np.round((np.sum(eigenv[0:i+1])/totalvar*100),2))+'%')
    print('===========================')
    print('Standard deviation of each princinal component:')
    print('[PC1, PC2, PC3, PC4]')
    print(new_X.std(axis=1).round(2))
    
    plt.figure()
    plot(range(len(eigenv)), eigenv, 'bo',linestyle='dashed')
    arr = eigenv[::-1]
    plt.xlabel('Principle components')
    plt.ylabel('variance (blue) + cumulative var (green)')
    plot(range(len(eigenv)), [np.sum(eigenv[0:i+1]) for i in range(len(eigenv))], 'go',linestyle='dashed')
	
print('================================')
print()
print('PCA applied to Covariance Matrix')
print()
print('================================')
A = (A - A.mean(axis=1)[:,np.newaxis]).round(2)
PCA_algo(A)


print('================================')
print()
print('PCA applied to Correlation Matrix')
print()
print('================================')
A = ((A - A.mean(axis=1)[:,np.newaxis])/A.std(axis=1)[:,np.newaxis]).round(2)
PCA_algo(A)


import numpy as np
import numpy.linalg as la
import csv
with open('RCsGoff.csv','r') as f:
    data = csv.reader(f,delimiter=',')
    Ab=np.array(list(data))
    
print(Ab[0,:])
print()

print(Ab[1,:])
print()

print(Ab.shape)
A = Ab[1:,1:].astype('float')

print(A.shape)



# The same function used in previous exercise
def PCA_algo(A):
    X = A
    Y = 1/np.sqrt(X.shape[1] - 1) * X.T
    U, sigma, Vt = la.svd(Y, full_matrices=False)
    
    eigenv = sigma**2    
    
    new_X = Vt.dot(X)
    return new_X, eigenv
	
	
	
# zero mean data
A0 = A - A.mean(axis=1)[:,np.newaxis]

pca, eigenv = PCA_algo(A0)

from matplotlib.pyplot import plot
from matplotlib import pyplot as plt
from matplotlib.ticker import EngFormatter
%matplotlib inline

print(pca.shape)
xlabel='PC1: '+str(np.round(eigenv[0]/eigenv.sum()*100,10))+' % variance'
ylabel='PC2: '+str(np.round(eigenv[1]/eigenv.sum()*100,10))+' % variance'
pyplot.xlabel(xlabel)
pyplot.ylabel(ylabel)

plot(pca[0,:], pca[1,:],'+')


sample = Ab[0,1:] # stands for day0_rep1,...,day18_rep3

# output matrix, colums = [sample, PC1, ..., PC20, Variance saved]
output = np.c_[sample, pca.T, eigenv[:,np.newaxis]/eigenv.sum()*100]

# writing output matrix to csv file
with open('RcsGoff_output.csv','w',newline='') as f:
    csvwriter = csv.writer(f, delimiter=',')
    csvwriter.writerows(output)
print('Done!')

print(output.shape)


plt.figure()
plot(range(len(eigenv)), eigenv, 'bo',linestyle='dashed')
arr = eigenv[::-1]
plt.xlabel('Principle components')
plt.ylabel('variance (blue) + cumulative var (green)')
plot(range(len(eigenv)), [np.sum(eigenv[0:i+1]) for i in range(len(eigenv))], 'go',linestyle='dashed')


