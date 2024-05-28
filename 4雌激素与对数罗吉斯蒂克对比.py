import numpy as np
import pylab as mpl
import matplotlib.pyplot as plt
import 拟合曲线
from scipy import optimize as op
from scipy import stats

#调整spacing 0.250和0.250

#修改全局字体，显示中文
from matplotlib import rcParams
rcParams['font.family'] = ['cmex10']
plt.rcParams['figure.dpi'] = 300
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'

# def gaussian(x, b, c):
#     return (1/np.sqrt(2*np.pi)/c) * np.exp(-((x - b) ** 2) / (2 * c**2))
def gaussian(x, a, b):
    return np.exp((x-a)/b)/(b*(1+np.exp((x-a)/b))*(1+np.exp((x-a)/b)))

def funU(x,xi,h):
    return (x-xi)/h

def funKu(u):
    pi = 3.14159
    e = 2.71828
    sqrt2pi=((2*pi)**0.5)
    ecifang = -u*u/2
    return (1/sqrt2pi)*(e**ecifang)

def fguji(n,h,x,xi):
    pi = 3.14159
    e = 2.71828
    a = 1/n/h
    b = 0
    for i in xi:
        b+=funKu(funU(x, i, h))
    return a*b

x = [1000,1690,640,5120,5120,10000,20924,1690,1280,5120,5120,6250,200,100,11400,25,50,
    200,100,460,600,830,2283,5707,6076,1,1.75,2,5,11.4,250.59,490,500,600,640,640,1000,
    1000,1000,3120,7760,2600,2200,5000,1400,20,15]
xlog = []
for i in x:
    xlog.append(np.log(i))
#画密度图
YYuan,bins,patches = plt.hist(xlog, density=True, bins=15,edgecolor='black',linewidth=1,color="#FFFFFF",label='实际频率\nActual frequence')

XYuan = []
for i in range(0,len(list(bins))-1):
    XYuan.append((bins[i+1]+bins[i])/2)

x1 = np.linspace(0,15,1000)
n = len(xlog)
#for h in [0.234,0.199,0.405,0.450]:
h = 0.395
y = []
for i in x1:
    y.append(fguji(n,h,i,xlog))
ywucha = []
for i in XYuan:
    ywucha.append(fguji(n,h,i,xlog))
plt.plot(x1,y,c='black',label='非参数核密度估计\nNon-parametric\nkernel density estimation',linewidth=1,linestyle='-.')

popt, pcov = op.curve_fit(gaussian, XYuan, list(YYuan))

sigema=popt[1]
mu=popt[0]

gss = []
for i in x1:
    gss.append(gaussian(i,mu,sigema))

plt.plot(x1,gss,c='black',label='Log-logistic',linewidth=2,linestyle=':')


plt.xlabel("对数化后的雌激素效应数据/(μg·L$^{{\scr -1}}$)\nLogarithmic data of estrogen effects/(μg·L$^{{\scr -1}}$)",size = 15) 
plt.ylabel("概率密度\nProbability density",size = 15)

plt.legend(frameon=False,loc="upper left")
plt.tight_layout()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.25, hspace=0.25)

# 取消顶部和右侧边框
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.show()

