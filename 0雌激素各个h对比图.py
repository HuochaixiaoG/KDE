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
#数据添加到这里
x = []

xlog = []
for i in x:
    xlog.append(np.log(i))
#画密度图
YYuan,bins,patches = plt.hist(xlog, density=True, bins=15,edgecolor='black',linewidth=1,color="#FFFFFF",label='  实际频率\nActual frequence')

XYuan = []
for i in range(0,len(list(bins))-1):
    XYuan.append((bins[i+1]+bins[i])/2)

x1 = np.linspace(0,15,1000)
n = len(xlog)
#for h in [0.234,0.199,0.405,0.450]:
h = [0.255,0.355,0.395,0.455]
lw = [1,1,2,1]
ls = ['-.',':','-','--']
for j in range(0,4):
    y = []
    for i in x1:
        y.append(fguji(n,h[j],i,xlog))
    ywucha = []
    for i in XYuan:
        ywucha.append(fguji(n,h[j],i,xlog))
    # print(h)
    # 拟合曲线.wuchajisuan(YYuan, ywucha)
    plt.plot(x1,y,c='black',label='   h='+str(h[j]),linewidth=lw[j],linestyle=ls[j])

plt.xlabel("对数化后的雌激素效应数据/(μg·L$^{{\scr -1}}$)\nLogarithmic data of estrogen effects/(μg·L$^{{\scr -1}}$)",size = 15) 
plt.ylabel("概率密度\nProbability density",size = 15)

plt.legend(frameon=False,loc="upper left")
plt.tight_layout()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.25, hspace=0.25)

# 取消顶部和右侧边框
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)



plt.show()

