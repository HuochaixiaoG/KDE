import numpy as np
import pylab as mpl
import matplotlib.pyplot as plt
from scipy.integrate import quad

#修改全局字体，显示中文
from matplotlib import rcParams
rcParams['font.family'] = ['cmex10']
plt.rcParams['figure.dpi'] = 300
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'

# for j in range(10,21,2):
#     n = j
#     p = []
#     print(n)
#     for i in range(1,n+1):
#         p.append((i-0.5)/n)
#     print(p)
# n = 106 
# p = []
# for i in range(1,n+1):
#     p.append((i-0.5)/n)
# print(p)
def leijigailv(n):
    p = []
    for i in range(1,n+1):
        p.append((i-0.5)/n)
    return p 
# 10
# [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95]
# 12
# [0.041666666666666664, 0.125, 0.20833333333333334, 0.2916666666666667, 0.375, 0.4583333333333333, 0.5416666666666666, 0.625, 0.7083333333333334, 0.7916666666666666, 0.875, 0.9583333333333334]
# 14
# [0.03571428571428571, 0.10714285714285714, 0.17857142857142858, 0.25, 0.32142857142857145, 0.39285714285714285, 0.4642857142857143, 0.5357142857142857, 0.6071428571428571, 0.6785714285714286, 0.75, 0.8214285714285714, 0.8928571428571429, 0.9642857142857143]
# 16
# [0.03125, 0.09375, 0.15625, 0.21875, 0.28125, 0.34375, 0.40625, 0.46875, 0.53125, 0.59375, 0.65625, 0.71875, 0.78125, 0.84375, 0.90625, 0.96875]
# 18
# [0.027777777777777776, 0.08333333333333333, 0.1388888888888889, 0.19444444444444445, 0.25, 0.3055555555555556, 0.3611111111111111, 0.4166666666666667, 0.4722222222222222, 0.5277777777777778, 0.5833333333333334, 0.6388888888888888, 0.6944444444444444, 0.75, 0.8055555555555556, 0.8611111111111112, 0.9166666666666666, 0.9722222222222222]
# 20
# [0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975]

def gaussian(x):
    return (1/np.sqrt(2*np.pi)/1.943) * np.exp(-((x - 7.371) ** 2) / (2 * 1.943**2))

def logistic(x):
    #return (1.0 / (1.0 + np.exp(-a * x))) - b
    try:
        return 0.139 / (1.0 + x**(-0.872)) 
    except:
        return 0.0
def loggaussian(x):
    return (1/(x*np.sqrt(2*np.pi)*0.249)) * np.exp(-((np.log(x) - 2.038) ** 2) / (2 * 0.249**2))

def logLogistic(x):
    return np.exp((x-7.346)/1.192)/(1.192*(1+np.exp((x-7.346)/1.192))*(1+np.exp((x-7.346)/1.192)))

def funU(x,xi,h):
    return (x-xi)/h

def funKu(u):
    pi = 3.14159
    e = 2.71828
    sqrt2pi=((2*pi)**0.5)
    ecifang = -u*u/2
    return (1/sqrt2pi)*(e**ecifang)

def fguji(x):
    pi = 3.14159
    e = 2.71828
    n = 47
    h = 0.395
    a = 1/n/h
    x0 = [1000,1690,640,5120,5120,10000,20924,1690,1280,5120,5120,6250,200,100,11400,25,50,
        200,100,460,600,830,2283,5707,6076,1,1.75,2,5,11.4,250.59,490,500,600,640,640,1000,
        1000,1000,3120,7760,2600,2200,5000,1400,20,15]
    h = 0.395
    xi = []
    for i in x0:
        xi.append(np.log(i))
    b = 0
    for i in xi:
        b+=funKu(funU(x, i, h))
    return a*b

if __name__ == '__main__':
    x = [1000,1690,640,5120,5120,10000,20924,1690,1280,5120,5120,6250,200,100,11400,25,50,
        200,100,460,600,830,2283,5707,6076,1,1.75,2,5,11.4,250.59,490,500,600,640,640,1000,
        1000,1000,3120,7760,2600,2200,5000,1400,20,15]
    h = 0.395
    xlog = []
    for i in x:
        xlog.append(np.log(i))

    xlog =sorted(xlog)
    x1 = np.linspace(0,12,100)
    #原始数据点
    y = leijigailv(len(x))
    plt.scatter(xlog, y,c='black', marker='.',lw=1,label='雌激素效应数据\nEstrogen effects')

    #非参数核密度累积概率
    # sum = 0
    # y2 = []
    # for i in x1:
    #     sum+=fguji(len(xlog),h,i,xlog)
    #     y2.append(sum)
    # plt.plot(x1,y2)
    #正态分布累积概率
    y1 = []
    for i in x1:
        result, error = quad(fguji, 0, i)
        y1.append(result)
    y2 = []
    for i in x1:
        result, error = quad(gaussian, 0, i)
        y2.append(result)
    y3 = []
    for i in x1:
        result, error = quad(logistic, 0, i)
        y3.append(result)
    y4 = []
    for i in x1:
        result, error = quad(loggaussian, 0, i)
        y4.append(result)    
    y5 = []
    for i in x1:
        result, error = quad(logLogistic, 0, i)
        y5.append(result)
    y6 = [0.05]*len(x1)
    
    plt.plot(x1,y1,c='black',label='非参数核密度估计\nNon-parametric kernel density estimation ',linewidth=2,linestyle='-')
    plt.plot(x1,y2,c='black',label='Normal',linewidth=1,linestyle='-.')
    #plt.plot(x1,y3,c='black',label='Logistic',linewidth=1,linestyle=':')
    plt.plot(x1,y4,c='black',label='Log-normal',linewidth=0.5,linestyle=':')
    plt.plot(x1,y5,c='black',label='Log-logistic',linewidth=1,marker=',')
    plt.plot(x1,y6,c='black',label='HC$_5$',linewidth=0.5,linestyle='--')

    plt.xlabel("对数化后的雌激素效应数据/(μg·L$^{{\scr -1}}$)\nLogarithmic data of estrogen effects/(μg·L$^{{\scr -1}}$)",size = 15) 
    plt.ylabel("累积概率密度\nCumulative probability density",size = 15)

    plt.legend(frameon=False,loc="upper left")
    plt.tight_layout()
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.25, hspace=0.25)

    # 取消顶部和右侧边框
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.show()
