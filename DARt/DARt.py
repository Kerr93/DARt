#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, gamma, uniform, poisson
import random
import pandas as pd
from tqdm import tqdm
import math

sigma_R = 0.3
MG = 1
# generation time
GT = [0,0,0.165720874545241,0.226350757019051,0.245007574714227,0.213515210247327,0.149405583474155]#1:7; non-zero index 3:7


D_s = np.array([0, 0,0, 0, 0,0, 0.0996906, 0.1130266, 0.1143032, 0.1069238,0.0937167999999999])
delay = list(D_s/sum(D_s))
dealy_start_day = 6
delay0 = delay[dealy_start_day:]+[0]* 2# shift two time steps


def readfile(datapath):
    # adjust the tail of onset curve using "delay_dist"
    # data input
    # sample_name: 'uk_report_1011.csv'
    dat0 = pd.read_csv(datapath)
    dat = dat0.iloc[::-1]
    Ct = list(dat.iloc[240:,1])
    date = list(dat.iloc[240:,0])
    #Ct = list(dat.iloc[:, 1])
    #date = list(dat.iloc[:, 0])
    np.array(Ct)
    Ct_mean = np.convolve(Ct, np.ones((7,)) / 7, mode='same')
    ct_diff = Ct_mean - Ct
    ct_diff2 = [ct ** 2 for ct in ct_diff]
    Ct_Var = np.convolve(ct_diff2, np.ones((7,)) / 7, mode='same')
    Ct_SD = [math.sqrt(ct) for ct in Ct_Var]
    return Ct, Ct_SD, date

filename = '/Users/kerr/Documents/GitHub/Rt/uk_report_1011.csv'
Ct, Ct_SD, date = readfile(filename)

# filtering ============================

# total number of time steps for Jt and RT; it is shorter than Ct due to observation delay
# N = len(Ct)-7
# print(N)
N = len(Ct)

# renewal process: pay attention to index
def renewal_fun(GT,Jt,i):
    inc_i = 0
    for j in range(len(GT)):
        inc_i = inc_i + Jt[i-j-1]*GT[j]     
    return inc_i

def independent_sample(fn_list):
    def sample_fn(n):
        return np.stack([fn(n) for fn in fn_list]).T

    return sample_fn

def cal_CI(particle,weight):
    order = np.argsort(particle)
    order_particle = particle[order]
    cdf = np.cumsum(weight[order])
    low_id = np.where(cdf>=0.05)[0][0]
    high_id = np.where(cdf>=0.95)[0][0]
    CI = [order_particle[low_id],order_particle[high_id]]
    return CI


Num = 200


def initialize_particle(Num):
# initialize particles
    prior_fn = independent_sample(
        [

            uniform(0,10).rvs,
            uniform(0,10).rvs,
            uniform(0,10).rvs,
            uniform(0,10).rvs,
            uniform(0,10).rvs,
            uniform(0,10).rvs,
            uniform(0,10).rvs,
            uniform(1,5).rvs,
            uniform(0,1).rvs
        ]
    )
    # particle location
    particle_previous_x = prior_fn(Num)#number of particles*8; the first 7 columns are for It; the last one is for Rt
    # particle weight
    particle_previous_w = 1/Num*np.ones(Num)

    for n in range(len(particle_previous_x)):
        if particle_previous_x[n][-1]<0.95:
    #     if particle_previous_x[n][-1]<2:
            particle_previous_x[n][-1] = 0
        else:
            particle_previous_x[n][-1] = 1
    return particle_previous_x, particle_previous_w



particle_previous_x, particle_previous_w = initialize_particle(Num)


def filtering(Num,dealy_start_day,N,particle_previous_w,particle_previous_x):
    # filtering
    all_particle = []
    all_weight = []
    It_est = np.zeros(N-dealy_start_day-1)
    Rt_est = np.zeros(N-dealy_start_day-1)
    It_CI = np.zeros((2,N-dealy_start_day-1))
    Rt_CI = np.zeros((2,N-dealy_start_day-1))
    M_est = np.zeros(N-dealy_start_day-1)
    # for t in range(7,11):
    for t in tqdm(range(7,N-dealy_start_day-1)):
        # resampling-------------
        particle_previous_selection = np.zeros(Num, dtype=int)
        for s in range(Num):
            u = random.uniform(0, 1)
            particle_previous_selection[s] = np.where(np.cumsum(particle_previous_w)>u)[0][0]
        particle_previous_selection = np.sort(particle_previous_selection)

        particle_previous_resampled = np.zeros(particle_previous_x.shape)
        for s in range(Num):
            particle_previous_resampled[s,] = particle_previous_x[particle_previous_selection[s],]

        # transition--------------
        particle_current_x = np.zeros(particle_previous_x.shape)
        for s in range(Num):
            It_s = particle_previous_resampled[s,:7]
            Rt_s = particle_previous_resampled[s,7]
            M_s = particle_previous_resampled[s,8]

            rdn = uniform(0,1).rvs(1)
            if rdn < 0.95:
                M_s_new = 0
            else:
                M_s_new = 1
    #         M_s_new = 0
            Rt_s_new = abs(norm(loc=Rt_s, scale=sigma_R*MG*0.1).rvs(size=1)[0])
    #         if t>20 and M_s_new == 1:
            if M_s_new == 1:
                Rt_s_new = uniform(0,4).rvs(1)

            It_end = max(Rt_s_new*sum(np.multiply(It_s,GT[::-1])),0)
    #         It_end_new = poisson.rvs(It_end, size=1)[0]
            It_end_new = It_end
            It_s_new = np.append(It_s[1:],It_end_new)

            particle_current_x[s,:7] = It_s_new
            particle_current_x[s,7] = Rt_s_new
            particle_current_x[s,8] = M_s_new

        # weight------------------
        particle_current_w = np.zeros(Num)
        for s in range(Num):
            It_s_new = particle_current_x[s,:7]
            Ct_s = sum(np.multiply(It_s_new,delay0[::-1]))
    #         particle_current_w[s] = poisson.pmf(Ct[t+3], Ct_s)
    #         particle_current_w[s] = pow(Ct_s,Ct[t+3])*math.exp(-Ct_s)
    #         particle_current_w[s] = Ct[t+3]*math.log(Ct_s)-Ct_s
    #         particle_current_w[s] = norm(Ct_s, math.sqrt(Ct_s)*2).pdf(Ct[t+3])
    #         particle_current_w[s] = norm(Ct[t+3], math.sqrt(Ct[t+3])*2).pdf(Ct_s)
            particle_current_w[s] = norm(Ct_s, Ct_SD[t+dealy_start_day+1]).pdf(Ct[t+dealy_start_day+1])

        # normalize------------------
    #     particle_current_w = particle_current_w-max(particle_current_w)
    #     particle_current_w = np.array([math.exp(w_i) for w_i in particle_current_w])
        particle_current_w = particle_current_w/sum(particle_current_w)

        particle_previous_x = particle_current_x
        particle_previous_w = particle_current_w

        # save mean
        It_est[t] = sum(np.multiply(particle_current_x[:,6],particle_current_w))
        Rt_est[t] = sum(np.multiply(particle_current_x[:,7],particle_current_w))
        M_est[t] = sum(np.multiply(particle_current_x[:,8],particle_current_w))
        # save confidence interval
        It_CI[:,t] = cal_CI(particle_current_x[:,6],particle_current_w)
        Rt_CI[:,t] = cal_CI(particle_current_x[:,7],particle_current_w)

        # save all paticles
        all_particle.append(particle_current_x)
        all_weight.append(particle_current_w)
    return  all_particle, all_weight,It_est,Rt_est,M_est,It_CI,Rt_CI

all_particle,all_weight,It_est,Rt_est,M_est,It_CI,Rt_CI = filtering(Num,dealy_start_day,N,particle_previous_w,particle_previous_x)



# smoothing ========================

def prob_fun(X_tPlus,X_t,SD_t):
    R_tPlus = X_tPlus[7]
    R_t = X_t[7]
    M_Plus_t = X_tPlus[8]
    if M_Plus_t == 1:
        Prob_Rt = 1
    else:
        Prob_Rt = norm.pdf(R_tPlus, loc=R_t, scale=sigma_R*MG)
    
    I_tPlus = X_tPlus[:7]
    I_t = X_t[:7]   

    mu = R_tPlus*sum(np.multiply(I_t,GT[::-1]))
#     Prob_It = poisson.pmf(I_tPlus[-1], mu)
    Prob_It = norm(mu, SD_t).pdf(I_tPlus[-1])
    prob = Prob_It*Prob_Rt
#     return Prob_It
    return prob


def smoothing(Num,dealy_start_day,N,all_particle,all_weight):
    particle_next = all_particle[-1]#N-4
    weight_next = all_weight[-1]
    Rt_back = np.zeros(N-dealy_start_day-1)
    It_back = np.zeros(N-dealy_start_day-1)
    M_back = np.zeros(N-dealy_start_day-1)
    Rt_back_CI = np.zeros((2,N-dealy_start_day-1))
    It_back_CI = np.zeros((2,N-dealy_start_day-1))

    #tail
    for i in range(7):
        It_back[-i-1] = sum(np.multiply(particle_next[:,6-i],weight_next))
        It_back_CI[:,-i-1] = cal_CI(particle_next[:,6-i],weight_next)

    #tail
    Rt_back_CI[:,-1] = Rt_CI[:,-1]
    Rt_back[-1] = Rt_est[-1]


    for t in tqdm(range(N-dealy_start_day-3,6,-1)):
        weight_t = all_weight[t-7]# all_weight's first element is save at t=7
        particle_t = all_particle[t-7]

        prob_mat = np.zeros((Num,Num))
        for k in range(Num):
            for n in range(Num):
                prob_mat[k,n] = prob_fun(particle_next[k,],particle_t[n,],Ct_SD[t+dealy_start_day+1+1])

        weight_now = np.zeros(Num)
        for i in range(Num):
            sum_update = 0
            for k in range(Num):
                fs = prob_mat[k,i]
                v = 0
                for n in range(Num):
                    v = v + weight_t[n]*prob_mat[k,n]
                sum_update = sum_update + weight_next[k]*fs/v

            weight_now[i] = weight_t[i]*sum_update
        weight_next = weight_now/sum(weight_now)

        # save smoothed mean
        Rt_back[t] = sum(np.multiply(particle_t[:,7],weight_next))
        It_back[t-6] = sum(np.multiply(particle_t[:,0],weight_next))
        M_back[t] = sum(np.multiply(particle_t[:,8],weight_next))
        # save smoothed CI
        Rt_back_CI[:,t] = cal_CI(particle_t[:,7],weight_next)
        It_back_CI[:,t-6] = cal_CI(particle_t[:,0],weight_next)

        particle_next = particle_t
    #     particle_next
    return Rt_back,It_back,M_back,Rt_back_CI,It_back_CI


Rt_back,It_back,M_back,Rt_back_CI,It_back_CI = smoothing(Num,dealy_start_day,N,all_particle,all_weight)




# In[92]:


PF_result = pd.DataFrame({'Date':date[7:(-1-dealy_start_day)],
                          'Rt_smooth':Rt_back[7:],'Rt_smooth_5%':Rt_back_CI[0,7:],'Rt_smooth_95%':Rt_back_CI[1,7:],
                      'Rt_filter':Rt_est[7:],'Rt_est_5%':Rt_CI[0,7:],'Rt_est_95%':Rt_CI[1,7:],
                      'Mt_smooth':M_back[7:], 'Mt_filter':M_est[7:],   
                     'Jt_smooth':It_back[7:],'Jt_smooth_5%':It_back_CI[0,7:],'Jt_smooth_95%':It_back_CI[1,7:],
                      'Jt_filter':It_est[7:],'Jt_est_5%':It_CI[0,7:],'Jt_est_95%':It_CI[1,7:], 'Ct':np.array(Ct)[7:(-1-dealy_start_day)]})



# In[97]:


plt.figure()
day = np.array(range(0, len(PF_result['Rt_smooth'])))

plt.plot(day, PF_result['Rt_smooth'], 'k')
plt.fill_between(day, PF_result['Rt_smooth_5%'], PF_result['Rt_smooth_95%'], color='k', alpha=0.1)
plt.legend(['Smoothing', 'CI'])
plt.xlabel('Date')
plt.ylabel('Rt')
plt.xticks(day[1::14], PF_result['Date'][1::14], rotation=60)

plt.grid('on')

# plt.xlim([1, 70])
plt.ylim([0, 3])


plt.figure()
plt.plot(day, PF_result['Ct'], color=u'#ff7f0e')
plt.plot(day, PF_result['Jt_smooth'],u'#1f77b4' )

plt.xticks(day[1::14], PF_result['Date'][1::14], rotation=60)
plt.grid('on')
# plt.xlim([1, 70])


plt.figure()
plt.plot(day, PF_result['Mt_smooth'], 'k')
plt.xticks(day[1::14], PF_result['Date'][1::14], rotation=60)
plt.grid('on')
# plt.xlim([1, 70])
plt.show()





