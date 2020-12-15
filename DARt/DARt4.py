#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, uniform, poisson
import random
import pandas as pd
from tqdm import tqdm
import math


class DARt:

    Num = 200
    sigma_R = 0.1
    MG = 1

    def __init__(self, GT, D_s, filename):
        if len(GT) > 10:
            GT_cut = GT[:10]
            self.GT = [g / sum(GT_cut) for g in GT_cut]
        else:
            for i in range(10 - len(GT)):
                GT.append(0)
            self.GT = GT

        delay_start = next((i for (i, x) in enumerate(D_s) if x), None)
        self.delay_start_day = delay_start
        nonzero_delay = D_s[delay_start:]
        if len(nonzero_delay) < 10:
            for i in range(10 - len(nonzero_delay)):
                nonzero_delay.append(0)
                self.delay0 = nonzero_delay
        else:
            delay = nonzero_delay[:10]
            self.delay0 = [d / sum(delay) for d in delay]
        self.len_delay0 = len(self.delay0)
        self.len_GT = len(self.GT)

        self.filename = filename

        self.Ct = [0,0,0,0,0,0,0,0,0,0]
        self.Ct_SD = []
        #10 days are added as placeholders ,because our method will start inference from 11 days.
        self.date = [1,2,3,4,5,6,7,8,9,10]
        self.N = 0

        self.Rt_back = []
        self.It_back = []
        self.M_back = []
        self.Rt_back_CI = []
        self.It_back_CI = []

        self.It_est = []
        self.Rt_est = []
        self.It_CI = []
        self.Rt_CI = []
        self.M_est = []

    def readfile(self):
        dat0 = pd.read_csv(self.filename)
        #dat = dat0.iloc[::-1]
        dat = dat0.iloc[::]
        self.Ct = self.Ct+(list(dat.iloc[:, 1]))
        self.date = self.date+(list(dat.iloc[:, 0]))
        np.array(self.Ct)
        Ct_mean = np.convolve(self.Ct, np.ones((7,)) / 7, mode='same')
        ct_diff = Ct_mean - self.Ct
        ct_diff2 = [ct ** 2 for ct in ct_diff]
        Ct_Var = np.convolve(ct_diff2, np.ones((7,)) / 7, mode='same')
        self.Ct_SD = [math.sqrt(ct) for ct in Ct_Var]
        self.N = len(self.Ct)
        self.Rt_back = np.zeros(self.N - self.delay_start_day - 1)
        self.It_back = np.zeros(self.N - self.delay_start_day - 1)
        self.M_back = np.zeros(self.N - self.delay_start_day - 1)
        self.Rt_back_CI = np.zeros((2, self.N - self.delay_start_day - 1))
        self.It_back_CI = np.zeros((2, self.N - self.delay_start_day - 1))

        self.It_est = np.zeros(self.N - self.delay_start_day - 1)
        self.Rt_est = np.zeros(self.N - self.delay_start_day - 1)
        self.It_CI = np.zeros((2, self.N - self.delay_start_day - 1))
        self.Rt_CI = np.zeros((2, self.N - self.delay_start_day - 1))
        self.M_est = np.zeros(self.N - self.delay_start_day - 1)

# filtering ============================

# total number of time steps for Jt and RT; it is shorter than Ct due to observation delay
# N = len(Ct)-7
# print(N)

    def renewal_fun(self, GT, Jt, i):
        inc_i = 0
        for j in range(len(GT)):
            inc_i = inc_i + Jt[i - j - 1] * GT[j]
        return inc_i

    def independent_sample(self, fn_list):
        def sample_fn(n):
            return np.stack([fn(n) for fn in fn_list]).T
        return sample_fn

    def cal_CI(self, particle, weight):
        order = np.argsort(particle)
        order_particle = particle[order]
        cdf = np.cumsum(weight[order])
        low_id = np.where(cdf >= 0.05)[0][0]
        high_id = np.where(cdf >= 0.95)[0][0]
        CI = [order_particle[low_id], order_particle[high_id]]
        return CI

    def initialize_particle(self, Num):
        independent_sample_list = [uniform(1, 5).rvs,uniform(0, 1).rvs]
        for i in range(self.len_GT):
            independent_sample_list.insert(0, uniform(1, max(self.Ct[0]*2,10)).rvs)
        prior_fn = self.independent_sample(independent_sample_list)
        # particle location
        particle_previous_x = prior_fn(Num)
        # particle weight
        particle_previous_w = 1 / Num * np.ones(Num)

        for n in range(len(particle_previous_x)):
            if particle_previous_x[n][-1] < 0.95:
                #     if particle_previous_x[n][-1]<2:
                particle_previous_x[n][-1] = 0
            else:
                particle_previous_x[n][-1] = 1
        return particle_previous_x, particle_previous_w

    def filtering(self, Num, delay_start_day, N, particle_previous_w, particle_previous_x):
        # filtering
        all_particle = []
        all_weight = []
        # for t in range(7,11):
        for t in tqdm(range(10, N - delay_start_day - 1)):
            # resampling-------------
            particle_previous_selection = np.zeros(Num, dtype=int)
            for s in range(Num):
                u = random.uniform(0, 1)
                particle_previous_selection[s] = np.where(np.cumsum(particle_previous_w) > u)[0][0]
            particle_previous_selection = np.sort(particle_previous_selection)

            particle_previous_resampled = np.zeros(particle_previous_x.shape)
            for s in range(Num):
                particle_previous_resampled[s, ] = particle_previous_x[particle_previous_selection[s],]

            # transition--------------
            particle_current_x = np.zeros(particle_previous_x.shape)
            for s in range(Num):
                It_s = particle_previous_resampled[s, :10]
                Rt_s = particle_previous_resampled[s, 10]
                M_s = particle_previous_resampled[s, 11]

                rdn = uniform(0, 1).rvs(1)
                if rdn < 0.95:
                    M_s_new = 0
                else:
                    M_s_new = 1
                Rt_s_new = abs(norm(loc=Rt_s, scale=DARt.sigma_R * DARt.MG ).rvs(size=1)[0])
                if M_s_new == 1:
                    Rt_s_new = uniform(0, Rt_s + 0.5).rvs(1)

                It_end = max(Rt_s_new * sum(np.multiply(It_s, self.GT[::-1])), 0)
                It_end_new = poisson.rvs(It_end, size=1)[0]
                It_s_new = np.append(It_s[1:], It_end_new)

                particle_current_x[s, :10] = It_s_new
                particle_current_x[s, 10] = Rt_s_new
                particle_current_x[s, 11] = M_s_new

            # weight------------------
            particle_current_w = np.zeros(Num)
            for s in range(Num):
                It_s_new = particle_current_x[s, :10]
                Ct_s = sum(np.multiply(It_s_new, self.delay0[::-1]))
                particle_current_w[s] = norm(Ct_s, self.Ct_SD[t + delay_start_day + 1]).pdf(self.Ct[t + delay_start_day + 1])

            # normalize------------------
            particle_current_w = particle_current_w / sum(particle_current_w)

            particle_previous_x = particle_current_x
            particle_previous_w = particle_current_w

            # save mean
            self.It_est[t] = sum(np.multiply(particle_current_x[:, 9], particle_current_w))
            self.Rt_est[t] = sum(np.multiply(particle_current_x[:, 10], particle_current_w))
            self.M_est[t] = sum(np.multiply(particle_current_x[:, 11], particle_current_w))
            # save confidence interval
            self.It_CI[:, t] = self.cal_CI(particle_current_x[:, 9], particle_current_w)
            self.Rt_CI[:, t] = self.cal_CI(particle_current_x[:, 10], particle_current_w)

            # save all paticles
            all_particle.append(particle_current_x)
            all_weight.append(particle_current_w)

        return all_particle, all_weight

# smoothing ========================
    def prob_fun(self, X_tPlus, X_t, SD_t):
        R_tPlus = X_tPlus[10]
        R_t = X_t[10]
        M_Plus_t = X_tPlus[11]
        if M_Plus_t == 1:
            Prob_Rt = 1/5
        else:
            Prob_Rt = norm.pdf(R_tPlus, loc=R_t, scale=DARt.sigma_R * DARt.MG)

        I_tPlus = X_tPlus[:10]
        I_t = X_t[:10]

        mu = R_tPlus * sum(np.multiply(I_t, self.GT[::-1]))
        Prob_It = poisson.pmf(I_tPlus[-1], mu)
        prob = Prob_It * Prob_Rt
        #     return Prob_It
        return prob

    def smoothing(self, Num, delay_start_day, N, all_particle, all_weight):
        particle_next = all_particle[-1]  # N-4
        weight_next = all_weight[-1]

        # tail
        for i in range(10):
            self.It_back[-i - 1] = sum(np.multiply(particle_next[:, 9 - i], weight_next))
            self.It_back_CI[:, -i - 1] = self.cal_CI(particle_next[:, 9 - i], weight_next)

        # tail
        self.Rt_back_CI[:, -1] = self.Rt_CI[:, -1]
        self.Rt_back[-1] = self.Rt_est[-1]

        for t in tqdm(range(N - delay_start_day - 3, 9, -1)):
            weight_t = all_weight[t - 10]  # all_weight's first element is save at t=7
            particle_t = all_particle[t - 10]

            prob_mat = np.zeros((Num, Num))
            for k in range(Num):
                for n in range(Num):
                    prob_mat[k, n] = self.prob_fun(particle_next[k,], particle_t[n,], self.Ct_SD[t + delay_start_day + 1 + 1])

            weight_now = np.zeros(Num)
            for i in range(Num):
                sum_update = 0
                for k in range(Num):
                    fs = prob_mat[k, i]
                    v = 0
                    for n in range(Num):
                        v = v + weight_t[n] * prob_mat[k, n]
                    sum_update = sum_update + weight_next[k] * fs / v

                weight_now[i] = weight_t[i] * sum_update
            weight_next = weight_now / sum(weight_now)

            # save smoothed mean
            self.Rt_back[t] = sum(np.multiply(particle_t[:, 10], weight_next))
            #self.It_back[t - 9] = sum(np.multiply(particle_t[:, 0], weight_next))
            self.It_back[t] = sum(np.multiply(particle_t[:, 9], weight_next))
            self.M_back[t] = sum(np.multiply(particle_t[:, 11], weight_next))
            # save smoothed CI
            self.Rt_back_CI[:, t] = self.cal_CI(particle_t[:, 10], weight_next)
            #self.It_back_CI[:, t - 9] = self.cal_CI(particle_t[:, 0], weight_next)
            self.It_back_CI[:, t] = self.cal_CI(particle_t[:, 9], weight_next)

            particle_next = particle_t
        #     particle_next

    def cal_r(self):
        self.readfile()
        particle_previous_x, particle_previous_w = self.initialize_particle(DARt.Num)
        all_particle, all_weight = self.filtering(DARt.Num, self.delay_start_day, self.N, particle_previous_w, particle_previous_x)
        self.smoothing(DARt.Num, self.delay_start_day, self.N, all_particle, all_weight)
        PF_result = pd.DataFrame({'Date': self.date[10:(-1 - self.delay_start_day)],
                                  'Rt_smooth': self.Rt_back[10:], 'Rt_smooth_5%': self.Rt_back_CI[0, 10:],
                                  'Rt_smooth_95%': self.Rt_back_CI[1, 10:],
                                  'Rt_filter': self.Rt_est[10:], 'Rt_est_5%': self.Rt_CI[0, 10:], 'Rt_est_95%': self.Rt_CI[1, 10:],
                                  'Mt_smooth': self.M_back[10:], 'Mt_filter': self.M_est[10:],
                                  'Jt_smooth': self.It_back[10:], 'Jt_smooth_5%': self.It_back_CI[0, 10:],
                                  'Jt_smooth_95%': self.It_back_CI[1, 10:],
                                  'Jt_filter': self.It_est[10:], 'Jt_est_5%': self.It_CI[0, 10:], 'Jt_est_95%': self.It_CI[1, 10:],
                                  'Ct': np.array(self.Ct)[10:(-1 - self.delay_start_day)]})

        groundtruth = pd.read_csv('../Simulation2_PF_data_Gaussian_noise-2.csv')
        plt.figure()
        day = np.array(range(0, len(PF_result['Rt_smooth'])))
        day2 = np.array(range(0, len(groundtruth['Rt'])))
        print(len(day))


        plt.plot(day, PF_result['Rt_smooth'], 'k')
        plt.plot(day2, groundtruth['Rt'], 'r')
        plt.fill_between(day, PF_result['Rt_smooth_5%'], PF_result['Rt_smooth_95%'], color='k', alpha=0.1)

        plt.legend(['Smoothing','groundtruth', 'CI'])
        plt.xlabel('Date')
        plt.ylabel('Rt')
        plt.xticks(day[0::10], PF_result['Date'][0::10], rotation=60)

        plt.grid('on')


        plt.ylim([0, 10])

        plt.figure()
        plt.plot(day, PF_result['Ct'], color=u'#ff7f0e')
        plt.plot(day, PF_result['Jt_smooth'], u'#1f77b4')

        plt.xticks(day[0::10], PF_result['Date'][0::10], rotation=60)
        plt.grid('on')

        plt.figure()
        plt.plot(day, PF_result['Mt_smooth'], 'k')
        plt.xticks(day[0::10], PF_result['Date'][0::10], rotation=60)
        plt.grid('on')
        plt.show()


if __name__ == "__main__":
    #GT = [0, 0, 0.165720874545241, 0.226350757019051, 0.245007574714227, 0.213515210247327,0.149405583474155]  # 1:7; non-zero index 3:7
    #D_s = [0, 0, 0, 0, 0, 0, 0.0996906, 0.1130266, 0.1143032, 0.1069238, 0.0937167999999999]
    #np.random.seed(12345)
    #random.seed(12345)
    D_s = [0, 0, 0.143259666443208, 0.254446153689323, 0.262700851779991, 0.204158756217722,0.135434571869756]  # 1:7; non-zero index 3:7 # attention index 0 is not considered here
    # generation time
    GT = [0, 0, 0.165720874545241, 0.226350757019051, 0.245007574714227, 0.213515210247327,0.149405583474155]  # 1:7; non-zero index 3:7
    #dart = DARt(GT=GT, D_s=D_s, filename='../uk_report_1011.csv')
    dart = DARt(GT=GT, D_s=D_s, filename='../simulation_data-2.csv')
    dart.cal_r()











