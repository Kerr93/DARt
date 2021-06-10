#!/usr/bin/env python
# coding: utf-8

import math
import random
import warnings
import matplotlib
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.stats import norm, uniform, poisson
from statsmodels.distributions.empirical_distribution import ECDF


def discrete_weibull(shape, scale, N2):
    x_pre_scale = np.random.weibull(shape, int(5e6))
    x = scale * x_pre_scale
    f = ECDF(x)
    h = np.zeros(N2)
    h[0] = f(1.5) - f(0)
    for i in range(1, N2):
        h[i] = (f(i+1.5) - f(i+0.5)) / (1-f(i+0.5))
    s = np.zeros(N2)
    s[0] = 1
    for i in range(1, N2):
        s[i] = s[i-1]*(1-h[i-1])
    SI0 = s * h
    SI1 = SI0[~np.isnan(SI0)]
    SI = np.zeros(N2)
    SI[0:len(SI1)] = SI1
    return SI


def discrete_lognormal(logmean, logsd, N2):
    x = np.random.lognormal(logmean, logsd, int(5e6))
    f = ECDF(x)
    h = np.zeros(N2)
    h[0] = f(1.5) - f(0)
    for i in range(1,N2):
        h[i] = (f(i+1.5) - f(i+0.5)) / (1-f(i+0.5))
    s = np.zeros(N2)
    s[0] = 1
    for i in range(1,N2):
        s[i] = s[i-1]*(1-h[i-1])
    SI0 = s * h
    SI1 = SI0[~np.isnan(SI0)]
    SI = np.zeros(N2)
    SI[0:len(SI1)] = SI1
    return SI


def cal_gt(shape, scale):
    dis_gt = discrete_weibull(shape, scale, 20)
    return dis_gt


def cal_inc(logmean, logsd):
    dis_inc = discrete_lognormal(logmean, logsd, 20)
    return dis_inc


def cal_rep(repmean, repsd, incmean, incsd):
    N2 = 30
    shape = (repmean**2) / (repsd**2)
    scale = repmean / shape
    x1 = np.random.gamma(shape, scale, int(5e6))
    x2 = np.random.lognormal(incmean, incsd, int(5e6))
    f = ECDF(x1+x2)
    h = np.zeros(N2)
    h[0] = f(1.5) - f(0)
    for i in range(1,N2):
        h[i] = (f(i+1.5) - f(i+0.5)) / (1-f(i+0.5))
    s = np.zeros(N2)
    s[0] = 1
    for i in range(1,N2):
        s[i] = s[i-1]*(1-h[i-1])
    SI0 = s * h
    SI1 = SI0[~np.isnan(SI0)]
    SI = np.zeros(N2)
    SI[0:len(SI1)] = SI1
    dis_rep = SI
    return dis_rep


def select_norm(SI, threshold):
    which = lambda lst: list(np.where(lst)[0])
    indexSI = which(SI > threshold)
    SI0 = np.zeros(indexSI[0])
    output = np.append(SI0, SI[indexSI])
    output_norm = output / sum(output)
    return output_norm

def cal_dis():
    gt_mean, gt_sd = map(eval, input(
        'Please input the shape and scale of the Weibull distribution for generation time, use blank as a separator \neg. 2.826 5.665, or press enter to use the default values:   ').split() or ['2.826', '5.665'])
    gt_dis = cal_gt(gt_mean, gt_sd)

    gt_cutoff_threshold = eval(input('\nPlease input the threshold of generation distribution \neg. 0.1, or press enter to use the default value:   ') or '0.1')
    gt_dis_cut = select_norm(gt_dis, gt_cutoff_threshold)

    observation_type = eval(
        input('\nPlease choose the type of input observation：press 0 for ONSET, press 1 for REPORT   ') or '1')

    if observation_type == 0:
        inc_mean, inc_sd = map(eval, input(
            '\nPlease input the logmean and logsd of the lognormal distribution for incubation time, use blank as a separator \neg.1.644 0.33, or press enter to use the default values:   ').split() or ['1.644', '0.33'])
        inc_dis = cal_inc(inc_mean, inc_sd)
        inc_cutoff_threshold = eval(input('\nPlease input the threshold of incubation distribution \neg. 0.1, or press enter to use the default value:   ') or '0.1')
        inc_dis_cut = select_norm(inc_dis, inc_cutoff_threshold)
        return gt_dis_cut, inc_dis_cut

    elif observation_type == 1:
        inc_mean, inc_sd = map(eval, input(
            '\nPlease input the logmean and logsd of the lognormal distribution for incubation time, use blank as a separator \neg.1.644 0.33, or press enter to use the default values:   ').split() or ['1.644', '0.33'])
        rep_mean, rep_sd = map(eval, input(
            '\nPlease input the mean and sd of the Gamma distribution for report time, use blank as a separator \neg.4.9 3.3, or press enter to use the default values:   ').split() or ['4.9', '3.3'])
        rep_dis = cal_rep(rep_mean, rep_sd, inc_mean, inc_sd)
        rep_cutoff_threshold = eval(input(
            '\nPlease input the threshold of report delay distribution \neg. 0.09, or press enter to use the default value:   ') or '0.09')
        rep_dis_cut = select_norm(rep_dis, rep_cutoff_threshold)
        return gt_dis_cut, rep_dis_cut

    else:
        print('Wrong input. Please rerun DARt again.')


class DARt:
    Num = 200
    sigma_R = 0.1
    MG = 1

    def __init__(self, filename):
        warnings.filterwarnings('ignore')
        GT, D_s = cal_dis()
        GT = list(GT)
        D_s = list(D_s)

        self.delay_start_day = next((i for (i, x) in enumerate(D_s) if x), None)
        nonzero_delay = D_s[self.delay_start_day:]
        if len(GT) >= len(nonzero_delay):
            for i in range(len(GT) - len(nonzero_delay)):
                nonzero_delay.append(0)
        else:
            for i in range(len(nonzero_delay) - len(GT)):
                GT.append(0)

        self.GT = GT
        self.len_GT = len(self.GT)
        self.delay0 = nonzero_delay

        self.filename = filename

        self.Ct = None
        self.Ct_SD = None
        self.date = None
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
        data = pd.read_csv(self.filename)
        Ct = data.iloc[:, 1]
        date = data.iloc[:, 0]

        first_nonzero = []
        for i in range(len(Ct)):
            if Ct.values[i] == 0:
                first_nonzero.append(i)
        if first_nonzero:
            start_index = first_nonzero[-1]
        else:
            start_index = 0

        self.Ct = list(Ct[start_index:])
        self.date = list(date[start_index:])

        Ct_mean = np.convolve(self.Ct, np.ones((self.len_GT,)) / self.len_GT, mode='same')
        ct_diff = Ct_mean - self.Ct
        ct_diff2 = [ct ** 2 for ct in ct_diff]
        Ct_Var = np.convolve(ct_diff2, np.ones((self.len_GT,)) / self.len_GT, mode='same')
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
        low_id = np.where(cdf >= 0.025)[0][0]
        high_id = np.where(cdf >= 0.975)[0][0]
        CI = [order_particle[low_id], order_particle[high_id]]
        return CI

    def initialize_particle(self, Num):
        independent_sample_list = [uniform(1, 5).rvs, uniform(0, 1).rvs]
        for i in range(self.len_GT):
            independent_sample_list.insert(0, uniform(1, max(self.Ct[0] * 2, 10)).rvs)
        prior_fn = self.independent_sample(independent_sample_list)
        # particle location
        particle_previous_x = prior_fn(Num)
        # particle weight
        particle_previous_w = 1 / Num * np.ones(Num)
        # mode
        for n in range(len(particle_previous_x)):
            if particle_previous_x[n][-1] < 0.95:
                particle_previous_x[n][-1] = 0
            else:
                particle_previous_x[n][-1] = 1
        return particle_previous_x, particle_previous_w

    def filtering(self, Num, delay_start_day, N, particle_previous_w, particle_previous_x):
        # filtering
        all_particle = []
        all_weight = []

        for t in tqdm(range(self.len_GT, N - delay_start_day - 1)):
            # resampling-------------
            particle_previous_selection = np.zeros(Num, dtype=int)
            for s in range(Num):
                u = random.uniform(0, 1)
                particle_previous_selection[s] = np.where(np.cumsum(particle_previous_w) > u)[0][0]
            particle_previous_selection = np.sort(particle_previous_selection)

            particle_previous_resampled = np.zeros(particle_previous_x.shape)
            for s in range(Num):
                particle_previous_resampled[s, ] = particle_previous_x[particle_previous_selection[s], ]

            # transition--------------
            particle_current_x = np.zeros(particle_previous_x.shape)
            for s in range(Num):
                It_s = particle_previous_resampled[s, :self.len_GT]
                Rt_s = particle_previous_resampled[s, self.len_GT]

                rdn = uniform(0, 1).rvs(1)
                if rdn < 0.95:
                    M_s_new = 0
                else:
                    M_s_new = 1
                Rt_s_new = abs(norm(loc=Rt_s, scale=DARt.sigma_R * DARt.MG).rvs(size=1)[0])
                if M_s_new == 1:
                    Rt_s_new = uniform(0, Rt_s + 0.5).rvs(1)

                It_end = max(Rt_s_new * sum(np.multiply(It_s, self.GT[::-1])), 0)
                It_end_new = poisson.rvs(It_end, size=1)[0]
                It_s_new = np.append(It_s[1:], It_end_new)

                particle_current_x[s, :self.len_GT] = It_s_new
                particle_current_x[s, self.len_GT] = Rt_s_new
                particle_current_x[s, self.len_GT+1] = M_s_new

            # weight------------------
            particle_current_w = np.zeros(Num)
            for s in range(Num):
                It_s_new = particle_current_x[s, :self.len_GT]
                Ct_s = sum(np.multiply(It_s_new, self.delay0[::-1]))
                particle_current_w[s] = norm(Ct_s, self.Ct_SD[t + delay_start_day + 1]).pdf(self.Ct[t + delay_start_day + 1])

            # normalize------------------
            particle_current_w = particle_current_w / sum(particle_current_w)
            particle_previous_x = particle_current_x
            particle_previous_w = particle_current_w

            # save mean
            self.It_est[t] = sum(np.multiply(particle_current_x[:, self.len_GT-1], particle_current_w))
            self.Rt_est[t] = sum(np.multiply(particle_current_x[:, self.len_GT], particle_current_w))
            self.M_est[t] = sum(np.multiply(particle_current_x[:, self.len_GT+1], particle_current_w))

            # save confidence interval
            self.It_CI[:, t] = self.cal_CI(particle_current_x[:, self.len_GT-1], particle_current_w)
            self.Rt_CI[:, t] = self.cal_CI(particle_current_x[:, self.len_GT], particle_current_w)

            # save all particles
            all_particle.append(particle_current_x)
            all_weight.append(particle_current_w)

        return all_particle, all_weight

    # smoothing ========================
    def prob_fun(self, X_tPlus, X_t):
        R_tPlus = X_tPlus[self.len_GT]
        R_t = X_t[self.len_GT]
        M_Plus_t = X_tPlus[self.len_GT+1]
        if M_Plus_t == 1:
            Prob_Rt = 1 / 5
        else:
            Prob_Rt = norm.pdf(R_tPlus, loc=R_t, scale=DARt.sigma_R * DARt.MG)

        I_tPlus = X_tPlus[:self.len_GT]
        I_t = X_t[:self.len_GT]

        mu = R_tPlus * sum(np.multiply(I_t, self.GT[::-1]))
        Prob_It = poisson.pmf(I_tPlus[-1], mu)
        prob = Prob_It * Prob_Rt
        return prob

    def smoothing(self, Num, N, all_particle, all_weight):
        particle_next = all_particle[-1]
        weight_next = all_weight[-1]

        for i in range(self.len_GT):
            self.It_back[-i-1] = sum(np.multiply(particle_next[:, self.len_GT-1-i], weight_next))
            self.It_back_CI[:, -i-1] = self.cal_CI(particle_next[:, self.len_GT-1-i], weight_next)

        self.Rt_back_CI[:, -1] = self.Rt_CI[:, -1]
        self.Rt_back[-1] = self.Rt_est[-1]

        for t in tqdm(range(N - self.delay_start_day-3, self.len_GT-1, -1)):
            weight_t = all_weight[t-self.len_GT]
            particle_t = all_particle[t-self.len_GT]

            prob_mat = np.zeros((Num, Num))
            for k in range(Num):
                for n in range(Num):
                    prob_mat[k, n] = self.prob_fun(particle_next[k, ], particle_t[n, ])

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

            self.Rt_back[t] = sum(np.multiply(particle_t[:, self.len_GT], weight_next))
            self.It_back[t] = sum(np.multiply(particle_t[:, self.len_GT-1], weight_next))
            self.M_back[t] = sum(np.multiply(particle_t[:, self.len_GT+1], weight_next))
            self.Rt_back_CI[:, t] = self.cal_CI(particle_t[:, self.len_GT], weight_next)
            self.It_back_CI[:, t] = self.cal_CI(particle_t[:, self.len_GT-1], weight_next)

            particle_next = particle_t

    def cal_r(self):
        self.readfile()
        particle_previous_x, particle_previous_w = self.initialize_particle(DARt.Num)
        all_particle, all_weight = self.filtering(DARt.Num, self.delay_start_day, self.N, particle_previous_w, particle_previous_x)
        self.smoothing(DARt.Num, self.N, all_particle, all_weight)
        PF_result = pd.DataFrame({'Date': self.date[self.len_GT:(-1 - self.delay_start_day)],
                                  'Rt_smooth': self.Rt_back[self.len_GT:], 'Rt_smooth_5%': self.Rt_back_CI[0, self.len_GT:],
                                  'Rt_smooth_95%': self.Rt_back_CI[1, self.len_GT:],
                                  'Rt_filter': self.Rt_est[self.len_GT:], 'Rt_est_5%': self.Rt_CI[0, self.len_GT:], 'Rt_est_95%': self.Rt_CI[1, self.len_GT:],
                                  'Mt_smooth': self.M_back[self.len_GT:], 'Mt_filter': self.M_est[self.len_GT:],
                                  'Jt_smooth': self.It_back[self.len_GT:], 'Jt_smooth_5%': self.It_back_CI[0, self.len_GT:],
                                  'Jt_smooth_95%': self.It_back_CI[1, self.len_GT:],
                                  'Jt_filter': self.It_est[self.len_GT:], 'Jt_est_5%': self.It_CI[0, self.len_GT:], 'Jt_est_95%': self.It_CI[1, self.len_GT:],
                                  'Ct': np.array(self.Ct)[self.len_GT:(-1 - self.delay_start_day)]})

        PF_result.to_csv('inference_result.csv')

    def plot(self):
        PF_result = pd.read_csv('inference_result.csv')
        
        plt.figure()
        day = np.array(range(0, len(PF_result['Rt_smooth'])))
        plt.plot(day, PF_result['Rt_smooth'], 'k')
        plt.fill_between(day, PF_result['Rt_smooth_5%'], PF_result['Rt_smooth_95%'], color='k', alpha=0.1)
        plt.legend(['Smoothing', 'CI'])
        plt.xlabel('Date')
        plt.ylabel('Rt')
        plt.xticks(day[0::10], PF_result['Date'][0::10], rotation=90)
        plt.grid('on')
        plt.ylim([0, 10])

        plt.gca().get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x))))
        plt.grid('on', axis='both', color='k', linestyle='--', linewidth=0.5, alpha=0.5)
        ax2 = plt.gca().twinx()
        day = np.array(range(0, len(PF_result['Mt_smooth'])))
        ax2.bar(day, PF_result['Mt_smooth'], color='g', alpha=0.5, width=1)

        plt.tight_layout()
        # plt.savefig('./Rt_estimate.svg')
        plt.savefig('./Rt_estimate.png', dpi=300)

        plt.figure()
        Jt0 = list(self.It_back) + [0] * (self.len_GT - 1)
        It = []
        It_05 = []
        It_95 = []
        Jt0_05 = list(self.It_back_CI[0]) + [0] * 6
        Jt0_95 = list(self.It_back_CI[1]) + [0] * 6
        for i in range(len(self.Ct)):
            It.append(Jt0[i:(i + self.len_GT)])
            It_05.append(Jt0_05[i:(i + self.len_GT)])
            It_95.append(Jt0_95[i:(i + self.len_GT)])

        pp = self.len_GT + self.delay_start_day
        Ct_infer = [0] * (self.len_GT + self.delay_start_day)
        Ct_infer_05 = [0] * (self.len_GT + self.delay_start_day)
        Ct_infer_95 = [0] * (self.len_GT + self.delay_start_day)
        for i in range(self.len_GT + self.delay_start_day, len(self.Ct)):
            Ct_i = round(sum(np.multiply(It[i - pp], self.delay0[::-1])))
            Ct_i_05 = round(sum(np.multiply(It_05[i - pp], self.delay0[::-1])))
            Ct_i_95 = round(sum(np.multiply(It_95[i - pp], self.delay0[::-1])))
            Ct_infer.append(Ct_i)
            Ct_infer_05.append(Ct_i_05)
            Ct_infer_95.append(Ct_i_95)
        plt.plot(self.Ct[self.len_GT:], label='Ct')
        plt.plot(Ct_infer[self.len_GT:], label='Recovered Ct')
        day = np.array(range(0, len(Ct_infer_05[self.len_GT:])))
        plt.fill_between(day, Ct_infer_05[self.len_GT:], Ct_infer_95[self.len_GT:], color='#1f77b4', alpha=0.2)
        day2 = np.array(range(0, len(PF_result['Jt_smooth'])))
        plt.plot(PF_result['Jt_smooth'], 'g')
        plt.fill_between(day2, PF_result['Jt_smooth_5%'], PF_result['Jt_smooth_95%'], color='g', alpha=0.2)
        plt.xticks(day2[0::10], PF_result['Date'][0::10], rotation=90)
        plt.legend(['ct', 'recovered ct', 'jt'])
        plt.grid('on')

        plt.tight_layout()
        # plt.savefig('./ct_estimate.svg')
        plt.savefig('./ct_estimate.png', dpi=300)
