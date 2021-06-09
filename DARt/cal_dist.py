import numpy as np
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
