# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 08:41:34 2018

@author: yiyuezhuo
"""

import numpy as np
import scipy.stats as stats

DEBUG = False

class TestPower:
    def __init__(self,alpha=0.05,power=0.95,N=None,reject_region=None, type='le'):
        #self.effect_size = effect_size
        self.alpha = alpha
        self.power = power
        self.N = N
        self.reject_region = reject_region
        self.type = type
    def get_H0_statistic_distribution(self):
        raise NotImplemented
    def get_H1_statistic_distribution(self):
        raise NotImplemented
    def get_reject_region(self):
        '''
        le:
        H0: \mu <= \mu_0
        H1: \mu >  \mu_1
        
        ge:
        H0: \mu >= \mu_0
        H1: \mu <  \mu_1
        
        eq:
        H0: \mu == \mu_0
        H1: \mu != \mu_1
        '''
        alpha,type = self.alpha,self.type
        dis = self.get_H0_statistic_distribution()
        
        if type == 'le':
            cut = dis.ppf(1-alpha)
        elif type == 'ge':
            cut = dis.ppf(alpha)
        elif type == 'eq':
            cut = (dis.ppf(alpha/2),1-dis.ppf(alpha/2))
        
        reject_region = dict(type=type,value=cut)
        return reject_region
    def get_power(self):
        reject_region = self.reject_region
        dis = self.get_H1_statistic_distribution()
        
        if reject_region['type'] == 'le':
            power = 1-dis.cdf(reject_region['value'])
        elif reject_region['type'] == 'ge':
            power = dis.cdf(reject_region['value'])
        elif reject_region['type'] == 'eq':
            power = dis.cdf(reject_region['value'][0]) + (1 - dis.cdf(reject_region['value'][1]))
        
        return power

    def get_N(self):
        power,_N,_reject_region = self.power,self.N,self.reject_region
        
        def N_to_power(N):
            self.N = N
            self.reject_region = self.get_reject_region()
            return self.get_power()
        
        N_min = 1
        N_max = 10
        while(N_to_power(N_max) < power):
            N_min,N_max = N_max,N_max*2
        while(N_min < N_max):
            N_center = (N_max+N_min)//2
            _power = N_to_power(N_center)
            if DEBUG:
                print(f'check {N_center} in [{N_min},{N_max}] {_power} ')
            if _power < power:
                N_min = N_center + 1
            elif _power > power:
                N_max = N_center - 1
            else:
                N = N_center # consider the power is a Continuous variable, it can only be realized for a minor probility
                break
        else:
            N = max(N_min,N_max)
        
        if N_to_power(N) < self.power:
            N += 1
        
        self.N = _N
        self.reject_region = _reject_region
        
        return N

class CohenTestPower(TestPower):
    # define Cohen style effect size distilled from H1 parameters
    def __init__(self,effect_size = 0.5,**kwargs):
        super().__init__(**kwargs)
        self.effect_size = effect_size
    def get_effect_size(self):
        raise NotImplemented
        
class TTestPower(CohenTestPower):
    def get_H0_statistic_distribution(self):
        N = self.N
        return stats.t(N-1)
    def get_H1_statistic_distribution(self):
        effect_size,N = self.effect_size,self.N
        ncp = self.get_ncp() #effect_size * np.sqrt(N)
        return stats.nct(N-1,ncp)
    def get_ncp(self):
        raise NotImplemented
    
class OneSampleTTestPower(TTestPower):
    '''
    H0: \mu - \mu_0 = 0
    H1: \mu - \mu_0 \neq 0
    
    t-statistic:
    t = \frac{\bar{X} - \mu_0}{\frac{s}{n}}
    '''
    def __init__(self,effect_size = 0.5,alpha=0.05,power=0.95,N = None,reject_region=None, type='le', 
                 mu=None,mu0=0.0,sigma=None):
        super().__init__(effect_size = effect_size, alpha = alpha, 
             power = power, N = N, reject_region = reject_region, type = type)
        
        self.mu0 = mu0
        self.mu = mu
        self.sigma = sigma
    def get_ncp(self):
        effect_size,N = self.effect_size,self.N
        return effect_size * np.sqrt(N)
    def get_effect_size(self):
        # Cohen'd \frac{\mu - \mu_0}{\sigma}
        mu,mu0,sigma = self.mu,self.mu0,self.sigma
        
        effect_size = (mu - mu0)/sigma
        return effect_size
    def get_mu(self):
        # ncp don't require exact H1, and it only require effect size, the  exact H1 require sigma
        effect_size,sigma,mu0 = self.effect_size,self.sigma,self.mu0
        
        mu = effect_size * sigma + mu0
        return mu

class TwoSampleTTestPower(TTestPower):
    '''
    H0: \mu_1 - \mu_2 = 0
    H1: \mu_1 - \mu_2 \neq 0
    
    t-statistic:
    t = \frac{\frac{1}{\sqrt{\frac{1}{N_1}+\frac{1}{N_2}}}(\bar{X}_1 - \bar{X_2})}{\frac{1}{N_1+N_2-2}((N_1 - 1)s_1^2 + (N_2 - 1)s_2^2)}
    '''
    def __init__(self,effect_size = 0.5,alpha=0.05,power=0.95,N = None,reject_region=None, type='le', 
                 mu1=0.0, mu2=1.0, sigma1=0.5, sigma2=0.5, allocation_p=(0.5,0.5)):
        super().__init__(effect_size = effect_size, alpha = alpha, power = power, N = N, reject_region = reject_region, type = type)
        
        self.mu1 = mu1
        self.mu2 = mu2
        self.sigma1 = sigma1
        self.sigma2 = sigma2
        self.allocation_p = allocation_p
    def get_ncp(self):
        effect_size,N,allocation_p = self.effect_size,self.N,self.allocation_p
        n1,n2 = int(self.N * allocation_p[0]) , int(self.N * allocation_p[1])
        return effect_size * np.sqrt((n1*n2)/(n1+n2))
    def get_effect_size(self):
        # Cohen'd \frac{\mu_1 - \mu_2}{\sigma} (for \sigma = \sigma_1 = \sigma_2, exact)
        # \sigma = \sqrt{\frac{\sigma_1^2 + \sigma_2^2}{2}} (for n_1 = n_2,approximated)
        # \sigma = \sqrt{\frac{n_1 \sigma_1^2 + n_2 \sigma_2}{n_1+n_2}} (relative general, but require n_1,n_2,approximated)
        # \sigma = \sqrt{p_1 \sigma_1^2 + p_2 \sigma_2} ( general, not require n_1,n_2,approximated)
        mu1,mu2,sigma1,sigma2,allocation_p = self.mu1, self.mu2, self.sigma1, self.sigma2, self.allocation_p
        
        sigma = np.sqrt(allocation_p[0] * sigma1 + allocation_p[1] * sigma2)
        effect_size = (mu1 - mu2)/sigma
        return effect_size
