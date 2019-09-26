import numpy as np
from scipy.stats import norm

mu = 0
sigma0 = 1


def d1(S,K,r,sigma,T):
    return (np.log(S/K) + (r + (sigma* sigma) / 2)* T)/(sigma* np.sqrt(T))

def d2(S,K,r,sigma,T):
    return d1(S,K,r,sigma,T) - (sigma *np.sqrt(T))

def call_price(S,K,r,sigma,T):
    return  S*norm.cdf(d1(S,K,r,sigma,T),mu,sigma0)- K * np.exp(-r*T) *norm.cdf(d2(S,K,r,sigma,T),mu,sigma0)

def put_price(S,K,r,sigma,T):
    return K * np.exp(-r*T) * norm.cdf(-d2(S,K,r,sigma,T),mu,sigma0)-S * norm.cdf(-d1(S,K,r,sigma,T),mu,sigma0)

def call_delta(S,K,r,sigma,T):
    return norm.cdf(d1(S,K,r,sigma,T),mu,sigma0)

def put_delta(S,K,r,sigma,T):
    return norm.cdf(d1(S,K,r,sigma,T),mu,sigma0)-1

def gamma(S,K,r,sigma,T):
    return norm.pdf(d1(S,K,r,sigma,T),mu,sigma0) /(S * sigma * np.sqrt(T))

def vega(S,K,r,sigma,T):
    return S * norm.pdf(d1(S,K,r,sigma,T),mu,sigma0) *np.sqrt(T)

def call_theta(S,K,r,sigma,T):
    return  -r* K* np.exp(-r* T)*norm. cdf(d2(S, K, r, sigma, T), mu, sigma0)-\
    (S* norm.pdf(d1(S, K, r, sigma, T), mu, sigma0)*sigma)/ (2 * np.sqrt(T))

def put_theta(S,K,r,sigma,T):
    return r * K * np.exp(-r*T) * norm.cdf(-d2(S,K,r,sigma,T),mu,sigma0) - \
                         (S * norm.pdf(d1(S,K,r,sigma,T),mu,sigma0) * sigma) / (2*np.sqrt(T))

def call_rho(S,K,r,sigma,T):
    return K * T * np.exp(-r*T) * norm.cdf(d2(S,K,r,sigma,T),mu,sigma0)

def put_rho(S,K,r,sigma,T):
    return -K * T * np.exp(-r*T) * norm.cdf(-d2(S,K,r,sigma,T),mu,sigma0)


if __name__ == '__main__':
    S = 10
    K = 10
    r = 0.05
    sigma = 0.5
    T = 1
    Cprice = call_price(S, K, r, sigma, T)
    Pprice = put_price(S, K, r, sigma, T)
    Cdelta = call_delta(S, K, r, sigma, T)
    Pdelta = put_delta(S, K, r, sigma, T)
    Cgamma = gamma(S, K, r, sigma, T)
    Cvega = vega(S, K, r, sigma, T)
    Ctheta = call_theta(S, K, r, sigma, T)
    Ptheta = put_theta(S, K, r, sigma, T)
    Crho = call_rho(S, K, r, sigma, T)
    Prho = put_rho(S, K, r, sigma, T)
    S1 = 10
    K1 = 30
    r1 = 0.04
    sigma1 = 0.5
    T1 = 1
    res = (call_price(S1,K1,r1,sigma1,T1),put_price(S1,K1,r1,sigma1,T1),call_delta(S1,K1,r1,sigma1,T1),put_delta(S1,K1,r1,sigma1,T1),gamma(S1,K1,r1,sigma1,T1),vega(S1,K1,r1,sigma1,T1),call_theta(S1,K1,r1,sigma1,T1),put_theta(S1,K1,r1,sigma1,T1),call_rho(S1,K1,r1,sigma1,T1),put_rho(S1,K1,r1,sigma1,T1))
    print(res)

