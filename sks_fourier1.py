#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 11:27:00 2018

@author: salonisaxena
"""
import numpy as np
import matplotlib.pyplot as plt
import time

def der1(uk_0):
    '''
    Parameters
    ----------
    uk_0 : FFT of current state.

    Returns
    -------
    der_k : FFT of the spatial derivative of u(x).
    '''
    der_k = 1j*coeff*uk_0
    der_k[N//2] = 0
    return der_k


def der2(uk_0):
    '''
    Parameters
    ----------
    uk_0 : FFT of current state.

    Returns
    -------
    der2_k : FFT of Second spatial derivative.

    '''
    der2_k = -coeff**2 * uk_0
    return der2_k

def der4(uk_0):
    '''
    Parameters
    ----------
    uk_0 : FFT of current state.

    Returns
    -------
    der4_k : FFT of fourth spatial derivative.

    '''
    der4_k = coeff**4 * uk_0
    return der4_k

def nonlinear(uk_0): 
    '''
    Parameters
    ----------
    uk_0 : FFT of current state. Converts it to position space, generates nonlinear term

    Returns
    -------
    FFT of nonlinear term

    '''
    der_k = der1(uk_0)
    
    temp = np.fft.ifft(der_k)      ###
    der_x = np.real(temp)
    return np.fft.fft(der_x * der_x)      ###
    
    
def real(N):
    '''
    Parameters
    ----------
    N : Number of lattice points.

    Returns
    -------
    a : array of random gaussian numbers.This is the real part of the FFT of the noise

    '''
    a_positive = np.random.normal(0, np.sqrt(0.5), (N-1)//2)
    a_negative = a_positive[::-1]
    a0 = np.random.normal(0,1)
    a_positive = np.insert(a_positive, 0, a0)
    a_mid = np.random.normal(0,1)
    a_negative = np.insert(a_negative, 0, a_mid)
    a = np.concatenate((a_positive, a_negative))
    return a

def imag(N):
    '''Parameters
    ----------
    N : Number of lattice points.

    Returns
    -------
    b : array of random gaussian numbers. Imaginary part of FFT of noise.

    '''
    b_positive = np.random.normal(0, np.sqrt(0.5), (N-1)//2)
    b_negative = -b_positive[::-1]
    b_positive = np.insert(b_positive, 0, 0)
    b_mid = 0 #np.random.normal(0,0.5)
    b_negative = np.insert(b_negative, 0, b_mid)
    b = np.concatenate((b_positive, b_negative))
    return b


def implicit1(uk_0): 
    '''
    Parameters
    ----------
    uk_0 : Current state.

    Returns
    -------
    FFT of solution at next time step

    '''
    noise_fourier = np.sqrt(2*eps*N*dt/h)*(real(N)+1j*imag(N))  ###white noise
    uk_1 = np.zeros(N, dtype=complex) ###
    nl_fourier = nonlinear(uk_0)
    uk_1 = (uk_0 + 1*dt*nl_fourier + noise_fourier)/(I - dt*(-alpha*I + 1.*coeff**2 - 1.0*coeff**4))    #first order euler
    return np.real(np.fft.ifft(uk_1)), noise_fourier ###
 
def growth_rate(Nc):
    '''
    Parameters
    ----------
    Nc : Index corresponding to some wavenumber.

    Returns
    -------
   Growth rate of that wavenumber.

    '''
    return -alpha + (2*np.pi*Nc/L)**2 - (2*np.pi*Nc/L)**4



#spatial grid, time step, etc.
start_time = time.time()
N = 4000 

half_N = int(N/2)
Nc = 225     #index corresponding to critical wavenumber
Nc_min = 200  # index of minimum allowed wavenumber
Nc_max = 240  # index of maximum allowed wavenumber
dt = 0.3  #time step

h = 0.5  #lattice spacing
L = N*h  # length of spatial domain 
alpha = 0.20  #control parameter that determines the wavenumber of the pattern                     
x = np.linspace(h,L,N) #lattice points

q = 2*np.pi*Nc/L  # wavenumber of initial state
qc = 2*np.pi*225/L

#beginning the integration

hist = np.zeros(half_N+1)  #histogram of wavenumbers i.e. the number of times a state with a given wavenumber occurs
freq = np.arange(0,half_N+1) #Fourier space frequencies
I = np.ones(N)   ###
k = N*h*np.fft.fftfreq(N,h)      ###
coeff = 2*np.pi*k/L

arr = []  #list to store histogram
arr1 = [] #list to store power spectrum

#initial conditions

u0 = 0.01*np.sin(q*x)        
uk_0 = np.fft.fft(u0)     

t = 0.0 # start time

n=0   # time step counter

eps = 0.004   # noise strength

T = 100000     # total time of integration

# plot the initial state                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
plt.figure()
plt.xlabel(r'$x$',fontsize=14)
plt.ylabel(r'$u(x,t)$',fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.plot(x,u0)
plt.show()

file1 = open('hist.txt', 'a')   #output files
file2 = open('PowerSpecrum.txt', 'a')

while t < T:
    arr.append(np.abs((growth_rate(Nc)*(uk_0[Nc])+(nonlinear(uk_0)[Nc]))))   # growth rate of critical wavenumber as a function of time
    u = implicit1(uk_0)[0] #new state in position space
    uk_0 = np.fft.fft(u)  # Fourier transform of new state
    
    n += 1   
    
    #calculating histogram  
    index = np.argmax(np.abs(uk_0[1:N//2])) + 1   #index of wavenumber with maximum amplitude
    if np.abs(uk_0[index])/np.abs(uk_0[index-1])>=2 and np.abs(uk_0[index])/np.abs(uk_0[index+1]) >= 2:   #if amplitude of this wavenumber is at least twice as great as all the other wavenumbers, increment appropriate histogram count by 1
            hist[index] += 1
    
    if n%6666 == 0:     # print output at regular time intervals to see if everything is working ok
            print(t, index, hist[Nc_min:Nc_max]/n) #
                    
    if n>225000 and n%100 == 0:  #save output to a file every 100 time steps
            arr1.append(t)
            arr.append(index)
            
            file1.write(str(t)+ str(hist[Nc_min:Nc_max+1]/n)+ '\n')
            file2.write(str(t) + str(np.abs(uk_0[Nc_min:Nc_max+1]))+'\n')
            
    
    t += dt

file1.close()  
file2.close()
    
plt.figure()
plt.xlabel(r'$q$',fontsize=14)
plt.ylabel(r'Fraction of time in state',fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.plot(coeff[190:240], hist[190:240]/n, 'o-')
plt.show()
#%%    
# plot power spectrum of current state
plt.figure() 
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)  
plt.xlabel('$q$', fontsize=14)
plt.ylabel(r'$|\widetilde{u}_{lq}|$', fontsize=14)
plt.plot(coeff[0:half_N], np.abs(uk_0)[0:half_N])
plt.show()



print("--- %s seconds ---" % (time.time() - start_time))

#plot the current state in position space
plt.figure()
plt.xlim(h,50)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('$x$', fontsize=14)
plt.ylabel('$u_{\\bar q,det}(x)$', fontsize=14)
plt.plot(x,u)
plt.show()

#selected wavenumber for different alpha
a1 = np.array([0.17, 0.20, 0.22, 0.24])  # alpha values used
qc = np.array([0.707, 0.707, 0.707, 0.707]) #critical wavenumber
q0 = np.array([0.6377, 0.6566, 0.6754, 0.6974])  #selected wavenumber
plt.figure()
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel(r'$\alpha$', fontsize=14)
plt.ylabel('Fourier mode with largest long-time amplitude', fontsize=14)
plt.plot(a1, q0, '-o', label='Stochastic')
plt.text(0.17,0.704,r'$q_c$',fontsize=14)
plt.plot(a1, qc, '--', label='Critical mode')
plt.legend(fontsize=14, loc=4)
plt.show()
