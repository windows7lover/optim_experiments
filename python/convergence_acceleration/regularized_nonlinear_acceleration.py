#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 17:45:57 2017

@author: dscieur
"""

import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import time


def rna(X,reg=0):
    
    # Regularized Nonlinear Acceleration
    # Take a matrix X of iterates, where X[:,i] is the ith iteration of the
    # fixed-point operation
    #   x_i = g(x_{i-1})
    #
    # reg is the regularization parameter used for solving the system
    #   (R'R + reg I)z = 1
    # where R is the matrix of residuals, i.e. R[:,i] = x_{i+1}-x_{i}
    
    
    # Recovers parameters, ensure X is a matrix
    (d,k) = np.shape(X);
    k = k-1;
    X = np.asmatrix(X);
    
    # Compute the matrix of residuals
    R = np.diff(X);
    
    # "Square" the matrix, and normalize it
    RR = np.dot(np.transpose(R),R);
    normRR = LA.norm(RR,2);
    RR = RR/normRR;
    
    # Solve (R'R + lambda I)z = 1
    reg_I = reg*np.eye(k);
    ones_k = np.ones(k);
    z = np.linalg.solve(RR+reg_I, ones_k);
    
    if(z.size != k):
        z = np.linalg.lstsq(RR+reg_I, ones_k, -1);
        z = z[0]
    
    # Recover weights c, where sum(c) = 1
    if( np.abs(np.sum(z)) < 1e-10):
        z = np.ones(k)
    else:
        c = np.asmatrix(z/np.sum(z)).T;
    
    # Compute the extrapolation / weigthed mean  "sum_i c_i x_i", and return
    extr = np.dot(X[:,0:k],c[:,0]);
    return np.array(extr),c





def rna_precomputed(X,RR,reg=0):
    
    # Regularized Nonlinear Acceleration, with RR precomputed
    # Same than rna, but RR is computed only once
    
    
    # Recovers parameters
    (d,k) = X.shape;
    k = k-1;
    
    # RR is already computed, we do not need this step anymore
    
    # # Compute the matrix of residuals
    # R = np.diff(X);
    
    # # "Square" the matrix, and normalize it
    # RR = np.dot(np.transpose(R),R);
    # normRR = LA.norm(RR,2);
    # RR = RR/normRR;
    
    # Solve (R'R + lambda I)z = 1
    reg_I = reg*np.eye(k);
    ones_k = np.ones(k);
    
    # In case of singular matrix, we solve using least squares instead
    try:
        z = np.linalg.solve(RR+reg_I, ones_k);
    except LA.linalg.LinAlgError:
        z = np.linalg.lstsq(RR+reg_I, ones_k, -1);
        z = z[0]
    
    # Recover weights c, where sum(c) = 1
    if( np.abs(np.sum(z)) < 1e-10):
        z = np.ones(k)
    
    c = np.asmatrix(z/np.sum(z)).T;
    
    # Compute the extrapolation / weigthed mean  "sum_i c_i x_i", and return
    extr = np.dot(X[:,0:k],c[:,0]);
    return np.array(extr),c

def grid_search(logmin,logmax,k,fun_obj):
    
    # Perform a logarithmic grid search between [10^logmin,10^logmax] using
    # k points. Return the best value found in the grid, i.e. the minimum value
    # of fun_obj. In other words, it returns the value satisfying
    #   argmin_{val in logspace(logmin,logmax)} fun_obj(val)
    
    #always test 0
    lambda_grid = np.append( [0] , np.logspace(logmin,logmax,k));
    k = k+1
    
    # pre-allocation    
    vec = np.zeros(k);
    
    # test all values in the grid
    for idx in range(len(lambda_grid)):
        vec[idx] = fun_obj(lambda_grid[idx]);
        
    # get the best value in the grid and return it
    idx = np.argmin(vec);
    return lambda_grid[idx];

def approx_line_search(obj_fun,x0,step):
    
    # Perform an approximate line search; i.e. find a good value t for the
    # problem
    #   min_t f(x0+t*step)
    
    # Define the anonymous function f(t), returning obj_fun(x0+t*step) for
    # fixed x0 and step
    d = len(x0);
    x0 = np.reshape(x0,(d,1))
    step = np.reshape(step,(d,1))
    objval_step = lambda t: obj_fun(x0+t*step);
    
    
    
    # We multiply the value of t at each iteration, then stop when the function
    # value increases.
    t = 1;
    oldval = objval_step(t);
    t = t*2;
    newval = objval_step(t);
    while(newval < oldval):
        t = 2*t
        oldval = newval;
        newval = objval_step(t)
    # Best value of t found
    t = t/2
    
    return x0+t*step, t;

def adaptive_rna(X,obj_fun,logrange=[-15,1]):
    
    # Adaptive regularized nonlinear acceleration
    #
    # Perform an adaptive search for lambda and the stepsize for the rna
    # algorithm. It automatically finds a good value of lambda wy using a grid
    # search, then find a good step size by an approximate line-search.
    #
    # X is the matrix of iterates, and obj_fun is the value of the objective
    # function that we want to minimize.
    
    d,k = X.shape;
    k = k-1;
    
    # Precompute the residual matrix
    R = np.diff(X);
    RR = np.dot(np.transpose(R),R);
    normRR = LA.norm(RR,2);
    RR = RR/normRR;
    
    #anonymous function, return the objective value of x_extrapolated(lambda)
    obj_extr = lambda lambda_val: obj_fun(rna_precomputed(X,RR,lambda_val)[0]);
    
    # grid search
    lambda_opt = grid_search(logrange[0],logrange[1],k,obj_extr);
    
    # Retrieve the best extrapolated point found in the grisd search
    x_extr,c = rna_precomputed(X,RR,lambda_opt)
    
    # Perform an approximate line-search
    step = x_extr[:,0]-X[:,0];
    (x_extr,t) = approx_line_search(obj_fun,X[:,0],step)
    
    # Return the extraplated point, the coefficients of the weigthed mean and
    # the step size.
    x_extr = np.reshape(x_extr,(d,1))
    return x_extr,c,t
   
def test():
    # Test of adaptive_rna on a quadratic function
    
    # function Parameters
    d = 100;
    L = 10.0;
    mu = 0.001;
    
    # Create a quadratic minimization problem: min 0.5*||A*(x-xstar)||
    # We build A so that mu*I <= A'*A <= L*I
    A = np.random.rand(d,d);
    AA = np.dot(np.transpose(A),A);
    AA = AA/LA.norm(AA);
    AA = (L-mu)*AA+mu*np.eye(d); # L control the smoothness parameter, 
                                 # and mu the strong convexity
    
    # Random initialization and solution 
    x0 = np.random.rand(d,1);
    xstar = np.random.rand(d,1);
    
    # Define the gradient of the function
    grad = lambda x: np.dot(AA,x-xstar);
    
    # We want to minimize the error, or the norm of the gradient
    err = lambda x: LA.norm(grad(x))
        
    
    # Total number of calls : k*nLoop 
    nLoop = 1000;
    k = 5;
    x_hist = np.zeros((d,k+1));
    
    # Record the error of accelerated gradient
    err_vec_rna = np.zeros((nLoop+1,1));
    err_vec_rna[0,0] = err(x0);
    
    # Record the error of gradient method
    err_vec_grad = np.zeros((nLoop+1,1));
    err_vec_grad[0,0] = err(x0);
    
    
    # accelerated gradient method
    start = time.time()
    x_new = x0;
    for i in range(0,nLoop):
        x_hist[:,:1] = x_new; # Restart at extrapolated point
        x = x_new;
        # Perform k gradient steps
        for j in range(1,k+1):
            x = x-(1/L)*grad(x); 
            x_hist[:,j] = x[:,0];
        # Adaptive regularization, using only k gradient steps
        (x_new,c,t) = adaptive_rna(x_hist,err)
        # Store the error of the extrapolation
        err_vec_rna[i+1,0] = err(x_new);
    stop = time.time()
    time_acc = stop-start;
    
    
    # classical gradient method
    start = time.time()
    x_plus = x0;
    for i in range(0,nLoop):
        for j in range(1,k+1):
            x_plus = x_plus-(1/L)*grad(x_plus);
        err_vec_grad[i+1,0] = err(x_plus);
    stop = time.time()
    time_grad = stop-start;
    
    
    iterates = range(0,nLoop+1);
    grad_plot, = plt.semilogy(iterates, err_vec_grad, 'b', label='Gradient descend')
    acc_grad_plot, = plt.semilogy(iterates, err_vec_rna, 'r', label='Reg. Nonlinear Acceleration on gradient')
    plt.legend()
    plt.show()
    print('Time taken (gradient / acceleration of gradient)')
    print((time_grad,time_acc))
    
test()