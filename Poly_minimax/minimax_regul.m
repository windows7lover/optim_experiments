% 

clear all
close all
clc

%% Just some plot to illustrate the impact of regularization
lambda = 1e-3;
degree = 7;
npoints = 10e1+1;
sigma = 1-1e-1;
x = 0:0.00001:sigma;

p = reg_minimax_U(sigma,lambda,degree,npoints);
% p = reg_minimax(sigma,lambda,degree,npoints);

maxval = max(polyval(p,x));
normp = norm(p);

display(maxval)
display(normp)

optCheby.spectrum = 'pos';
c_cheby = find_cheby_poly(sigma,degree+1,optCheby);


figure
hold on
plot(x,polyval(p,x))
plot(x,polyval(c_cheby,x),'r-.')

%% Norm in function of lambda

lambdaVec = logspace(-15,5,35);

normPvec = zeros(length(lambdaVec),1);
maxvalVec = zeros(length(lambdaVec),1);

for i=1:length(lambdaVec)
    i
    [p,maxval,normp] = reg_minimax(sigma,lambdaVec(i),degree,npoints);
    normPvec(i) = normp;
    maxvalVec(i) = maxval;
end
figure
loglog(lambdaVec, normPvec)
figure
loglog(lambdaVec, maxvalVec)

%% Norm in function of k

kVec = 3:25;
lambdaVec = kVec.*lambda;

normPvec = zeros(length(kVec),1);
maxvalVec = zeros(length(kVec),1);

x = 0:0.00001:sigma;

for i=1:length(kVec)
    i
%     [p,maxval,normp] = reg_minimax(sigma,lambda,kVec(i),npoints);
    [p,maxval,normp] = reg_minimax_U(sigma,lambdaVec(i),kVec(i),npoints);
    normPvec(i) = normp;
    maxvalVec(i) = max(polyval(p,x));
    display(maxval - maxvalVec(i))
end
figure
plot(kVec, normPvec)
figure
hold on
plot(kVec, maxvalVec)
plot(kVec, maxvalVec + lambdaVec'.*normPvec,'r')

%% Create some mapping with fixed sigma

% map
