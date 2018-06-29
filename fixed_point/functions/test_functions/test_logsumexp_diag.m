% test_logsumexp_diag

clear all
close all
clc

n = 1;
dir_vec = rand(n,1);
dir_vec = dir_vec/norm(dir_vec);
x = -10:0.001:10;
x=x';
npoints = length(x);

[f,fp,fpp,L,xstar] = logsumexp_diag(n);

fx = zeros(npoints,1);
fpx = zeros(npoints,1);
fppx = zeros(npoints,1);

for i=1:length(x)
    fx(i) = f(x(i)*dir_vec);
    fpx(i) = fp(x(i)*dir_vec)'*dir_vec;
    fppx(i) = dir_vec'*fpp(x(i)*dir_vec)*dir_vec;
end
figure
subplot(4,1,1)
plot(x,fx)
subplot(4,1,2)
plot(x,fpx)
hold on
plot(x(1:end-1),diff(fx)./diff(x),'r-.')
subplot(4,1,3)
plot(x,fppx)
hold on
plot(x,L*ones(npoints,1),'g')
plot(x(1:end-1),diff(fpx)./diff(x),'r-.')
subplot(4,1,4)
plot(x(1:end-1),diff(fppx)./diff(x),'r-.')