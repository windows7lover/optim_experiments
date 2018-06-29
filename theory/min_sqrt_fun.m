clear all
close all
clc

a = 1;
b = 1;
c = 10;
x = linspace(0,a/b,100);

f = @(y) sqrt(a-b*y) +c*y;

figure
plot(x,f(x))