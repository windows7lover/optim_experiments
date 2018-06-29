
addpath(genpath('/home/dscieur/Projects/toolboxes/cvx_free'))
clear all
close all
clc

alpha = 0.01;
kappa = 0.1;
npoints = 100;
npointstest = 10000;

degree_vec = 3:2:21;
n_degree = length(degree_vec);

max_val_poly_vec = zeros(n_degree,1);
max_val_var_vec = zeros(n_degree,1);
max_val_problem_vec = zeros(n_degree,1);

max_val_avepoly_vec = zeros(n_degree,1);
max_val_avevar_vec = zeros(n_degree,1);
max_val_aveproblem_vec = zeros(n_degree,1);

for i = 1:n_degree
    display(i)
    degree = degree_vec(i);
    warning off
    [p, max_val_poly_vec(i), max_val_var_vec(i), max_val_problem_vec(i)] = stoch_cheby_poly(degree,alpha,kappa,npoints,npointstest);
    warning on
    
    avepoly = ones(degree+1,1)/(degree+1);
    x = linspace(0,1-kappa,npointstest)';
    [max_val_avepoly_vec(i),max_val_avevar_vec(i),max_val_aveproblem_vec(i)] = return_val_poly(avepoly,x,alpha);
    
    display('---------------------------------------------')
    display('[max_val_poly, max_val_var, max_val_problem]')
    display([ max_val_poly_vec(i), max_val_var_vec(i), max_val_problem_vec(i)])
    display([ max_val_avepoly_vec(i), max_val_avevar_vec(i), max_val_aveproblem_vec(i)])
    display('---------------------------------------------')
end


%%
x = linspace(0,1-kappa,10*npointstest);
figure
hold on
for i=2:length(p)
    plot(x,polyval(p(end-i+1:end),x))
end

figure
hold on
legend_cell = {};
plot(degree_vec,max_val_poly_vec,'r'); legend_cell = [legend_cell 'max poly'];
plot(degree_vec,max_val_var_vec,'g'); legend_cell = [legend_cell 'max var'];
plot(degree_vec,max_val_problem_vec,'c'); legend_cell = [legend_cell 'max problem'];
plot(degree_vec,max_val_avepoly_vec,'-.r'); legend_cell = [legend_cell 'max ave poly'];
plot(degree_vec,max_val_avevar_vec,'-.g'); legend_cell = [legend_cell 'max ave var'];
plot(degree_vec,max_val_aveproblem_vec,'-.c'); legend_cell = [legend_cell 'max ave problem'];
legend(legend_cell);
