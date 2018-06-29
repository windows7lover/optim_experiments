% 

clear all
close all
clc

norm_x0_xstar = 1e-4;
M = 0.1;
L = 1;
mu = 0.1;

%lambda_vec = logspace(-9,0,10);
% k_vec = 3:2:19;
k_vec = [(1:1:9) , (10:2:20) ];

kappa = L/mu;
sigma = 1-mu/L;

nite = length(k_vec);

valbound = zeros(nite,1);

npoints_cheb = 1e2;
npointstest_cheb = 1e4;

for i = 1:nite
    display(i);
    k = k_vec(i);
    
    norm_R = geom(sigma,k)*norm_x0_xstar;
    gammavec = zeros(i,1);
    gammavec(1) = (M/(2*L)) * (sigma^2) * norm_x0_xstar^2; % "gamma2vec(0)" = 0 by definition
    for j=2:(i+1)
        gammavec(j) = sigma*gammavec(j-1) + (M/(2*L)) * (sigma^2) * norm_x0_xstar^2;
    end
    gamma = sum(gammavec);
    
    sigmavec_var = sigma.^(0:k);
    coeff_bary = sqrt(k+1)*std(sigmavec_var);
    
    norm_X = geom(sigma,k+1)*norm_x0_xstar/10;
    
%     [norm_X sqrt(coeff_bary)*norm_x0_xstar]
    
    norm_P = 4* (norm_R*gamma + gamma^2);
    
    lambda = 10*(norm_X^2*norm_P^2)/(kappa^2);
    lambda = lambda^(1/3);
%     lambda = lambda / k;
    lambdabar = lambda/norm_x0_xstar;
    
    [val_cheby,p,~,optval] = reg_minimax(sigma,lambdabar,k,npoints_cheb,npointstest_cheb);
    val_cheby = sqrt(val_cheby^2 + lambdabar*norm(p)^2);
    
    subterm1 = kappa^2;
    subterm2 = norm_X^2*norm_P^2/(lambda^3);
    term1 = sqrt(subterm1 + subterm2) * val_cheby * norm_x0_xstar ;
    term2 =  (gamma/sqrt(k+1)) * sqrt((1+norm_R^2/lambda));
%     [subterm1 subterm2]
%     [term1 term2]
    valbound(i) = term1+term2;
end

%%
bound_nest = kappa*((1-sqrt(1/kappa)).^k_vec')*norm_x0_xstar;
bound_grad = kappa*((1-(1/kappa)).^k_vec')*norm_x0_xstar;
bound_ref = bound_nest;

figure
hold on
plot(k_vec,bound_ref./bound_nest,'o-k','linewidth',3)
plot(k_vec,bound_ref./valbound,'o-r','linewidth',3)
plot(k_vec,bound_ref./bound_grad,'o-b','linewidth',3)
legend({'Nesterov-reference','normalized-rate-extrapolation','Gradient method'},'Location','east')
xlabel('degree','fontsize',16)
ylabel('lambda','fontsize',16)
set(gca,'FontSize',16);
set(gca,'PlotBoxAspectRatio',[1 0.85 1]);
xlabel('degree')
ylabel('speedup')
