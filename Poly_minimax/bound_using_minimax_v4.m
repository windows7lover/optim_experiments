% minimax bound

clear all
close all
clc

normX0 = 1e-4;

mu  = 1e1;
L   = 1e2;
M   = 1e-1;

r = sqrt((L-mu)/(L+mu));
kappa = L/mu;

% npoints1 = 5000;
npoints2 = 1e4;
kvec = [2:5:30,30];
% kvec = kvec(2:end);
optValVec = zeros(size(kvec));
lambdaVec = zeros(size(kvec));
minLambdaVec = zeros(size(kvec));
boundVec = zeros(size(kvec));

npoints1 = 401;

for i=1:length(kvec)
    k = kvec(i)

    normU = kappa*(1-(1-1/kappa)^k)*normX0;
    normX = ((1-r^k)/(1-r))*normX0;
%     normXXtilde = (1+kappa)^2*(M/(2*L))*(0.5-(1-1/kappa)^k+0.5*((1-1/kappa)/(1+1/kappa))^k)*normX0^2;
    normXXtilde = compute_normXXtilde(r,L,mu,M,normX0,k);
    normE = 2*normXXtilde;
    normPert = normE*(2*normU+normE);
    normUTilde = sqrt(normU^2 + normE^2);
    
    beta = 1;
    lambda = beta*normPert;
    
    normC = sqrt((lambda+normU^2)/(k*lambda));
    normCtilde = sqrt((lambda+normU^2)/(k*lambda));
    
    [p,maxval,normp,optValVec(i)] = reg_minimax_U(1-1/kappa,lambda/normX0^2,k-1,npoints1,npoints2);
    
    lambdaVec(i) = lambda;
    
    boundVec(i) =   (sqrt(optValVec(i))*normX0)*sqrt(...
                        kappa^2 + (1/lambda)*...
                        (1+normPert/lambda)^2 * ...
                        (normXXtilde + kappa*normPert/(2*sqrt(lambda)))^2 ...
                    );
                
    [
    sqrt(optValVec(i))*normX0 ...
    kappa^2 ...
    ((1+normPert/lambda)*normXXtilde/sqrt(lambda))^2 ...
    ((1+normPert/lambda)*(kappa*normPert/(2*sqrt(lambda)))/sqrt(lambda))^2 ...
    ]
                
                
    minLambdaVec(i) = normPert;
    
end
%%


colors = [[224,69,71]; ...
[250,126,63]; ...
[115,123,13]; ...
[69,180,235]; ...
[235,111,217]];
colors = colors/255;


racc = sqrt(1-sqrt(1/kappa));
figure
rate_grad = normX0*(r.^kvec);
plot(kvec,rate_grad./rate_grad,'-x','Color',colors(1,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
hold on
plot(kvec,rate_grad./(normX0*kappa*racc.^kvec),'-*','Color',colors(2,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
plot(kvec,rate_grad./(boundVec),'-s','Color',colors(3,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
xlabel('xlabel','fontsize',16)
ylabel('ylabel','fontsize',16)
set(gca,'FontSize',16);
set(gca,'PlotBoxAspectRatio',[1 0.85 1]);
legend({'Grad','Nest','RMPE'},'location','NW')


figure
semilogy(kvec,minLambdaVec/normX0^2,'-s','Color',colors(3,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
set(gca,'FontSize',16);
xlabel('xlabel','fontsize',16)
ylabel('ylabel','fontsize',16)
set(gca,'PlotBoxAspectRatio',[1 0.85 1]);
% legend({})


