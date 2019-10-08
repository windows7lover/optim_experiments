% regularized chebyshev


condition_number = 1e-3;
sigma = 1-condition_number;

% loadfile = true;
% filename = ['cheby_reg_sigma_1e' , num2str(num2str(round(log10(condition_number))))];

ndegree = 10;
ntau = 10;
npoints = 1e2;
npointstest = 1e4;

degree = 3:(3+ndegree-1);
tau = [0 logspace(-0,5,ntau)];
ntau = ntau+1;


val_constrained_cheby = zeros(ndegree,ntau);
poly = cell(ndegree,ntau);
val_cheby = zeros(ndegree,1); % lambda = 0

points = linspace(0,sigma,10*npoints);

if(~loadfile)
    for i=1:ndegree
        for j=1:ntau
            display([i,j])
            [val_constrained_cheby(i,j), poly{i,j}] = const_minimax(sigma,tau(j),degree(i),npoints,npointstest);
        end
        beta = (1-sqrt(1-sigma))/(1+sqrt(1+sigma));
        val_cheby(i) = 2*(beta^degree(i))/(1+beta^(2*degree(i)));
    end
    save([filename , 'temp']);
else
%     load(filename)
end

%% plot
val_constrained_cheby = min(val_constrained_cheby,1); % not very accurate :'(
% Compute ratio between unregularized version and regularized version

val_constrained_cheby_ratio = val_constrained_cheby;
for i=1:ndegree
    beta = (1-sqrt(1-sigma))/(1+sqrt(1-sigma));
    val_cheby(i) = 2*(beta^degree(i))/(1+beta^(2*degree(i)));
    for j=1:ntau
        val_constrained_cheby_ratio(i,j) = log(val_cheby(i,1))/log(val_constrained_cheby(i,j)) ;
    end
end

[S,L] = meshgrid(degree,tau);
fig = figure;
contourf(S,L,val_constrained_cheby_ratio')
set(gca,'yscale','log');
%set(gca,'xscale','log');
colorbar()
xlabel('degree','fontsize',16)
ylabel('rho','fontsize',16)
set(gca,'FontSize',16);
set(gca,'PlotBoxAspectRatio',[1 0.85 1]);
%title(['Regularized Chebyshev for sigma = 1-1e' num2str(round(log10(condition_number)))])
hold on


%%


figure
degree2plot = 10;
tau2plot = [2 5 10];
pgrad = 0*poly{degree2plot, 1};
pgrad(1) = 1;
legendcell = {};
for i = tau2plot
    p = poly{degree2plot, i};
    plot(points, polyval(p,points))
    hold on
    taustr = round(tau(i),1-floor(log(tau(i))));
    legendcell{end+1} = ['Const. Cheb. (\tau=' num2str(taustr) ')'];
end
plot(points, polyval(pgrad,points))
legendcell{end+1} = 'Gradient polynomial';
legend(legendcell)

