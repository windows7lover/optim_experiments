% regularized chebyshev


condition_number = 1e-3;
sigma = 1-condition_number;

loadfile = true;
filename = ['cheby_reg_sigma_1e' , num2str(num2str(round(log10(condition_number))))];

ndegree = 12;
nlambda = 15;
npoints = 1e2;
npointstest = 1e4;

degree = 3:(3+ndegree-1);
lambda = logspace(-5,0,nlambda);


val_regularized_cheby = zeros(ndegree,nlambda);
val_cheby = zeros(ndegree,1); % lambda = 0

if(~loadfile)
    for i=1:ndegree
        for j=1:nlambda
            display([i,j])
            val_regularized_cheby(i,j) = reg_minimax(sigma,lambda(j),degree(i),npoints,npointstest);
        end
        beta = (1-sqrt(1-sigma))/(1+sqrt(1+sigma));
        val_cheby(i) = 2*(beta^degree(i))/(1+beta^(2*degree(i)));
    end
    save([filename , 'temp']);
else
    load(filename)
end

%% plot
val_regularized_cheby = min(val_regularized_cheby,1); % not very accurate :'(
% Compute ratio between unregularized version and regularized version

val_regularized_cheby_ratio = val_regularized_cheby;
for i=1:ndegree
    beta = (1-sqrt(1-sigma))/(1+sqrt(1-sigma));
    val_cheby(i) = 2*(beta^degree(i))/(1+beta^(2*degree(i)));
    for j=1:nlambda
        val_regularized_cheby_ratio(i,j) = log(val_cheby(i,1))/log(val_regularized_cheby(i,j)) ;
    end
end

[S,L] = meshgrid(degree,lambda);
fig = figure;
contourf(S,L,val_regularized_cheby_ratio')
set(gca,'yscale','log');
%set(gca,'xscale','log');
colorbar()
xlabel('degree','fontsize',16)
ylabel('lambda','fontsize',16)
set(gca,'FontSize',16);
set(gca,'PlotBoxAspectRatio',[1 0.85 1]);
%title(['Regularized Chebyshev for sigma = 1-1e' num2str(round(log10(condition_number)))])
