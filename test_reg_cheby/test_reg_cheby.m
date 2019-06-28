
kappa = 0.01;
alpha = 0.0;
tx = @(x) (2*x/(1-kappa)) - 1;

points = 0:kappa*0.1:(1-kappa);
points_tx = tx(points);

N = 1;

figure
hold on
for i=1:N
    chebpoly = reg_cheb(i,alpha);
    plot(points,polyval(chebpoly,points_tx)/polyval(chebpoly,tx(1)));
    ones_poly = ones(size(chebpoly));
    plot(points,polyval(ones_poly,points)/polyval(ones_poly,1),'-.');
end