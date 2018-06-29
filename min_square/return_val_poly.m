function [poly_val,var_val,problem_val] = return_val_poly(c,points,alpha)

degree = length(c)-1;
npt = length(points);


poly_val_points = zeros(npt,1);
poly_val_points(:,1) = polyval(flipud(c),points).^2;


var_val_points = zeros(npt,degree+1);
for i=1:degree+1
    vec_idx = (degree+2-i):degree+1;
    poly_i = flipud(c(vec_idx));
    var_val_points(:,i) = polyval(poly_i,points).^2;
end
max_val_var = sum(var_val_points,2);

[problem_val, idx] = max( poly_val_points + alpha * max_val_var);

poly_val = poly_val_points(idx);
var_val = max_val_var(idx);