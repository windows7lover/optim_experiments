function [X,y,xstar] = gen_problem_svm(npoints,nfeatures)

% points1 = rand(npoints/2,nfeatures);
% points2 = rand(npoints/2,nfeatures);
% 
% X = [points1 ; points2 ];
% y = [ones(npoints/2,1); -ones(npoints/2,1)];

load classif/Data/Madelon/madelon_data.mat
X = [100*A/norm(A,'fro'), ones(size(A,1),1)];

[Q, p, A, b] = transform_svm_dual(1,X,y);

options = optimoptions('quadprog','MaxIterations',max(3*length(p),500),'OptimalityTolerance',1e-10);
 [xstar] = quadprog(Q,p,A,b,[],[],[],[],[],options);