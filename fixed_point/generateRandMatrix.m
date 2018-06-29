function A = generateRandMatrix(n,type,opt)

L = opt.L;

if(L >= 1 || L<= 0)
    warning('L should sometimes be in ]0,1[')
end

A = rand(n);
A = A'*A;
    
if strcmpi(type,'Random')
    [U,S,V] = svd(A);
    A = U*S*V';
%     A = eye(n)-A/eigs(A,1); % I-A has more eigenvalues close to one
    A = L*A/norm(A);
end

if strcmpi(type,'RandomStrongConvex')
    mu = opt.mu;
    [U,S,V] = svd(A);
    A = U*S*V';
%     A = eye(n)-A/eigs(A,1); % I-A has more eigenvalues close to one
    [U,S,V] = svd(A);
    s = diag(S);
    s = s-min(s);
    s = s/max(s);
    s = (L-mu)*s+mu;
    A = U*diag(s)*V';
end

if strcmpi(type,'UniformSpectrum')
    [U,S,V] = svd(A);
    S = diag([0 ; L ; L*rand(size(S,1)-2,1)]);
    A = U*S*V';
%     A = A*(L/norm(A));
end

if strcmpi(type,'UniformLogSpectrum')
    [U,S,V] = svd(A);
    mu = opt.mu;
%     s = (1:n).*(1+0.1*(0.5-rand(1,n)));
%     s = exp(opt.alpha*s);
%     s = s/max(s);
%     s = s-min(s);
%     s = mu+(L-mu)*s;
    s = logspace(log10(mu),log10(L),n);
    S = diag(s);
    A = U*S*V';
%     A = A*(L/norm(A));
end

if strcmpi(type,'StrongConvexUniformSpectrum')
    mu = opt.mu;
    [U,S,V] = svd(A);
    S = diag([mu ; L ; mu+(L-mu)*rand(size(S,1)-2,1)]);
    A = U*S*V';
%     A = A*(L/norm(A));
end


if strcmpi(type,'StrongConvexDecaySpectrum')
    mu = opt.mu;
    alpha = opt.alpha;
    [U,S,V] = svd(A);
    S = diag((1./linspace((1/mu)^(1/alpha),(1/L)^(1/alpha),n)).^alpha);
    A = U*S*V';
%     A = A*(L/norm(A));
end

if strcmpi(type,'StrongConvexClusteredSpectrum')
    nclustmu = floor(n/2);
    nclustL = floor(n/2);
    nclustMid = n-nclustmu-nclustL;
    
    mu = opt.mu;
    
    intervalSize = L-mu;
    clustSizeRatio = 0.05;
    clustSize = clustSizeRatio*intervalSize;
    
    clustmu = mu+[rand(nclustmu-1,1)*clustSize ; 0];
    clustL = L-[rand(nclustL-1,1)*clustSize ; 0];
    clustMid = ((L+mu)/2)+(1-rand(nclustMid,1))*clustSize;
    
    [U,S,V] = svd(A);
    S = diag([clustmu; clustL; clustMid]);
    A = U*S*V';
end


if strcmpi(type,'ExpOneSpectrum')
    [U,S,V] = svd(A);
    S = diag(1-([rand(size(S,1)-2,1);0;1]).^15);
    A = U*S*V';
    A = A*L; % eigs(A,1) = 1
end

if strcmpi(type,'ExpZeroSpectrum')
    [U,S,V] = svd(A);
    S = diag(([rand(size(S,1)-2,1);0;1]).^15);
    A = U*S*V';
    A = A*L; % eigs(A,1) = 1
end

if strcmpi(type,'OnePositiveEigenVal')
    [U,S,V] = svd(A);
    S = diag([zeros(size(S,1)-1,1);L]);
    A = U*S*V';
end

if strcmpi(type,'OneZeroEigenVal')
    [U,S,V] = svd(A);
    S = diag([L*ones(size(S,1)-1,1);0]);
    A = U*S*V';
end

if strcmpi(type,'SymmetricEigenval')
    [U,S,V] = svd(A);
    S = diag([0*ones(size(S,1)-2,1);L;-L]);
    A = U*S*V';
end