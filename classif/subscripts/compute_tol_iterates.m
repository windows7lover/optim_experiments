function [epoch,tolF,tolX,time,accuracy] = compute_tol_iterates(algoinfo,problem,nIterTol)

if(isnan(nIterTol))
    nIterTol = size(algoinfo.w,2);
end
iterates = algoinfo.w;
epochvec = round(algoinfo.grad_calc_count/problem.samples());
timevec = algoinfo.time;

nIte = size(iterates,2);

temp = linspace(1,nIte,nIterTol+1);
itertol = floor(temp);

if(any(temp~=itertol))
    warning('nIteMax/nIterTol not an integer, rounding...')
    display('nIteMax/nIterTol not an integer, rounding...')
end

sizevec = length(itertol);
epoch = zeros(sizevec,1);
tolF = zeros(sizevec,1);
tolX = zeros(sizevec,1);
time = zeros(sizevec,1);
accuracy = zeros(sizevec,1);

for i=1:length(itertol)
    idx = itertol(i);
    
    epoch(i) = epochvec(idx);
    time(i) = timevec(idx);
    
    point = iterates(:,idx) ;
    tolF(i) = problem.cost(point)-problem.fstar;
    tolX(i) = norm(point-problem.xstar);
    if(isfield(problem,'accuracy'))
        accuracy(i) = problem.accuracy(problem.prediction(point));
    end
end