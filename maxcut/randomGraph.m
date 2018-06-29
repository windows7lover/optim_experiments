function G = randomGraph(nNodes,nEdges)

n= nNodes;
E=nEdges;
adj = spalloc(n, n, E);
idx = randperm(n * n, E+n);
idx(ismember(idx, 1:n+1:n*n)) = [];
idx = idx(1:E);
adj(idx) = rand(size(idx));
A = adj + adj';
G = graph(A);