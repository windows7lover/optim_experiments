function [legend,out] = plot_algo(algo,lw,legend,scaling)


if(isfield(algo,'iter'))
    x = algo.iter;
else
    x = (1:length(algo.fval))-1;
end

if(nargin == 4 && scaling)
    x = x/max(x);
    x = x*algo.time;
end
legend = {legend{:},algo.name};
out = semilogy(x,algo.fval,'color',algo.color,'linestyle',algo.linestyle,'linewidth',lw);
