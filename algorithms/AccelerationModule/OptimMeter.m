classdef OptimMeter
    properties
        fval;
        gradval;
        normx;
        finfo;
        nIter;
        iter;
        innerCounter;
        
        color;
        linewidth;
        linestyle;
        marker;
    end
    methods
        function obj = OptimMeter(nIter,finfo,xinit)
            if(nIter > 0)
                obj.nIter = nIter;
            else
                error('nIter should be > 0')
            end
            if(nargin < 3)
                x0 = finfo.x0;
            else
                x0 = xinit;
            end
            
            obj.finfo = finfo;
            obj.fval = zeros(nIter+1,1);
            obj.gradval = zeros(nIter+1,1);
            obj.normx = zeros(nIter+1,1);
            obj.iter = 0:nIter;
            obj.innerCounter = 0;
            
            obj = obj.plotParam();
            obj = obj.store(x0);
        end
        
        function obj = plotParam(obj,color, linewidth, linestyle, marker)
            if(nargin < 2)
                color = 'r';
            end
            if(nargin < 3)
                linewidth = 1;
            end
            if(nargin < 4)
                linestyle = '-';
            end
            if(nargin < 5)
                marker = 'none';
            end
            
            obj.color = color;
            obj.linewidth = linewidth;
            obj.linestyle = linestyle;
            obj.marker = marker;
            
            
        end
        
        function obj = store(obj,x)
            if(obj.innerCounter > obj.nIter)
                error('Too much iteration')
            end
            
            obj.innerCounter = obj.innerCounter+1;
            obj.fval(obj.innerCounter) = obj.finfo.f(x)-obj.finfo.f(obj.finfo.xstar);
            obj.gradval(obj.innerCounter) = norm(obj.finfo.fp(x));
            obj.normx(obj.innerCounter) = norm(x-obj.finfo.xstar);
        end
        
        
        function plot(obj,dataName,typePlot)
            if(nargin < 2)
                dataName = 'valF';
            end
            if(nargin < 3)
                typePlot = 'semilogy';
            end
            
            if strcmpi(dataName, 'valf')
                data = obj.fval;
            elseif strcmpi(dataName, 'normx')
                data = obj.gradval;
            elseif strcmpi(dataName, 'normgrad')
                data = obj.normx;
            else
                error('Data name not recognized')
            end
            
            if strcmpi(typePlot,'plot')
                plotfun = @plot;
            end
            if strcmpi(typePlot,'semilogx')
                plotfun = @semilogx;
            end
            if strcmpi(typePlot,'semilogy')
                plotfun = @semilogy;
            end
            
            iter_data = obj.iter(1:length(data));
            plotfun(iter_data,data,'color',obj.color,'linewidth',obj.linewidth,'linestyle',obj.linestyle,'marker',obj.marker);
        end
        
    end
end