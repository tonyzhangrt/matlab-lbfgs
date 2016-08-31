classdef optimLineSearch < handle
    % optimLineSearch Summary of this class goes here
    %   optimization with line search
    
    properties
        % basic properties
        iter;
        feval;
        method;
        lsmethod;
        maxIter;
        Display;
        TolFun;
        TolX;
        flag;
        
        %
        curX;
        func;
        curfval;
        curgrad;
        curDir;
        curstepsize;
        inv2Precond;
    end
    
    methods
        function [this] = optimLineSearch() % default optimization method
            this.method = lbfgs(); 
%             this.lsmethod = simpleLineSearch();
            this.lsmethod = mtLineSearch();
        end
        
        function [] = optionInit(this,options)
            if (isfield(options, 'maxIter'))
                this.maxIter = options.maxIter;
            else
                this.maxIter = 1000;
            end
            if (isfield(options, 'TolFun'))
                this.TolFun = options.TolFun;
            else
                this.TolFun = 1e-6;
            end
            if (isfield(options, 'TolX'))
                this.TolX = options.TolX;
            else
                this.TolX = 1e-6;
            end
            if (isfield(options, 'Display'))
                switch options.Display
                    case 'iter-detailed'
                        this.Display = 'iter-detailed';
                    case 'off'
                        this.Display = 'off';
                    otherwise
                        this.Display = 'iter-detailed';
                end
            end      
        
        end
        
        function [X, fval] = runOptimLineSearch(this, func, X0, options,inv2Precond,RESCALE)
            this.func = func;
            this.curX = X0;
            this.optionInit(options);
            this.iter = 0;
            this.feval = 0;
            
            if (nargin == 5)
                this.method.Init(size(this.curX,1), inv2Precond);
            elseif (nargin == 6)
                this.method.Init(size(this.curX,1), inv2Precond, RESCALE);
            else
                this.method.Init(size(this.curX,1));
            end
            
            this.flag = 1;
            [this.curfval, this.curgrad] = func(this.curX);
            this.feval = this.feval + 1;
            this.optimdisplay();
            while (this.iter<=this.maxIter && this.flag~=0 )
                this.iter = this.iter + 1;
                this.method.newIteration(this.curgrad);                
                this.curDir = this.method.computeNewDir();
%                 linefunc = @(stepsize) func(this.curX+stepsize*this.curDir);
                [this.curstepsize, tmpX, tmpfval, tmpgrad ]= this.lsmethod.doLineSearch(this.func, this.curX, this.curfval, this.curgrad, this.curDir);
                this.feval = this.feval + this.lsmethod.iter;
                if (this.curstepsize == 0)
                    this.flag = 0;
                else
                    oldX = this.curX;
                    oldfval = this.curfval;
%                     oldgrad = this.curgrad;
                    this.curfval = tmpfval;
                    this.curgrad = tmpgrad;
                    this.curX = tmpX;                  
                    
%                     tmpX = this.curX + this.curstepsize * this.curDir;
                    
                    this.method.finishIteration(this.curstepsize, this.curgrad);
                    
                    % display and check end conditions
                    this.optimdisplay();
                    this.checkTolFun(oldfval);
                    this.checkTolX(oldX);
                    
                
                end    
                
            end
            X = this.curX;
            fval = this.curfval;            
        end
        
        function [] = optimdisplay(this)
            switch this.Display
                case 'iter-detailed'
                    if (this.iter == 0)
                        fprintf('%s\n','Iteration        feval        f(x)          step          optimality');
                        fprintf('%6d%14d%16g%16g%16g\n',this.iter, this.feval, this.curfval, 0, max(this.curgrad));
                    else
                        fprintf('%6d%14d%16g%16g%16g\n',this.iter, this.feval, this.curfval, this.curstepsize, max(this.curgrad));
                    end
                case 'off'    
            end
        end
        
        function [] = checkTolFun(this, oldfval)
            fdiff = abs(oldfval-this.curfval);
            if ( fdiff < this.TolFun )
               this.flag = 0;
               switch this.Display
                   case 'off'
                  
                   otherwise
                       fprintf('TolFun reached: %g.\n', fdiff);
               end
            end                         
        end
        
        function [] = checkTolX(this, oldX)                       
            Xdiff = norm(oldX-this.curX);
            if ( Xdiff < this.TolX)
                this.flag = 0;
                switch this.Display
                    case 'off'
                        
                    otherwise
                        fprintf('TolX reached: %g.\n', Xdiff);
                end
            end
        end
    end
    
end

