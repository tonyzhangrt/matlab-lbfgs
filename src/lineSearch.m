classdef lineSearch < handle
    %  lineSearch Summary of this class goes here
    %  line search Wolfe condtions
    
    properties        
        % Wolfe conditions
        c1wolfe;
        c2wolfe;
        tau1;
        tau2;
        maxIter; % maximum line search steps
        % fucntion handle to evaluate fval and gradient and only depends on step size
        fold;
        gradold;
        dir;
        desc;
    end
    
    methods
        function this = lineSearch(c1,t1,c2,t2,maxIter)
            this.c1wolfe = 1e-4;
            this.c2wolfe = 0.9;
            this.tau1 = 0.5;
            this.tau2 = 4;
            this.maxIter = 20;
            switch nargin
                case 0
                case 2
                    this.c1wolfe = c1;
                    this.tau1 = t1;
                case 4
                    this.c1wolfe = c1;
                    this.tau1 = t1;
                    this.c2wolfe = c2;
                    this.tau2 = t2;
                case 5
                    this.c1wolfe = c1;
                    this.tau1 = t1;
                    this.c2wolfe = c2;
                    this.tau2 = t2;
                    this.maxIter = maxIter;
                otherwise
                    error('LineSearch Initiation: Wrong Input !')
            end
            
        end
        
        function [stepsize, fval, grad] = doLineSearch(this,linefunc,fold, gradold, dir, stepsize0)
            if (nargin == 5)
                stepsize = 1;
            elseif (nargin == 6)
                stepsize = stepsize0;
            else
                error('doLineSearch: Wrong Input !')
            end
            
            
            this.fold = fold;
            this.gradold = gradold;
            this.dir = dir;
            this.desc = dot(gradold, dir);
            
            iter = 0;
            flag = 2;
            while (iter<=this.maxIter && flag~=0 )
                iter  = iter +1;
                [fval, grad] = linefunc(stepsize);
                flag = this.checkWolfe(fval, grad);
                if (flag == -1)
                    stepsize = this.tau1 * stepsize;
                elseif (flag == 1)
                    stepsize = this.tau2 * stepsize;
                else
                    return;
                end
            end
            
            fprintf('doLineSearch: maxIter %d reached, flag %d !\n', iter, flag);
            stepsize = 0;
            return;
            
            
        end
        
        
        function [flag] = checkWolfe(this, fval, grad)
            if ( fval-this.fold > this.c1wolfe * this.desc)
                flag = -1;
%                 fprintf('Stepsize too large !\n');
            elseif ( dot(grad, this.dir) < this.c2wolfe * this.desc)
                flag = 1;
%                 fprintf('Stepsize too small !\n');
            else
                flag = 0;
            end                       
        end
    end
    
end

