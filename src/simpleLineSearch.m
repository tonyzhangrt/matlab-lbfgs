classdef simpleLineSearch < handle
    %  simpleLineSearch Summary of this class goes here
    %  simple line search with Wolfe condtions
    
    properties        
        % Wolfe conditions
        c1wolfe;
        c2wolfe;
        tau1;
        tau2;
        maxIter; % maximum line search steps
        iter;
        % fucntion handle to evaluate fval and gradient and only depends on step size
        Xold;
        fold;
        gradold;
        dir;
        desc;
    end
    
    methods
        function this = simpleLineSearch(c1,t1,c2,t2,maxIter)
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
        
        function [stepsize, X, fval, grad] = doLineSearch(this,func, Xold, fold, gradold, dir, stepsize0)
            if (nargin == 6)
                stepsize = 1;
            elseif (nargin == 7)
                stepsize = stepsize0;
            else
                error('doLineSearch: Wrong Input !')
            end
            
            this.Xold = Xold;
            this.fold = fold;
            this.gradold = gradold;
            this.dir = dir;
            this.desc = dot(gradold, dir);
            
            if (this.desc >= 0.0)
                error('doLineSearch: wrong descent direction !')
            end

            
            this.iter = 0;
            flag = 2;
            while (this.iter<=this.maxIter && flag~=0 )
                this.iter  = this.iter +1;
                X = this.Xold + stepsize*this.dir;
                [fval, grad] = func(X);
                flag = this.checkWolfe(fval, grad);
                if (flag == -1)
                    stepsize = this.tau1 * stepsize;
                elseif (flag == 1)
                    stepsize = this.tau2 * stepsize;
                else
                    return;
                end
            end
            
%             fprintf('doLineSearch: maxIter %d reached, flag %d !\n', this.iter, flag);
%             stepsize = 0;
%             return;

            if (fval > this.fold)
                fprintf('Function cannot be decresed in current direction !\n');
                stepsize = 0;
            end
            
            
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

