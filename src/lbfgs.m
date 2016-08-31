classdef lbfgs < handle
    %LBFGS Summary of this class goes here
    %   L-BFGS update of optimization
    %   Update of direction and required storage are implemented here
    
    properties
        n; % problem size
        m; % stored iterations
        iter; % number of iterations
        lbfgsDir; % search direction
        buffS; % difference of position
        buffY; % difference of gradients
        buffRho; % 1/(s.'*y)
%         fold; % buffered value of f
        buffGrad; % buffered Gradient
        Hk0;
        
        % precondioner information
        isPrecond; % if precond
        inv2Precond; % C^(-1) * C^(-T), can be matrix or function handle 
        
        % RESCALING
        RESCALE;
        
    end
    
    methods
        function this = lbfgs(m)
            this.m = 10;
            if (nargin == 1)
                this.m = m;
            end
        end        
        
        function [] = Init( this, np, P, RESCALE)
            this.iter = 0;
            this.n = np;
            this.Hk0 = 1;
            this.RESCALE = 1;
            if (nargin == 2)
                this.isPrecond = 0;
                this.buffS = zeros(this.n, this.m);
                this.buffY = zeros(this.n, this.m);
                this.buffGrad = zeros(this.n, this.m);
            elseif (nargin == 3)
                this.isPrecond = 1;
                this.inv2Precond = P;
                this.buffS = zeros(this.n, 2*this.m);
                this.buffY = zeros(this.n, 2*this.m);
                this.buffGrad = zeros(this.n, 2*this.m);
            elseif (nargin == 4)
                this.isPrecond = 1;
                this.inv2Precond = P;
                this.RESCALE = RESCALE;
                this.buffS = zeros(this.n, 2*this.m);
                this.buffY = zeros(this.n, 2*this.m);
                this.buffGrad = zeros(this.n, 2*this.m);
            else
                error ('LBFGS initiation: Wrong input !')
            end
        end
        
        function vectorW = applyInv2Precond(this, vector)
            if isa(this.inv2Precond, 'function_handle')
                vectorW = this.inv2Precond(vector);
            else
                vectorW = this.inv2Precond * vector;
            end                
        end
        
        function []= newIteration (this, grad)
            if (this.RESCALE ~= 1)
                this.newIterationRESCALE(grad);
                return;
            end
            
            
            if ( size(grad,1) ~= this.n )
                error ('LBFGS update: Wrong size of grad !')
            end     
                       
            this.iter = this.iter+1;
            if (this.isPrecond == 0)
                this.buffGrad = grad;
            else
                this.buffGrad(:,1) = grad;
                this.buffGrad(:,2) = this.applyInv2Precond(grad);
            end
            
        end
        
        function [] = finishIteration (this, stepsize, grad)
            if (this.RESCALE ~= 1)
                this.finishIterationRESCALE(stepsize, grad);
                return;
            end
            
            pos = mod(this.iter-1, this.m)+1;
            
            if (this.isPrecond == 0)
                this.buffS(:,pos) = stepsize * this.lbfgsDir;
                this.buffY(:,pos) = grad - this.buffGrad;
                this.buffRho(pos) = dot(this.buffS(:,pos), this.buffY(:,pos));
                this.Hk0 = dot(this.buffY(:,pos), this.buffY(:,pos));
                this.Hk0 = this.buffRho(pos) / this.Hk0;
                this.buffRho(pos) = 1 / this.buffRho(pos);
            else
                this.buffS(:,pos) = stepsize * this.lbfgsDir(:,1);
                this.buffS(:,this.m+pos) = stepsize * this.lbfgsDir(:,2);
                gradW = this.applyInv2Precond(grad);
                this.buffY(:,pos) = grad - this.buffGrad(:,1);
                this.buffY(:,this.m+pos) = gradW - this.buffGrad(:,2);
                this.buffRho(pos) = dot(this.buffS(:,this.m+pos), this.buffY(:,pos));
                this.Hk0 = dot(this.buffY(:,pos), this.buffY(:,this.m+pos));
                if (this.Hk0 < eps || this.buffRho(pos) < eps)
                    warning('L-BFGS : Hessian inverse calculation, near singular Hessian !');
                end                
                this.Hk0 = this.buffRho(pos) / this.Hk0;
                this.buffRho(pos) = 1 / this.buffRho(pos);                
            end
            
        end
        
        
        function [searchDir] = computeNewDir (this) 
            if (this.RESCALE ~=1 )
                searchDir = this.computeNewDirRESCALE();
                return;
            end
            
            this.lbfgsDir = this.buffGrad;        
            alpha = zeros(this.m,1); 
            
            if (this.iter > 1)
                k = min(this.iter-1, this.m);
                
                for i = 1:k
                    pos = mod(this.iter-i-1, this.m)+1;
                    if (this.isPrecond == 0)                    
                        alpha(pos) = dot(this.buffS(:,pos), this.lbfgsDir(:,1));
                    else
                        alpha(pos) = dot(this.buffS(:,this.m+pos), this.lbfgsDir(:,1));
                    end
                    alpha(pos) = alpha(pos) * this.buffRho(pos);
                    
                    if (this.isPrecond == 0)
                        this.lbfgsDir = this.lbfgsDir - alpha(pos) * this.buffY(:,pos);
                    else
                        this.lbfgsDir(:,1) = this.lbfgsDir(:,1) - alpha(pos) * this.buffY(:,pos);
                        this.lbfgsDir(:,2) = this.lbfgsDir(:,2) - alpha(pos) * this.buffY(:,this.m+pos);
                    end                    
                end
                
                this.lbfgsDir = this.Hk0 * this.lbfgsDir;
                
                for i = (this.iter-k):(this.iter-1)
                    pos = mod(i-1, this.m)+1;
                    
                    if (this.isPrecond == 0)
                        beta = dot(this.buffY(:,pos), this.lbfgsDir);
                    else
                        beta = dot(this.buffY(:,pos), this.lbfgsDir(:,2));                        
                    end
                    beta = beta * this.buffRho(pos);
                    
                    if (this.isPrecond == 0)
                        this.lbfgsDir = this.lbfgsDir + (alpha(pos)-beta) * this.buffS(:,pos);
                    else
                        this.lbfgsDir(:,1) = this.lbfgsDir(:,1) + (alpha(pos)-beta) * this.buffS(:,pos);
                        this.lbfgsDir(:,2) = this.lbfgsDir(:,2) + (alpha(pos)-beta) * this.buffS(:,this.m+pos);
                    end
                                                       
                end
                
            end
            
            this.lbfgsDir = -this.lbfgsDir;
                        
            if (this.isPrecond == 0)
                searchDir = this.lbfgsDir;
            else
                searchDir = this.lbfgsDir(:,2);
            end
            
        end
        
        
        function []= newIterationRESCALE (this, grad)
            if ( size(grad,1) ~= this.n )
                error ('LBFGS update: Wrong size of grad !')
            end     
                       
            this.iter = this.iter+1;
            if (this.isPrecond == 0)
                this.buffGrad = grad  * this.RESCALE; % multiplied by RESCALE^1
            else
                this.buffGrad(:,1) = grad  * this.RESCALE; % multiplied by RESCALE^1
                this.buffGrad(:,2) = this.applyInv2Precond(grad) * this.RESCALE; % multiplied by RESCALE^1
            end
            
        end
        
        function [] = finishIterationRESCALE (this, stepsize, grad)
            grad = grad * this.RESCALE;
            
            pos = mod(this.iter-1, this.m)+1;
            
            if (this.isPrecond == 0)
                this.buffS(:,pos) = stepsize * this.lbfgsDir / this.RESCALE; % S not rescaled
                this.buffY(:,pos) = grad - this.buffGrad; % multiplied by RESCALE^1
                this.buffRho(pos) = dot(this.buffS(:,pos), this.buffY(:,pos)); % mutiplied by RESCALE^1
                this.Hk0 = dot(this.buffY(:,pos), this.buffY(:,pos)); % multiplied by RESCALE^2
                if (this.Hk0 < eps || this.buffRho(pos) < eps)
                    warning('L-BFGS : Hessian inverse calculation, near singular Hessian !');
                end 
                this.Hk0 = this.buffRho(pos) / this.Hk0 * this.RESCALE; % not RESCALE
                this.buffRho(pos) = 1 / this.buffRho(pos); % RESCALE^(-1)
            else
                this.buffS(:,pos) = stepsize * this.lbfgsDir(:,1) / this.RESCALE;
                this.buffS(:,this.m+pos) = stepsize * this.lbfgsDir(:,2) / this.RESCALE;
                gradW = this.applyInv2Precond(grad); % multiplied by RESCALE^1 
                this.buffY(:,pos) = grad - this.buffGrad(:,1); % multiplied by RESCALE^1
                this.buffY(:,this.m+pos) = gradW - this.buffGrad(:,2); % multiplied by RESCALE^1
                this.buffRho(pos) = dot(this.buffS(:,this.m+pos), this.buffY(:,pos)); % multiplied by RESCALE^1
                this.Hk0 = dot(this.buffY(:,pos), this.buffY(:,this.m+pos)); % multiplied by RESCALE^2
                if (this.Hk0 < eps || this.buffRho(pos) < eps)
                    warning('L-BFGS : Hessian inverse calculation, near singular Hessian !');
                end                
                this.Hk0 = this.buffRho(pos) / this.Hk0 * this.RESCALE; % not RESCALE
                this.buffRho(pos) = 1 / this.buffRho(pos); % multiplied by RESCALE^(-1)               
            end
            
        end        
        
        
        function [searchDir] = computeNewDirRESCALE (this) 
            this.lbfgsDir = this.buffGrad;        
            alpha = zeros(this.m,1); 
            
            if (this.iter > 1)
                k = min(this.iter-1, this.m);
                
                for i = 1:k
                    pos = mod(this.iter-i-1, this.m)+1;
                    if (this.isPrecond == 0)                    
                        alpha(pos) = dot(this.buffS(:,pos), this.lbfgsDir(:,1)); % mutiplied by RESCALE^1
                    else
                        alpha(pos) = dot(this.buffS(:,this.m+pos), this.lbfgsDir(:,1)); % mutiplied by RESCALE^1
                    end
                    alpha(pos) = alpha(pos) * this.buffRho(pos); % not RESCALE
                    
                    if (this.isPrecond == 0)
                        this.lbfgsDir = this.lbfgsDir - alpha(pos) * this.buffY(:,pos); % mutiplied by RESCALE^1
                    else
                        this.lbfgsDir(:,1) = this.lbfgsDir(:,1) - alpha(pos) * this.buffY(:,pos); % mutiplied by RESCALE^1
                        this.lbfgsDir(:,2) = this.lbfgsDir(:,2) - alpha(pos) * this.buffY(:,this.m+pos); % mutiplied by RESCALE^1
                    end                    
                end
                
                this.lbfgsDir = this.Hk0 * this.lbfgsDir; % mutiplied by RESCALE^1
                
                for i = (this.iter-k):(this.iter-1)
                    pos = mod(i-1, this.m)+1;
                    
                    if (this.isPrecond == 0)
                        beta = dot(this.buffY(:,pos), this.lbfgsDir); % mutiplied by RESCALE^2
                    else
                        beta = dot(this.buffY(:,pos), this.lbfgsDir(:,2));  % mutiplied by RESCALE^2                       
                    end
                    beta = beta * this.buffRho(pos) / this.RESCALE; % not RESCALE
                    
                    if (this.isPrecond == 0)
                        this.lbfgsDir = this.lbfgsDir + (alpha(pos) * this.RESCALE - beta * this.RESCALE) * this.buffS(:,pos); % mutiplied by RESCALE^1
                    else
                        this.lbfgsDir(:,1) = this.lbfgsDir(:,1) + (alpha(pos) * this.RESCALE - beta * this.RESCALE) * this.buffS(:,pos); % mutiplied by RESCALE^1
                        this.lbfgsDir(:,2) = this.lbfgsDir(:,2) + (alpha(pos) * this.RESCALE - beta * this.RESCALE) * this.buffS(:,this.m+pos); % mutiplied by RESCALE^1
                    end
                                                       
                end
                
            end
            
            this.lbfgsDir = -this.lbfgsDir;
                        
%             if (this.isPrecond == 0)
%                 searchDir = this.lbfgsDir;
%             else
%                 searchDir = this.lbfgsDir(:,2);
%             end
            
            if (this.isPrecond == 0)
                searchDir = this.lbfgsDir / this.RESCALE;
            else
                searchDir = this.lbfgsDir(:,2) / this.RESCALE;
            end
            
        end
        
        
    end
    
end

