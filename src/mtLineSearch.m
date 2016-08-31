classdef mtLineSearch < handle
    %  MTLineSearch Summary of this class goes here
    %  More and Thuente line search with Wolfe condtions
    
    properties        
        % Parameters of MT Line Search
        c1wolfe;
        c2wolfe;
        iter;
        maxIter; % maximum line search steps i.e. f evluation
        stepmin; 
        stepmax;
        tolX; % tolerance on x value, will be of sqrt of wolfe's condition around a local minimum
        
        % fucntion handle to evaluate fval and gradient and only depends on step size
        func;
        Xold;
        fold;
        gradold;
        dir;
        desc;
    end
    
    methods
        function this = mtLineSearch(c1,c2,maxIter)
            this.c1wolfe = 1e-4;
            this.c2wolfe = 0.9;
            this.maxIter = 20;
            this.stepmin = 1e-20;
            this.stepmax = 1e20;
            this.tolX = 1e-10;
            switch nargin
                case 0
                case 1
                    this.c1wolfe = c1;
                case 2
                    this.c1wolfe = c1;
                    this.c2wolfe = c2;
                case 3
                    this.c1wolfe = c1;
                    this.c2wolfe = c2;
                    this.maxIter = maxIter;
                otherwise
                    error('LineSearch Initiation: Wrong Input !')
            end
            
        end
        
        function [stepsize, X, fval, grad] = doLineSearch(this, func, Xold, fold, gradold, dir, stepsize0)
            if (nargin == 6)
                stepsize = 1;
            elseif (nargin == 7)
                stepsize = stepsize0;
            else
                error('doLineSearch: Wrong Input !')
            end
            
            this.func = func;
            this.Xold = Xold;
            this.fold = fold;
            this.gradold = gradold;
            this.dir = dir;
            this.desc = dot(gradold, dir);
            
            if (this.desc >= 0.0)
                fprintf('doLineSearch: wrong descent direction !\n');
		stepsize = 0;
		X = Xold;
		fval = fold;
		grad = gradold;
		return ;
            end
                
                
            [X, fval, grad, stepsize, info, this.iter] = this.mt_cvsrch( stepsize);
            
            
            
%             if ( info ~= 1)
%                 fprintf('Wolfe conditions not satisfied !\n');
%                 stepsize = 0;
%             end
            if (fval > this.fold)
                fprintf('Function cannot be decresed in current direction !\n');
                stepsize = 0;
            end
            
                
            
            
        end
        
        function [x,f,g,stp,info,nfev] = mt_cvsrch(this, stp)
%     **********
      p5 = .5;
      p66 = .66;
      xtrapf = 4;
      info = 0;
      infoc = 1;
        
      
      x = this.Xold;
      f = this.fold;
      g = this.gradold;
      s = this.dir;
      
      fcn = this.func;
      ftol = this.c1wolfe;
      gtol = this.c2wolfe;
      xtol = this.tolX;
      stpmin = this.stepmin;
      stpmax = this.stepmax;
      maxfev = this.maxIter;
      
      dginit = this.desc;
%     Initialize local variables.
%
      brackt = 0;
      stage1 = 1;
      nfev = 0;
      finit = f;
      dgtest = ftol*dginit;
      width = stpmax - stpmin;
      width1 = 2*width;
      wa = x;
%
%     The variables stx, fx, dgx contain the values of the step, 
%     function, and directional derivative at the best step.
%     The variables sty, fy, dgy contain the value of the step,
%     function, and derivative at the other endpoint of
%     the interval of uncertainty.
%     The variables stp, f, dg contain the values of the step,
%     function, and derivative at the current step.
%
      stx = 0.0;
      fx = finit;
      dgx = dginit;
      sty = 0.0;
      fy = finit;
      dgy = dginit;
%
%     Start of iteration.
%
   while (1)   
%
%        Set the minimum and maximum steps to correspond
%        to the present interval of uncertainty.
%
         if (brackt) 
            stmin = min(stx,sty);
            stmax = max(stx,sty);
         else
            stmin = stx;
            stmax = stp + xtrapf*(stp - stx);
         end 
%
%        Force the step to be within the bounds stpmax and stpmin.
%
         stp = max(stp,stpmin);
         stp = min(stp,stpmax);
%
%        If an unusual termination is to occur then let 
%        stp be the lowest point obtained so far.
%
         if ((brackt & (stp <= stmin | stp >= stmax)) ...
            | nfev >= maxfev-1 | infoc == 0 ...
            | (brackt & stmax-stmin <= xtol*stmax)) 
            stp = stx;
         end
%
%        Evaluate the function and gradient at stp
%        and compute the directional derivative.
%
         x = wa + stp * s;
         [f,g] = fcn(x);
         nfev = nfev + 1;
         dg = g' * s;
         ftest1 = finit + stp*dgtest;
%
%        Test for convergence.
%
         if ((brackt & (stp <= stmin | stp >= stmax)) | infoc == 0) 
                  info = 6;
         end
         if (stp == stpmax & f <= ftest1 & dg <= dgtest) 
                  info = 5;
         end
         if (stp == stpmin & (f > ftest1 | dg >= dgtest)) 
                  info = 4;
         end
         if (nfev >= maxfev) 
                  info = 3;
         end
         if (brackt & stmax-stmin <= xtol*stmax) 
                  info = 2;
         end
         if (f <= ftest1 & abs(dg) <= gtol*(-dginit)) 
                  info = 1;
         end
%
%        Check for termination.
%
         if (info ~= 0) 
                  return
         end
%
%        In the first stage we seek a step for which the modified
%        function has a nonpositive value and nonnegative derivative.
%
         if (stage1 & f <= ftest1 & dg >= min(ftol,gtol)*dginit) 
                stage1 = 0;
         end
%
%        A modified function is used to predict the step only if
%        we have not obtained a step for which the modified
%        function has a nonpositive function value and nonnegative 
%        derivative, and if a lower function value has been  
%        obtained but the decrease is not sufficient.
%
         if (stage1 & f <= fx & f > ftest1) 
%
%           Define the modified function and derivative values.
%
            fm = f - stp*dgtest;
            fxm = fx - stx*dgtest;
            fym = fy - sty*dgtest;
            dgm = dg - dgtest;
            dgxm = dgx - dgtest;
            dgym = dgy - dgtest;
% 
%           Call cstep to update the interval of uncertainty 
%           and to compute the new step.
%
            [stx,fxm,dgxm,sty,fym,dgym,stp,fm,dgm,brackt,infoc] ...
             = this.cstep(stx,fxm,dgxm,sty,fym,dgym,stp,fm,dgm, ...
                     brackt,stmin,stmax);
%
%           Reset the function and gradient values for f.
%
            fx = fxm + stx*dgtest;
            fy = fym + sty*dgtest;
            dgx = dgxm + dgtest;
            dgy = dgym + dgtest;
         else
% 
%           Call cstep to update the interval of uncertainty 
%           and to compute the new step.
%
            [stx,fx,dgx,sty,fy,dgy,stp,f,dg,brackt,infoc] ...
             = this.cstep(stx,fx,dgx,sty,fy,dgy,stp,f,dg, ...
                     brackt,stmin,stmax);
         end
%
%        Force a sufficient decrease in the size of the
%        interval of uncertainty.
%
         if (brackt) 
            if (abs(sty-stx) >= p66*width1) 
              stp = stx + p5*(sty - stx);
            end
            width1 = width;
            width = abs(sty-stx);
         end
%
%        End of iteration.            
            
            
            
            
       

   end
        
        
        

        
        
        
        end
    
        
             function  [stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,info]  = cstep(this,stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,stpmin,stpmax)
%   Translation of minpack subroutine cstep 
%   Dianne O'Leary   July 1991
%     **********
%
%     Subroutine cstep
%
%     The purpose of cstep is to compute a safeguarded step for
%     a linesearch and to update an interval of uncertainty for
%     a minimizer of the function.
%
%     The parameter stx contains the step with the least function
%     value. The parameter stp contains the current step. It is
%     assumed that the derivative at stx is negative in the
%     direction of the step. If brackt is set true then a
%     minimizer has been bracketed in an interval of uncertainty
%     with endpoints stx and sty.
%
%     The subroutine statement is
%
%       subroutine cstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
%                        stpmin,stpmax,info)
% 
%     where
%
%       stx, fx, and dx are variables which specify the step,
%         the function, and the derivative at the best step obtained
%         so far. The derivative must be negative in the direction
%         of the step, that is, dx and stp-stx must have opposite 
%         signs. On output these parameters are updated appropriately.
%
%       sty, fy, and dy are variables which specify the step,
%         the function, and the derivative at the other endpoint of
%         the interval of uncertainty. On output these parameters are 
%         updated appropriately.
%
%       stp, fp, and dp are variables which specify the step,
%         the function, and the derivative at the current step.
%         If brackt is set true then on input stp must be
%         between stx and sty. On output stp is set to the new step.
%
%       brackt is a logical variable which specifies if a minimizer
%         has been bracketed. If the minimizer has not been bracketed
%         then on input brackt must be set false. If the minimizer
%         is bracketed then on output brackt is set true.
%
%       stpmin and stpmax are input variables which specify lower 
%         and upper bounds for the step.
%
%       info is an integer output variable set as follows:
%         If info = 1,2,3,4,5, then the step has been computed
%         according to one of the five cases below. Otherwise
%         info = 0, and this indicates improper input parameters.
%
%     Subprograms called
%
%       FORTRAN-supplied ... abs,max,min,sqrt
%                        ... dble
%
%     Argonne National Laboratory. MINPACK Project. June 1983
%     Jorge J. More', David J. Thuente
%
%     **********
      p66 = 0.66;
      info = 0;
%
%     Check the input parameters for errors.
%
      if ((brackt & (stp <= min(stx,sty) | ...
          stp >= max(stx,sty))) | ...
          dx*(stp-stx) >= 0.0 | stpmax < stpmin) 
         return
      end
%
%     Determine if the derivatives have opposite sign.
%
      sgnd = dp*(dx/abs(dx));
%
%     First case. A higher function value.
%     The minimum is bracketed. If the cubic step is closer
%     to stx than the quadratic step, the cubic step is taken,
%     else the average of the cubic and quadratic steps is taken.
%
      if (fp > fx) 
         info = 1;
         bound = 1;
         theta = 3*(fx - fp)/(stp - stx) + dx + dp;
         s = norm([theta,dx,dp],inf);
         gamma = s*sqrt((theta/s)^2 - (dx/s)*(dp/s));
         if (stp < stx) 
             gamma = -gamma;
         end
         p = (gamma - dx) + theta;
         q = ((gamma - dx) + gamma) + dp;
         r = p/q;
         stpc = stx + r*(stp - stx);
         stpq = stx + ((dx/((fx-fp)/(stp-stx)+dx))/2)*(stp - stx);
         if (abs(stpc-stx) < abs(stpq-stx)) 
            stpf = stpc;
         else
           stpf = stpc + (stpq - stpc)/2;
         end 
         brackt = 1;
%
%     Second case. A lower function value and derivatives of
%     opposite sign. The minimum is bracketed. If the cubic
%     step is closer to stx than the quadratic (secant) step, 
%     the cubic step is taken, else the quadratic step is taken.
%
      elseif (sgnd < 0.0) 
         info = 2;
         bound = 0;
         theta = 3*(fx - fp)/(stp - stx) + dx + dp;
         s = norm([theta,dx,dp],inf);
         gamma = s*sqrt((theta/s)^2 - (dx/s)*(dp/s));
         if (stp > stx) 
            gamma = -gamma;
         end
         p = (gamma - dp) + theta;
         q = ((gamma - dp) + gamma) + dx;
         r = p/q;
         stpc = stp + r*(stx - stp);
         stpq = stp + (dp/(dp-dx))*(stx - stp);
         if (abs(stpc-stp) > abs(stpq-stp))
            stpf = stpc;
         else
            stpf = stpq;
         end 
         brackt = 1;
%
%     Third case. A lower function value, derivatives of the
%     same sign, and the magnitude of the derivative decreases.
%     The cubic step is only used if the cubic tends to infinity 
%     in the direction of the step or if the minimum of the cubic
%     is beyond stp. Otherwise the cubic step is defined to be 
%     either stpmin or stpmax. The quadratic (secant) step is also 
%     computed and if the minimum is bracketed then the the step 
%     closest to stx is taken, else the step farthest away is taken.
%
      elseif (abs(dp) < abs(dx)) 
         info = 3;
         bound = 1;
         theta = 3*(fx - fp)/(stp - stx) + dx + dp;
         s = norm([theta,dx,dp],inf);
%
%        The case gamma = 0 only arises if the cubic does not tend
%        to infinity in the direction of the step.
%
         gamma = s*sqrt(max(0.,(theta/s)^2 - (dx/s)*(dp/s)));
         if (stp > stx) 
             gamma = -gamma;
         end
         p = (gamma - dp) + theta;
         q = (gamma + (dx - dp)) + gamma;
         r = p/q;
         if (r < 0.0 & gamma ~= 0.0)
            stpc = stp + r*(stx - stp);
         elseif (stp > stx)
            stpc = stpmax;
         else
            stpc = stpmin;
         end 
         stpq = stp + (dp/(dp-dx))*(stx - stp);
         if (brackt) 
            if (abs(stp-stpc) < abs(stp-stpq)) 
               stpf = stpc;
            else
               stpf = stpq;
            end
         else
            if (abs(stp-stpc) > abs(stp-stpq)) 
               stpf = stpc;
            else
               stpf = stpq;
            end 
         end 
%
%     Fourth case. A lower function value, derivatives of the
%     same sign, and the magnitude of the derivative does
%     not decrease. If the minimum is not bracketed, the step
%     is either stpmin or stpmax, else the cubic step is taken.
%
      else
         info = 4;
         bound = 0;
         if (brackt) 
            theta = 3*(fp - fy)/(sty - stp) + dy + dp;
            s = norm([theta,dy,dp],inf);
            gamma = s*sqrt((theta/s)^2 - (dy/s)*(dp/s));
            if (stp > sty) 
                gamma = -gamma;
            end
            p = (gamma - dp) + theta;
            q = ((gamma - dp) + gamma) + dy;
            r = p/q;
            stpc = stp + r*(sty - stp);
            stpf = stpc;
         elseif (stp > stx)
            stpf = stpmax;
         else
            stpf = stpmin;
         end 
      end 
%
%     Update the interval of uncertainty. This update does not
%     depend on the new step or the case analysis above.
%
      if (fp > fx) 
         sty = stp;
         fy = fp;
         dy = dp;
      else
         if (sgnd < 0.0)
            sty = stx;
            fy = fx;
            dy = dx;
         end 
         stx = stp;
         fx = fp;
         dx = dp;
      end
%
%     Compute the new step and safeguard it.
%
      stpf = min(stpmax,stpf);
      stpf = max(stpmin,stpf);
      stp = stpf;
      if (brackt & bound)
         if (sty > stx) 
            stp = min(stx+p66*(sty-stx),stp);
         else
            stp = max(stx+p66*(sty-stx),stp);
         end
      end
%       return
%
%     Last card of subroutine cstep.
%
             end       
        
        
    
    end

end

