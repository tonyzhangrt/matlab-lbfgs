% L-BFGS with preconditioner
clear optim
optim = optimLineSearch();
options = struct('maxIter',1000,'Display','iter-detailed','TolX',1e-12,'TolFun',1e-12);


% MATLAB supplied fminunc with BFGS
options1 = struct('Diagnostic','on','GradObj','on','TolX',1e-12,'TolFun',1e-12,'Display','iter-detailed','Algorithm','quasi-newton','HessUpdate','bfgs');
% [X,fval] = fminunc(@testfunc, [2;3;4;5], options1)
% 
[X,fval] = optim.runOptimLineSearch(@testfunc, [2;3;4;5], options )
[X,fval] = optim.runOptimLineSearch(@testfunc, [2;3;4;5], options, diag([0.01;1;1;1]) )
[X,fval] = optim.runOptimLineSearch(@testfunc, [2;3;4;5], options, diag([0.01;1;1;1]), 100 )

% % [X,fval] = optim.runOptimLineSearch(@testfunc, [4;5], options, diag([0.01;1;1;1]) )
% [X,fval] = fminunc(@rosenbrockfunc, [-1;2], options1)
% [X,fval] = optim.runOptimLineSearch(@rosenbrockfunc, [-1;2], options )
