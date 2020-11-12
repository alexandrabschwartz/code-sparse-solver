%% test problem for SPARSE solver

%   min ||Ax-b||_2^2 + weight * sum(|sign(x)|)

% solution for weight in [0, 9]:
    % x = [-4; 3; 3] 
    % f_opt = 0 
    % support_opt = [1; 1; 1]
% solution for weight in [9, 18]:
    % x = [1; 0; 0]
    % f_opt = 18
    % support_opt = [1; 0; 0]
% solution for weight in [18, inf]:
    % x = [0; 0; 0]
    % f_opt = 36
    % support_opt = [0; 0; 0]
    
% x = [0; 06923; 06923] with f_opt = 11.0769 and support_opt = [0; 1; 1]
% minimizes f over all feasible points with cardinality 2, but is not optimal
% for the sparse problem for any weight.

weight = 25;
problem.objective = @objective_SPARSE;
problem.weights = weight * ones(1,3);
problem.xl = [];
problem.xu = [];
problem.A = [];
problem.bl = [];
problem.bu = [];
problem.nlcons = [];
problem.cl = [];
problem.cu = [];
problem.x_start = [-4; 3; 3];
problem.dimension = 3; 

options.objectiveGradient = true;
options.constraintsJacobian = true;
options.NLPsolver = 'fmincon' ;
options.yl = 0;
options.y_start = 1;


% call of the different solvers for sparse problems

disp('fmincon ============')
options.algorithm = 'direct';
[x_opt, f_opt, support, information] = solveSPARSE(problem,options);
x_opt
f_opt
% support
% information.message

disp('scholtes ============')
options.algorithm = 'relaxation';
options.relaxation = 'scholtes';
[x_opt, f_opt, support, information] = solveSPARSE(problem,options);
x_opt
f_opt
% support
% information.message

disp('steffensen ============')
options.algorithm = 'relaxation';
options.relaxation = 'steffensen';
[x_opt, f_opt, support, information] = solveSPARSE(problem,options);
x_opt
f_opt
% support
% information.message

disp('schwartz ============')
options.algorithm = 'relaxation';
options.relaxation = 'schwartz';
[x_opt, f_opt, support, information] = solveSPARSE(problem,options);
x_opt
f_opt
% support
% information.message

disp('kadrani ============')
options.algorithm = 'relaxation';
options.relaxation = 'kadrani';
[x_opt, f_opt, support, information] = solveSPARSE(problem,options);
x_opt
f_opt
% support
% information.message


%% objective functions and nonlinear constraints

function [f, Df] = objective_SPARSE(x)
    A = [0  3 -3; ...
         3  2  2; ...
         3  3  3];
    b = [0; 0; 6];

    f = sum((A*x-b).^2);
    
    if nargout > 1
        % derivative should be a row vector
        Df = (2*A'*(A*x-b))';
    end
end