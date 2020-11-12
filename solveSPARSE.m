function [x_opt, f_opt, support, information] = solveSPARSE(problem, options)

% This function is given an optimization problem with a weighted sparsity
% term in the objective function of the form
%    min f(x) + weights*|sign(x)|  s.t. xl <=   x  <= xu
%                                       bl <=  A*x <= bu
%                                       cl <= c(x) <= cu
% and solves it using different algorithms.

% The problem should be provided as a struct with the following fields: 
    % problem.objective = @objective, may be empty
    % problem.weights = (1 x n_x) weight vector, where entry w_i > 0,
    %                    if x_i affects the sparsity term, and w_i = 0 else
    % problem.xl = xl
    % problem.xu = xu
    % problem.A = A
    % problem.bl = bl
    % problem.bu = bu
    % problem.nlcons = @nlcons
    % problem.cl = cl
    % problem.cu = cu
    % problem.x_start = x_start
    % problem.dimension = n_x 
% For the objective function and the nonlinear constraints the respective  
% functions can either only return the function value or additionally the
% gradients (oriented row-wise). The default assumption is that no gradients
% are provided. 

% If you want  to use gradient information, specify the type of solution 
% algorithm, the used relaxation function or the NLP solver or decide on the
% use of a lower bound or an initial value for  the auxiliary variable y, 
% additionally provide an options struct with the fields
    % options.objectiveGradient = true or false
    % options.constraintsJacobian = true or false
    % options.NLPsolver = 'fmincon' or 'snopt'
    % options.algorithm = 'direct' or 'relaxation'
    % options.relaxation = 'scholtes' or 'steffensen' or 'schwartz' or 'kadrani'
    % options.yl = 0, negative value or -inf
    % options.y_start = y_start, should be in [yl, 1]

% The function returns
    % x_opt                    computed solution
    % f_opt                    objective function f(x_opt)
    % support_opt              support of x_opt
    % information.message      exit message of the solver
    % information.maxVio_box   maximum violation of box constraints
    % information.maxVio_lin   maximum violation of linear constraints
    % information.maxVio_nln   maximum violation of nonlinear constraints
    % information.iterations   number of NLPs solved


%% set up missing options using default values

if nargin == 1
    options = [];
end
options = setupSPARSE_defaultOptions(options);


%% call the specified solution algorithm for SPARSE problems

switch lower(options.algorithm)
    case 'direct'
        [x_opt, f_opt, support, information] = solveSPARSE_direct(problem, options);
    case 'relaxation'
        [x_opt, f_opt, support, information] = solveSPARSE_relaxation(problem, options);
    otherwise
        disp('Unknown SPARSE algorithm, will use direct NLP reformulation instead')
        [x_opt, f_opt, support, information] = solveSPARSE_direct(problem, options);
end
        
        

