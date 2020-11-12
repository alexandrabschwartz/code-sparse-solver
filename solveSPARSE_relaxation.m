function [x_opt, f_opt, support_opt, information] = solveSPARSE_relaxation(problem, options)

% This function is given an optimization problem with a weighted sparsity
% term in the objective function of the form
%    min f(x) + weights*|sign(x)|  s.t. xl <=   x  <= xu
%                                       bl <=  A*x <= bu
%                                       cl <= c(x) <= cu
% and solves it by reformulating it as as the relaxed problem
%    min_{x,y} f(x) + weights*(1-y)  s.t. xl <=   x  <= xu
%                                         bl <=  A*x <= bu
%                                         cl <= c(x) <= cu
%                                         yl <=   y  <= e
%                                         phi(x,y;t) <= 0
%                                         phi(-x,y;t) <= 0
% and solving the resulting NLP for decreasing values of t > 0.

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

% If you want  to use gradient information, specify the NLP solver or the 
% relaxation function, decide on the use of a lower bound or an initial 
% value for the auxiliary variable y, additionally provide an options struct
% with the fields
    % options.objectiveGradient = true or false
    % options.constraintsJacobian = true or false
    % options.NLPsolver = 'fmincon' or 'snopt'
    % options.relaxation = 'scholtes' or 'steffensen' or 'schwartz' or 'kadrani'
    % options.yl = 0, negative value or -inf
    % options.y_start = y_start in [yl, 1]
    
% The function returns
    % x_opt                    computed solution
    % f_opt                    objective function f(x_opt)
    % support_opt              support of x_opt
    % information.message      final exit message of the NLP solver
    % information.maxVio_box   maximum violation of box constraints
    % information.maxVio_lin   maximum violation of linear constraints
    % information.maxVio_nln   maximum violation of nonlinear constraints
    % information.iterations   number of NLPs solved


%% parameters

t_start = 1; % initial regularization parameter
sigma = 0.1; % reducing factor for the regularization parameter in each iteration
t_min = 10^-8; % lower bound for the regularization parameter
    
complementarity_tol = 10^-6; % tolerance for x .* y <= tol
x_tol = 10^-6; % tolerance to decide if x_i = 0 or not

iterations = 0; % counter for the number of iterations
iterations_max = 100;  % maximum number of iterations
    

%% set up missing options using default values

if nargin == 1
    options = [];
end
options = setupSPARSE_defaultOptions(options);


%% check problem data for completeness and set up missing entries using default values

[problem, n_x, n_lin, n_nln] = setupSPARSE_missingData(problem);


%% gather data for the auxiliary variable y

% The reformulation maximizes y as part of the objective function. Thus, a
% lower bound yl on y is not needed. But a lower bound yl <= 0 can be 
% provided to obtain a compact feasible set for y. 

% box constraints on y
yl = options.yl * ones(n_x,1);
yu = ones(n_x,1);

% initial value for auxiliary variable y
y_start = options.y_start * ones(n_x,1);


%% bounds on the constraints
% phi(-x,y;t) <= cl_xy
% phi( x,y;t) <= cu_xy

% If w_i = 0, we do not need to enforce x_i*y_i = 0, because y_i has no
% effect on feasibility or the objective function.

% If xl_i >= 0, we know that all feasible pairs (x_i, y_i) with maximal y_i
% are in the first quadrant, i.e. satisfy x_i*y_i >= 0. Thus we do not need
% a lower bound on x_i*y_i.

% If xu_i <= 0, we know that all feasible pairs (x_i, y_i) with maximal y_i
% are in the fourth quadrant, i.e. satisfy x_i*y_i <= 0. Thus we do not need
% an upper bound on x_i*y_i.

cl_xy = zeros(n_x,1);
cl_xy(problem.weights == 0) = inf; % free x-variables
cl_xy(problem.xl >= 0) = inf; % nonnegative x-variables

cu_xy = zeros(n_x,1);
cu_xy(problem.weights == 0) = inf; % free x-variables
cu_xy(problem.xu <= 0) = inf; % nonpositive x-variables


%% define nonlinear reformulation

NLPproblem.objective = @objective_with_auxiliary;
NLPproblem.xl =        [problem.xl; yl];
NLPproblem.xu =        [problem.xu; yu];
NLPproblem.A =         [problem.A zeros(n_lin, n_x)];
NLPproblem.bl =        problem.bl;
NLPproblem.bu =        problem.bu;
NLPproblem.nlcons =    @nlcons_with_auxiliary; % [c(x); phi(x,y;t); phi(-x,y;t)]
NLPproblem.cl =        [problem.cl; -inf(n_x,1); -inf(n_x,1)];
NLPproblem.cu =        [problem.cu; cu_xy; cl_xy];
NLPproblem.x_start =   [problem.x_start; y_start];
NLPproblem.dimension = 2*n_x;


%% solve the nonlinear reformulation for decreasing relaxation parameters

t = t_start;

while (iterations == 0) || ((iterations < iterations_max) && (t > t_min) &&...
        any(abs(x_opt(problem.weights > 0) .* y_opt(problem.weights > 0)) > complementarity_tol))
    
    % solve the relaxed problem
    [X_opt, ~, NLPinformation] = solveNLP(NLPproblem, options);
    
    % compute return values
    X_opt = X_opt(:);
    x_opt = X_opt(1:n_x);
    y_opt = X_opt(n_x+1:end);
    
    % update the initial point
    NLPproblem.x_start = X_opt;
    
    % decrease the relaxation parameter
    t = t*sigma;
    
    % update the iteration counter
    iterations = iterations + 1;   
end



%% compute return values

f_opt = problem.objective(x_opt);

support_opt = (abs(x_opt) > x_tol);

information.message = NLPinformation.message;
information.maxVio_box = max([max(x_opt-problem.xu, 0);...
                               max(problem.xl-x_opt, 0)]);
information.maxVio_lin = max([max(problem.A*x_opt-problem.bu, 0);...
                               max(problem.bl-problem.A*x_opt, 0)]);
information.maxVio_nln = max([max(problem.nlcons(x_opt)-problem.cu, 0);...
                               max(problem.cl-problem.nlcons(x_opt), 0)]);


%% auxiliary functions

function [F, DF] = objective_with_auxiliary(X)
    % rewritten objective with auxiliary variables y: 
    % F(x,y) = f(x) + weights*(1-y)
    % X = (x,y)
    X = X(:);
    x = X(1:n_x);
    y = X(n_x+1:end);
    
    if nargout == 1
        F = problem.objective(x) + problem.weights*(1-y);
        
    elseif nargout > 1
        [f, Df] = problem.objective(x);
        F = f + problem.weights*(1-y);
        % gradient of the objective is a row vector
        DF = [Df -problem.weights];
    end
end

function [C,DC] = nlcons_with_auxiliary(X)
    % rewritten nonlinear constraints with auxiliary variables y: [c(x); phi(x,y;t); phi(-x,y;t)]
    % X = (x,y)
    X = X(:);
    x = X(1:n_x);
    y = X(n_x+1:end);

    if nargout == 1
        c = problem.nlcons(x);
        phi_upper = relaxationSPARSE(x, y, t, options.relaxation);
        phi_lower = relaxationSPARSE(-x, y, t, options.relaxation);
        C = [c; phi_upper; phi_lower];
        
    elseif nargout > 1
        [c, Dc] = problem.nlcons(x);
        [phi_upper, Dphi_upper] = relaxationSPARSE( x, y, t, options.relaxation);
        [phi_lower, Dphi_lower] = relaxationSPARSE(-x, y, t, options.relaxation);
        C = [c; phi_upper; phi_lower];
        % gradients of the constraints are row vectors
        DC = [ Dc       zeros(n_nln, n_x);...
               diag(Dphi_upper(:,1))  diag(Dphi_upper(:,2))...
              -diag(Dphi_lower(:,1))  diag(Dphi_lower(:,2))];
    end
end            

end