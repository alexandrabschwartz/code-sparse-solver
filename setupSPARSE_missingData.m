function [problem, n_x, n_lin, n_nln] = setupSPARSE_missingData(problem)

% This function is given an optimization problem with a weighted sparsity term
% in the objective function of the form
%    min f(x) + weights*|sign(x)|  s.t. xl <=   x  <= xu
%                                       bl <=  A*x <= bu
%                                       cl <= c(x) <= cu.

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
    
% The function sets up missing data using default values and determines
    % n_x = number of variables
    % n_lin = number of linear constraints
    % n_nln = number of nonlinear constraints
% In case linear or nonlinear constraints are missing, a dummy constraint of
% this type with n_... = 0 is inserted to simplify the later reformulation of 
% the sparse problem.


%% parameters

x_start_default = 0;
weight_default = 1;


%% dimension of the optimization variable

if ~isfield(problem, 'dimension') || isempty(problem.dimension)
    % try to determine problem dimension from initial value or box constraints
    if ~isfield(problem, 'x_start') || isempty(problem.x_start)
        if ~isfield(problem, 'xl') || isempty(problem.xl)
            if ~isfield(problem, 'xu') || isempty(problem.xu)
                error('Problem dimension is missing')
            else
                problem.dimension = length(problem.xu);
            end
        else
            problem.dimension = length(problem.xl);
        end
    else
        problem.dimension = length(problem.x_start);
    end   
end
n_x = problem.dimension;


%% initial value for x

if ~isfield(problem, 'x_start') || isempty(problem.x_start)
    % use default initial value for x
    problem.x_start = x_start_default * ones(n_x,1);
elseif length(problem.x_start) == 1
    % use given initial value for all x_i
    disp('Only one initial value for x provided, using it for all variables')
    problem.x_start = problem.x_start * ones(n_x,1);
end


%% objective function

if ~isfield(problem, 'objective') || isempty(problem.objective)
    % no objective function f is present, use dummy function instead
    disp('no objective functions is provided, minimizing weighted support only')
    problem.objective = @(x) noObjective(x);
end


%% weight matrix

if ~isfield(problem, 'weights') || isempty(problem.weights)
    disp(['Weight matrix is missing. I will assume that all variables have the same weight' num2str(weight_default) '.'] )
    problem.weights = weight_default*ones(1,n_x);
end


%% box constraints on x

if ~isfield(problem, 'xl') || isempty(problem.xl)
    problem.xl = -inf(n_x,1);
end

if ~isfield(problem, 'xu') || isempty(problem.xu)
    problem.xu = inf(n_x,1);
end


%% linear constraints

if ~isfield(problem, 'A') || isempty(problem.A)
    problem.A = zeros(0,n_x);
end
n_lin = size(problem.A,1);

if ~isfield(problem, 'bl') || isempty(problem.bl)
    problem.bl = -inf(n_lin,1);
end

if ~isfield(problem, 'bu') || isempty(problem.bu)
    problem.bu = inf(n_lin,1);
end


%% nonlinear constraints

if ~isfield(problem, 'nlcons') || isempty(problem.nlcons)
    problem.nlcons = @(x) noNlcons(x);
end
n_nln = length(problem.nlcons(problem.x_start));

if ~isfield(problem, 'cl') || isempty(problem.cl)
    problem.cl = -inf(n_nln,1);
end

if ~isfield(problem, 'cu') || isempty(problem.cu)
    problem.cu = inf(n_nln,1);
end


%% auxiliary functions

function [c, Dc] = noNlcons(~)
    % dummy function in case no nonlinear constraints are present
    c = zeros(0,1);
    if nargout > 1
        % gradient is a row vector
        Dc = zeros(0,n_x);
    end
end

function [f, Df] = noObjective(~)
    % dummy function in case no objective function besides the sparsity term is present
    f = zeros(0,1);
    if nargout > 1
        Df = zeros(0,n_x);
    end
end

end