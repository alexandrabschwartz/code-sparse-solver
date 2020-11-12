function options = setupSPARSE_defaultOptions(options)

% This function is given an empty argument or an options struct with some
% or all of the following fields:
    % options.objectiveGradient = true or false
    % options.constraintsJacobian = true or false
    % options.NLPsolver = 'fmincon' or 'snopt'
    % options.algorithm = 'direct' or 'relaxation'
    % options.relaxation = 'scholtes' or 'steffensen' or 'schwartz' or 'kadrani'
    % options.yl = 0, negative value or -inf
    % options.y_start = y_start, should be in [yl; 1]
    
% It sets up missing or empty fields using the following default values:
    % options.objectiveGradient = false
    % options.constraintsJacobian = false
    % options.NLPsolver = 'fmincon' 
    % options.algorithm = 'direct' 
    % options.relaxation = 'scholtes'
    % options.yl = 0
    % options.y_start = 1

%% parameters

y_lower_default = 0;
y_start_default = 1;


%% set up default options

if isempty(options)
    options.objectivegradient = [];
    options.constraintsjacobian = [];
    options.NLPsolver = [];
    options.algorithm = [];
    options.relaxation = [];
    options.yl = [];
    options.y_start = [];
end

if ~isfield(options, 'objectiveGradient') || isempty(options.objectiveGradient)
    % default value is no gradient information
    options.objectiveGradient = false;
end

if ~isfield(options, 'constraintsJacobian') || isempty(options.constraintsJacobian)
    % default value is no gradient information
    options.constraintsJacobian = false;
end

if ~isfield(options, 'NLPsolver') || isempty(options.NLPsolver)
    % default NLP solver is fmincon
    options.NLPsolver = 'fmincon';
end

if ~isfield(options, 'algorithm') || isempty(options.algorithm)
    % default algorithm is a direct NLP reformulation
    options.algorithm = 'direct';
end

if ~isfield(options, 'relaxation') || isempty(options.relaxation)
    % default relaxation is scholtes
    options.relaxation = 'scholtes';
end

if ~isfield(options, 'yl') || isempty(options.yl)
    % default lower bound on the auxiliary variables y
    options.yl = y_lower_default;
end

if ~isfield(options, 'y_start') || isempty(options.y_start)
    % default initial value for auxiliary variable y
    options.y_start = y_start_default;
end