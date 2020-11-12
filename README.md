# code-sparse-solver
Some solvers to tackle optimization problems with a sparsity term in the objective function. The included functions are described below.

I did test this code, but it may still contain bugs. If you find errors or have suggestions for improvements, please let me know.

## solveSPARSE

This function takes a nonlinear optimization problem with a sparsity term in the objective function as well as optional options as input. Within the options, you can specify which solver for sparse problems you want to use. At them moment, DIRECT and RELAXATION are possible. The function then passes the problem on to the specified solver.

## solveSPARSE_direct

This function takes a nonlinear optimization problem with a sparsity term in the objective function as well as optional options as input. It reformulates the problem into an equivalent fully nonlinear optimization problem using continuous (i.e. not binary) auxiliary variables. The reformulation is the solved using an NLP solver of your choice.

## solveSPARSE_relaxation

This function takes a nonlinear optimization problem with a sparsity term in the objective function as well as optional options as input. It reformulates the problem into an equivalent fully nonlinear optimization problem using continuous (i.e. not binary) auxiliary variables. The reformulation contains complementarity-type constraints, which are relaxed using a relaxation function of your choice and a sequence of relaxed problems is solved using an NLP solver of your choice.

## setupSPARSE_missingData

This function takes a nonlinear optimization problem with a sparsity term in the objective function as input. It checks thee problem data for completeness and inserts missing data -- if possible -- using default values. E.g. if you did not specify box constraints on the variable it inserts -inf/inf as lower/upper bounds.

## setupSPARSE_defaultOptions

This function takes an options struct as input and sets up missing options using default values. E.g. if you did not specify the solver for spare problems, it chooses DIRECT.

## relaxationSPARSE

This function takes two vectors, a relaxation parameter and a type of relaxation functions as input and passes this data on to the specified relaxation function. at the momen, SCHOLTES, STEFFENSEN, SCHWARTZ and KADRANI are possible. If you want to add more types, you have to register them here.

## relaxationSPARSE_scholtes, relaxationSPARSE_steffensen, relaxationSPARSE_schwartz, relaxationSPARSE_kadrani

This function takes two vectors and a relaxation parameter as input an evaulates the relaxation function. The relaxation functions are called after one of the authors of the initial papers:
* S. Scholtes: *Convergence properties of a regularization scheme for mathematical programs with complementarity constraints*, SIAM J. Optim. 11, 918–936, 2001
* S. Steffensen and M. Ulbrich: *A new regularization scheme for mathematical programs with equilibrium constraints*, SIAM J. Optim. 20, 2504–2539, 2010
* C. Kanzow and A. Schwartz: *A new regularization method for mathematical programs with complementarity constraints with strong convergence properties*, SIAM J. Optim. 23, 770–798, 2013
* A. Kadrani et al: *A new regularization scheme for mathematical programs with complementarity constraints*, SIAM J. Optim. 20, 78–103, 2009

## testSPARSE

This script contains a toy problem to illustrate how the functions are called.
