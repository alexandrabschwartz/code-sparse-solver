# code-sparse-solver
Some solvers to tackle optimization problems with a sparsity term in the objective function. The included functions are described below.

I did test this code, but it may still contain bugs. If you find errors or have suggestions for improvements, please let me know.

## solveSPARSE

This function takes a nonlinear optimization problem with a sparsity term in the objective function as well as optional options as input. Within the options, you can specify, which solver for sparse problems you want to use. At the moment, DIRECT and RELAXATION are possible. The function then passes the problem on to the specified solver. If you want to add more solvers, you have to register them here.

## solveSPARSE_direct

This function takes a nonlinear optimization problem with a sparsity term in the objective function as well as optional options as input. It reformulates the problem into an equivalent fully nonlinear optimization problem using continuous (i.e. not binary) auxiliary variables. The reformulation is then solved using an NLP solver of your choice.

## solveSPARSE_relaxation

This function takes a nonlinear optimization problem with a sparsity term in the objective function as well as optional options as input. It reformulates the problem into an equivalent fully nonlinear optimization problem using continuous (i.e. not binary) auxiliary variables. The reformulation contains complementarity-type constraints, which are relaxed using a relaxation function of your choice and a sequence of relaxed problems is solved using an NLP solver of your choice.

Some of these algorithms -- partially in a version for cardinality-constrained problems -- are described here:
* O. Burdakov, C. Kanzow and A. Schwartz: *Mathematical Programs with Cardinality Constraints: Reformulation by Complementarity-type Constraints and a Regularization Method*, SIAM J. Optim. 26, 397–425, 2016
* M. Branda, M. Bucher, M. Cervinka, Michal and A. Schwartz: *Convergence of a Scholtes-type Regularization Method for Cardinality-Constrained Optimization Problems with an Application in Sparse Robust Portfolio Optimization*, Computational Optimization and Applications 70, 2017
* M. Feng, M, J.E. Mitchell, J.-S. Pang, X. Shen, A. Wächter: *Complementarity formulation of l0-norm optimization problems*, Technical Report, Industrial Engineering and Management Sciences, Northwestern University, Evanston, IL, USA, September 2013

## setupSPARSE_missingData

This function takes a nonlinear optimization problem with a sparsity term in the objective function as input. It checks the problem data for completeness and inserts missing data -- if possible -- using default values. E.g. if you did not specify box constraints on the variable, it inserts -inf/inf as lower/upper bounds.

## setupSPARSE_defaultOptions

This function takes an options struct as input and sets up missing options using default values. E.g. if you did not specify the solver for sparse problems, it chooses DIRECT.

## relaxationSPARSE

This function takes two vectors, a relaxation parameter and a type of relaxation functions as input and passes this data on to the specified relaxation function. At the moment, SCHOLTES, STEFFENSEN, SCHWARTZ and KADRANI are possible. If you want to add more types, you have to register them here.

## relaxationSPARSE_scholtes, relaxationSPARSE_steffensen, relaxationSPARSE_schwartz, relaxationSPARSE_kadrani

These functions take two vectors and a relaxation parameter as input and evaulate the relaxation function. The relaxation functions are called after one of the authors of the respective initial papers:
* S. Scholtes: *Convergence properties of a regularization scheme for mathematical programs with complementarity constraints*, SIAM J. Optim. 11, 918–936, 2001
* S. Steffensen and M. Ulbrich: *A new regularization scheme for mathematical programs with equilibrium constraints*, SIAM J. Optim. 20, 2504–2539, 2010
* C. Kanzow and A. Schwartz: *A new regularization method for mathematical programs with complementarity constraints with strong convergence properties*, SIAM J. Optim. 23, 770–798, 2013
* A. Kadrani et al: *A new regularization scheme for mathematical programs with complementarity constraints*, SIAM J. Optim. 20, 78–103, 2009

## testSPARSE

This script contains a toy problem to illustrate how the functions are called.
