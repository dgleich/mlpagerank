Multilinear PageRank
====================

#### David F. Gleich (@dgleich), Lek-Heng Lim, Yongyang Yu (@yuyongyang800)

The third order multilinear PageRank vector is the solution of:

    x = alpha*R*kron(x,x) + (1-alpha)*v
    
where *R* is an *n x n^2* column stochastic Matrix, *alpha* is a probability and
*v* is a column stochastic vector. The solution *x* is a probability distribution
that models the stationary distribution of a spacey random surfer.

The multilinear PageRank vector is also a "tensor PageRank" vector 

Synopsis
--------

    load tensors/mtm3.mat
    x = tensorpr3(R3_1,0.45).solve() % solve a problem with alpha = 0.45
    tpr = tensorpr3(R3_1,0.45); % create a class to 
    x = tpr.newton() % solve with Newton's iteration
    x = tpr.inverseiter() % solve with an Inverse Iteration
    tpr.residual(x)
    x = tpr.newton('tol',1e-15) % high accuracy solution
    tpr.residual(x)

Description:
------------

The code repository computes the eigenvector of the high-order Markov chains.
Formally, the problem is to solve the non-linear system of the following form,

    alpha R (x kron x) + (1-alpha) v = x, or
    alpha * R * kron(x,x) + (1-alpha) * v == x (in Matlab syntax)

where x is the solution of the high-order eigen-problem, R the n x n^2
transition matrix (derived from 3 dimensional tensor transition), alpha is 
a damping factor 0 \< alpha \< 1, and v is a the localization vector.

We include five algorithms to solve this problem:

1.  fixed point iterative method
2.  shifted fixed point iterative method
3.  an inner-outer iteration
4.  an inverse iteration
5.  a Newton's method

We have shown that the problem has the unique solution when alpha \< 1/2. For
this case, all of the algorithms produce the same results. However, the problem
becomes complicated when alpha approaches 1. The convergence of Newton's method
depends on the initial choice of x. We also found that problems that do not
converge for non-shift iterative method will not converge for shifted iterative
method either.

Errata
------

The original code contains a bug in the convergence detection which 
causes it detect convergence even when the iterate did not converge to
a sufficiently small residual. Thanks to Federico Poloni for pointing
this out. 

Details
--------    

`tensorpr3` : a class to work with 3rd order multilinear
              PageRank problem. This problem class
              takes a stochastic tensor represented
              by R or a 3rd order tensor P, alpha,
              and v.

`tensorpr3/solve` : our recommended solver for
                    tensorpr3 problems for use in 
                    other codes

`tensorpr3/power` : a fixed point
`tensorpr3/shifted` : a shifted fixed point method
`tensorpr3/inverse_iter` : an inverse iteration
`tensorpr3/newton` : a Newton iteration
`tensorpr3/innout` : an inner-outer iteration

`tensorpr3/residual` : check the residual

`li_gamma` : compute the value of gamma from Li and Ng 
             to determine if a solution of a tensor
             problem is unique.

Paper
-----

We have a preprint associated with this code. In order to reproduce
the figures from the preprint, please use the following 

### Example 3.1 solution

    >> setup
    >> cd examples
    >> higherorder_pagerank_examples
    
### Example 3.1 multilinear PageRank solution

    >> setup % can skip again
    >> cd examples
    >> tensor_pagerank_examples
    
### Examples 4.2 multiple solutions
        
    >> setup % can skip again
    >> cd examples
    >> multiple_solutions
    
### Figures 2, 3, 4, 5, 7

    >> setup
    >> cd examples
    >> illustration_figures
    
## Figure 6

    >> setup
    >> cd experiments/newton-variations/
    >> newton_difference_figs
    
## Figure 8

    >> setup
    >> cd experiments/shiftstudy/
    >> shiftstudy_R4_19_final
    
## Table 1, Table 2    

    >> setup
    >> cd experiments/solver_table
    >> generate_solver_table_shifting % Table 1
    >> generate_solver_table; generate_solver_table_long; write_final_table % Table 2
