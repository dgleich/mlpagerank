###Yongyang Yu (yu163@purdue.edu)

###David F. Gleich

==================================


Description:
------------

The code repository computes the eigenvector of the high-order Markov chains.
Formally, the problem is to solve the non-linear system of the following form,

x = alpha\*R\*kron(x,x) + (1-alpha)\*v

where x is the solution of the high-order eigen-problem, R the n x n^2
transition matrix (derived from 3 dimensional tensor transition), alpha damping
factor 0 \< alpha \< 1, and v the personalized vector.

In order to solve the problem efficiently, the following 3 algorithms are
considered.

1.  shifted iterative method

2.  non-shift iterative method

3.  Newton's method

We have shown that the problem has the unique solution when alpha \< 1/2. For
this case, all of the algorithms produce the same results. However, the problem
becomes complicated when alpha approaches 1. The convergence of Newton's method
depends on the initial choice of x. We also found that problems that do not
converge for non-shift iterative method will not converge for shifted iterative
method either.

### Shifted iterative method

An intuitive iterative method follows the step,

x_{k+1} = alpha\*R\*kron(x_k, x_k) + (1-alpha)\*v.

However, a bunch of counterexamples can show that this intuitive iterative
method works poorly when alpha -\> 1. A modified shifted version of the
iterative method has the step,

x_{k+1} = 1 / (1 + gamma)(alpha\*R\*(x_k, x_k) + gamma\*x_k + (1-alpha)\*v).

By choose gamma = 0.5, the shifted version of the iterative method works well
for most of the cases when alpha -\> 1. The shifted iterative method is not
perfect. We have examples to show that it takes a long time to converge.

### Non-shift iterative method

By manipulating the original system,

x = alpha\*R\*(kron(I, x)\*x / 2 + kron(x, I)\*x / 2) + (1-alpha)\*v

i.e.,

(I - alpha/2 \* R\*(kron(I, x) + kron(x, I)))\*x = (1-alpha)\*v,

a new version of the iterative method can be obtained.

(I - alpha/2 \* R\*(kron(I, x_k) + kron(x_k, I)))\*x_{k+1} = (1-alpha)\*v

During each iteration, a linear system of x_{k+1} is computed based on the
previous solution x_k.

### Newton's method

By viewing the original system as

F(x) = x - alpha\*R\*kron(x,x) - (1-alpha)\*v

Newton's method can be formed as

J(x)(x_{k+1} - x_k) = -F(x_k)

where J(x) is the Jacobian matrix of F(x).

To be concrete,

(I - alpha\*R\*(kron(x_k, I) + kron(I, x_k)))\*x_{k+1} = -alpha\*R\*kron(x_k,
x_k) + (1-alpha)\*v

Since Newton's method is a second order method, it is expected that the
convergence should be fast if the initial point is well chosen.

===================================

### Code

The code is implemented in Matlab. The core functions are listed below,

tensorRank.m

This function implements the shifted iterative method with default tolerance of
1e-12. The shift factor gamma can be specified by the user or 0.5 by default.

tensorRankNonShift.m

This function implements the non-shift iterative method with default tolerance
of 1e-12.

tensorNewton.m

This function implements the Newton's method with default tolerance of 1e-12.

valSol.m

This function check if the solution satisfies the original equation with
tolerance of 1e-10.

bad_example_generator.m

This function tries to generates "bad examples" that do not converge for the
algorithms.

gen6x6x6_bad_example.m

We use this script to generate the non-converge examples in the test case
directory.

genSymScript.m

This function takes R as the argument and assumes the personalized vector v is
uniformly distributed.

A new file called 'mySym.m' is generated. The new file can be used to computed
the symbolic solution of the original problem.

NOTE: This script only works with Matlab 2012a and subsequent versions.

====================================

oneTwoStepChain.m

This function computes the eigenvector X of the following one-step chain and
two-step chain.

(alpha\*P + (1-alpha)\*V)\*X1 = X1

(alpha\*P + (1-alpha)\*V)(alpha\*P + (1-alpha)\*V)\*X2 = X2

twoStepChainSimplified.m

This function assumes X can be written in the form X = kron(x, x) and computes
x. x can be viewed as the solution of the following equation.

alpha\*R\*(alpha\*P + (1-alpha)\*V)\*kron(x, x) + (1-alpha)\*v = x

convertR2P.m

P is the n^2-by-n^2 transition matrix derived from the transition tensor. R is
the n-by-n^2 transition matrix used in the tensorPageRank problem. This function
provides a way to convert R to P.

====================================

plot_fig:

This directory contain all the Matlab files to generate the triangle plots.

Two kinds of triangle plots are produced by the code:

1) ratio of two consecutive steps (kappa)

2) second largest eigenvalue of the Jacobian

====================================

### Test Cases

The test_cases directory contains special examples we have discovered. The
examples are collected by testing all 3-by-3-by-3 transition tensors whose
elements are 0 or 1 and randomly sampling of 6-by-6-by-6 tensors.

example.mat

This file contains 4 example 3-by-3-by-3 tensors, R1, R2, R3, and R4.

These examples are used to generate the triangle plots for the kappa values and
second largest eigenvalues of Jacobians.

counterExample.mat

This file contains a 6-by-6-by-6 tensor that has 3 different eigenvectors for
alpha=0.99. This example has multiple solutions when alpha \> 0.5.

6x6x6_slow_converge_shift.mat

This file contains two 6-by-6-by-6 tensors. The first does not converge for
shifted iterative method and the second converge very slowly for gamma=0.5; y1
and y2 are solutions to the two problem.

6x6x6_not_converge_non_shift.mat

This file contains two 6-by-6-by-6 tensors. Both do not converge for non-shift
method and shifted method. The solution x1 and x2 are computed by Newton's
method with randomly choice of the initial starting point.

  
  

