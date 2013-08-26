Authors: Yongyang Yu (yu163@purdue.edu)

David F. Gleich

=============================

Description:
------------

The code repository computes the eigenvector of the high-order Markov chains.
Formally, the problem is to solve the non-linear system of the following form,

x = alpha\*R\*kron(x,x) + (1-alpha)\*v

where x is the solution of the high-order eigen-problem, R the nxn^2 transition
matrix (derived from 3 dimensional tensor transition), alpha damping factor 0 \<
alpha \< 1, and v the personalized vector.

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

===========================


