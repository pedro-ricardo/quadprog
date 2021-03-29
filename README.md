# QuadProg Solver
Simple Quadratic programming solver. The solver credit goes to a `R` package.

This repo simply translate that solver to use a **Matlab** like interface.

* Both C and Fortran usage examples
* Free from memmory leaks

## Tested With:
| Library  | Version  |
| -------- | -------- |
| GCC      | 10.2.1   |
| Valgrind | 3.15.0   |

## Usage
```
                                              /
     /                         \             |      [A].{x} <= {b}
    |  1                        |            |
min | --- {x}'[H]{x} + {f}'.{x} | such that <|    [Aeq].{x} = {beq}
 x  |  2                        |            |
     \                         /             |  {lb} <= {x} <= {ub}
                                              \
```
### Arguments:
* `n` (int): Number of variables
* `q` (int): Total number of contraints (equality + inequality)
* `qeq` (int): Number of equality constrains qeq (between 0 and q)
* `H` (double) [n,n]: Matrix [H] as in the problem description  
* `f` (double) [n]: Vector {f} as in the problem description
* `Aeq` (double) [q,n]: Joined matrices [A] and [Aeq] like:
```
            --     --
            | [Aeq] | [qeq,n]
Aeq [q,n] = |  [A]  | [q-qeq,n]
            --     --
```
* `beq` (double) [q]: Joined vector [ [beq] [b] ]
* `lb` (double) [n] Optional: Lower limits for x variable
* `ub` (double) [n] Optional: Upper limits for x variable
* `tol` (double): Solver tolerance

* `x` (double) [n]: Optimun values
* `res` (double): Solver residue
* `err` (int): Error code, can be one of the following:
  * 0 - Success
  * 1 - Minimization problem has no solution
  * 2 - Problem in matrix [H] decomposition
  * 3 - Minimization problem solved but violates equality contrains
  * 4 - Minimization problem solved but violates inequality contrains

### Fortran callings :
```fortran
    call quadprog(n,q,qeq,H,f,Aeq,beq,lb,ub, x,tol,res,err)
```
or
```fortran
    call quadprog(n,q,qeq,H,f,Aeq,beq, x,tol,res,err)
```

### C callings :
```c
    quadprog(n,q,qeq,H,f,Aeq,beq,lb,ub, x,tol,res,err)
```
or
```c
    quadprog(n,q,qeq,H,f,Aeq,beq,NULL,NULL, x,tol,res,err)
```
