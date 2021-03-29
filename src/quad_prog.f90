! #############################################################################
! This module has interfaces to call the 'R' quadprog routine, inplemented 
! in F77 code.
! #Uses a MATLAB sintax for quadprog.
!                                               /
!      /                         \             |      [A].{x} <= {b}
!     |  1                        |            |
! min | --- {x}'[H]{x} + {f}'.{x} | such that <|    [Aeq].{x} = {beq}
!  x  |  2                        |            |
!      \                         /             |  {lb} <= {x} <= {ub}
!                                               \
!
! #Fortran callings :
!   quadprog(n,q,qeq,H,f,Aeq,beq,lb,ub, x,tol,res,err)
!       OR
!   quadprog(n,q,qeq,H,f,Aeq,beq, x,tol,res,err)
!
! #C callings :
!   quadprog(n,q,qeq,H,f,Aeq,beq,lb,ub, x,tol,res,err)
!       OR
!   quadprog(n,q,qeq,H,f,Aeq,beq,NULL,NULL, x,tol,res,err)
!
! #Arguments:
! ##Inputs
!   `n` (int): Number of variables
!   `q` (int): Total number of contraints (equality + inequality)
!   `qeq` (int): Number of equality constrains qeq [0,q]
!   `H` (dble) [n,n]: Matrix [H] as in the problem description  
!   `f` (dble) [n]: Vector {f} as in the problem description
!   `Aeq` (dble) [q,n]: Joined matrices [A] and [Aeq] like:
!                        --     --
!                        | [Aeq] | [qeq,n]
!            Aeq [q,n] = |  [A]  | [q-qeq,n]
!                        --     --
!   `beq` (dble) [q]: Joined vector [ [beq] [b] ]
!   `lb` (dble) [n] Optional: Lower limits for x variable
!   `ub` (dble) [n] Optional: Upper limits for x variable
!   `tol` (dble): Solver tolerance
! ##Outputs
!   `x` (dble) [n]: Optimun values
!   `res` (dble): Solver residue
!   `err` (int): Error code, can be one of the following:
!        0 - Success
!        1 - Minimization problem has no solution
!        2 - Problem in matrix [H] decomposition
!        3 - Minimization problem solved but violates equality contrains
!        4 - Minimization problem solved but violates inequality contrains
!
! #############################################################################

module quad_prog
    use solve_qp
    use iso_c_binding
    
    implicit none

    interface quadprog
        module procedure:: quadprogF77_w_lim
        module procedure:: quadprogF77
    end interface

    !#####################################
    !            Module Usage:
    !#####################################
    private

    !The following subroutines/variables are the only
    !ones accessible outside this module
    public:: quadprog, qp_res_correction
    !-------------------------------------

contains

! #############################################################################
! C interface for Fortran quadprog solver
subroutine quadprog_C(n,q,qeq,vH,f,vAeq,beq,plb,pub, x,tol,res,err) bind(c,name='quadprog')
    implicit none
    !Inputs:
    integer(c_int), value, intent(in):: n, q, qeq
    type(c_ptr), intent(in):: plb, pub
    real(c_double), dimension(n), intent(in):: f
    real(c_double), dimension(q), intent(in):: beq
    real(c_double), dimension(n*n), intent(in):: vH
    real(c_double), dimension(q*n), intent(in):: vAeq
    real(c_double), intent(in):: tol
    !Outputs:
    real(c_double), dimension(n), intent(out):: x
    real(c_double), intent(out):: res
    integer(c_int), intent(inout):: err
    !Local:
    integer:: i, j
    double precision, dimension(n,n):: H
    double precision, pointer:: lb(:), ub(:)
    double precision, dimension(q,n):: Aeq

    ! --------------
    ! Transform row_based memory from C to Fortran
    ! --------------
    !{vH} to [H]
    H = 0.d0
    do i=1,n
        j = 1+n*(i-1)
        H(i,:) = vH(j:j+n-1)
    end do
    !{vAeq} to [Aeq]
    Aeq = 0.d0
    do i=1,qeq
        j = 1+n*(i-1)
        Aeq(i,:) = vAeq(j:j+n-1)
    end do

    if (c_associated(plb).and.c_associated(pub)) then
        call c_f_pointer(plb, lb, shape=[n])
        call c_f_pointer(pub, ub, shape=[n])
        call quadprog(n,q,qeq,H,f,Aeq,beq,lb,ub, x,tol,res,err)
    else
        call quadprog(n,q,qeq,H,f,Aeq,beq, x,tol,res,err)
    end if
    
end subroutine quadprog_C
! #############################################################################

! #############################################################################
! Prepare matrices to F77 solver with separated limit array
subroutine quadprogF77_w_lim(n,q,qeq,H,f,Aeq,beq,lb,ub, x,tol,res,err)
    implicit none
    !Inputs:
    integer, intent(in):: n, q, qeq
    double precision, dimension(n,n), intent(in):: H
    double precision, dimension(n), intent(in):: f, lb, ub
    double precision, dimension(q,n), intent(in):: Aeq
    double precision, dimension(q), intent(in):: beq
    double precision, intent(in):: tol
    !Outputs:
    double precision, dimension(n), intent(out):: x
    double precision, intent(out):: res
    integer, intent(inout):: err
    !Local:
    integer:: i, j, q2
    double precision, allocatable, dimension(:,:):: Aeq_loc
    double precision, allocatable, dimension(:):: beq_loc
    
    ! Number of total constrains
    q2 = q+2*n
    allocate(Aeq_loc(q2,n))
    allocate(beq_loc(q2))
    Aeq_loc = 0.d0
    beq_loc= 0.d0
    
    ! Insert limits on arrays
    if (q>0) then
        Aeq_loc(1:q,1:n) = Aeq
        beq_loc = [beq, ub, -lb]
    else
        beq_loc = [ub, -lb]
    end if

    ! Fill constrain matrix with limits
    do i=1,n
        j = q + i
        Aeq_loc(j,i) = 1.d0
    end do
    do i=1,n
        j = q + n + i
        Aeq_loc(j,i) = -1.d0
    end do

    call quadprogF77(n,q2,qeq,H,f,Aeq_loc,beq_loc, x,tol,res,err)
    
end subroutine quadprogF77_w_lim
! #############################################################################

! #############################################################################
! Prepare WorkSpace and call the quadprog solver
subroutine quadprogF77(n,q,qeq,H,f,Aeq,beq, x,tol,res,err)
    implicit none
    !Inputs:
    integer, intent(in):: n, q, qeq
    double precision, dimension(n,n), intent(in):: H
    double precision, dimension(n), intent(in):: f
    double precision, dimension(q,n), intent(in):: Aeq
    double precision, dimension(q), intent(in):: beq
    double precision, intent(in):: tol
    !Outputs:
    double precision, dimension(n), intent(out):: x
    double precision, intent(out):: res
    integer, intent(inout):: err
    !Local:
    double precision, allocatable, dimension(:):: f_loc, beq_loc
    double precision, allocatable, dimension(:,:):: H_loc, Aeq_loc
    !Específicas para o Solver
    integer:: nact, iter(2,1)=0, iw, r
    integer, allocatable, dimension(:):: iact
    double precision, allocatable, dimension(:):: lagr, work
    
    ! --------------
    ! Duplicate arrays to pass to F77 routine
    ! Signals changed to adequate MATLAB quadprog sintax to 'qpgen2' solver
    ! --------------
    allocate(f_loc,source=-f)
    allocate(H_loc,source=H)
    allocate(beq_loc,source=-beq)
    allocate(Aeq_loc,source=transpose(-Aeq))
    ! Check positive definite H matrix
    call H_check(H_loc)

    ! --------------
    ! Create needed workspace for F77 function
    ! --------------
    ! Work Space
    r = min(n,q)
    iw = 2*n+r*(r+5)/2 + 2*q +1
    allocate(work(iw))
    ! Other variables passed to F77
    allocate(lagr(q))
    allocate(iact(q))
    
    ! --------------
    ! F77 quadprog call
    ! --------------
    
    ! Initialize solutiol {x} values
    x = 1.d40
    call qpgen2(H_loc, f_loc, n, n, x, lagr, res, Aeq_loc, beq_loc, n, q, qeq, iact, nact, iter, work, err)

    call constrain_check(x,Aeq_loc,beq_loc,qeq,q,tol,err)

    ! --------------
    ! Free local memory
    ! --------------
    deallocate(H_loc, Aeq_loc)
    deallocate(f_loc, beq_loc)
    deallocate(work, lagr, iact)
    
    
end subroutine quadprogF77
! #############################################################################

! #############################################################################
! H matrix need to be decomposed by quadprog algorithm.
! Needs better and more reliable solution than this.
subroutine H_check(H)
    implicit none
    !Inputs
    double precision, dimension(:,:), intent(inout):: H
    !Local:
    integer:: i
    double precision, parameter:: eps = 1.d-10
    
    do i=1,minval(shape(H))
        ! Check for zero value in diagonal
        if (H(i,i)>-eps .and. H(i,i)<eps) then
            ! Insert small value on it
            H(i,i) = eps
        end if
    end do
    
end subroutine H_check
! #############################################################################

! #############################################################################
subroutine constrain_check(x,Aeql,beql,qeq,q,tol,errcode)
    implicit none
    !Inputs:
    integer, intent(in):: qeq, q
    double precision, dimension(:), intent(in):: x, beql
    double precision, dimension(:,:), intent(in):: Aeql
    double precision, intent(in):: tol
    !Outputs:
    integer, intent(inout):: errcode
    !Local:
    integer:: n
    double precision, allocatable, dimension(:,:):: x_loc, Res_eq, Res_ineq
    
    ! Get sizes
    n = size(x)
    
    if (errcode==0) then
        !Allocate temporary variable
        allocate(x_loc(n,1))
        !Copy Solution
        x_loc(:,1)=x
        
        if(qeq>0)then
            allocate(Res_eq(qeq,1))
            !Calculate Residue in equalities [Aeq]^T * {x}  - {beq}
            Res_eq  = matmul( transpose(Aeql(:,1:qeq)), x_loc )
            Res_eq(:,1) = Res_eq(:,1) - beql(1:qeq)
            !Error check
            if (.not. all(Res_eq(:,1)>=-tol) .and. all(Res_eq(:,1)<=tol) ) then
                errcode=3
                deallocate(x_loc, Res_eq)
                return
            end if
        end if
        
        !Calculate Residue in equalities [Aineq]^T * {x}  - {bineq}
        allocate(Res_ineq(q-qeq,1))
        Res_ineq = matmul( transpose(Aeql(:,qeq+1:)), x_loc )
        Res_ineq(:,1) = Res_ineq(:,1) - beql(qeq+1:)
        !Error check
        if (.not. all(Res_ineq(:,1)>=-tol) ) then
            errcode=4
            deallocate(x_loc, Res_ineq)
            if (allocated(Res_eq)) deallocate(Res_eq)
            return
        end if
        
        !Free memory of temporary variables
        deallocate(x_loc, Res_ineq)
        if (allocated(Res_eq)) deallocate(Res_eq)
        
    end if !errcode
    
end subroutine constrain_check
! #############################################################################

! #############################################################################
! Export the real equation system residue
function qp_res_correction(n,f,res) result(res_corr) bind(c)
    implicit none
    !Inputs:
    integer(c_int), value, intent(in):: n
    real(c_double), intent(in):: res
    real(c_double), dimension(n), intent(in):: f
    !Outputs:
    real(c_double):: res_corr
    !Local:
    double precision:: corr(1,1), f_loc(n,1)
    
    !Copy f variable to column matrix
    f_loc(:,1) = f
    !Calculate residue
    corr = res*(2.d0**2) + matmul(transpose(f_loc), f_loc)
    !Export to variable
    res_corr = corr(1,1)
    
end function
! #############################################################################

end module quad_prog