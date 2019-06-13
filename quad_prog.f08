module quad_prog
    use solve_qp
    use iso_c_binding
    
implicit none

private

public:: quadprog, qp_res_correction

contains

! #############################################################################
! Import 'Matlab like' structures and convert to use F77 quadprog algorithm
! Error codes are:
! 0 - Success
! 1 - Minimization problem has no solution
! 2 - Problem in matrix [H] decomposition
! 3 - Minimization problem solved but violates equality contrains
! 4 - Minimization problem solved but violates inequality contrains

subroutine quadprog(n,qeq,vH,f,vAeq,beq,lb,ub, x,tol,res,err) bind(c)
    implicit none
    !Entrada
    real(c_double), dimension(n), intent(in):: f, lb, ub
    real(c_double), dimension(qeq), intent(in):: beq
    real(c_double), dimension(n*n), intent(in):: vH
    real(c_double), dimension(qeq*n), intent(in):: vAeq
    real(c_double), intent(in):: tol
    integer(c_int), intent(in):: n, qeq
    !Saida
    real(c_double), dimension(n), intent(out):: x
    real(c_double), intent(out):: res
    integer(c_int), intent(out):: err
    !Local:
    double precision, dimension(n,n):: H
    double precision, dimension(qeq,n):: Aeq
    double precision, allocatable, dimension(:):: f_loc, beq_loc
    double precision, allocatable, dimension(:,:):: H_loc, Aeq_loc
    integer:: q, i ,j
    ! Solver specific
    integer:: nact, iter(2,1)=0, iw, r
    integer, allocatable, dimension(:):: iact
    double precision, allocatable, dimension(:):: lagr, work
    
    ! --------------
    !Transform row_based memory from C to Fortran
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
    
    ! Number of total constrains
    q = qeq+2*n
    
    ! --------------
    ! Allocate local variables
    ! --------------
    allocate(f_loc(n))
    allocate(H_loc(n,n))
    allocate(beq_loc(q))
    allocate(Aeq_loc(n,q))
    ! Unused variables passed to F77
    allocate(lagr(q))
    allocate(iact(q))
    
    ! --------------
    ! Translate Matlab syntax to F77 function syntax
    ! --------------
    ! {f} vector
    f_loc = -f
    ! [H] matrix
    H_loc = H
    call H_check(H_loc)
    ! {beq} vector
    beq_loc = [beq, lb, -ub]
    ! [Aeq] matrix
    Aeq_loc = 0.d0
    Aeq_loc(1:n,1:qeq) = transpose(Aeq)
    do i=1,n
        j = qeq + i
        Aeq_loc(i,j) = 1.d0
    end do
    do i=1,n
        j = qeq + n + i
        Aeq_loc(i,j) = -1.d0
    end do
    ! Inform [H] status (not decomposed)
    err = 0
    
    
    ! --------------
    ! Create needed workspace for F77 function
    ! --------------
    ! Work Space
    r = min(n,q)
    iw = 2*n+r*(r+5)/2 + 2*q +1
    allocate(work(iw))
    
    
    ! --------------
    ! F77 quadprog call
    ! --------------
    ! Initialize solutiol {x} values
    x = 1.d40
    call qpgen2(H_loc, f_loc, n, n, x, lagr, res, Aeq_loc, beq_loc, n, q, qeq, iact, nact, iter, work, err)

    call constrain_check(x,Aeq_loc,beq_loc,qeq,tol,err)

    ! --------------
    ! Free local memory
    ! --------------
    deallocate(H_loc, Aeq_loc)
    deallocate(f_loc, beq_loc)
    deallocate(work, lagr, iact)
    
end subroutine quadprog
! #############################################################################

! #############################################################################
! H matrix need to be decomposed by quadprog algorithm.
! Needs better and more reliable solution than this.
subroutine H_check(H)
    implicit none
    !Entrada
    double precision, dimension(:,:), intent(inout):: H
    !Local:
    integer:: i
    double precision, parameter:: eps = 1.d-10
    
    
    do i=1,minval(shape(H))
        ! Check for zero value in diagonal
        if (H(i,i)==0.d0) then
            ! Insert small value on it
            H(i,i) = eps
        end if
    end do
    
end subroutine H_check
! #############################################################################

! #############################################################################
subroutine constrain_check(x,Aeql,beql,qeq,tol,errcode)
    implicit none
    !Entrada:
    integer, intent(in):: qeq
    double precision, dimension(:), intent(in):: x, beql
    double precision, dimension(:,:), intent(in):: Aeql
    double precision, intent(in):: tol
    integer, intent(inout):: errcode
    !Local:
    integer:: n, q
    double precision, allocatable, dimension(:,:):: x_loc, Res_eq, Res_ineq
    
    ! Get sizes
    n = size(x)
    q = 2*n + qeq
    
    if (errcode==0) then
        !Allocate temporary variable
        allocate(x_loc(n,1))
        allocate(Res_eq(qeq,1), Res_ineq(2*n,1))
        !Copy Solution
        x_loc(:,1)=x
        
        !Calculate Residue in equalities [Aeq]^T * {x}  - {beq}
        Res_eq  = matmul( transpose(Aeql(:,1:qeq)), x_loc )
        Res_eq(:,1) = Res_eq(:,1) - beql(1:qeq)
        !Error check
        if (.not. all(Res_eq(:,1)>=-tol) .and. all(Res_eq(:,1)<=tol) ) then
            errcode=3
            deallocate(x_loc, Res_eq, Res_ineq)
            return
        end if
        
        !Calculate Residue in equalities [Aineq]^T * {x}  - {bineq}
        Res_ineq = matmul( transpose(Aeql(:,qeq+1:)), x_loc )
        Res_ineq(:,1) = Res_ineq(:,1) - beql(qeq+1:)
        !Error check
        if (.not. all(Res_ineq(:,1)>=-tol) ) then
            errcode=4
            deallocate(x_loc, Res_eq, Res_ineq)
            return
        end if
        
        !Free memory of temporary variables
        deallocate(x_loc, Res_eq, Res_ineq)
        
    end if !errcode
    
end subroutine constrain_check
! #############################################################################

! #############################################################################
function qp_res_correction(n,f,res) result(res_corr) bind(c)
    implicit none
    !Entrada
    integer(c_int), intent(in):: n
    real(c_double), intent(in):: res
    real(c_double), dimension(n), intent(in):: f
    !Saida
    real(c_double):: res_corr
    !Local
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