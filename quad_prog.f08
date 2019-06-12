module quad_prog
    use solve_qp
    
implicit none

private

public:: quadprog

contains

! #############################################################################
subroutine quadprog(H,f,Aeq,beq,lb,ub, x,res,err)
    implicit none
    !Entrada
    double precision, dimension(:), intent(in):: f, beq, lb, ub
    double precision, dimension(:,:), intent(in):: H, Aeq
    !Saida
    double precision, dimension(:), intent(out):: x
    double precision, intent(out):: res
    integer, intent(out):: err
    !Local:
    double precision, allocatable, dimension(:):: f_loc, beq_loc
    double precision, allocatable, dimension(:,:):: H_loc, Aeq_loc, x_loc, res_loc, fm_loc
    integer:: n, q, qeq, i ,j
    ! Solver specific
    integer:: nact, iter(2,1)=0, iw, r
    integer, allocatable, dimension(:):: iact
    double precision, allocatable, dimension(:):: lagr, work
    double precision:: crval
    
    
    ! --------------
    ! Get sizes
    ! --------------
    ! Number of variables
    n = size(x) 
    ! Number of total constrains
    q = size(beq)+ size(lb)+ size(ub)
    ! Number of equality contrains
    qeq = size(beq)
    
    ! --------------
    ! Allocate local variables
    ! --------------
    allocate(f_loc(n))
    allocate(H_loc(n,n))
    allocate(beq_loc(q))
    allocate(Aeq_loc(n,q))
    allocate(x_loc(n,1))
    allocate(fm_loc(n,1))
    allocate(res_loc(1,1))
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
    call qpgen2(H_loc, f_loc, n, n, x, lagr, crval, Aeq_loc, beq_loc, n, q, qeq, iact, nact, iter, work, err)
    
    ! --------------
    ! Resirue Calculation
    ! --------------
    x_loc(:,1) = x
    fm_loc(:,1) = f
    res_loc = (0.5d0)*matmul( transpose(x_loc), matmul(H, x_loc)  ) + matmul( transpose(fm_loc), x_loc)
    res = res_loc(1,1)
        
    ! --------------
    ! Free local memory
    ! --------------
    deallocate(H_loc, Aeq_loc, x_loc, fm_loc)
    deallocate(f_loc, beq_loc, res_loc)
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

end module quad_prog