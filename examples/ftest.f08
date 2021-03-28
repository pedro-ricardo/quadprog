program Ftest
    use quad_prog
    implicit none

    ! Sizes
    integer, parameter:: n=3, qeq=2, q=2
    ! Matrix and Vector declarations
    double precision, dimension(n,n):: H
    double precision, dimension(n):: f, lb, ub
    double precision, dimension(qeq,n):: Aeq
    double precision, dimension(qeq)::beq
    ! Output
    integer:: err, i
    double precision, dimension(n):: x
    double precision:: res
    double precision:: tol=1.d-3

    ! Matrix Set
    H(1,:) = [2.d0, 0.d0, 0.d0]
    H(2,:) = [0.d0, 0.d0, 0.d0]
    H(3,:) = [0.d0, 0.d0, 0.d0]

    f(:) = [-68.72d0, 0.d0, 0.d0]   

    Aeq(1,:) = [1.d0, -1.d0, -1.d0]  
    Aeq(2,:) = [3.4144d3, -3.0951d3, -2.7316d3]
    beq(:) = [0.d0, 218.d2]

    lb = [1.944d0, 0.d0, 1.944d0]
    ub = [73.16d0, 49.30d0, 33.33d0]
    x = [0.d0, 0.d0, 0.d0]
    
    ! Calculate
    call quadprog(n,q,qeq,H,f,Aeq,beq,lb,ub, x,tol,res,err)

    ! Print Results
    write(*,*)
    write(*,*)'Calculation Done'
    write(*,'(a,i2)')' Error code :', err
    write(*,'(a,es10.3)')' Eq. Residue:', qp_res_correction(n,f,res)
    write(*,*)
    write(*,*)'Solution:'
    do i=1,n
        write(*,'(8(1x,f8.3))')x(i)
    end do
    

end program Ftest