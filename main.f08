program main
    use quad_prog
    
    implicit none
    ! Setup Variables
    integer:: n, qeq
    parameter( n = 3 )
    parameter( qeq = 2 )
    ! Matrix and Vector declarations
    double precision, dimension(n,n):: H
    double precision, dimension(n):: f, tg1_min, tg1_max, lb, ub
    double precision, dimension(qeq,n):: Aeq
    double precision, dimension(qeq)::beq
    ! Constant inputs
    double precision:: tg1_flow_rate, tg1_power, efficiency
    double precision:: h1, h2, h3, contrain_tol
    ! Output
    integer:: err
    double precision, dimension(n):: x
    double precision:: res

    ! -------------------
    ! Real input
    ! -------------------
    tg1_flow_rate = 123.7d0 !117.7d0 !t/h
    tg1_power = 20.71d0!21.5d0 !MW
    efficiency = 0.95d0

    ! Enthalpy
    h1 = 3.4144d3 !3.4112d3 !kJ/kg
    h2 = 3.0951d3 !3.1105d3 !kJ/kg
    h3 = 2.7316d3 !2.7425d3 !kJ/kg

    ! Limits
    tg1_min = [7.d0, 0.d0, 7.d0] !t/h
    tg1_max = [263.4d0, 177.5d0, 120.d0] !t/h
    contrain_tol = 5*1.d3/3.6d3

    ! -------------------
    ! Problem Setup
    ! -------------------

    ! Fill Matrix [H]
    H = 0.d0
    H(1,1)=2.d0
    ! Fill {f} vector
    f = [-2.d0*tg1_flow_rate*1.d3/3.6d3, 0.d0, 0.d0]
    ! Fill [Aeq]
    Aeq (1,:) = [1.d0, -1.d0, -1.d0]
    Aeq (2,:) = [h1  , -h2  , -h3  ]
    ! Fill {beq}
    beq = [0.d0, tg1_power*1.d3/efficiency]
    ! Fill {lb}
    lb = tg1_min*1.d3/3.6d3
    ! Fill {ub}
    ub = tg1_max*1.d3/3.6d3

    ! -------------------
    ! Resolve
    ! -------------------
    call quadprog(H,f,Aeq,beq,lb,ub, x,contrain_tol,res,err)

    ! -------------------
    ! Show results
    ! -------------------
    write(*,*)
    write(*,*)'Calculation Done'
    write(*,*)'Error code:', err
    write(*,*)'Residue   :', res
    write(*,*)
    write(*,*)'Solution:'
    write(*,*)'TG1 entry     ', x(1)*3.6d3/1.d3
    write(*,*)'TG1 extraction', x(2)*3.6d3/1.d3
    write(*,*)'TG1 exaustion ', x(3)*3.6d3/1.d3

end program main