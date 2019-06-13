program Ftest
    use quad_prog
    
    implicit none
    ! Setup Variables
    integer:: n, qeq
    parameter( n = 8 )
    parameter( qeq = 5 )
    ! Matrix and Vector declarations
    double precision, dimension(n*n):: H
    double precision, dimension(n):: f, tg1_min, tg1_max, lb, ub
    double precision, dimension(qeq*n):: Aeq
    double precision, dimension(qeq)::beq
    ! Constant inputs
    double precision, parameter:: rho_water = 997.d0
    double precision:: contrain_tol
    double precision:: boiler_steam_flow_rate, tg2_condens_flow_rate
    double precision:: tg1_flow_rate, tg1_power, tg1_efficiency
    double precision:: tg2_flow_rate, tg2_power, tg2_efficiency
    double precision:: tg1_h1, tg1_h2, tg1_h3
    double precision:: tg2_h1, tg2_h2, tg2_h3
    double precision, dimension(n):: var_min, var_max
    ! Output
    integer:: err, i
    double precision, dimension(n):: x
    double precision:: res

    ! -------------------
    ! Real input
    ! -------------------
    boiler_steam_flow_rate = 294.8 !t/h
    
    tg1_flow_rate = 129.4d0 !t/h
    tg1_power = 23.05d0 !MW
    tg1_efficiency = 0.95d0
    
    tg2_condens_flow_rate = 72.06d0 !m³/h
    tg2_flow_rate = rho_water*tg2_condens_flow_rate/1.d3 ! t/h
    tg2_power = 35.28d0 !MW
    tg2_efficiency = 0.95d0

    ! Enthalpy
    tg1_h1 = 3.4095d3
    tg1_h2 = 3.1065d3
    tg1_h3 = 2.7244d3
    
    tg2_h1 = 3.4212d3
    tg2_h2 = 2.7235d3
    tg2_h3 = 2.3760d3

    ! Limits  !  Fin1     Fext1   Fexh1    Fin2   Fext2   Fexh2   Ftot      Rprv
    var_min = [   7.d0,    0.d0,   7.d0,  40.d0,   0.d0,   0.d0,  0.d0  ,  0.d0 ]
    var_max = [263.4d0, 177.5d0, 120.d0, 160.d0, 100.d0, 100.d0,  400.d0, 80.d0 ]
    
    contrain_tol = 1.d-5!5*1.d3/3.6d3

    ! -------------------
    ! Problem Setup
    ! -------------------

    ! Fill Matrix [H]
    H = [2, 0, 0, 0, 0, 0, 0, 0,&
         0, 0, 0, 0, 0, 0, 0, 0,&
         0, 0, 0, 0, 0, 0, 0, 0,&
         0, 0, 0, 0, 0, 0, 0, 0,&
         0, 0, 0, 0, 0, 0, 0, 0,&
         0, 0, 0, 0, 0, 2, 0, 0,&
         0, 0, 0, 0, 0, 0, 2, 0,&
         0, 0, 0, 0, 0, 0, 0, 0 ]
    
    ! Fill {f} vector
    f = 0.d0
    f(1) = -2.d0*tg1_flow_rate*1.d3/3.6d3
    f(6) = -2.d0*tg2_flow_rate*1.d3/3.6d3
    f(7) = -2.d0*boiler_steam_flow_rate*1.d3/3.6d3
    
    ! Fill [Aeq]
    Aeq = [1d0, -1d0, -1d0, 0d0, 0d0, 0d0, 0d0, 0d0,&
           tg1_h1, -tg1_h2, -tg1_h3, 0d0, 0d0, 0d0, 0d0, 0d0,&
           0d0, 0d0, 0d0, 1d0, -1d0, -1d0, 0d0, 0d0,&
           0d0, 0d0, 0d0, tg2_h1, -tg2_h2, -tg2_h3, 0d0, 0d0,&
           -1d0, 0d0, 0d0, -1d0, 0d0, 0d0, 1d0, -1d0 ]
           
    ! Fill {beq}
    beq = 0.d0
    beq(2) = 1.d3*tg1_power/tg1_efficiency
    beq(4) = 1.d3*tg2_power/tg2_efficiency
    
    ! Fill {lb}
    lb = var_min*1.d3/3.6d3
    ! Fill {ub}
    ub = var_max*1.d3/3.6d3
    
    write(*,*)'H'
    do i=1,n
    write(*,'(8(1x,es10.3))')H(1+(i-1)*n:(i)*n)
    end do
    write(*,*)'f'
    write(*,'(8(1x,es10.3))')f
    write(*,*)'Aeq'
    do i=1,qeq
        write(*,'(8(1x,es10.3))')Aeq(1+(i-1)*n:(i)*n)
    end do
    write(*,*)'beq'
    write(*,'(8(1x,es10.3))')beq
    write(*,*)'-----'
    ! -------------------
    ! Resolve
    ! -------------------
    call quadprog(n,qeq,H,f,Aeq,beq,lb,ub, x,contrain_tol,res,err)

    ! -------------------
    ! Show results
    ! -------------------
    write(*,*)
    write(*,*)'Calculation Done'
    write(*,*)'Error code:', err
    write(*,*)'Residue   :', res
    write(*,*)'Residue_C :', qp_res_correction(n,f,res)
    write(*,*)
    write(*,*)'Solution:'
    write(*,*)'TG1 entry     ', x(1)*3.6d3/1.d3
    write(*,*)'TG1 extraction', x(2)*3.6d3/1.d3
    write(*,*)'TG1 exaustion ', x(3)*3.6d3/1.d3
    write(*,*)'TG2 entry     ', x(4)*3.6d3/1.d3
    write(*,*)'TG2 extraction', x(5)*3.6d3/1.d3
    write(*,*)'TG2 exaustion ', x(6)*3.6d3/1.d3
    write(*,*)'Total         ', x(7)*3.6d3/1.d3
    write(*,*)'Valvula       ', x(8)*3.6d3/1.d3

end program Ftest