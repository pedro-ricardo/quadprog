program main
    use solve_qp
    
implicit none
! Setup Variables
integer:: n, q, qeq, ierr, i
double precision, allocatable, dimension(:,:):: H, Aeq
double precision, allocatable, dimension(:):: f, x, beq
double precision, allocatable, dimension(:):: tg1_min, tg1_max, lb, ub
double precision:: Fm, h1, h2, h3, eta, Pot
double precision:: tg1_flow_rate, tg1_power
! Output Variables
double precision, allocatable, dimension(:):: lagr, work
integer, allocatable, dimension(:):: iact
double precision:: crval
integer:: nact, iter(2,1)=0, r, iw

! -------------------
! Basic input
! -------------------

! Number of variables
n = 3
! Number of constrains
q = 2+6
! Number of equality constrains
qeq = 2
! H matrix status 
ierr = 0

! -------------------
! Allocation
! -------------------

! Allocate matrices
allocate(H(n,n))
H=0.d0
allocate(Aeq(q,n))
Aeq=0.d0
! Allocate vectors
allocate(f(n), x(n), lb(n), ub(n), tg1_min(n), tg1_max(n))
f=0.d0; x=0.d0
allocate(beq(q), lagr(q)); 
beq=0.d0; lagr=0.d0
allocate(iact(q))
iact=0

! -------------------
! Real input
! -------------------
tg1_flow_rate = 117.7d0 !t/h
tg1_power = 21.5d0 !MW

! Enthalpy
h1 = 3.4112d3 !kJ/kg
h2 = 3.1105d3 !kJ/kg
h3 = 2.7425d3 !kJ/kg

! Limits
tg1_min = [7.d0, 0.d0, 7.d0] !t/h
tg1_max = [263.4d0, 177.5d0, 120.d0] !t/h

! -------------------
! Conversions
! -------------------

! Mass flux measured
Fm = tg1_flow_rate*1.d3/3.6d3
! Power
eta = 0.95d0
Pot = tg1_power*1.d3
! Limits
lb = tg1_min*1.d3/3.6d3
ub = tg1_max*1.d3/3.6d3

! Work Space
r = min(n,q)
iw = 2*n+r*(r+5)/2 + 2*q +1
allocate(work(iw)); work=0.d0

! -------------------
! Problem Setup
! -------------------

! Fill Matrix H
H(1,:)=[2.d0, 0.d0, 0.d0]
H(2,:)=[0.d0, 1.d-10, 0.d0]
H(3,:)=[0.d0, 0.d0, 1.d-10]


! Fill f vector
f = [-2.d0*Fm, 0.d0, 0.d0]

! Fill Aeq
Aeq (1,:) = [1.d0, -1.d0, -1.d0]
Aeq (2,:) = [h1, -h2, -h3]

Aeq (3,1) = 1.d0
Aeq (4,2) = 1.d0
Aeq (5,3) = 1.d0

Aeq (6,1) = -1.d0
Aeq (7,2) = -1.d0
Aeq (8,3) = -1.d0


beq = [0.d0, Pot/eta, lb, -ub ]

do i=1,n
    write(*,*)H(i,:)
end do
write(*,*)'----'

! -------------------
! Resolve
! -------------------
call qpgen2(H, -f, n, n, x, lagr, crval, transpose(Aeq), beq, n, q, qeq, iact, nact, iter, work, ierr)


write(*,*)'Done', ierr
write(*,*) x*3.6d3/1.d3

deallocate(H, Aeq, f, x, beq, lagr, work, iact, lb, ub, tg1_min, tg1_max)
end program main