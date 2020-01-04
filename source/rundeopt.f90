program main
    use derived
    use svd, only: readin,coeff,pt_coeff,fini
    implicit none
    type(s_struct) :: struct
    integer :: i,j,k
    type(s_mse) :: bestval
    real(8) :: Vr_bestmem(ndim)
    integer :: nfeval
    external objfun
    ! valueto: "Value to Reach" (stop when ofunc < valueto)
    struct%valueto=0.0001d0
    ! Vr_minbound,Vr_maxbound: vector of lower and upper bounds of initial population.
    !     The algorithm seems to work especially well if [Vr_minbound,Vr_maxbound]
    !     covers the region where the global minimum is expected.
    !struct%Vr_minbound=-100.d0
    !struct%Vr_maxbound=100.0d0
    ! ibnd_constr: 1: use bounds as bound constraints, 0: no bound constraints.
    struct%ibnd_constr=1
!    ! itermax: maximum number of iterations (generations)
!    struct%itermax=1000000
!    ! weight: DE-stepsize weight ex [0,2]
!    struct%weight=0.85d0
!    ! crossover: crossover probability constant ex [0,1]
!    struct%crossover=1.d0
    open(unit=11,file='input.de')
    read(11,*) struct%itermax
    read(11,*) struct%weight
    read(11,*) struct%crossover
    read(11,*) np
    !read(11,*) rati1
    read(11,*)
    do i=1,ndim
        read(11,*) struct%Vr_minbound(i),struct%Vr_maxbound(i)
    end do
    close(11)
    ! irefresh: intermediate output will be produced after "irefresh"
    !           iterations. No intermediate output will be produced 
    !           if irefresh is < 1.
    struct%irefresh=10
    
!***************************************************************************
! Problem dependent but constant values. For speed reasons these values are 
! defined here. Otherwise we have to redefine them again and again in the
! cost function or pass a large amount of parameters values.
!***************************************************************************

    call readin
    call coeff
    call pt_coeff

!***************************************************************************
! Start of optimization.
!***************************************************************************
    call deopt(objfun,struct,Vr_bestmem,bestval,nfeval)

!***************************************************************************
! Sorting.
!***************************************************************************
    write(*,*) 
    write(*,'(" agamma          afactor")')
    do j=1,ndim
        write(*,'(e18.7)') Vr_bestmem(j)
    end do
    
    ! from the module 'svd'
    call fini
end program main
