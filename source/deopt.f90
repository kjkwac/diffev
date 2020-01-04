subroutine deopt(objfun,struct,Vr_bestmem,bestval,nfeval)
    use derived
    use dummy
    implicit none
    type(s_mse) :: objfun,bestval
    type(s_struct) :: struct
    real(8) :: Vr_bestmem(ndim)
    integer :: nfeval
    external objfun

    ! local variables.
    integer :: ibnd_constr,itermax,irefresh
    real(8) :: weight,crossover,Vr_minbound(ndim),Vr_maxbound(ndim),valueto
    real(8) :: M_pop(np,ndim),tmat1(1,ndim),tmat2(np,1),tmat(np,ndim),tran
    real(8) :: ranm_np(np,1),ranm_ndim(1,ndim)
    real(8),external :: ran2
    real(8) :: f1,fnp(np,1),np1(np,1)
    real(8) :: M_popold(np,ndim)
    real(8) :: M_origin(np,ndim)
    real(8) :: Vr_bestmemit(ndim)
    integer :: k,iter,j,n
    integer :: best_index
    type(s_mse) :: val(np),bestvalit,tempval

    ! population matrix 1 - 5
    real(8) :: M_pm1(np,ndim),M_pm2(np,ndim),M_pm3(np,ndim),M_pm4(np,ndim),M_pm5(np,ndim)     
    real(8) :: M_bm(np,ndim)    ! Vr_bestmember matrix.
    real(8) :: M_ui(np,ndim)    ! intermediate population of perturbed vectors.
    real(8) :: M_mui(np,ndim)   ! mask for intermediate population.
    real(8) :: M_mpo(np,ndim)   ! mask for old popoulation.
    integer :: Vr_rot(np)     ! rotating index array (size np).
    integer :: Vr_rotd(ndim)  ! rotating index array (size ndim).
    integer :: Vr_rt(np)        ! another rotating index array.
    integer :: Vr_rtd(ndim)     ! rotating index array for exponential crossover.
    integer :: Vr_a1(np)        ! index array
    integer :: Vr_a2(np)        ! index array
    integer :: Vr_a3(np)        ! index array
    integer :: Vr_a4(np)        ! index array
    integer :: Vr_a5(np)        ! index array
    integer :: Vr_ind(4)

    interface
        function randperm(num)
            use dummy
            implicit none
            integer, intent(in) :: num
            integer :: randperm(num)
        end function randperm
        integer function left_win(S_x,S_y)
            use derived
            implicit none
            type(s_mse) :: S_x,S_y
        end function left_win
    end interface

    idum=-1   ! seed for the random number generator.
    !-----This is just for notational convenience and to keep the code uncluttered.----
    weight=struct%weight
    crossover=struct%crossover
    Vr_minbound=struct%Vr_minbound
    Vr_maxbound=struct%Vr_maxbound
    ibnd_constr=struct%ibnd_constr
    itermax=struct%itermax
    valueto=struct%valueto
    irefresh=struct%irefresh
    
    !-----Check input variables-----------------------------------------
    if ((crossover.LT.0.d0).OR.(crossover.GT.1.d0)) then
        crossover=0.5d0
        write(*,*) "'crossover' should be from interval [0,1]; set to default value 0.5"
    end if
    if (itermax.LE.0) then
        itermax=200
        write(*,*) "'itermax' should be > 0; set to default value 200"
    end if

    !----M_pop is a matrix of size I_NPx(I_D+1). It will be initialized------
    !----with random values between the min and max values of the-------------
    !----parameters-----------------------------------------------------------
    do k=1,np
       call random_mat(tmat1,1,ndim)
       M_pop(k,:)=Vr_minbound+tmat1(1,:)*(Vr_maxbound-Vr_minbound)
    end do 

    nfeval=0                            ! number of function evaluation.

    !------Evaluate the best member after initialization----------------------

    best_index=1                        ! start with first popyulation member.
    val(1)=objfun(M_pop(best_index,:),struct)
    
    bestval=val(1)                      ! best objective function value so far.
    nfeval=nfeval+1
    do k=2,np
        val(k)=objfun(M_pop(k,:),struct)
        nfeval=nfeval+1
        if (left_win(val(k),bestval).EQ.1) then
            best_index=k
            bestval=val(k)
        end if
    end do  
    Vr_bestmemit=M_pop(best_index,:)    ! best member of current iteration.
    bestvalit=bestval                   ! best value of current iteration.

    Vr_bestmem=Vr_bestmemit             ! best member ever.

!------DE-Minimization---------------------------------------------
!------FM_popold is the population which has to compete. It is--------
!------static through one iteration. FM_pop is the newly--------------
!------emerging population.----------------------------------------

    Vr_rot=(/(j,j=0,np-1)/)
    Vr_rotd=(/(j,j=0,ndim-1)/)
    np1=1.d0

    iter=1
    do
        if (iter.GT.itermax.OR.bestval%Vr_oa(1).LT.valueto) exit
        M_popold=M_pop

        Vr_ind=randperm(4)              ! index pointer array.

        Vr_a1=randperm(np)                  ! shuffle locations of vectors.
        Vr_rt=mod(Vr_rot+Vr_ind(1),np)      ! rotate indices by ind(1) positions.
        Vr_a2=Vr_a1(Vr_rt+1)                ! rotate vector locations.
        Vr_rt=mod(Vr_rot+Vr_ind(2),np)
        Vr_a3=Vr_a2(Vr_rt+1)
        Vr_rt=mod(Vr_rot+Vr_ind(3),np)
        Vr_a4=Vr_a3(Vr_rt+1)
        Vr_rt=mod(Vr_rot+Vr_ind(4),np)
        Vr_a5=Vr_a4(Vr_rt+1)
        
        M_pm1=M_popold(Vr_a1,:)         ! suffled population 1. 
        M_pm2=M_popold(Vr_a2,:)         ! suffled population 2. 
        M_pm3=M_popold(Vr_a3,:)         ! suffled population 3. 
        !M_pm4=M_popold(Vr_a4,:)         ! suffled population 4. 
        !M_pm5=M_popold(Vr_a5,:)         ! suffled population 5. 

        do k=1,np                       ! population filled with the best member
            M_bm(k,:)=Vr_bestmemit      ! of the last iteration.
        end do

        ! M_mui: all random numbers < crossover are 1, 0 otherwise.
        ! M_mpo: inverse mask to M_mui
        call make_mask(M_mui,M_mpo,np,ndim,crossover)

        M_ui=M_pm3+weight*(M_pm1-M_pm2)      ! differential variation
        M_ui=M_popold*M_mpo+M_ui*M_mui       ! crossover
        M_origin=M_pm3

!-----Optional parent+child selection-----------------------------------------
  
!-----Select which vectors are allowed to enter the new population------------
        do k=1,np

            !=====Only use this if boundary constraints are needed==================
            if (ibnd_constr.EQ.1) then
                do j=1,ndim !----boundary constraints via bounce back-------
                    if (M_ui(k,j).GT.Vr_maxbound(j)) then
                        tran=ran2(idum)
                        M_ui(k,j)=Vr_maxbound(j)+tran*(M_origin(k,j)-Vr_maxbound(j))
                    end if
                    if (M_ui(k,j).LT.Vr_minbound(j)) then
                        tran=ran2(idum)
                        M_ui(k,j)=Vr_minbound(j)+tran*(M_origin(k,j)-Vr_minbound(j))
                    end if
                end do
            end if
            !=====End boundary constraints==========================================

            tempval=objfun(M_ui(k,:),struct)    ! check cost of competitor.
            nfeval=nfeval+1
            if (left_win(tempval,val(k)).EQ.1) then
                M_pop(k,:)=M_ui(k,:)            ! replace old vector with new one (for new iteration).
                val(k)=tempval                  ! save value in "cost array".

                !----we update S_bestval only in case of success to save time-----------
                if (left_win(tempval,bestval).EQ.1) then
                    bestval=tempval             ! new best value
                    Vr_bestmem=M_ui(k,:)        ! new best parameter vector ever
                end if
            end if
        end do ! k=1,np                

        Vr_bestmemit=Vr_bestmem    ! freeze the best member of this iteration for the coming
                                   ! iteration. This is needed for some of the strategies.

        !----Output section----------------------------------------------------------

        if (irefresh.GT.0) then
            if ((mod(iter,irefresh).EQ.0).OR.iter.EQ.1) then
                write(*, '("Iteration: ",i7,",   Best:  ",f12.6,",   F_weight: ",f12.6,",   F_CR: ",f12.6,",   NP: ",i7)') &
                    iter,bestval%Vr_oa(1),weight,crossover,np
            end if
        end if
    
        iter=iter+1
    end do

end subroutine deopt


subroutine make_mask(M_mui,M_mpo,np,ndim,cr)
    use dummy
    implicit none
    real(8) :: M_mui(np,ndim),M_mpo(np,ndim),cr
    integer :: np,ndim
    real(8) :: r,ran2
    integer :: i,j
    external ran2

    do i=1,np
        do j=1,ndim
            r=ran2(idum)
            if (r.LT.cr) then
                M_mui(i,j)=1.d0
                M_mpo(i,j)=0.d0
            else
                M_mui(i,j)=0.d0
                M_mpo(i,j)=1.d0
            end if
        end do
    end do
end subroutine make_mask
