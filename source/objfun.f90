!=================================================================================
! Description:      Implements the cost function to be minimized.
! Parameters:       Vr_temp       (I)    Parameter vector
!                   struct        (I)    Contains a variety of parameters.
!                                       For details see Rundeopt.m
! Return value:     objfun.I_nc   (O)    Number of constraints
!                   objfun.FVr_ca (O)    Constraint values. 0 means the constraints
!                                       are met. Values > 0 measure the distance
!                                       to a particular constraint.
!                   objfun.I_no   (O)    Number of objectives.
!                   objfun.FVr_oa (O)    Objective function values.
!=================================================================================
type(s_mse) function objfun(Vr_temp,struct)
    use derived
    use svd, only: rati,rati_svd,training,pt_rmse
    implicit none
    real(8) :: Vr_temp(ndim)
    type(s_struct) :: struct
    real(8) :: rmse,Vr_ca,rtot
    integer :: i

    rtot=sum(Vr_temp(1:ndim))
    do i=1,ndim
        rati(i)=Vr_temp(i)/rtot
    end do
    call rati_svd
    call training
    
    Vr_ca=0.d0 ! No constraint.
    
    objfun%Vr_ca=Vr_ca  ! constraint array
    objfun%Vr_oa(1)=pt_rmse
end function objfun
