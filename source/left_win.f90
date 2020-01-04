!================================================================================
! Function:         left_win(S_x,S_y)
! Description:      left_win(S_x,S_y) takes structures S_x and S_y as an argument.
!                   The function returns 1 if the left structure of the input structures,
!                   i.e. S_x, wins. If the right structure, S_y, wins, the result is 0.
! Parameters:       
!                   S_x.Vr_ca    (I)    Constraint array containing the constraint violation values.
!                                       If the value is 0 the constraint is met. If it is > 0 it is
!                                       still violated.
!                   S_x.Vr_oa    (I)    Objective array containing cost values which are supposed to be
!                                       minimized.
! Return value:     left_win     (O)    If S_x wins over S_y then 1 else 0.
!================================================================================
integer function left_win(S_x,S_y)
    use derived
    implicit none
    integer :: z,k
    type(s_mse) :: S_x,S_y

    z=1 ! start with z=1.
    
    !----deal with the constraints first. If constraints are not met------
    !----S_x can't win.---------------------------------------------------
    if (nc.GT.0) then
        do k=1,nc
            if (S_x%Vr_ca(k).GT.0) then ! if constaint is not yet met
                if (S_x%Vr_ca(k).GT.S_y%Vr_ca(k)) then ! if just one constraint of S_x is not improved
                    z=0
                end if
            end if
        end do
    end if

    if (no.GT.0) then
        do k=1,no
            if (S_x%Vr_oa(k).GT.S_y%Vr_oa(k)) then ! if just one objective of S_x is less
                z=0
            end if
        end do
    end if
    
    left_win=z
end function left_win
