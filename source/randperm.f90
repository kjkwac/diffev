module dummy
    implicit none
    save
    integer :: idum
end module dummy
!program main
!    implicit none
!    integer,parameter :: num=4
!    integer :: a1(num)
!    integer :: i,j
!    real(8) :: tmat(4,3)
!    interface
!        function randperm(num)
!            use dummy
!            implicit none
!            integer, intent(in) :: num
!            integer :: randperm(num)
!        end function randperm
!        function ranmat(n,m)
!            use dummy
!            implicit none
!            integer,intent(in) :: n,m
!            real(8) :: ranmat(n,m)
!        end function ranmat
!    end interface
!    do i=1,10
!    a1=randperm(4)
!    write(*,*) a1
!    end do
!    tmat=ranmat(3,4)
!    do i=1,4
!        write(*,*) (tmat(i,j),j=1,3)
!    end do
!    write(*,*) 
!    tmat=ranmat(3,4)
!    do i=1,4
!        write(*,*) (tmat(i,j),j=1,3)
!    end do
!    write(*,*) 
!    tmat=ranmat(3,4)
!    do i=1,4
!        write(*,*) (tmat(i,j),j=1,3)
!    end do
!end program main


subroutine random_mat(ranmat,n,m)
    use dummy
    implicit none
    integer,intent(in) :: n,m
    real(8) :: ranmat(n,m)
    integer :: i,j
    real(8),external :: ran2
    do i=1,n
        do j=1,m
            ranmat(i,j)=ran2(idum)
        end do
    end do
end subroutine random_mat


function randperm(num)
    use dummy
    implicit none
    integer,intent(in) :: num
    integer :: inum,i,j,k
    integer :: randperm(num)
    real(8) :: mrand(num)
    real(8) :: ran2
    external ran2
    do i=1,num
        mrand(i)=ran2(idum)
    end do
    do i=1,num
        inum=1
        do j=1,num
            if (mrand(i) > mrand(j)) then
                inum=inum+1
            end if
        end do
        do k=1,i-1
            if (mrand(i).LE.mrand(k).AND.mrand(i).GE.mrand(k)) then
                inum=inum+1
            end if
        end do
        randperm(i)=inum
    end do
end function randperm


FUNCTION ran2(idum)
    IMPLICIT NONE
    INTEGER, PARAMETER :: K4B=selected_int_kind(9)
    INTEGER(K4B), INTENT(INOUT) :: idum
    REAL(8) :: ran2
    INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
    REAL(8), SAVE :: am
    INTEGER(K4B), SAVE :: ix=-1,iy=-1,k
    if (idum <= 0 .or. iy < 0) then
        am=nearest(1.0,-1.0)/IM
        iy=ior(ieor(888889999,abs(idum)),1)
        ix=ieor(777755555,abs(idum))
        idum=abs(idum)+1
    end if
    ix=ieor(ix,ishft(ix,13))
    ix=ieor(ix,ishft(ix,-17))
    ix=ieor(ix,ishft(ix,5))
    k=iy/IQ
    iy=IA*(iy-k*IQ)-IR*k
    if (iy < 0) iy=iy+IM
    ran2=am*ior(iand(IM,ieor(ix,iy)),1)
END FUNCTION ran2

