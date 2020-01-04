!  ifort -o singular_partition singular_partition.f90  -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include  ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
! 2019-01-09: singular_ratition.f90: modified from singular_step.f90 to ratition a frequency shift.
module svd
    implicit none
    integer,parameter :: f_in=10,f_sh=11,f_tr=12,f_co=13,f_pt=14,f_ou=15
    character(len=256) :: shotfile,pt_file
    integer :: nshot,maxatom
    integer,allocatable :: natom(:)
    character(len=4),allocatable :: atname(:)
    real(8),allocatable :: shift(:),shift0(:),coord(:,:,:),wat(:,:,:)
    real(8),allocatable :: rshift(:,:)

    integer :: nanma,nawat
    integer :: nwater
    real(8),allocatable :: nmacrd(:,:)
    integer :: nsite
    integer,parameter :: nwsite=2,norder=14
    integer,allocatable :: site(:)
    real(8) :: rati(norder)
    real(8),allocatable :: rinver(:,:,:,:),rmatrix(:,:),pmatrix(:,:)
    real(8),allocatable :: avector(:),params(:,:,:)

    real(8) :: rmse,pt_rmse
    integer :: pt_nshot,pt_maxatom
    integer,allocatable :: pt_natom(:)
    character(len=4),allocatable :: pt_atname(:)
    real(8),allocatable :: pt_shift(:),pt_coord(:,:,:),pt_rinver(:,:,:,:)

contains

    subroutine init
    implicit none
    integer :: maxat
    maxat=max(maxatom,pt_maxatom)
    allocate(nmacrd(3,nanma))
    allocate(site(nsite))
    allocate(natom(nshot))
    allocate(shift(nshot))
    allocate(shift0(nshot))
    allocate(atname(maxatom))
    allocate(coord(3,maxatom,nshot))
    allocate(wat(3,3,maxatom))
    allocate(rinver(nsite,nwsite,norder,nshot))
    !allocate(rmatrix(nshot,nsite*nwsite*norder))
    !allocate(pmatrix(nsite*nwsite*norder,nshot))
    !allocate(avector(nsite*nwsite*norder))
    allocate(rmatrix(nshot,nsite*nwsite))
    allocate(pmatrix(nsite*nwsite,nshot))
    allocate(avector(nsite*nwsite))
    allocate(params(nsite,nwsite,norder))
    allocate(rshift(norder,nshot))
    !allocate(rati(norder))
    allocate(pt_natom(pt_nshot))
    allocate(pt_atname(maxat))
    allocate(pt_shift(pt_nshot))
    allocate(pt_coord(3,maxat,pt_nshot))
    allocate(pt_rinver(nsite,nwsite,norder,pt_nshot))
    open(unit=f_ou,file='output')
    end subroutine init

    subroutine fini
    implicit none
    deallocate(pt_rinver)
    deallocate(pt_coord)
    deallocate(pt_shift)
    deallocate(pt_atname)
    deallocate(pt_natom)
    !deallocate(rati)
    deallocate(rshift)
    deallocate(params)
    deallocate(avector)
    deallocate(pmatrix)
    deallocate(rmatrix)
    deallocate(rinver)
    deallocate(wat)
    deallocate(coord)
    deallocate(atname)
    deallocate(shift0)
    deallocate(shift)
    deallocate(natom)
    deallocate(site)
    deallocate(nmacrd)
    close(f_ou)
    end subroutine fini

    subroutine readin
    implicit none
    integer :: ishot,ia,num
    character(len=4) :: col

    !site=(/1,5,6,7,8,9/)

    open(unit=f_in,file="svd.in")
    read(f_in,'(a)') shotfile
    read(f_in,'(a)') pt_file
    read(f_in,*) nanma,nawat
    read(f_in,*) nsite

    open(unit=f_sh,file=shotfile)
    read(f_sh,*) nshot,maxatom
    open(unit=f_pt,file=pt_file)
    read(f_pt,*) pt_nshot,pt_maxatom
    call init
    read(f_in,*) (site(ia),ia=1,nsite)
    close(f_in)
    do ishot=1,nshot
        read(f_sh,*) natom(ishot),shift0(ishot)
        do ia=1,natom(ishot)
            read(f_sh,*) atname(ia),coord(1:3,ia,ishot)
        end do
    end do
    do ishot=1,pt_nshot
        read(f_pt,*) pt_natom(ishot),pt_shift(ishot)
        do ia=1,pt_natom(ishot)
            read(f_pt,*) pt_atname(ia),pt_coord(1:3,ia,ishot)
        end do
    end do
    close(f_sh)
    end subroutine readin

    subroutine coeff
    implicit none
    integer :: ishot,ia,iwater,ibegin,ins,ja,iwt,iod
    integer :: nwater
    real(8) :: a(3),b(3),rij(3),rijab,rijinv

    rinver=0.d0
    do ishot=1,nshot
        do ia=1,nanma
            nmacrd(1:3,ia)=coord(1:3,ia,ishot)
        end do
        nwater=(natom(ishot)-nanma)/nawat
        do iwater=1,nwater
            ibegin=nanma+(iwater-1)*3
            do ja=1,3
                wat(1:3,ja,iwater)=coord(1:3,ibegin+ja,ishot)
            end do
        end do

        do ins=1,nsite
            do iwater=1,nwater
                do ja=1,3
                    a=nmacrd(:,site(ins))
                    b=wat(:,ja,iwater)
                    rij=a-b
                    rijab=sqrt(dot_product(rij,rij))
                    rijinv=1.d0/rijab
                    if (ja.EQ.1) then
                        if (atname(nanma+(iwater-1)*3+ja).NE.'O') stop 'Should be O atom!'
                        iwt=1  ! Water O atom.
                    else
                        if (atname(nanma+(iwater-1)*3+ja).NE.'H') stop 'Should be H atom!'
                        iwt=2  ! Water H atom.
                    end if
                    do iod=1,norder
                        rinver(ins,iwt,iod,ishot)=rinver(ins,iwt,iod,ishot)+rijinv**iod
                    end do
                end do 
            end do
        end do
    end do
    end subroutine coeff

    subroutine rati_svd
    implicit none
    integer :: ishot,i,iod,ip,iwt,ins

    do ishot=1,nshot
        do i=1,norder
            rshift(i,ishot)=rati(i)*shift0(ishot)
        end do
    end do

    do iod=1,norder
        do ishot=1,nshot
            ip=0
            do iwt=1,nwsite
                do ins=1,nsite
                    ip=ip+1
                    rmatrix(ishot,ip)=rinver(ins,iwt,iod,ishot)
                end do
            end do
        end do
        call decomp(iod) ! params(ins,iwt,iod) is obtained.
    end do

    open(unit=f_co,file='params.dat')
    ip=0
    do iod=1,norder
        do iwt=1,nwsite
            do ins=1,nsite
                ip=ip+1
                write(f_co,'(i4,3i4,1x,e22.15)') ip,ins,iwt,iod,params(ins,iwt,iod)
            end do
        end do
    end do
    close(f_co)
    end subroutine rati_svd

    subroutine pt_coeff
    implicit none
    integer :: ishot,ia,nwater,iwater,ibegin,ja,ins,iod,iwt
    real(8) :: a(3),b(3),rij(3),rijab,rijinv

    pt_rinver=0.d0
    do ishot=1,pt_nshot
        do ia=1,nanma
            nmacrd(1:3,ia)=pt_coord(1:3,ia,ishot)
        end do
        nwater=(pt_natom(ishot)-nanma)/nawat
        do iwater=1,nwater
            ibegin=nanma+(iwater-1)*3
            do ja=1,3
                wat(1:3,ja,iwater)=pt_coord(1:3,ibegin+ja,ishot)
            end do
        end do
        do ins=1,nsite
            do iwater=1,nwater
                do ja=1,3
                    a=nmacrd(:,site(ins))
                    b=wat(:,ja,iwater)
                    rij=a-b
                    rijab=sqrt(dot_product(rij,rij))
                    rijinv=1.d0/rijab
! fortest
    !if (ishot.EQ.1) write(f_ou,*) 'rijinv',ins,iwater,ja,rijinv
    !if (ishot.EQ.1) write(f_ou,*) 'rijinv',ins,iwater,ja,nwater
                    if (ja.EQ.1) then
                        if (pt_atname(nanma+(iwater-1)*3+ja).NE.'O') stop 'Should be O atom!'
                        iwt=1  ! Water O atom.
                    else
                        if (pt_atname(nanma+(iwater-1)*3+ja).NE.'H') stop 'Should be H atom!'
                        iwt=2  ! Water H atom.
                    end if
                    do iod=1,norder
                        pt_rinver(ins,iwt,iod,ishot)=pt_rinver(ins,iwt,iod,ishot)+rijinv**iod
                    end do
                end do
            end do
        end do
    end do
    end subroutine pt_coeff

    subroutine decomp(iod)
    USE lapack95, ONLY: GESVD
    implicit none
    integer,intent(in) :: iod
    integer :: m,n,info,i,ip,iwt,ins
    integer :: j
    real(8),allocatable :: a(:,:),aa(:,:),s(:),u(:,:),v(:,:),vt(:,:),ww(:)
    real(8),allocatable :: smat(:,:),s1mat(:,:)
    real(8),allocatable :: rtest(:,:)

    m=nshot 
    n=nwsite*nsite
    allocate(a(m,n),aa(m,n),s(min(m,n)),u(m,m),v(n,n),vt(n,n),ww(min(m,n)-1))
    allocate(smat(m,n),s1mat(n,m))
    allocate(rtest(n,n))

    a=rmatrix
    aa=a
    call GESVD(aa,s,u=u,vt=vt,ww=ww,job='U',info=info)

    smat=0.d0
    s1mat=0.d0
    do i=1,n
        smat(i,i)=s(i)
        if (abs(s(i)).GT.1.d-15) s1mat(i,i)=1.d0/s(i)
    end do
    pmatrix=matmul(transpose(vt),matmul(s1mat,transpose(u)))
! fortest
   !write(f_ou,*) 'shape(rmatrix)',shape(rmatrix)
   !write(f_ou,*) 'shape(pmatrix)',shape(pmatrix)
   rtest=matmul(pmatrix,rmatrix)
   !do i=1,n
   !    write(f_ou,*) 'rtest',i,rtest(i,i),rtest(i,1)
   !end do
   avector=matmul(pmatrix,rshift(iod,:))

    ip=0
    do iwt=1,nwsite
        do ins=1,nsite
            ip=ip+1
            params(ins,iwt,iod)=avector(ip)
        end do
    end do

    deallocate(rtest)
    deallocate(smat,s1mat)
    deallocate(a,aa,s,u,v,vt,ww)
    end subroutine decomp

    subroutine training
!   Compare with the ab initio frequency shifts and calculated ones for the training set.
    implicit none
    integer :: ishot,iod,iwt,ins
    real(8) :: dfrqcalc,fac

    open(unit=f_tr,file='training.dat')
    fac=0.d0
    do ishot=1,nshot
        dfrqcalc=0.d0
        do iod=1,norder
            do iwt=1,nwsite
                do ins=1,nsite
                    dfrqcalc=dfrqcalc+params(ins,iwt,iod)*rinver(ins,iwt,iod,ishot)
                end do
            end do
        end do
        write(f_tr,*) ishot,shift0(ishot),dfrqcalc
        if (shift(ishot).GT.300.d0) write(f_ou,*) 'large:'
        fac=fac+(shift0(ishot)-dfrqcalc)**2
    end do
    rmse=sqrt(fac/dble(nshot))

    fac=0.d0
    do ishot=1,pt_nshot
        dfrqcalc=0.d0
        do iod=1,norder
            do iwt=1,nwsite
                do ins=1,nsite
                    dfrqcalc=dfrqcalc+params(ins,iwt,iod)*pt_rinver(ins,iwt,iod,ishot)
                end do
            end do
        end do
        fac=fac+(pt_shift(ishot)-dfrqcalc)**2
    end do
    pt_rmse=sqrt(fac/dble(pt_nshot))
    !write(*,*) pt_rmse
    close(f_tr)
    end subroutine training
end module svd

!program main
!    use svd
!    implicit none
!    call readin
!    call pt_coeff
!
!    rati=0.5d0
!    call coeff
!    call training
!
!    call fini
!end program main
