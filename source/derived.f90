module derived
    implicit none
    integer,parameter :: nc=1 ! number of constraints.
    integer,parameter :: no=1 ! number of objectives (costs).
    integer,parameter :: ndim=14 ! number of parameters of the objective function.
    !integer,parameter :: np=130 ! number of population members (minimum value is 5).
    integer :: np ! number of population members (minimum value is 5).
    character(len=256) :: command,fileout
    real(8) :: rati1
    type s_mse
        real(8) :: Vr_ca(nc)
        real(8) :: Vr_oa(no)
    end type s_mse
    type s_struct
        real(8) :: weight
        real(8) :: crossover
        real(8) :: Vr_minbound(ndim),Vr_maxbound(ndim)
        integer :: ibnd_constr
        integer :: itermax
        real(8) :: valueto
        integer :: irefresh
    end type s_struct
end module derived
