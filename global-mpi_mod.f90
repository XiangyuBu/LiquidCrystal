MODULE global_parameters
IMPLICIT NONE

DOUBLE PRECISION, PARAMETER :: pi=3.1415926535897932384D0
DOUBLE PRECISION, PARAMETER :: eps = 1d-10

type node
	DOUBLE PRECISION :: x
	DOUBLE PRECISION :: y
	DOUBLE PRECISION :: z
end type
type(node),allocatable :: polymer(:,:)     !polymer in the ball 
type(node),allocatable :: bond_vector(:)
type(node),allocatable :: azo(:,:)         !azo-polymer in the substrate

Integer :: Nm, N_azo, Nm_chain, N_chain  
Integer :: Nz, Nr, i_interface
Integer :: N_theta, N_phi
Integer :: Movestep    !each MC attempt move
Integer :: i_azo, isomer
Integer :: num
Double Precision :: dr, Lr, Lz, dz, Lz_2,dtheta,dphi
Double Precision :: r_sphere_2,r_sphere, roL, rho_0 
Double Precision :: move_max, rotate, rotate_s,loa_polymer,csoL

Integer, DIMENSION(:,:), ALLOCATABLE :: ir, iz, ir_azo, iz_azo, itheta_azo, iphi_azo
Integer, DIMENSION(:), ALLOCATABLE :: itheta, iphi
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xx, zz
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: v_tide

DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: w, w_new

Double PRECISION :: r_dr,r_dz
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: r_a, z_i, phi_rtot, phi_r, rr_r, zz_r
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: trial_move_rate

DOUBLE PRECISION :: nu, tau, epsilon, epsilon_azo, Lbox
DOUBLE PRECISION ::  Loa, deltaS,lambda
Integer :: n_iter, Max_iter, N_pre, Npre, Nmove, moves, NMCs,NMCstot, MCS, ncount
Integer :: seed
logical ::keepruning   !if keep MC simulating in the w field ,keepruning=.true.
!while the Npre will times 10 to keep a longer MC simulation time ,you can times
!a much longer time if you want


END MODULE global_parameters