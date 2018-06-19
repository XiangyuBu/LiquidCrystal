subroutine initialize()
USE global_parameters
USE utility_routines 
IMPLICIT NONE

Integer*8 npts,mynpts,change_2,change_3
Integer :: flag_c, index_con, n_graft_point, x_flag, x_stat
Integer :: i, j, k, jp,ip
integer :: n, m, l, length
integer :: change, change_1  ! change is the flag of MC, 1 for accepte the move.
integer :: ddt(8) 

DOUBLE PRECISION :: axis(3)
DOUBLE PRECISION :: r, cos_t, sin_t, phi
DOUBLE PRECISION :: yy, yyp
DOUBLE PRECISION :: unew(3), uold(3) 
DOUBLE PRECISION :: alpha, beta, angle, dotp
DOUBLE PRECISION :: ratio_sigma, a_plane,loa_azo 
DOUBLE PRECISION :: azo_position, phi_azo
DOUBLE PRECISION, PARAMETER :: TOL = 1.0D-5	
DOUBLE PRECISION :: r_radius,r_radius_1, rr, temp, x_r, y_r, z_r


type(node) :: graft_point(0:401)  
character*7 res,resres,res_s
character res0,res1,res2
logical alive,check

call date_and_time(values=ddt)
seed=(ddt(8)-500)*54321 + 11223344
!open(unit=15,file='log.txt')
open(unit=10,file='input.txt')

read(10,*) azo_position    ! azo position / L from the free end
read(10,*) ratio_sigma     ! grafting density ratio between substrate and sphere sigma_p/sigma_s
read(10,*) nu              ! interaction parameter
read(10,*) Loa             ! L/a
read(10,*) roL             ! radius of sphere/L
read(10,*) N_azo           ! number of grafted chains
read(10,*) lambda          ! mixing parameter
read(10,*) Nm              ! number of bonds on azo polymer
read(10,*) index_con       ! number of conformations  
read(10,*) Nr              ! number of points in r direction
read(10,*) Nz              ! number of points in z direction
read(10,*) Npre            ! number of prerotated   
read(10,*) Nmove           ! number of move in a MC step 
read(10,*) Max_iter        ! Max of iterations
read(10,*) num             ! num==1 ,pivot move with no small pivot move 
read(10,*) rotate          ! the limitation of the random move angle 
read(10,*) rotate_s        ! the limitation of the rotate of sphere angle
read(10,*) move_max        ! max of polymer trial move
read(10,*) loa_polymer           ! the real loa      !!!!!It will be change!!!!!
read(10,*) loa_azo         ! loa_azo
read(10,*) N_chain         ! number of polymers
read(10,*) Nm_chain        ! number of bonds on liquid crystal polymers
read(10,*) N_theta         ! lattice of theta in [0,pi/2] by discretlizing delta_z in [0, 1]
read(10,*) N_phi           ! lattice of phi in [-pi,pi] by discretlizing delta_x in [-1, 1]
read(10,*) isomer          !  0 for trans (straight line) ; 1 for sis (folding)  
read(10,*) length          !test initial MC move
close(10)

allocate(polymer(1:N_chain,0:Nm_chain))
allocate(w(1:Nr+1,1:Nz,0:N_theta,0:N_phi))
allocate(ir(1:N_chain,0:Nm_chain), iz(1:N_chain,0:Nm_chain), bond_vector(1:N_chain))
allocate(itheta(1:N_chain), iphi(1:N_chain))
allocate(r_a(1:Nr+1), rr_r(1:Nr+1), zz_r(1:Nz))
allocate(xx(1:N_phi),zz(1:N_theta))
allocate(v_tide(1:N_theta,1:N_phi,1:N_theta,1:N_phi))
allocate(trial_move_rate(1:4))

! trial move rates
trial_move_rate(1) = 0.4d0 ! rate of azo move 
trial_move_rate(2) = 0.08d0 + trial_move_rate(1) ! rate of sphere rotate 
trial_move_rate(3) = 0.3d0 + trial_move_rate(2) ! rate of polymer transit
trial_move_rate(4) = 0.22d0 + trial_move_rate(3) ! rate of polymer rotate


NMCs = 5*10**index_con

!!!!!!!!!!!!!!!!!!!!
! angle position 
!!!!!!!!!!!!!!!!!!!!

dtheta = 1.0d0/N_theta
dphi = 2.0d0/N_phi
do i=1,N_theta
    zz(i) = dtheta * i    
end do
do i=1,N_phi
    xx(i) = 1 - dphi * i
end do

v_tide = 0
change_2 = 1
change_3 = 1
do jp = 1,N_theta 
    do ip=1,N_phi
        temp = xx(ip)*xx(ip) + zz(jp)*zz(jp)
        if (temp < 1.0d0) then
            change_2 = 1
            yyp = dsqrt(1.0d0-xx(ip)*xx(ip)-zz(jp)*zz(jp))
        else if (temp == 1.0d0) then
            change_2 = 1
            yyp = 0
        else
            change_2 = 0
        end if
		do i=1,N_phi
        	do j=1,N_theta
                temp = xx(i)*xx(i) + zz(j)*zz(j)
                if (temp < 1.0d0) then
                    change_3 = 1
                    yy  = dsqrt(1.0d0-xx(i)*xx(i)-zz(j)*zz(j))
                else if (temp == 1.0d0) then
                    change_3 = 1
                    yy = 0
                else
                    change_3 = 0
                end if                     
                if (change_2 .and. change_3 == 1) then
                    v_tide(j,i,jp,ip) = dsqrt( 1.0d0 - 0.999999999d0*(xx(i)*xx(ip) + yy*yyp + zz(j)*zz(jp))**2 ) &
                                      + dsqrt( 1.0d0 - 0.999999999d0*(xx(i)*xx(ip) - yy*yyp + zz(j)*zz(jp))**2 )
 !                   print*,xx(i)*xx(ip) + yy*yyp + zz(j)*zz(jp),xx(i)*xx(ip) - yy*yyp + zz(j)*zz(jp)
                    if (v_tide(j,i,jp,ip) /= 0) then
                    end if 
                end if
            end do
        end do
    end do
end do 
!!!!!!!!!!!!!!!!!!
! calculate rho_0
!!!!!!!!!!!!!!!!!!!!!
!rho_0 = 1.0d0*Nm * (4*n_azo + N_chain) / (4*Nz*Nr*Nr )
rho_0 = 1.0d0*Nm * (n_azo + N_chain) / (Nz*Nr*Nr ) 
!write(15,*) "rho_0", rho_0

!!!!!!!!!!!!!!!!!!!!!
! the i_azo
!!!!!!!!!!!!!!!!!!!!!
i_azo = floor(azo_position*Nm)     


! depend on Loa , hahah depend on the real loa
! compute bending coefficent of chain
epsilon_azo = 1.0d0*Nm/(4.0d0*loa_azo)
epsilon = 1.0d0*Nm/(4.0d0*loa_polymer)    
deltaS = 1.0d0*Loa/Nm


! depend on roL
! r used in MC
r_sphere = roL*Nm                                   
r_sphere_2 = r_sphere*r_sphere

  
Lz = Nm*(2.0d0*roL+4.0d0) 
Lz_2 = 0.5d0*Lz
Lr = Lz_2
Lbox = Lr
dr = 1.0d0*Lr/Nr
dz = 1.0d0*Lz/Nz


r_a(1)= 0.5d0*dr
do i= 2, nr + 1
    r_a(i) = r_a(i-1) + dr    
end do

! r used in physical world
r_dr = 1.0d0*Loa*Lr/Nr/Nm
r_dz = 1.0d0*Loa*Lz/Nz/Nm

rr_r(1)= 0.5d0*r_dr
do i= 2, nr + 1
    rr_r(i) = rr_r(i-1) + r_dr    
end do
zz_r(1)= r_dz
do i= 2, nz
    zz_r(i) = zz_r(i-1) + r_dz    
end do

!!!!!!!!!!!!
! initialize the omega
!!!!!!!!!!!!                 
w = 0
w_new = 0


!!!!!!!!!!!!!initialize azo-polymer on the sphere

allocate(azo(1:n_azo,0:Nm),ir_azo(1:n_azo,0:Nm),iz_azo(1:n_azo,0:Nm))
allocate(itheta_azo(1:n_azo,0:Nm),iphi_azo(1:n_azo,0:Nm))
open(unit=43,file='grafting_points.txt')
if ( isomer == 0 ) then
	
	do j=1, N_azo
		read(43,*) bond_vector(j)%x,bond_vector(j)%y,bond_vector(j)%z
		r_radius = dsqrt( bond_vector(j)%x**2 + bond_vector(j)%y**2 + bond_vector(j)%z**2 )
    	bond_vector(j)%x = bond_vector(j)%x / r_radius
    	bond_vector(j)%y = bond_vector(j)%y / r_radius
    	bond_vector(j)%Z = bond_vector(j)%z / r_radius
    	azo(j,0)%x = r_sphere * bond_vector(j)%x 
    	azo(j,0)%y = r_sphere * bond_vector(j)%y
    	azo(j,0)%z = r_sphere * bond_vector(j)%z
    
		do i=1, Nm        
    		unew(1) = azo(j,i-1)%x + bond_vector(j)%x 
			unew(2) = azo(j,i-1)%y + bond_vector(j)%y
    		unew(3) = azo(j,i-1)%z + bond_vector(j)%z 
        	azo(j,i)%x = unew(1)
            azo(j,i)%y = unew(2)
        	azo(j,i)%z = unew(3)             	        
		end do!!do i
	end do!! do j    

else ! isomer

	do j=1, N_azo
		read(43,*) bond_vector(j)%x,bond_vector(j)%y, bond_vector(j)%z
		r_radius = dsqrt( bond_vector(j)%x**2 + bond_vector(j)%y**2 + bond_vector(j)%z**2 )
    	bond_vector(j)%x = bond_vector(j)%x / r_radius
    	bond_vector(j)%y = bond_vector(j)%y / r_radius
    	bond_vector(j)%Z = bond_vector(j)%z / r_radius
    	azo(j,0)%x = r_sphere * bond_vector(j)%x 
    	azo(j,0)%y = r_sphere * bond_vector(j)%y
    	azo(j,0)%z = r_sphere * bond_vector(j)%z
    
		do i=1, Nm        	
			if ( i == i_azo)then  
        		azo(j,i)%x = azo(j,i-1)%x + 1.0d0
        		azo(j,i)%y = azo(j,i-1)%y                 
        		azo(j,i)%z = azo(j,i-1)%z
        	else if (i == i_azo+1) then
        		azo(j,i)%x = azo(j,i-1)%x - 0.5d0      
        		azo(j,i)%y = azo(j,i-1)%y + 0.5d0*dsqrt(3.0d0)  
        		azo(j,i)%z = azo(j,i-1)%z 
			else
    			unew(1) = azo(j,i-1)%x + bond_vector(j)%x 
				unew(2) = azo(j,i-1)%y + bond_vector(j)%y
    			unew(3) = azo(j,i-1)%z + bond_vector(j)%z 
        		azo(j,i)%x = unew(1)
            	azo(j,i)%y = unew(2)
        		azo(j,i)%z = unew(3)     
			end if
		end do!!do i
	end do!! do j    

end if
close(43)
! azo is initialized

! initialize the liquid crystal polymer 

do j=1,N_chain
    r = ran2(seed)
    if (r<=0.5) then
        temp = ran2(seed)*0.999999999d0*(Lz_2-roL*nm) + roL*nm
    else
        temp = -ran2(seed)*0.999999999d0*(Lz_2-roL*nm) - roL*nm
    end if
    cos_t = (2*ran2(seed)-1)*0.99999999d0
    sin_t = dsqrt(1.0d0-cos_t**2) 
    phi = (2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)

    polymer(j,0)%x = sin_t*dcos(phi)*temp 
    polymer(j,0)%y = sin_t*dsin(phi)*temp
    polymer(j,0)%z = cos_t*temp
     
   change = 0 
   do while(change == 0)                 ! make sure wall is impentrate
        change = 1
        cos_t=(2*ran2(seed)-1)*0.99999999d0
        sin_t=dsqrt(1.0d0-cos_t**2)
        phi = (2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
        axis(1) = sin_t*dcos(phi)
        axis(2) = sin_t*dsin(phi)
        axis(3) = cos_t    
        do i=1,Nm_chain
       
          unew(1) = polymer(j,0)%x + i*axis(1) 
          unew(2) = polymer(j,0)%y + i*axis(2)
          unew(3) = polymer(j,0)%z + i*axis(3)                    
 
          r_radius = unew(1)*unew(1) + unew(2)*unew(2) + unew(3)*unew(3)
          if (r_sphere_2 > r_radius ) then     ! jiao die le
              change = 0 
          end if
      enddo                             
  end do !! do while
  do i = 1,Nm_chain
    polymer(j,i)%x = polymer(j,0)%x + i*axis(1) 
    polymer(j,i)%y = polymer(j,0)%y + i*axis(2)
    polymer(j,i)%z = polymer(j,0)%z + i*axis(3)     
  end do
end do



open(22,file='iriz_azo.dat')
do j=1,N_azo
	r_radius = dsqrt( azo(j,0)%x*azo(j,0)%x + azo(j,0)%y*azo(j,0)%y )
    ir_azo(j,0) = floor( r_radius / dr ) + 1    
    iz_azo(j,0) = floor( (Lz_2 + azo(j,0)%z) / dz ) + 1
    do i=1,Nm   
	    r_radius = dsqrt( azo(j,i)%x*azo(j,i)%x + azo(j,i)%y*azo(j,i)%y )
    	ir_azo(j,i) = floor( r_radius / dr ) + 1    
    	iz_azo(j,i) = floor( (Lz_2 + azo(j,i)%z) / dz ) + 1

        itheta_azo(j,i) = floor( 0.99999999d0*abs(azo(j,i)%z - azo(j,i-1)%z)/dtheta ) + 1
        phi_azo = 1.0d0 - 0.99999999d0*((azo(j,i)%x - azo(j,i-1)%x)*azo(j,i)%x   &
                                       +(azo(j,i)%y - azo(j,i-1)%y)*azo(j,i)%y)  &
                                  /dsqrt(azo(j,i)%x*azo(j,i)%x + azo(j,i)%y*azo(j,i)%y)
        iphi_azo(j,i) = floor( 0.99999999d0*phi_azo/dphi ) + 1 
    	write(22,*) j,i,ir_azo(j,i),iz_azo(j,i)
    end do   
end do
close(22)

open(22,file='iriz.dat')
do j=1, N_chain
	do i=0, Nm_chain
        if ( polymer(j,i)%x>Lbox) then
            x_r = polymer(j,i)%x - Lz
        else if ( polymer(j,i)%x<(-Lbox) ) then
            x_r = polymer(j,i)%x + Lz
        else
            x_r = polymer(j,i)%x
        end if
        
        if ( polymer(j,i)%y>Lbox ) then
            y_r = polymer(j,i)%y - Lz
        else if ( polymer(j,i)%y<(-Lbox)  ) then
            y_r = polymer(j,i)%y + Lz
        else
            y_r = polymer(j,i)%y
        end if
                    
        if ( polymer(j,i)%z>Lbox ) then
            z_r = polymer(j,i)%z - Lz
        else if ( polymer(j,i)%z<(-Lbox) ) then
            z_r = polymer(j,i)%x + Lz
        else 
            z_r = polymer(j,i)%z
        end if    
        
        r_radius = dsqrt( x_r*x_r + y_r*y_r )
        
        if(r_radius<Lr)then        
			ir(j,i) = floor( r_radius / dr ) + 1
        else
            ir(j,i) = Nr + 1
		end if
    	iz(j,i) = floor( (Lz_2 + z_r) / dz ) + 1	       
        	
		write(22,*) j,i, ir(j,i), iz(j,i)
	end do   
    itheta(j) = floor( 0.99999999d0*abs(polymer(j,1)%z - polymer(j,0)%z)/dtheta ) + 1
    phi_azo = 1.0d0- 0.99999999d0*((polymer(j,1)%x - polymer(j,0)%x)*polymer(j,1)%x   &
                                  +(polymer(j,1)%y - polymer(j,0)%y)*polymer(j,1)%y)  &
                             /dsqrt(polymer(j,1)%x*polymer(j,1)%x+polymer(j,1)%y*polymer(j,1)%y)
    iphi(j) = floor( 0.99999999d0*phi_azo/dphi ) + 1
end do
close(22)
open(unit=30, file='polymer.txt')
    do j = 1, N_chain
        do i = 0,Nm_chain
            write(30,*) polymer(j,i)%x, polymer(j,i)%y, polymer(j,i)%z
        end do
    end do 
close(30)
open(unit=30, file='azo.txt')
    do j = 1, N_azo
        do i = 0,Nm
            write(30,*) azo(j,i)%x, azo(j,i)%y, azo(j,i)%z
        end do
    end do 
close(30)
call checkpolymer (flag_c)
print*,"initial OK"


j=0
do i=1,length
 
    call pivot_azo(change)
    if (change==1) then
        j = j + 1
    end if
end do
call checkpolymer (flag_c)
print*, 1.0d0*j/length,"pivot azo"

j=0
do i=1,length 
    call polymer_move(change)
    if (change==1) then
        j = j + 1
    end if
end do
call checkpolymer (flag_c)
print*, 1.0d0*j/length,"polymer_move azo"

 
j=0
do i=1,length 
    call pivot(change)
    if (change==1) then
        j = j + 1
    end if
end do
call checkpolymer (flag_c)
print*, 1.0d0*j/length,"pivot move" 


j = 0
do i=1,length
 
    call rotate_sphere(change)
    if (change==1)then
    	j = j + 1
    end if 

end do
print*, 1.0d0 * j/length,"rotate move" 
call checkpolymer (flag_c)

Deallocate(bond_vector)    

write(*,*) "CREATED OK"

open(unit=30, file='polymer_ini.txt')
    do j = 1, N_chain
        do i = 0,Nm_chain
            write(30,*) polymer(j,i)%x, polymer(j,i)%y, polymer(j,i)%z
        end do
    end do 
close(30)
open(unit=30, file='azo_ini.txt')
    do j = 1, N_azo
        do i = 0,Nm
            write(30,*) azo(j,i)%x, azo(j,i)%y, azo(j,i)%z
        end do
    end do 
close(30)
end subroutine initialize          
