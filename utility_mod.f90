MODULE utility_routines
IMPLICIT NONE
CONTAINS

subroutine pivot_azo (change)
    USE global_parameters
    IMPLICIT NONE
    INTEGER i, j, length, jj
    INTEGER :: change
    type(node) :: new(0:401)
    INTEGER :: ir_temp(0:401), iz_temp(0:501), itheta_new(0:501), iphi_new(0:401)
    DOUBLE PRECISION :: axis(3)
    DOUBLE PRECISION :: r_radius, r, cos_t, sin_t, phi
    DOUBLE PRECISION :: unew(3), uold(3)  
    DOUBLE PRECISION :: alpha, beta, angle, dotp
    DOUBLE PRECISION :: x_r, y_r
    DOUBLE PRECISION :: DE1, DE2
    DOUBLE PRECISION :: phi_azo 
    
    change = 1
    jj = floor(ran2(seed)*0.9999999d0*(N_azo)) + 1! random pickup the chain in [1,N_azo] to be rotated 
    i = floor(ran2(seed)*0.9999999d0*(Nm-1)) + 1  ! random pickup the monomer in [0,Nm-1] to be rotated         
    if (isomer/=0) then
    	if (i/=i_azo) then
        	cos_t=(2*ran2(seed) - 1)*0.99999999d0
	        sin_t=dsqrt(1.0d0 - cos_t**2) 
    	    phi = (2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
        	axis(1) = sin_t*dcos(phi)
	        axis(2) = sin_t*dsin(phi)
    	    axis(3) = cos_t
	     else if (i==i_azo) then
    	    axis(1) = azo(jj,i)%x - azo(jj,i-1)%x
        	axis(2) = azo(jj,i)%y - azo(jj,i-1)%y
        	axis(3) = azo(jj,i)%z - azo(jj,i-1)%z      
        end if
    else ! isomer ==0
        	cos_t=(2*ran2(seed) - 1)*0.99999999d0
	        sin_t=dsqrt(1.0d0 - cos_t**2) 
    	    phi = (2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
        	axis(1) = sin_t*dcos(phi)
	        axis(2) = sin_t*dsin(phi)
    	    axis(3) = cos_t    	         
    end if
    angle = rotate*(2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
    alpha = dcos(angle)
    beta = dsin(angle)

    do j = i+1, Nm     
        uold(1) = azo(jj,j)%x - azo(jj,i)%x
        uold(2) = azo(jj,j)%y - azo(jj,i)%y
        uold(3) = azo(jj,j)%z - azo(jj,i)%z
                        
        dotp = axis(1)*uold(1) + axis(2)*uold(2) + axis(3)*uold(3)
        unew(1) = uold(1)*alpha + axis(1)*dotp*(1-alpha) + ( uold(2)*axis(3) - uold(3)*axis(2) )*beta
        unew(2) = uold(2)*alpha + axis(2)*dotp*(1-alpha) + ( uold(3)*axis(1) - uold(1)*axis(3) )*beta
        unew(3) = uold(3)*alpha + axis(3)*dotp*(1-alpha) + ( uold(1)*axis(2) - uold(2)*axis(1) )*beta

        new(j)%x = azo(jj,i)%x + unew(1)
        new(j)%y = azo(jj,i)%y + unew(2)
        new(j)%z = azo(jj,i)%z + unew(3)
        r_radius = new(j)%x*new(j)%x + new(j)%y*new(j)%y + new(j)%z*new(j)%z

        if (r_sphere_2 > r_radius ) then
            change = 0 
            exit
        end if
    end do         
if (change == 1) then

    if (i==i_azo) then
        DE1 = 0
    else
        DE1 = (new(i+1)%x - 2*azo(jj,i)%x + azo(jj,i-1)%x)**2   &
            + (new(i+1)%y - 2*azo(jj,i)%y + azo(jj,i-1)%y)**2   &  
            + (new(i+1)%z - 2*azo(jj,i)%z + azo(jj,i-1)%z)**2   &
            - (azo(jj,i+1)%x - 2*azo(jj,i)%x + azo(jj,i-1)%x)**2   &
            - (azo(jj,i+1)%y - 2*azo(jj,i)%y + azo(jj,i-1)%y)**2   &  
            - (azo(jj,i+1)%z - 2*azo(jj,i)%z + azo(jj,i-1)%z)**2           
    end if   !endif i
    DE2 = 0.0d0    
    itheta_new(i+1) = floor( 0.99999999d0*abs(new(i+1)%z - azo(jj,i)%z)/dtheta ) + 1
    phi_azo = 1 - 0.99999999d0*((new(i+1)%x - azo(jj,i)%x)*new(i+1)%x   & 
                               +(new(i+1)%y - azo(jj,i)%y)*new(i+1)%y) &
                               /dsqrt(new(i+1)%x*new(i+1)%x + new(i+1)%y*new(i+1)%y)
    iphi_new(i+1) = floor( 0.99999999d0*phi_azo/dphi ) + 1
    do j = i+2, Nm
        itheta_new(j) = floor( 0.99999999d0*abs(new(j)%z - new(j-1)%z)/dtheta ) + 1
        phi_azo = 1 - 0.99999999d0*((new(j)%x - new(j-1)%x)*new(j)%x   & 
                                   +(new(j)%y - new(j-1)%y)*new(j)%y)  &
                              /dsqrt(new(j)%x*new(j)%x + new(j)%y*new(j)%y)
        iphi_new(j) = floor( 0.99999999d0*phi_azo/dphi ) + 1 
!        print*,j,ir_azo(jj,j), iz_azo(jj,j), itheta_azo(jj,j), iphi_azo(jj,j)        
    end do    
    do j = i+1, Nm
        r_radius = dsqrt( new(j)%x*new(j)%x + new(j)%y*new(j)%y )
        ir_temp(j) = floor( r_radius / dr ) + 1
        iz_temp(j) = floor( (new(j)%z + Lz_2) / dz ) + 1
!        print*,ir_temp(j), iz_temp(j), itheta_new(j), iphi_new(j)
        DE2 = DE2 + w(ir_temp(j), iz_temp(j), itheta_new(j), iphi_new(j))  &
                  - w(ir_azo(jj,j), iz_azo(jj,j), itheta_azo(jj,j), iphi_azo(jj,j) )          
    end do
    r = ran2(seed)
   
    if ( r < dexp ( -epsilon_azo*DE1 - deltaS*DE2 ))then
        do j = i+1,Nm          
            azo(jj,j)%x = new(j)%x
            azo(jj,j)%y = new(j)%y
            azo(jj,j)%z = new(j)%z
            iz_azo(jj,j) = iz_temp(j)
            ir_azo(jj,j) = ir_temp(j)
            itheta_azo(jj,j) = itheta_new(j)
            iphi_azo(jj,j) = iphi_new(j)
        end do
    else
        change = 0
    end if
end if       
end subroutine pivot_azo


subroutine pivot (change)
    USE global_parameters
    IMPLICIT NONE
    INTEGER i, j, length,jj
    INTEGER :: change
	
    type(node) :: new(0:401)
    INTEGER :: ir_temp(0:401), iz_temp(0:401)
    INTEGER :: itheta_new, iphi_new
    DOUBLE PRECISION :: axis(3)
    DOUBLE PRECISION :: r_radius,r_radius_1, r, cos_t, sin_t, phi
    DOUBLE PRECISION :: unew(3), uold(3)  
    DOUBLE PRECISION :: alpha, beta, angle, dotp    
    DOUBLE PRECISION :: DE1, DE2
    DOUBLE PRECISION :: phi_azo
    DOUBLE PRECISION :: x_r, y_r, z_r  
   
    change = 1
    jj =  1 + floor( ran2(seed)*0.9999999d0*(N_chain) ) ! random pickup the chain in [1,N_chain] to be rotated        
 
    cos_t=(2*ran2(seed)-1)*0.99999999d0
    sin_t=dsqrt(1.0d0-cos_t**2) 
    phi = (2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)

    axis(1) = sin_t*dcos(phi)
    axis(2) = sin_t*dsin(phi)
    axis(3) = cos_t

    angle = rotate*(2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
    alpha = dcos(angle)
    beta = dsin(angle)
 
    do j = 1, Nm_chain     
        uold(1) = polymer(jj,j)%x - polymer(jj,0)%x
        uold(2) = polymer(jj,j)%y - polymer(jj,0)%y
        uold(3) = polymer(jj,j)%z - polymer(jj,0)%z
                        
        dotp = axis(1)*uold(1) + axis(2)*uold(2) + axis(3)*uold(3)
        unew(1) = uold(1)*alpha + axis(1)*dotp*(1-alpha) + ( uold(2)*axis(3) - uold(3)*axis(2) )*beta
        unew(2) = uold(2)*alpha + axis(2)*dotp*(1-alpha) + ( uold(3)*axis(1) - uold(1)*axis(3) )*beta
        unew(3) = uold(3)*alpha + axis(3)*dotp*(1-alpha) + ( uold(1)*axis(2) - uold(2)*axis(1) )*beta

        new(j)%x = polymer(jj,0)%x + unew(1)
        new(j)%y = polymer(jj,0)%y + unew(2)
        new(j)%z = polymer(jj,0)%z + unew(3)
        r_radius = new(j)%x*new(j)%x + new(j)%y*new(j)%y + new(j)%z*new(j)%z

        if (r_sphere_2 > r_radius) then
            change = 0 
            exit
        end if
    end do
    if (change == 1) then     
        itheta_new = floor( 0.99999999d0*abs(new(1)%z - polymer(jj,0)%z)/dtheta ) + 1
        phi_azo = 1 - 0.99999999d0*((new(1)%x - polymer(jj,0)%x)*new(1)%x     &
                                  + (new(1)%y - polymer(jj,0)%y)*new(1)%y)  &
                              /dsqrt(new(1)%x*new(1)%x + new(1)%y*new(1)%y)
        iphi_new = floor( 0.99999999d0*phi_azo/dphi ) + 1
        do j = 1, Nm_chain  
            if ( new(j)%x>Lbox) then
                x_r = new(j)%x - Lz
            else if ( new(j)%x<(-Lbox) ) then
                x_r = new(j)%x + Lz
            else
                x_r = new(j)%x
            end if
        
            if ( new(j)%y>Lbox ) then
                y_r = new(j)%y - Lz
            else if ( new(j)%y<(-Lbox)  ) then
                y_r = new(j)%y + Lz
            else
                y_r = new(j)%y
            end if
                    
            if ( new(j)%z>Lbox ) then
                z_r = new(j)%z - Lz
            else if ( new(j)%z<(-Lbox) ) then
                z_r = new(j)%x + Lz
            else 
                z_r = new(j)%z
            end if    
            r_radius = dsqrt( x_r*x_r + y_r*y_r )
            if(r_radius<Lr)then
        	    ir_temp(j) = floor( r_radius / dr ) + 1
            else 
                ir_temp(j) = Nr + 1
            end if            
            iz_temp(j) = floor( (Lz_2 + z_r ) / dz ) + 1
!        print*,ir_temp(j),iz_temp(j),itheta_new,iphi_new      
            DE2 = DE2 + w(ir_temp(j),iz_temp(j),itheta_new,iphi_new)  &
                      - w(ir(jj,j), iz(jj,j), itheta(jj), iphi(jj))
        end do

    r = ran2(seed)
   
    if ( r < dexp ( - deltaS*DE2 ))then

        do j = 1, Nm_chain          
            polymer(jj,j)%x = new(j)%x
            polymer(jj,j)%y = new(j)%y
            polymer(jj,j)%z = new(j)%z
            iz(jj,j) = iz_temp(j)
            ir(jj,j) = ir_temp(j)
        end do  
        itheta(jj) = itheta_new
        iphi(jj) = iphi_new
    else
        change = 0
    end if
 
end if       

end subroutine pivot 



subroutine polymer_move (change)
    use global_parameters
    implicit none
    integer i, j, length, jj
    integer :: change
	
    type(node) :: new(0:401)
    Integer :: ir_temp(0:401), iz_temp(0:401)

    DOUBLE PRECISION :: axis(3)
    DOUBLE PRECISION :: r_radius,r_radius_1, r, cos_t, sin_t, phi
    DOUBLE PRECISION :: unew(3), uold(3)  
    DOUBLE PRECISION :: alpha, beta, angle, dotp    
    DOUBLE PRECISION :: DE1, DE2, DE3  
    DOUBLE PRECISION :: x_r, y_r, z_r
    change = 1
    !this is the big MC move
    jj =  1 + floor( ran2(seed)*0.9999999d0*(N_chain) ) ! random pickup the chain in [1,N_chain] to be rotated        
    
    unew(1) = rotate*(2*ran2(seed) - 1)*move_max
    unew(2) = rotate*(2*ran2(seed) - 1)*move_max
    unew(3) = rotate*(2*ran2(seed) - 1)*move_max        
    do j = 0, Nm_chain     

        new(j)%x = polymer(jj,j)%x + unew(1)
        new(j)%y = polymer(jj,j)%y + unew(2)
        new(j)%z = polymer(jj,j)%z + unew(3)
        r_radius = new(j)%x*new(j)%x + new(j)%y*new(j)%y + new(j)%z*new(j)%z

        if (r_sphere_2 > r_radius ) then
            change = 0 
            exit
        end if
    end do 
      
    if (change == 1) then   

        if ( new(0)%x>Lbox) then   
            do j = 0,Nm_chain
                new(j)%x = new(j)%x - Lz
            end do
        else if (new(0)%x<(-Lbox)) then
            do j = 0,Nm_chain
                new(j)%x = new(j)%x + Lz
            end do
        end if

        if ( new(0)%y>Lbox) then   
            do j = 0,Nm_chain
                new(j)%y = new(j)%y - Lz
            end do
        else if (new(0)%y<(-Lbox)) then
            do j = 0,Nm_chain
                new(j)%y = new(j)%y + Lz
            end do
        end if

        if ( new(0)%z>Lbox) then   
            do j = 0,Nm_chain
                new(j)%z = new(j)%z - Lz
            end do
        else if (new(0)%z<(-Lbox)) then
            do j = 0,Nm_chain
                new(j)%z = new(j)%z + Lz
            end do
        end if

        do j = 0, Nm_chain

            if ( new(j)%x>Lbox) then
                x_r = new(j)%x - Lz
            else if ( new(j)%x<(-Lbox) ) then
                x_r = new(j)%x + Lz
            else
                x_r = new(j)%x
            end if
        
            if ( new(j)%y>Lbox ) then
                y_r = new(j)%y - Lz
            else if ( new(j)%y<(-Lbox)  ) then
                y_r = new(j)%y + Lz
            else
                y_r = new(j)%y
            end if
                    
            if ( new(j)%z>Lbox ) then
                z_r = new(j)%z - Lz
            else if ( new(j)%z<(-Lbox) ) then
                z_r = new(j)%x + Lz
            else 
                z_r = new(j)%z
            end if    

            r_radius = dsqrt( x_r*x_r + y_r*y_r )
            if(r_radius<Lr)then
        	    ir_temp(j) = floor( r_radius / dr ) + 1
            else 
                ir_temp(j) = Nr + 1
            end if            
            iz_temp(j) = floor( (Lz_2 + z_r ) / dz ) + 1
            DE2 = DE2 + w(ir_temp(j), iz_temp(j), itheta(jj), iphi(jj) )  &
                      - w(ir_azo(jj,j), iz_azo(jj,j), itheta(jj), iphi(jj) )
        end do


    r = ran2(seed)
   
    if ( r < dexp ( - deltaS*DE2 ))then
        
        do j = 0, Nm_chain          
            polymer(jj,j)%x = new(j)%x
            polymer(jj,j)%y = new(j)%y
            polymer(jj,j)%z = new(j)%z
            iz(jj,j) = iz_temp(j)
            ir(jj,j) = ir_temp(j)
        end do  
 
    else     
        change = 0
    end if
 
end if       

end subroutine polymer_move



subroutine rotate_sphere (change)
    use global_parameters
    implicit none
    integer i, j, k, l 
    integer :: change 
    type(node) :: new(1:N_azo,0:Nm)
    Integer :: ir_temp(1:N_azo,0:Nm), iz_temp(1:N_azo,0:Nm)
    Integer :: itheta_new(1:N_azo,0:Nm), iphi_new(1:N_azo,0:Nm)
    DOUBLE PRECISION :: axis(3)
    DOUBLE PRECISION :: r_radius, r, cos_t, sin_t, phi 
    DOUBLE PRECISION :: uold(3)  
    DOUBLE PRECISION :: alpha, beta, angle, dotp
    DOUBLE PRECISION :: DE2
    DOUBLE PRECISION :: phi_azo
   
    change = 1       
    cos_t=(2*ran2(seed)-1)*0.99999999d0
    sin_t=dsqrt(1.0d0-cos_t**2) 
    phi = (2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)

    axis(1) = sin_t*dcos(phi)
    axis(2) = sin_t*dsin(phi)
    axis(3) = cos_t

    angle = rotate_s*(2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
    alpha = dcos(angle)
    beta = dsin(angle)
 
    do j = 1,N_azo
        do i = 0, Nm     
            uold(1) = azo(j,i)%x 
            uold(2) = azo(j,i)%y 
            uold(3) = azo(j,i)%z 
                        
            dotp = axis(1)*uold(1) + axis(2)*uold(2) + axis(3)*uold(3)
            new(j,i)%x = uold(1)*alpha + axis(1)*dotp*(1-alpha) + ( uold(2)*axis(3) - uold(3)*axis(2) )*beta
            new(j,i)%y = uold(2)*alpha + axis(2)*dotp*(1-alpha) + ( uold(3)*axis(1) - uold(1)*axis(3) )*beta
            new(j,i)%z = uold(3)*alpha + axis(3)*dotp*(1-alpha) + ( uold(1)*axis(2) - uold(2)*axis(1) )*beta                
        end do
    end do


    DE2 = 0
    do j = 1, N_azo
        i = 0
            r_radius = dsqrt(new(j,i)%x*new(j,i)%x + new(j,i)%y*new(j,i)%y)
            ir_temp(j,i) = floor( r_radius / dr ) + 1 
            iz_temp(j,i) = floor( ( new(j,i)%z + Lz_2 ) / dz ) + 1
            DE2 = DE2 + w(ir_temp(j,i), iz_temp(j,i), 0, 0)  &
                      - w(ir_azo(j,i), iz_azo(j,i), 0, 0)              
        do i= 1, Nm 
            r_radius = dsqrt(new(j,i)%x*new(j,i)%x + new(j,i)%y*new(j,i)%y)
        	ir_temp(j,i) = floor( r_radius / dr ) + 1 
            iz_temp(j,i) = floor( ( new(j,i)%z + Lz_2 ) / dz ) + 1  
            itheta_new(j,i) = floor( 0.99999999d0*abs(new(j,i)%z - new(j,i-1)%z)/dtheta ) + 1
            phi_azo = 1 - 0.99999999d0*((new(j,i)%x - new(j,i-1)%x)*new(j,i)%x    &
                                       +(new(j,i)%y - new(j,i-1)%y)*new(j,i)%y)   &
                                  /dsqrt(new(j,i)%x*new(j,i)%x + new(j,i)%y*new(j,i)%y)
            iphi_new(j,i) = floor( 0.99999999d0*phi_azo/dphi ) + 1 

            DE2 = DE2 + w(ir_temp(j,i), iz_temp(j,i), itheta_new(j,i), iphi_new(j,i))  &
                      - w(ir_azo(j,i), iz_azo(j,i), itheta_azo(j,i), iphi_azo(j,i) )            
            
        end do
    end do

    r = ran2(seed)
    if ( r < dexp ( - deltaS*DE2 ) )then
        do j = 1, N_azo
            do i = 0, Nm          
                azo(j,i)%x = new(j,i)%x
                azo(j,i)%y = new(j,i)%y
                azo(j,i)%z = new(j,i)%z
                iz_azo(j,i) = iz_temp(j,i)
                ir_azo(j,i) = ir_temp(j,i)
            end do
            do i=1, nm
                itheta_azo(j,j) = itheta_new(j,i)
                iphi_azo(j,i) = iphi_new(j,i)
            end do
        end do            

    else
        change = 0
    end if
           
End subroutine rotate_sphere 

!!!!!!!!!!This is the function that creat a random number of [0-1]!!!!!!!!!!!!!!
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
        IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,	 &
        NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
        end do
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=AM*iy
      if(ran2.gt.RNMX) ran2=RNMX
      return
      END function  


subroutine checkpolymer (flag_c)
use global_parameters
implicit none
integer :: i, j, flag_c,phi,theta
double precision :: r
if(isomer==1)then
    do j=1,N_azo
        r = (azo(j,i_azo)%x - azo(j,i_azo-1)%x)*(azo(j,i_azo+1)%x - azo(j,i_azo)%x)    &
           +(azo(j,i_azo)%y - azo(j,i_azo-1)%y)*(azo(j,i_azo+1)%y - azo(j,i_azo)%y)    &
           +(azo(j,i_azo)%z - azo(j,i_azo-1)%z)*(azo(j,i_azo+1)%z - azo(j,i_azo)%z) 
        if( abs(r-dcos(2.0d0*pi/3.0d0))>1.0d-6 )then  
            print*, i_azo,"i_azo is broken,error 1"    
            print*,azo(j,i_azo-1)%x,azo(j,i_azo-1)%y,azo(j,i_azo-1)%z
            print*,azo(j,i_azo)%x,azo(j,i_azo)%y,azo(j,i_azo)%z
            print*,azo(j,i_azo+1)%x,azo(j,i_azo+1)%y,azo(j,i_azo+1)%z
            stop        
        end if
    end do
end if

do j = 1, N_azo
	do i=1, Nm
	   r = (azo(j,i)%x - azo(j,i-1)%x)**2    &
		 + (azo(j,i)%y - azo(j,i-1)%y)**2    &
         + (azo(j,i)%z - azo(j,i-1)%z)**2         
      if (abs(r-1)>1.d-5) then
          print*, "azobond is broken,error 2", j,i, "length",abs(r-1)
          print*,azo(j,i)%x,azo(j,i-1)%x
          print*,azo(j,i)%y,azo(j,i-1)%y
          print*,azo(j,i)%z,azo(j,i-1)%z
          flag_c = 1
          stop
      end if
  end do
end do

do j= 1, N_azo
    do i=1, Nm
        r = dsqrt(azo(j,i)%z*azo(j,i)%z &
    			+ azo(j,i)%y*azo(j,i)%y &
    			+ azo(j,i)%x*azo(j,i)%x )
        if ( (r-r_sphere) <= -1.0d-5) then
           print*, "azomonoer",i,"on chain", j,"overlap sphere_1"
           print*,azo(j,i)%x
           print*,azo(j,i)%y
           print*,azo(j,i)%z
           flag_c = 1
           stop
        end if        
    end do
end do

do j=1,N_azo
    do i=1,Nm   
        if (ir_azo(j,i)>Nr+1) then
            print*,"azo ir is beyond the box",j,i,ir_azo(j,i)
        else if (iz_azo(j,i)>Nz) then
            print*,"azo iz is beyond the box",j,i,iz_azo(j,i)
        else if (itheta_azo(j,i)>N_theta) then
            print*,"azo itheta is beyond the box",j,i,itheta_azo(j,i)
        else if (iphi_azo(j,i)>N_phi) then
            print*,"azo iphi is beyond the box",j,i,iphi_azo(j,i)
        flag_c = 1
        stop
        end if
    end do   
end do

do j = 1, N_chain
    do i=1, Nm_chain
		    r = (polymer(j,i)%x - polymer(j,i-1)%x)**2    &
			  + (polymer(j,i)%y - polymer(j,i-1)%y)**2    &
       		  + (polymer(j,i)%z - polymer(j,i-1)%z)**2         
        if (abs(r-1)>1.d-5) then
            print*, "polymerbond is broken", j,i, "length",abs(r-1)
            print*,polymer(j,i)%x,polymer(j,i-1)%x
            print*,polymer(j,i)%y,polymer(j,i-1)%y
            print*,polymer(j,i)%z,polymer(j,i-1)%z
            flag_c = 1
            stop
        end if
    enddo
end do

do j= 1, N_chain
    do i=1, Nm_chain
        r = dsqrt(polymer(j,i)%z*polymer(j,i)%z &
    		    + polymer(j,i)%y*polymer(j,i)%y &
    			+ polymer(j,i)%x*polymer(j,i)%x )
        if ( (r-r_sphere) <= -1.0d-5) then
            print*, "monoer",i,"on chain", j,"overlap sphere_1"
            print*,polymer(j,i)%x
            print*,polymer(j,i)%y
            print*,polymer(j,i)%z
            flag_c = 1
            stop
        end if      
    end do
end do

do j=1, N_chain
    do i=0, Nm_chain
        if (ir(j,i)>Nr+1) then
            print*,"polymer ir is beyond the box",j,i,ir(j,i)
        else if (iz(j,i)>Nz) then
            print*,"polymer iz is beyond the box",j,i,iz(j,i)
        flag_c = 1
        stop
        end if       
    end do   
    if (itheta(j)>N_theta) then
        print*,"polymer itheta is beyond the box",j,i,itheta(j)
    else if (iphi(j)>N_phi) then
        print*,"polymer iphi is beyond the box",j,i,iphi(j)
        flag_c = 1
        stop
    end if
end do

end subroutine checkpolymer
        

             
  SUBROUTINE simpson(g,h,rr)
  IMPLICIT NONE
  Double precision, INTENT(IN) :: h
  Double precision, INTENT(OUT) :: rr
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: g
  integer :: n,i
  rr=0.0d0
  n=size(g)

  do i=4,n-3
    rr=rr+g(i)
  end do
  rr=3.0d0*(g(1)+g(n))/8.0d0 &
     +7.0d0*(g(2)+g(n-1))/6.0d0 &
      +23.0d0*(g(3)+g(n-2))/24.0d0 &
       + rr
  rr=rr*h
  END SUBROUTINE simpson 
  

SUBROUTINE comformation_write()
use global_parameters
implicit none
integer :: i,j
open(unit=27,file='azo_ini.txt')
do j=1,n_azo
    do i=0,nm
        write(27,*) azo(j,i)%x, azo(j,i)%y,azo(j,i)%z
    end do
end do
close(27)
open(unit=27,file='polymer_ini.txt')
do j=1,n_azo
    do i=0,nm
        write(27,*) polymer(j,i)%x, polymer(j,i)%y,polymer(j,i)%z
    end do
end do
close(27)
END SUBROUTINE comformation_write

SUBROUTINE comformation_out()
use global_parameters
implicit none
integer :: i, flag_c,j
character*7 res
character res0,res1,res2

!call checkpolymer (flag_c)

do j=1,N_chain
        
    res0=achar(48+mod(j,10))
    res1=achar(48+mod(int(j/10),10))
    res2=achar(48+int(j/100))
    res=res2 // res1 // res0 // '.dat'
    open(4,file=res,status='new')
 
	do i=0,nm
    	!write(4,"(A4,27X,F7.3,1X,F7.3,1X,F7.3)") "ATOM",polymer(j,i)%x,polymer(j,i)%y,polymer(j,i)%z
    	write(4,*) polymer(j,i)%x,polymer(j,i)%y,polymer(j,i)%z
	end do
	close(4) 
end do
END SUBROUTINE comformation_out
SUBROUTINE comformation_azo_out()
use global_parameters
implicit none
integer :: i, flag_c,j
character*7 res
character res0,res1,res2

!call checkpolymer (flag_c)

do j=1,N_azo
        
    res0=achar(48+mod(j,10))
    res1=achar(48+mod(int(j/10),10))
    res2=achar(48+int(j/100))
    res=res2 // res1 // res0 // '.dat'
    open(4,file=res,status='new')
 
	do i=0,nm
    	!write(4,"(A4,27X,F7.3,1X,F7.3,1X,F7.3)") "ATOM",azo(j,i)%x,azo(j,i)%y,azo(j,i)%z
    	write(4,*) azo(j,i)%x,azo(j,i)%y,azo(j,i)%z
	end do
	close(4) 
end do
END SUBROUTINE comformation_azo_out

END MODULE utility_routines
