subroutine SCMFT()
USE global_parameters
USE utility_routines    

implicit none

Integer :: i, j, k,jp,ip
integer :: change ! change is the flag of MC, 1 for accepte the move.
Integer :: i_move,i_rotate, i_pivot, i_small
integer :: n_move,n_rotate, n_pivot, n_small 
integer :: flag_c
DOUBLE PRECISION, PARAMETER :: TOL = 1.5D-3
DOUBLE PRECISION ::  w_erro, r
DOUBLE PRECISION :: derro(1:70),erro(1:70),errodown
DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: density


character*7 res,resres,res_s,res_o
character res0,res1,res2

allocate( density(1:Nr+1,1:Nz,1:N_theta,1:N_phi) )
allocate( w_new(1:Nr+1,1:Nz,0:N_theta,0:N_phi) )
!!! iteration
 
open(unit=15, file='Npre.txt')
open(unit=61, file='w.txt')

do n_iter = 1, Max_iter
print*, "start", n_iter,"iteration"

!    res0=achar(48+mod(n_iter,10))
!    res1=achar(48+mod(int(n_iter/10),10))
!    res2=achar(48+int(n_iter/100))
!    res_s =res2 // res1 // res0 // '.ome'
!    res_o =res2 // res1 // res0 // '.pol'
!    res =res2 // res1 // res0 // '.den'
!    resres = res2 // res1 // res0 // '.azo'
!    open(61,file=res_s,status='new')
!    open(62,file=res_o,status='new')
!    open(63,file=res,status='new')
!    open(64,file=resres,status='new')         
    ! pre-move
        N_pre = 0
        n_move = 0
        n_rotate = 0 
        n_pivot = 0
        n_small = 0
        i_move = 0
        i_rotate = 0 
        i_pivot = 0
        i_small = 0
       
        do while(N_pre < Npre)              
            r = ran2(seed)
            if ( r<trial_move_rate(1) ) then
                n_pivot = n_pivot + 1
                call checkpolymer (flag_c)
                call pivot_azo(change)
                if (change == 1) then
                    i_pivot = i_pivot + 1
                end if 
            else if (r>=trial_move_rate(1) .and. r<trial_move_rate(2) ) then
                n_rotate = n_rotate + 1
                call checkpolymer (flag_c)
                call rotate_sphere(change)
                if (change == 1) then
                    i_rotate = i_rotate + 1
                end if 

            else if ( r>=trial_move_rate(2) .and. r<trial_move_rate(3)  ) then
                n_small = n_small + 1
                call checkpolymer (flag_c)
                call pivot(change)
                if (change == 1) then
                    i_small = i_small + 1
                end if 
                 
            else
                n_move = n_move + 1
                call checkpolymer (flag_c)
                call polymer_move(change)
                if (change == 1) then
                    i_move = i_move + 1
                end if    
            end if 
            
            if(change == 1)then
                N_pre = 1 + N_pre         
            end if
        end do  
        write(15,*) n_iter
        write(15,*) i_move, 1.0d0 * i_move/n_move,"polymer move"
        write(15,*) i_rotate, 1.0d0 * i_rotate/n_rotate,"rotate move"
        write(15,*) i_pivot, 1.0d0 * i_pivot/n_pivot,"azo pivot"    
        write(15,*) i_small, 1.0d0 * i_small/n_small,"polymer pivot"
!! find out w_new
    
    MCS = 0
    density = 0
    do while(MCS < NMCs)
               
        MCS = MCS + 1
  
        moves = 0

        do while(moves < Nmove)              
        
            r = ran2(seed)

            if ( r<trial_move_rate(1) ) then
                call checkpolymer (flag_c)
                call pivot_azo(change)

            else if (r>=trial_move_rate(1) .and. r<trial_move_rate(2) ) then
                call checkpolymer (flag_c)
                call rotate_sphere(change)

            else if ( r>=trial_move_rate(2) .and. r<trial_move_rate(3)  ) then
                call checkpolymer (flag_c)
                call pivot(change)
                 
            else 
                call checkpolymer (flag_c)
                call polymer_move(change)
    
            end if 
     
            if(change == 1)then
                moves = 1 + moves         
            end if  
        end do   
    
        do j=1,N_chain
            density(ir(j,0),iz(j,0),0,0) = density(ir(j,0),iz(j,0),0,0) + 1
            do i=1,Nm_chain
                density(ir(j,i),iz(j,i),itheta(j),iphi(j)) &
                = density( ir(j,i),iz(j,i),itheta(j),iphi(j) ) + 1
!                print*,density(ir(j,i),iz(j,i),itheta(j),iphi(j))                
            end do
        end do

        do j=1,N_azo
            density(ir_azo(j,0),iz_azo(j,0),0,0) = density(ir_azo(j,0),iz_azo(j,0),0,0) + 1
            do i=1,Nm
                density(ir_azo(j,i),iz_azo(j,i),itheta_azo(j,i),iphi_azo(j,i)) &
                = density(ir_azo(j,i),iz_azo(j,i),itheta_azo(j,i),iphi_azo(j,i)) + 1                
            end do
        end do
    
    end do   ! MCS
    density = 0.5d0*deltaS*density / MCS  ! here 0.5 is consider the symetry of f(varphi) = f(-varphi) 
    do i=1,nr+1  
        density (i,:,:,:) = density (i,:,:,:)/( r_a(i)*r_dr*r_dz )
    end do
    w_new = 0
    do j = 1, N_theta
        do i = 1, N_phi
            do jp = 1, N_theta
                do ip = 1, N_phi
                    w_new(:,:,j,i) = w_new(:,:,j,i) + nu*density(:,:,jp,ip) * v_tide(j,i,jp,ip)
                end do
            end do
        end do
    end do
   ! compute erros
    w_erro = 0

    do j = 1, Nz
        do i = 1, Nr+1
            do jp =1, N_theta
                do ip = 1, N_phi
                    w_erro = w_erro + abs(w_new(i,j,jp,ip) - w(i,j,jp,ip))  
                end do
            end do
       end do
    end do
 
    w_erro = 1.d0*w_erro/Nz/Nr/N_theta/N_phi    
    
    print*, "SCMFT", n_iter, w_erro
   
    erro(n_iter) = w_erro
    if (w_erro<TOL .and. n_iter>3) then
        exit
    end if
    
    if (n_iter>5) then
        
        derro(n_iter) = erro(n_iter) - erro(n_iter-1)
        if (n_iter>10) then
            errodown = derro(n_iter) + derro(n_iter-1) + derro(n_iter-2) + derro(n_iter-3) + derro(n_iter-4) 
        else 
            errodown = -1.0d0                   
        end if
    end if
    if (n_iter==10) then
        do j = 1, Nr+1 
            do i = 1, Nz
                do jp =1, N_theta
                    do ip = 1, N_phi
                        write(61,*) w_new(j,i,jp,ip)  
                    end do
                end do
            end do
        end do        
    end if

    if (errodown>0.0d0 .or. n_iter==20) then
        print*,"errodown=",errodown
        do j = 1, Nr+1 
            do i = 1, Nz
                do jp =1, N_theta
                    do ip = 1, N_phi
                        write(61,*) w_new(j,i,jp,ip)  
                    end do
                end do
            end do
        end do
        exit
    end if    
    
    !simple mixing scheme
    w = lambda*w_new + (1-lambda)*w

    ! boundary condition
    do j = 1, Nz
        do jp =1, N_theta
            do ip = 1, N_phi
              w(nr+1,j,jp,ip) = w(nr,j,jp,ip)
            end do
        end do 
    end do         
    call checkpolymer (flag_c)             
end do  ! enddo n_iter
close(15)
deallocate(w_new)
stop"SCMFT is ok"
end subroutine SCMFT
