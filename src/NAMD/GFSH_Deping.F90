!
!...Surface hopping simulation of model system, like Tully.
!...Deping Hu,10/13/2022
!---------------------------------------------------------------------
!...initialize the simulation
subroutine initial
use global
implicit none
integer id
write(*,*)"initial begin"
write(*,*)"start_state=",start_state

call allocate_global

!set initial phase
phase=1.d0

!read positions and velocities
call read_xv

if(trim(adjustl(start_rep))=="adiabatic")then
   cs_ad=(0.d0,0.d0)
   cs_ad(start_state)=(1.d0,0.d0)
else if(trim(adjustl(start_rep))=="diabatic")then
   cs_d=(0.d0,0.d0)
   cs_d(start_state)=(1.d0,0.d0)
else
   stop"start_rep we don't know."
endif

active_state=start_state

!open files
call open_files

write(*,*)"initial end"
end subroutine initial
!---------------------------------------------------------------------
!...calculate the probability of current active state
subroutine calc_prob_active_state
use global
implicit none

prob_active_state=conjg(cs_ad(active_state))*cs_ad(active_state)
!write(*,*)"current occupation prob=",prob_active_state

end subroutine calc_prob_active_state
!---------------------------------------------------------------------
!...calculate density matrix and other properties
subroutine property
use global
implicit none
integer::is,js,id
real*8::dentemp
complex*16,allocatable::U_temp(:,:)
complex*16::trace_rho

!nuclear kintetic energy
Ekin=0.d0
do id=1,ndof
   Ekin=Ekin+0.5d0*mass(id)*vel(id)*vel(id)
enddo

!density matrix in adiabatic
do is=1,ncs
   do js=1,ncs
      denM_ad(is,js)=conjg(cs_ad(is))*cs_ad(js)
   enddo
enddo

!density matrix in diabatic
call Ad2Di(ncs,active_state,U_ss,denM_ad,denM)

!electric potential energy
Epot=H_ad(active_state,active_state)
Etot=Ekin+Epot

trace_rho=trace(denM_ad,ncs)
!write(*,*)"trace(denMad)=",trace_rho

call write_files

end subroutine property
!---------------------------------------------------------------------
!check the nuclear position whether or not in the [xmin,xmax]
subroutine check_xrange(outside)
use global,only:pos,xmin,xmax
implicit none
logical,intent(out)::outside

outside=.true.
if(pos(1)>xmin.and.pos(1)<xmax)outside=.false.

end subroutine check_xrange
!---------------------------------------------------------------------
!...close files and deallocate arraies
subroutine finish
use global
implicit none

call close_files

call deallocate_global

end subroutine finish
!---------------------------------------------------------------------
!...calculate hopping probility
subroutine calc_hop_prob
use global
implicit none
integer::is,js,kd
real*8::temp,hop_sign,rho_decress_all

hop_prob=0.d0

!write(*,*)"prob_active_state_old=",prob_active_state_old
!write(*,*)"prob_active_state=",prob_active_state
if(prob_active_state>prob_active_state_old)then
  !the population of active state increase, so there is no hopping
  return
else
  hop_sign=1.d0
endif

!density matrix
do is=1,ncs
   do js=1,ncs
      denM_ad(is,js)=conjg(cs_ad(is))*cs_ad(js)
      denM_ad_old(is,js)=conjg(cs_ad_old(is))*cs_ad_old(js)
   enddo
enddo

!calculate the jump prob from start_state to any other states
!do is=1,ncs
!   if(is==active_state)cycle
!   temp=real(denM_ad(active_state,is))*tdc(active_state,is)
!   !write(*,*)"temp=",temp
!   call pass_sign_A2B(hop_sign,temp)
!   hop_prob(is)=2.d0*timestep/real(denM_ad(active_state,active_state))*temp
!   hop_prob(is)=max(hop_prob(is),0.d0)
!enddo

if (real(denM_ad(active_state, active_state)) .lt. real(denM_ad_old(active_state, active_state))) then
   rho_decress_all = 0.d0
   do is=1,ncs
      if (real(denM_ad(is, is)) .lt. real(denM_ad_old(is, is))) then
         rho_decress_all = rho_decress_all + (real(denM_ad_old(is, is)) - real(denM_ad(is, is)))
      end if
   end do
   do is=1,ncs
      if(is==active_state)cycle
      hop_prob(is) = (real(denM_ad(is,is)) - real(denM_ad_old(is,is))) &
                   /real(denM_ad_old(active_state, active_state)) &
                   *(real(denM_ad_old(active_state, active_state)) - real(denM_ad(active_state, active_state))) &
                   /rho_decress_all
      hop_prob(is)=max(hop_prob(is),0.d0)
   end do
endif

end subroutine calc_hop_prob
!...hopping between different electronic states
subroutine hopping
use global
implicit none
real*8::rand,sum_prob1,sum_prob2,aIJ,bIJ,gammaIJ,delta
integer::is,js,kd,jump2state,old_state,new_state
real*8::pro_sum(ncs)

hop=.false.

call random_number(rand)
!write(*,*)"rand=",rand

pro_sum = 0.0
do is = 1, ncs
  do js = 1, is
     pro_sum(is) = pro_sum(is) + hop_prob(js)
  end do
end do

old_state = active_state
new_state = active_state

if ((0 .le. rand) .and. (rand .lt. pro_sum(1)) ) then
   new_state = 1
else
   do is = 2, ncs
      if ((pro_sum(is - 1) .le. rand) .and. (rand .lt. pro_sum(is))) then
         new_state = is
      end if
   end do
end if

if (new_state .ne. old_state) then
   hop=.true.
   jump2state=new_state
endif

if(hop.eqv..true.)then
  aIJ=0.d0 ; bIJ=0.d0
  do kd=1,ndof
     aIJ=aIJ+NACv(active_state,jump2state,kd)**2/mass(kd)
     bIJ=bIJ+NACv(active_state,jump2state,kd)*vel(kd)
  enddo
  aIJ=aIJ*0.5d0
  !write(*,*)"aIJ=",aIJ
  !write(*,*)"bIJ=",bIJ

  !calculate available kinetic energy
  delta=bIJ**2+4.d0*aIJ*(H_ad(active_state,active_state)-H_ad(jump2state,jump2state))
  !write(*,*)"delta=",delta
  if(delta<0.d0)then
     !not enough kinetic energy; reject the hop
     hop=.false.
  else if(bIJ<0.d0)then
     gammaIJ=(bIJ+sqrt(delta))/aIJ*0.5d0
  else
     gammaIJ=(bIJ-sqrt(delta))/aIJ*0.5d0
  endif
endif
write(*,*) "5"

if(hop.eqv..true.)then
   !write(*,*)"hop=",hop
   !write(*,*)"gamma",(bIJ+sqrt(delta))/aIJ*0.5d0,(bIJ-sqrt(delta))/aIJ*0.5d0
   write(*,*)"current step=",cstep
   !write(*,*)"gammaIj=",gammaIJ
   !write(*,*)"velocity=",vel(:)
   !write(*,*)"NACv=",NACv(active_state,jump2state,1)

   !rescale the velocities
   do kd=1,ndof
      vel(kd)=vel(kd)-gammaIJ*NACv(active_state,jump2state,kd)/mass(kd)
   enddo
   !write(*,*)"velocity=",vel(:)

   !new time derivative coupling
   call calc_tdc

   !calculate new force
   active_state=jump2state
   call calc_force

endif
write(*,*) "6"

end subroutine hopping
!...energy based decoherence correction
subroutine decoherence
use global
implicit none
integer::is,id
real*8::tao,temp,sumc

!no decoherence correction
if(decoh_alpha<=0.d0)return

!nuclear kintetic energy
Ekin=0.d0
do id=1,ndof
   Ekin=Ekin+0.5d0*mass(id)*vel(id)*vel(id)
enddo

sumc=0.d0
do is=1,ncs
   if(is==active_state)cycle
   temp=abs(H_ad(is,is)-H_ad(active_state,active_state))
   !write(*,*)"temp=",temp
   tao=(1.d0+decoh_alpha/Ekin)/temp
   cs_ad(is)=cs_ad(is)*exp(-timestep/tao)
   sumc=sumc+abs(cs_ad(is))**2
enddo

cs_ad(active_state)=cs_ad(active_state)*sqrt(1.d0-sumc)/abs(cs_ad(active_state))

!write(*,*)"coherence2=",conjg(cs_ad(1))*cs_ad(2)
end subroutine decoherence
!---------------------------------------------------------------------
!change the denM to denMad
subroutine Di2Ad(ncs,U_ss,denM,denM_ad)
use linalg
implicit none
integer,intent(in)::ncs
real*8,intent(in)::U_ss(ncs,ncs)
complex*16,intent(in)::denM(ncs,ncs)
complex*16,intent(out)::denM_ad(ncs,ncs)

denM_ad=denM
call Simil_Trans("UtAU",ncs,U_ss,denM_ad)

end subroutine Di2Ad
!---------------------------------------------------------------------
subroutine Ad2Di(ncs,active_state,U_ss,denM_ad,denM)
use linalg
implicit none
integer,intent(in)::ncs,active_state
real*8,intent(in)::U_ss(ncs,ncs)
complex*16,intent(in)::denM_ad(ncs,ncs)
complex*16,intent(out)::denM(ncs,ncs)
integer::is,js

denM=denM_ad

do is=1,ncs
   denM(is,is) = 0
enddo

denM(active_state,active_state) = 1

call Simil_Trans("UAUt",ncs,U_ss,denM)

end subroutine Ad2Di
!---------------------------------------------------------------------
!transform cs_d to cs_ad
subroutine cs_d2_ad(ncs,U_ss,cs_d,cs_ad)
use linalg
implicit none
integer,intent(in)::ncs
real*8,intent(in)::U_ss(ncs,ncs)
complex*16,intent(in)::cs_d(ncs)
complex*16,intent(out)::cs_ad(ncs)
real*8,allocatable::Ut(:,:)
integer::is

allocate(Ut(ncs,ncs))
Ut=transpose(U_ss)

cs_ad=0.d0
do is=1,ncs
   cs_ad(is)=cs_ad(is)+dot_product(Ut(is,:),cs_d(:))
enddo

deallocate(Ut)
end subroutine cs_d2_ad

!determine the initial index of the adiabatic state accrording 
! the initial density
subroutine initial_state(ncs,cs_ad,active_state)
use linalg
implicit none
integer,intent(in)::ncs
complex*16,intent(in)::cs_ad(ncs)
real*8::cs_ad_pro(ncs),inte_cs_ad(ncs)
real*8::tmp, ran
integer::is,js,active_state

do is=1,ncs
   cs_ad_pro(is) = conjg(cs_ad(is))*cs_ad(is)
enddo

tmp = 0.d0
do is=1,ncs
   tmp = tmp + cs_ad_pro(is)
   inte_cs_ad(is) = tmp
enddo

call random_number(ran)

do is=1,ncs
   if (ran .lt. inte_cs_ad(is)) then
      active_state = is
      exit
   endif
enddo

write(*,*) "active_state", active_state

end subroutine initial_state


!---------------------------------------------------------------------
!...single trajectory Surface hopping dynamics in cavity
subroutine sh
use global
implicit none
integer::istep,iter,is,js
logical::outside

cstep=0

!initialize the dynamics
call initial
call get_H_d   !Hamiltonian and its gradient of diabatic model
call get_H_ad  !Hamiltonian and its gradient in adiabatic basis

!call energy_shift(nstate,H_ad)
call NACvec
call calc_tdc

!transform denM.out to denMad.out
if(trim(adjustl(start_rep))=="diabatic")then

  call cs_d2_ad(ncs,U_ss,cs_d,cs_ad)
  call initial_state(ncs,cs_ad,active_state)

endif

call calc_force
call calc_prob_active_state
call save_old

!call energy_shift_back(nstate,H_ad)
!density matrix and other properties
call property

do istep=1,Nstep
   cstep=istep
   write(*,*)"cstep=",cstep

   !propagate nuclear positions
   call VVerlet_xstep
!   call check_xrange(outside)
!   if(outside)exit

   !get new Hamiltonian
   call get_H_d   !Hamiltonian and its gradient of diabatic model
   call get_H_ad  !Hamiltonian and its gradient in adiabatic basis
   call get_overlap
   call get_phase
   call correct_phase
   !call energy_shift(nstate,H_ad)
   call NACvec
   call calc_tdc

   !new force
   call calc_force

   !update nuclear velocities
   call VVerlet_vstep

   !Solving TDSE
   cs_ad_old = cs_ad
   call propa_wf
!   write(59,'(i6,1x,100e23.15)')istep,cs_ad(:)
   call calc_prob_active_state
   !calculate hopping probility
   call calc_hop_prob

   !hopping
   call hopping

   !decoherence correction
   call decoherence

   !save old quantities for the future
   call save_old

   !call energy_shift_back(nstate,H_ad)
   !density matrix and other properties
   if(mod(istep,dump_freq)==0)then
      call property
   endif

enddo

call finish

end subroutine sh
!
!...main program starts here
!
program sh_ensemble
use global
implicit none
integer::i,ini_seed,clock

open(unit=88,file="seed",action="read")
read(88,*) ini_seed
close(88)

call system_clock(count=clock)
ini_seed = ini_seed + clock
write(*,*) "ini_seed", ini_seed
call init_random_seed(ini_seed)

call read_input
if(np>1)then
  ncs=np*nstate !# of cavity states
else
  ncs=nstate
endif

!#if(ndof/=1)then
!#  stop"We only implement ndof=1"
!#endif

write(*,*)"Ntraj=",Ntraj
if(Ntraj==1)then
  workdir='.'
  call sh
else if(Ntraj>1)then
   do i=start_traj,start_traj+Ntraj-1
      call construct_name("Trajectory/TRAJ_",i,workdir,5)
      write(*,*)"workdir=",workdir
      call sh
   enddo
else
  write(*,*)"Ntraj=",Ntraj," is wrong!"
  stop
endif

end program sh_ensemble
