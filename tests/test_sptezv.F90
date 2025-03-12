! This is a test from the NCEPLIBS-sp project.
!
! This test tests the sptez() and sptezv() subrroutines.
!
! Kyle Gerheiser
program test_sptezv
  use sp_mod
  use iso_fortran_env, only: real64
  implicit none

#if(LSIZE==d)
  integer,parameter:: iromb=0,maxwv=7
  integer,parameter:: idrtg=4,idrte=0,imax=16,jmaxg=8,jmaxe=17
  real(real64) :: MAX_DIFF = 1e-9
  
  call test_scalar(iromb,maxwv,idrtg,imax,jmaxg)
  call test_scalar(iromb,maxwv,idrte,imax,jmaxe)
  call test_vector(iromb,maxwv,idrtg,imax,jmaxg)
  call test_vector(iromb,maxwv,idrte,imax,jmaxe)
  
  call test_scalar(0,126,4,256,128)
  call test_scalar(0,126,0,256,257)
  call test_vector(0,126,4,256,128)
  call test_vector(0,126,0,256,257)

contains
  
  subroutine test_scalar(iromb,maxwv,idrt,imax,jmax)
    implicit none
    integer,intent(in):: iromb,maxwv,idrt,imax,jmax
    real(real64) :: wave((maxwv+1)*((iromb+1)*maxwv+2)/2*2)
    real(real64) :: wave2((maxwv+1)*((iromb+1)*maxwv+2)/2*2)
    real(real64) :: grid(imax,jmax)
    real(real64) :: avg_diff
    wave=1d0
    wave(2:2*maxwv+2:2)=0d0
    call sptez(iromb,maxwv,idrt,imax,jmax,wave,grid,+1)
    call sptez(iromb,maxwv,idrt,imax,jmax,wave2,grid,-1)
    avg_diff = sqrt(sum((wave2-wave)**2)/size(wave))

    print *, "avg_diff = ", avg_diff

    if (avg_diff > MAX_DIFF) then
       print *, "average difference > MAX_DIFF: ", avg_diff, " > ", MAX_DIFF
       error stop 
    endif
 
  end subroutine test_scalar
  
  subroutine test_vector(iromb,maxwv,idrt,imax,jmax)
    implicit none
    integer,intent(in):: iromb,maxwv,idrt,imax,jmax
    real(real64) :: waved((maxwv+1)*((iromb+1)*maxwv+2)/2*2)
    real(real64) :: wavez((maxwv+1)*((iromb+1)*maxwv+2)/2*2)
    real(real64) :: waved2((maxwv+1)*((iromb+1)*maxwv+2)/2*2)
    real(real64) :: wavez2((maxwv+1)*((iromb+1)*maxwv+2)/2*2)
    real(real64) :: gridu(imax,jmax)
    real(real64) :: gridv(imax,jmax)
    real(real64) :: avg_diff
    waved=1d0
    waved(2:2*maxwv+2:2)=0d0
    waved(1)=0d0
    wavez=1d0
    wavez(2:2*maxwv+2:2)=0d0
    wavez(1)=0d0
    call sptezv(iromb,maxwv,idrt,imax,jmax,waved,wavez,gridu,gridv,+1)
    call sptezv(iromb,maxwv,idrt,imax,jmax,waved2,wavez2,gridu,gridv,-1)
    avg_diff = sqrt((sum((waved2-waved)**2)+sum((wavez2-wavez)**2))/(2*size(waved)))

    print *, "avg_diff = ", avg_diff
    
    if (avg_diff > MAX_DIFF) then
       print *, "average difference > MAX_DIFF: ", avg_diff, " > ", MAX_DIFF
       error stop 
    endif
    
  end subroutine test_vector
#endif
  
end program test_sptezv
