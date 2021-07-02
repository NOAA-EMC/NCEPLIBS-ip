 program ipolates_driver

!-----------------------------------------------------------------
! Interpolate a global lat/lon grid of scalars to several
! grids of various projections using all ipolates 
! interpolation options.
!
! The ipolates series of routines is threaded, so it is compiled
! and run with threads.
!-----------------------------------------------------------------

 use omp_lib
 use get_input_data

 implicit none

 integer :: tid

!$OMP PARALLEL PRIVATE(TID)
 tid=omp_get_thread_num()
 print*,'- HELLO WORLD FROM THREAD: ',tid
!$OMP END PARALLEL

 call degrib_input_data

 call interp

 print*,"- NORMAL TERMINATION"

 stop
 end program ipolates_driver
