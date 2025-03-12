C> @file
C> Set number of cpus.
C>
C> ### Program History Log
C> Date | Programmer | Comments
C> -----|------------|---------
C> 94-08-19 | Iredell | Initial.
C> 98-11-09 | Vuong | Add doc>block and remove cray references.
C> 1998-12-18 | Iredell | IBM SMP version.
C> 2010-11-16 | Slovacek | Linux must have different call.
C> 2012-11-01 | Mirvis | Multi-threading on LINUX-IBM/TIDE.
C>
C> @author Iredell @date 94-08-19

C> Set number of CPUs - the number of processors over which
C> to parallelize.
C>
C> @return Number of CPUs assigned.
C>
C> @author Iredell @date 94-08-19
      FUNCTION NCPUS()
      INTEGER TID
C    Obtain thread number
#ifdef OPENMP
      INTEGER OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
!$OMP PARALLEL PRIVATE(TID)
      TID = OMP_GET_THREAD_NUM()
!      PRINT *, '...............thread # ', TID
      if (TID. eq. 0) then
        NCPUS=OMP_GET_NUM_THREADS()
!        PRINT *, 'totaly #------------------- of threads = ',NCPUS
      endif
!$OMP END PARALLEL
#else
      TID = 0
      NCPUS = 1
#endif
      RETURN
      END
