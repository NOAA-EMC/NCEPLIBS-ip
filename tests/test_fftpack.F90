PROGRAM test_fftpack
#if(LSIZE==4)
  use test_fftpack_mod_4
#elif(LSIZE==D)
  use test_fftpack_mod_d
#elif(LSIZE==8)
  use test_fftpack_mod_8
#endif
  implicit none

  call run_fftpack_tests()

END PROGRAM test_fftpack