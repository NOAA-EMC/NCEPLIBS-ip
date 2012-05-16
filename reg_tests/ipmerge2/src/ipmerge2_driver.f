 program ipmerge2_driver

! test routine ipmerge2 which does not seem to be called
! by any iplib routines.  is this routine even used anymore?

 implicit none

 integer :: no, nf, m1, m2, mo

 logical*1, allocatable :: l1(:,:), l2(:,:), lo(:,:)
 
 real, allocatable :: f1(:,:), f2(:,:), fo(:,:)

 nf = 2  ! # fields
 no = 4  ! # points in each output field
 mo = no ! # first dimension of output field arrays

 allocate (fo(mo,nf))
 allocate (lo(mo,nf))

 m1 = no

 allocate (l1(m1,nf))
 allocate (f1(m1,nf))
 l1(1,1) = .true. 
 l1(2,1) = .true.
 l1(3,1) = .false.
 l1(4,1) = .false.
 f1(:,1) = 1.0 
 l1(1,2) = .false.
 l1(2,2) = .false.
 l1(3,2) = .true.
 l1(4,2) = .true.
 f1(:,2) = 3.0 

 m2 = no

 allocate (l2(m2,nf))
 allocate (f2(m2,nf))
 l2(1,:) = .true. 
 l2(2,:) = .false.
 l2(3,:) = .true.
 l2(4,:) = .false.
 f2(:,1) = 2.0
 f2(:,2) = 4.0

 print*,'CALL ROUTINE IPMERGE'
 call ipmerge2(no,nf,m1,l1,f1,m2,l2,f2,mo,lo,fo)

 print*,'LO FIELD1 ',lo(:,1)
 print*,'FO FIELD1 ',fo(:,1)
 print*,'LO FIELD2 ',lo(:,2)
 print*,'FO FIELD2 ',fo(:,2)

 print*,'NORMAL TERMINATION'

 stop

 end program ipmerge2_driver
