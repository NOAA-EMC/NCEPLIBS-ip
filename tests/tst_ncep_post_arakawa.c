/**************************************************************
  Test to ensure the C interfaces created via the BIND(C)
  attributes for handling the ncep_post_arakawa flag are
  working correctly from here.
**************************************************************/

#include <stdio.h>
#include <stdbool.h>

extern bool ncep_post_arakawa; // C interface to the Fortran logical public save scalar, ncep_post_arakawa
extern void use_ncep_post_arakawa(); // C interfaces to Fortran subroutine use_ncep_post_arakawa
extern void unuse_ncep_post_arakawa(); // C interfaces to Fortran subroutine unuse_ncep_post_arakawa

int main()
{

   int testint = 10;
   printf("Default value of ncep_post_arakawa = %d\n", ncep_post_arakawa);
   printf("Initial value of testint = %d\n", testint);
   testint = (int)ncep_post_arakawa;
   printf("Value of testint, now should be Fortran logical ncep_post_arakawa = %d ... ", testint);
   if (testint == 10)
      return -1;
   printf("SUCCESS!\n");

   printf("Calling use_ncep_post_arakawa()\n");
   use_ncep_post_arakawa();
   printf("\tAfter call to use_ncep_post_arakawa(), value of ncep_post_arakawa = %d ... ", ncep_post_arakawa);
   if (ncep_post_arakawa == 0)
      return -1;
   printf("SUCCESS!\n");

   printf("Calling unuse_ncep_post_arakawa()\n");
   unuse_ncep_post_arakawa();
   printf("\tAfter call to unuse_ncep_post_arakawa(), value of ncep_post_arakawa = %d ... ", ncep_post_arakawa);
   if (ncep_post_arakawa == 1)
      return -1;
   printf("SUCCESS!\n");

   printf("Testing C interfaces to ncep_post_arakawa flags ... SUCCESS!\n");

   return 0;
}
