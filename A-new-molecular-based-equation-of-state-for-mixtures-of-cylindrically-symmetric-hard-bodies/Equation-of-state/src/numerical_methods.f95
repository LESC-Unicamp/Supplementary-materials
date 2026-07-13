! ************************************************************************************************ !
!                                        NUMERICAL ROUTINES                                        !
! ************************************************************************************************ !
!            This subroutine contains all numerical routines used in the main program.             !
! ************************************************************************************************ !
! => AUTHOR:     Nathan Barros de Souza                                                            !
! => E-MAIL:     n264179@dac.unicamp.br                                                            !
! => SUPERVISOR: Luís Fernando Mercier Franco                                                      !
! ************************************************************************************************ !
!                       University of Campinas, Campinas - São Paulo, Brazil                       !
!                                         March 7th, 2024                                          !
! ************************************************************************************************ !
! Version 1.0.0                                                                                    !
! ************************************************************************************************ !

! ************************************************************************************************ !
! Calculates the factorial of a number                                                             !
! ************************************************************************************************ !
DOUBLE PRECISION FUNCTION Factorial( Number )

! Uses one module: global variables
USE GlobalVariables

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: Number  ! Number to calculate the factorial
INTEGER( Kind= Int64 ) :: Counter ! Counter

! Factorial of the number
Factorial = 1.D0
DO Counter = 1, Number
  Factorial = Factorial * DBLE( Counter )
END DO

RETURN

END FUNCTION Factorial