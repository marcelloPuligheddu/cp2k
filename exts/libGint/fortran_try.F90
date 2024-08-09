program test1

use libGint, only: libgint_init
implicit none

write (*,*) "start"

call libgint_init()

write (*,*) "end"

end program test1
