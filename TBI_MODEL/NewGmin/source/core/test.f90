program test
  use basin_hopping
  implicit none 
  
  call allocate_arrays()
  ! Okay, so this is fairly ghetto hack. 
  ! Some functions require a subroutine to be called at the start of the run,
  ! others do not. To avoid putting a init() wrapper in all functions we add 
  ! a logical parameter to the data file. 
  !
  ! if this is false then the optimised compiler will ignore this subroutine
  ! call and so it need not exist.
  !
print *, 'test here'
  if (run_init) then 
    call init() 
    write(*,*) "IO finished"
  end if

  write(*,*) "IO Opened"
  call open_io()
  write(*,*) "IO Closed"
  call basin_hop()
  X(1:(N/2))=ANGLES
  X((N/2)+1:N)=LENGTHS
OPEN(1,FILE='coords.perturb')
WRITE(1,'(F20.10)') X
CLOSE(1)
 print *, 'closing io now'
  call close_io()
end program
