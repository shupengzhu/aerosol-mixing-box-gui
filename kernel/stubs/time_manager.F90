!> Standalone replacement for CAM's time_manager, used by the SCRAM box kernel.
!> Only get_nstep is needed (RDB diagnostics); the box driver sets the step.
module time_manager

  implicit none
  private
  public :: get_nstep, box_set_nstep

  integer, save :: box_nstep = 0

contains

  integer function get_nstep()
    get_nstep = box_nstep
  end function get_nstep

  subroutine box_set_nstep(n)
    integer, intent(in) :: n
    box_nstep = n
  end subroutine box_set_nstep

end module time_manager
