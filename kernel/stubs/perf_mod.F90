!> Standalone replacement for CESM's perf_mod timers, used by the SCRAM box kernel.
!> Timers are no-ops; the box driver measures wall time itself.
module perf_mod

  implicit none
  private
  public :: t_startf, t_stopf

contains

  subroutine t_startf(name)
    character(len=*), intent(in) :: name
  end subroutine t_startf

  subroutine t_stopf(name)
    character(len=*), intent(in) :: name
  end subroutine t_stopf

end module perf_mod
