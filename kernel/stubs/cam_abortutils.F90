!> Standalone replacement for CAM's cam_abortutils, used by the SCRAM box kernel.
module cam_abortutils

  implicit none
  private
  public :: endrun

contains

  !> Print the message and stop with a non-zero status.
  subroutine endrun(msg)

    character(len=*), intent(in), optional :: msg

    if (present(msg)) then
       write(0,'(A,A)') 'SCRAM ENDRUN: ', trim(msg)
       write(*,'(A,A)') 'SCRAM ENDRUN: ', trim(msg)
    else
       write(0,'(A)') 'SCRAM ENDRUN'
    end if
    flush(0)
    flush(6)
    error stop 2

  end subroutine endrun

end module cam_abortutils
