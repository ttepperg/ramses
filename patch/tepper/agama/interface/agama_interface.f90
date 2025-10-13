! Author:
!   Thor Tepper-García
! About:
!   Provides a Fortran interface to the functions defined in agama_wrapper.f90
! Created:
!   12 OCT 2025
!
!   Last Modified:
!
!
! Notes:
!   This file is intended to be used in conjunction with agama_wrapper.cpp. !   These two files demonstrate how to add functionality to an existing AGAMA !   installation without the need to modify its source and recompile.

module agama_interface
  use iso_c_binding
  implicit none
  private
  public :: agama_initfromfile_wrapper, agama_delete

  interface

    subroutine agama_initfromfile_wrapper(verbose, c_obj, filename)
      use iso_c_binding
      integer(c_int), value :: verbose
      character(len=8) :: c_obj
      character(len=*)::filename
    end subroutine agama_initfromfile_wrapper

    subroutine agama_delete(verbose, c_obj)
      use iso_c_binding
      integer(c_int), value :: verbose
      character(len=8) :: c_obj
    end subroutine agama_delete

  end interface

end module agama_interface
