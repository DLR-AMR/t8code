! This file is part of t8code.
! t8code is a C library to manage a collection (a forest) of multiple
! connected adaptive space-trees of general element classes in parallel.
!
! Copyright (C) 2026 the developers
!
! t8code is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! t8code is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with t8code; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

! This file contains a module with an example callback of the Fortran interface
! that refines and coarsenes some elements based on their centroids' coordinates.
module t8_fortran_example_adapt_mod
    use t8_fortran_interface_mod

contains

      ! Example callback for adapting based on the centroid coordinates x, y, z.
      function example_fortran_adapt_by_coordinates_callback(x, y, z, is_family) &
            & result(ret) bind(c, name="example_fortran_adapt_by_coordinates_callback")
            use, intrinsic :: ISO_C_BINDING, only: c_int, c_double
            real(c_double), value :: x
            real(c_double), value :: y
            real(c_double), value :: z
            integer(c_int), value :: is_family
            integer(c_int) :: ret
            ret = 0
            if(x < 0.5d0 .and. y < 0.5d0 .and. z < 0.5d0) then
                  ret = 1
            elseif(x > 0.5d0 .and. y > 0.5d0) then
                  if(is_family.ne.0) then
                        ret = -1
                  endif
            endif

      end function example_fortran_adapt_by_coordinates_callback

end module t8_fortran_example_adapt_mod