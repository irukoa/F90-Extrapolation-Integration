program example_extrapolation_integration

  use extrapolation_integration

  integer, parameter :: dp = 8

  !Arrays in memory layout.
  real(kind=dp), allocatable :: rescalar(:), revector(:, :)
  complex(kind=dp), allocatable :: cescalar(:), cevector(:, :)

  real(kind=dp), allocatable :: rresultv(:)
  complex(kind=dp), allocatable :: cresultv(:)

  !dim = 1 + 1 vector arrays.
  !1st dim stores the point index in memory layout.
  !2nd dim represents vector freedoms in array layout.
  !This type of array would be a prototype of Wannier90-s standard working array.
  real(kind=dp), allocatable :: rev(:, :)
  complex(kind=dp), allocatable :: cev(:, :)

  integer :: l1 = 33, & !Discretization in dim = 1.
             l2 = 33, & !Discretization in dim = 2.
             l3 = 33, & !Discretization in dim = 3.
             v1 = 1 !Vector freedoms.

  real(kind=dp) :: int_bounds(6)

  integer :: a1, a2, a3, u1, count, info
  real(kind=dp) :: x1, x2, x3

  !Integral bounds in sequential order.
  int_bounds = (/0.0_dp, 2.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 2.0_dp/)

  !Allocation of array layout arrays.
  allocate (rev(l1*l2*l3, -v1:v1))
  allocate (cev(l1*l2*l3, -v1:v1))

  !Data obtention in the specified mesh.
  count = 1
  do a1 = 1, l1
    x1 = int_bounds(1) + (int_bounds(2) - int_bounds(1))*real(a1 - 1, dp)/real(l1 - 1, dp)
    do a2 = 1, l2
      x2 = int_bounds(3) + (int_bounds(4) - int_bounds(3))*real(a2 - 1, dp)/real(l2 - 1, dp)
      do a3 = 1, l3
        x3 = int_bounds(5) + (int_bounds(6) - int_bounds(5))*real(a3 - 1, dp)/real(l3 - 1, dp)

        do u1 = lbound(rev, 2), ubound(rev, 2)
          cev(count, u1) = z(x1, real(u1, dp))*z(x2, real(u1, dp))*z(x3, real(u1, dp))
        enddo
        count = count + 1

      enddo
    enddo
  enddo
  rev = real(cev, dp)

  !Contract to memory layout.
  allocate (rescalar(l1*l2*l3), revector(l1*l2*l3, size(rev(1, :))))
  allocate (cescalar(l1*l2*l3), cevector(l1*l2*l3, size(rev(1, :))))

  do count = 1, l1*l2*l3
    call shrink_array(cev(count, :), cevector(count, :), info)
  enddo
  revector = real(cevector, dp)
  cescalar = cevector(:, 1)
  rescalar = revector(:, 1)

  !At this point, cevector(:, 1) contains the values of the integrand:
  !cos(x_1)*exp(sin(-x_1)) + i*cos(x_1)*exp(sin(-2*x_1))*
  !cos(x_2)*exp(sin(-x_2)) + i*cos(x_2)*exp(sin(-2*x_3))*
  !cos(x_3)*exp(sin(-x_3)) + i*cos(x_2)*exp(sin(-2*x_3)).
  !We will integrate this function from 0 to 2 in all dimensions, wich gives
  !-0.0481480 +0.352825 i.

  allocate (rresultv(size(revector(1, :))))
  allocate (cresultv(size(cevector(1, :))))
  call integral_extrapolation(cevector, (/l1, l2, l3/), int_bounds, cresultv, info)
  rresultv = real(cresultv, dp)
  if (info .eq. 1) then
    print *, "Extrapolation integration in", l1, "x", l2, "x", l3, "mesh."
    print *, "Result = ", cresultv(1)
    print *, "True result = -0.0481480 +0.352825 i"
  elseif (info .eq. 0) then
    print *, "Regular rectangle integration in", l1, "x", l2, "x", l3, "mesh."
    print *, "Result = ", cresultv(1)
    print *, "True result = -0.0481480 +0.352825 i"
  endif

contains

  !Test function.
  function z(x, v) result(u)

    real(kind=dp), intent(in) :: x, v

    complex(kind=dp) :: u

    u = cmplx(cos(x)*exp(sin(v*x)), cos(x)*exp(sin(2*v*x)), dp)

  end function z

end program example_extrapolation_integration
