program basic1 
 implicit none
 real :: r,area
 real, parameter :: pi = 4.*atan(1.)

 r = 2.0
 area = pi*r**2

 print*, "Value of pi is",pi
 print*, "Circle whose radius is",r,"have area =",area 

end program basic1