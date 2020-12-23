program wave_function

    implicit none
    
    ! Calclation is done in Atomic units

    real(kind=8):: E, hbar, omega, m_e, epsilon, k, h, psi_0, psi_1, x_end
    integer :: n, grid, i, j 
    real, ALLOCATABLE:: x(:) , y(:), f(:) , psi(:), zeta(:)
    grid=100
    x_end = 10
    h = x_end/grid
    

   

    n = 0 ! Ground State
    hbar  = 1
    m_e = 1 
    k = 2 ! Here Energy is of th ground state is set to be 1 a.u.

    omega = (k/m_e)**(1/2)
    E = (n+1/2)*hbar*omega

    ! Transformation 
    epsilon = E/(hbar*omega) 
    ALLOCATE( x(1:grid) , y(1:grid), f(1:grid) , psi(1:grid), zeta(1:grid))

    psi_0 = 0
    psi_1 = 0.001
    
    psi(1) = psi_0
    psi(2) = psi_1

    x(1) = 0

    do i = 1, grid

            x(i+1) = x(i) + h

            zeta(i) = ((m_e*omega)/hbar)**(1/2)*x(i)
            f(i) = -2*(epsilon - zeta(i)**2/2)
            y(i) = 1 - (h)**2/(12*f(i))

    end do 

    do j = 1, grid
        
        psi(j+2) = (psi(j+1)*(12-10*y(j+1)) - y(j)*psi(j))/y(j+2)
        print '(3e16.8, 8X, 3e16.16)',x(j), psi(j)
        print *, " "

    end do




end program wave_function
