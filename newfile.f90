program wave_function

    implicit none
    
    ! Calclation is done in Atomic units

    real(kind=8):: E, hbar, omega, m_e, epsilon, zeta, k, h, psi_0, psi_1, x_end
    integer :: n, grid, i 
     real, ALLOCATABLE:: x(:) , y(:), f(:) , psi(:)
    !real, dimension(1:grid):: x,y,f,psi
    x_end = 10
    h = x_end/grid
    

   

    n = 0 ! Ground State
    grid=100
    hbar  = 1
    m_e = 1 
    k = 2 ! Here Energy is of th ground state is set to be 1 a.u.

    omega = (k/m_e)**(1/2)
    E = (n+1/2)*hbar*omega

    ! Transformation 
    epsilon = E/(hbar*omega) 
    ALLOCATE( x(1:grid) , y(1:grid), f(1:grid) , psi(1:grid))

    psi_0 = 0
    psi_1 = 0.001
    
    psi(1) = psi_0
    psi(2) = psi_1

    x(1) = 0.0

    do i = 2, grid
        
        x(i) = x(i-1) + h

        !zeta = ((m*omega)/hbar)**(1/2)*x(i-1)
        !f(i-1) = -2*(epsilon - zeta**2/2)
        !y(i-1) = 1 - (h)**2/(12*f(i-1))

        !zeta = ((m*omega)/hbar)**(1/2)*x(i)
        !f(i) = -2*(epsilon - zeta**2/2)
        !y(i) = 1 - (h)**2/(12*f(i))

        !zeta = ((m*omega)/hbar)**(1/2)*x(i+1)
        !f(i+1) = -2*(epsilon - zeta**2/2)
        !y(i+1) = 1 - (h)**2/(12*f(i+1))

        !psi(i+1) = (psi(i)(12-10*y(i)) - y(i-1)psi(i-1))/y(i+1)

        print *, x 
    end do 





end program wave_function
