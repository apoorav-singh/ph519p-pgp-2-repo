program wave_function

    implicit none
    
    ! Calclation is done in Atomic units

    real(kind=8):: E, hbar, omega, m_e, epsilon, k, h, psi_0, psi_1, x_end, x_beg,psi_l, psi_r, c
    integer :: n, grid, i1, i2, j1, j2, l 
    real, ALLOCATABLE :: x(:) , y(:), f(:) , psi(:)
    grid = 10000
    x_beg = -10
    x_end = 10
    h = (x_end-x_beg)/grid
    

   

    n = 100 ! Ground State
    hbar  = 1
    m_e = 1 
    k = 2 ! Here Energy is of the ground state is set to be 1 a.u.

    omega = (k/m_e)**(1/2)
    E = (n+1/2)*hbar*omega

    ! Transformation 
    epsilon = E/(hbar*omega) 
    ALLOCATE( x(1:grid) , y(1:grid), f(1:grid) , psi(1:grid))

    psi_0 = 0
    psi_1 = 0.001
    
    psi(1) = psi_0
    psi(2) = psi_1
    psi(grid) = psi_0
    psi(grid-1) = psi_1

    x(1) = x_beg
    x(grid) = x_end

    do i1 = 1, grid

            x(i1+1) = x(i1) + h

            f(i1) = -(2*m_e/hbar**2)*(E - 0.5*k*x(i1)**2)
            y(i1) = 1 - (h**2)/(12)*f(i1)
    
    end do 

    do j1 = 1, (grid/2-4)
        
        psi(j1+2) = (psi(j1+1)*(12-(10*y(j1+1)))-y(j1)*psi(j1))/y(j1+2)

    end do
    !psi_r = psi(grid/2)
    do j2 = grid, (grid/2), -1

        psi(j2-2) = (psi(j2-1)*(12-(10*y(j2-1)))-y(j2)*psi(j2))/y(j2-2)

        !psi(grid-1-j2) = (psi(grid-j2)*(12-(10*y(grid-j2)))-y(grid-j2+1)*psi(grid-j2+1))/y(grid-1-j2)
    end do
    !psi_l = psi(grid/2)
    !c = psi_r/psi_l
    !psi(1:grid) = psi(1:grid)*c
    do  l = 0, grid

        print '(3e16.8, 8X, 3e16.8)',x(l), psi(l)

    end do

!print *, h



end program wave_function
