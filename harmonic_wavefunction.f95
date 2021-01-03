program wave_function

    implicit none
    
    ! Calclation is done in Atomic units

    integer, parameter :: dp = selected_real_kind(14,200)
    real(dp):: E, hbar, omega, m_e, epsilon, k, h, psi_0, psi_1, x_end, x_beg,psi_l, psi_r, c, norm
    integer(dp) :: n, grid, i1, i2, j1, j2, l 
    real(dp), ALLOCATABLE :: x(:) , y(:), f(:) , psi(:), a(:)
    grid = 10000
    x_beg = -10.0_dp
    x_end = 10_dp
    h = (x_end-x_beg)/grid
    

   

    n = 10 ! Ground State
    hbar  = 1.0_dp
    m_e = 1.0_dp
    k = 2.0_dp ! Here Energy is of the ground state is set to be 1 a.u.

    omega = (k/m_e)**(0.5_dp)
    E = (n+0.5_dp)*hbar*omega

    ! Transformation 
    ALLOCATE( x(1:grid) , y(1:grid), f(1:grid) , psi(1:grid), a(1:grid))

    psi_0 = 1d-21
    
    
    

    x(0) = x_beg
    x(grid) = x_end

    do i1 = 0, grid

            
            f(i1) = -(2*m_e/hbar**2)*(E - 0.5*k*x(i1)**2)
            y(i1) = 1 - (h**2)/(12)*f(i1)
            x(i1+1) = x(i1) + h

    end do 

   
    psi(0) = psi_0
    psi(1) = (12.0_dp-10.0_dp*y(0))*psi(0)/(2.0_dp*y(1)) 
    norm = psi(0) 
    do j1 = 1, (grid/2)
        psi(j1+1) = (psi(j1)*(12.0_dp-(10.0_dp*y(j1)))-y(j1-1)*psi(j1-1))/y(j1+1)
        norm = norm + psi(j1)
    end do

    i1 = 1
    do j2 = grid, (grid/2+1), -1

        psi(j2) = (-1_dp)**n*psi(i1)
        i1 = i1 + 1

    end do

    do j2 = 0, grid

        a(j2) = psi(j2)/sqrt(2*norm)

    end do


    do  l = 0, grid

        print '(3e16.8, 4X, 3e16.8)',x(l), a(l)

    end do





end program wave_function