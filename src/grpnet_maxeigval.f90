!   grpnet_maxeigval.f90 - maximum eigenvalue of real symmetric matrix
!   Nathaniel E. Helwig (helwig@umn.edu)
!   Department of Psychology and School of Statistics
!   University of Minnesota
!   Date: 2023-07-10
!
! INPUTS
!   A = real symmetric matrix (N x N)
!   N = number of rows/columns of A
! OUTPUT:
!   MAXEV = maximum eigenvalue of A

SUBROUTINE grpnet_maxeigval(A, N, MAXEV)
    IMPLICIT NONE
    INTEGER N
    DOUBLE PRECISION A(N,N), MAXEV, X(N), OLDEV
    X = 1.0D0 / SQRT(DBLE(N))
    MAXEV = 1.0D0
    OLDEV = 0.0D0
    DO WHILE (ABS(MAXEV - OLDEV) >= 1E-8)
        OLDEV = MAXEV
        X = MATMUL(A, X)
        MAXEV = SQRT(SUM(X**2))
        IF (MAXEV > 0.0D0) X = X / MAXEV
    END DO
END SUBROUTINE
