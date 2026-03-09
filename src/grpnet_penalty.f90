!   grpnet_penalty.f90 - group elastic net shrinkage penalty
!   Nathaniel E. Helwig (helwig@umn.edu)
!   Department of Psychology and School of Statistics
!   University of Minnesota
!   Date: 2023-06-20
!
! INPUTS
!   znorm = Euclidean norm of least squares estimate
!   penid = penalty id: 1 = LASSO, 2 = MCP, 3 = SCAD
!   penone = L1 penalty weight
!   pentwo = L2 penalty weight
!   gamma = extra tuning parameter (for MCP and SCAD)
! OUTPUT:
!   shrink = shrinkage factor in interval (0,1)

SUBROUTINE grpnet_penalty(znorm, penid, penone, pentwo, gamma, shrink)
    IMPLICIT NONE
    INTEGER penid
    DOUBLE PRECISION znorm, penone, pentwo, gamma, shrink
    IF (penid == 1) THEN
        shrink = MAX(0.0D0, 1.0D0 - penone / znorm)
        shrink = shrink / (1.0D0 + pentwo)
    ELSE IF (penid == 2) THEN
        IF (znorm <= gamma * penone * (1.0D0 + pentwo)) THEN
            shrink = MAX(0.0D0, 1.0D0 - penone / znorm)
            shrink = shrink / (1.0D0 + pentwo - 1.0D0 / gamma)
        ELSE
            shrink = 1.0D0 / (1.0D0 + pentwo)
        END IF
    ELSE IF (penid == 3) THEN
        IF (znorm <= (penone * (1.0D0 + pentwo) + penone)) THEN
            shrink = MAX(0.0D0, 1.0D0 - penone / znorm) / (1.0D0 + pentwo)
        ELSE IF (znorm <= gamma * penone * (1.0D0 + pentwo)) THEN
            shrink = MAX(0.0D0, 1.0D0 - (penone / znorm) * (gamma / (gamma - 1.0D0)) )
            shrink = shrink / (1.0D0 + pentwo - 1.0D0 / (gamma - 1.0D0) )
        ELSE
            shrink = 1.0D0 / (1.0D0 + pentwo)
        END IF
    END IF
END SUBROUTINE

