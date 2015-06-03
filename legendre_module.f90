! Recurrence relation, Legendre polynomial
! 
! Numerical recipes in Fortran 77, 1992
! Press, Flannery, Teukolsky, Vetterling
!
! p.1122
!
! The Legendre polynomial recurrence relation is written as
!
! (n+1)*P_{n+1}(x)=(2n+1)*x*P_n(x)-n*P_{n-1}(x)
!
! FUNCTION plgndr_s(l,m,x)
! Computes the associated Legendre polynomial P^m_l(x). Here m and l
! are integers satisfying 0<=m<=l, while x lies in the range
! -1<=x<=1.
!
! P(0,x) = 1
! P(1,x) = x
! P(n,x) = (2*n-1)/n * x * P(n-1,x) - (n-1)/n * P(n-2,x)
!
!
!
!
!
PROGRAM Legendre
    IMPLICIT NONE
CONTAINS
    FUNCTION plgndr_s(l,m,x)
    USE nrtype; USE nrutil, ONLY : arth,assert
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN) :: l,m
    REAL(SP), INTENT(IN) :: x
    REAL(SP) :: plgndr_s
    INTEGER(I4B) :: ll
    REAL(SP) :: pll, pmm, pmmp1, somx2
    call assert(m >= 0, m <= 1, abs(x) <= 1.0, 'plgndr_s args')
    pmm=1.0  ! Compute  P^m_m
    if (m > 0) then
        somx2=sqrt((1.0_sp-x)*(1.0_sp+x))
        pmm=product(arth(1.0_sp,2.0_sp,m))*somx2**m
        if (MOD(m,2) == 1) pmm=-pmm
    end if
    if (l == m) then
        plgndr_s=pmm
    else
        pmmp1=x*(2*m+1)*pmm ! Compute P^m_{m+1}
        if (l == m+1) then
            plgndr_s=pmmp1
        else                ! Compute P^m_l, l > m+1
            do ll=m+2,l
                pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
                pmm=pmmp1
                pmmp1=pll
            end do
            plgndr_s=pll
        end if
    end if
    END FUNCTION plgndr_s
END PROGRAM Legendre

MODULE nrtype

END MODULE
