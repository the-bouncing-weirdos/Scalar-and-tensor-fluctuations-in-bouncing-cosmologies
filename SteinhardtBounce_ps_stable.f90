program evol; implicit none;
logical break
real y(2)
real, parameter :: dt = 1.d1
integer  j, i



do i = 1,6
    y(1) = 31.9
    y(2) = 0.175 - 0.01*(i-1)
    break= .true.
    do while (break)
        call gl8_background(y,dt)
        write (*,*) y(1), y(2)
        if (y(2)< 0.12 .or. isNaN(y(1)) .or. y(1)>63.0 .or. y(2)>0.3) then
            break = .false.
        end if
    end do
write (*,*) ''; write (*,*) ''
end do

do i = 1,10
    y(1) = 31.9
    y(2) = 0.175 + 0.015*i
    break= .true.
    do while (break)
        call gl8_background(y,dt)
        write (*,*) y(1), y(2)
        if (y(2)< 0.12 .or. isNaN(y(1)) .or. y(1)>63.0 .or. y(2)>0.3) then
            break = .false.
        end if
    end do
write (*,*) ''; write (*,*) ''
end do

do i = 1,10
    y(1) = 32.0 + 0.5*i
    y(2) = 0.125
    break= .true.
    do while (break)
        call gl8_background(y,dt)
        write (*,*) y(1), y(2)
        if (y(2)< 0.12 .or. isNaN(y(1)) .or. y(1)>63.0 .or. y(2)>0.3) then
            break = .false.
        end if
    end do
write (*,*) ''; write (*,*) ''
end do

do i = 1,12
y(1) = 66.0 - 0.5*i
y(2) = 0.1467
break= .true.
do while (break)
call gl8_background(y,dt)
write (*,*) y(1), y(2)
if (y(2)< 0.12 .or. isNaN(y(1)) .or. y(2)>0.3) then
break = .false.
end if
end do
write (*,*) ''; write (*,*) ''
end do


do i = 1,20
    y(1) = 63.3102
    y(2) = 0.18 + 0.005*i
    break= .true.
    do while (break)
        call gl8_background(y,dt)
        write (*,*) y(1), y(2)
        if (y(2)< 0.12 .or. isNaN(y(1)) .or. y(1)>63.32 .or. y(2)>0.3) then
            break = .false.
        end if
    end do
write (*,*) ''; write (*,*) ''
end do

!end do


contains


subroutine evalf(y, dydx)
! y = (/ phi_0, mu_0, h_0, a_0, L_0, L_0, 1/R /)
! dydx = (/ phi dot, phi dotdot, h dot, a dot, L dot, L dot, S_h dot /)
real y(2), dydx(2)

! time
dydx(1) = y(2)*den(y)

! scale factor
dydx(2) = (-kprime(y)*y(2)**2/2.d0-75.d-2*qprime(y)*y(2)**4-75.d-2*q(y)*y(2)**6-3.d0*(k(y)+q(y)*y(2)**2+&
15.d-1*y(2)**4)*hubble(y)*y(2))

end subroutine evalf


! -------------------------------------------------------------------------
! 8th order implicit Gauss-Legendre integrator
subroutine gl8_background(y, dt)
integer, parameter :: s = 4, n = 2
real y(n), g(n,s), dt; integer i, k

! Butcher tableau for 8th order Gauss-Legendre method
real, parameter :: a(s,s) = reshape((/ &
0.869637112843634643432659873054998518Q-1, -0.266041800849987933133851304769531093Q-1, &
0.126274626894047245150568805746180936Q-1, -0.355514968579568315691098184956958860Q-2, &
0.188118117499868071650685545087171160Q0,   0.163036288715636535656734012694500148Q0,  &
-0.278804286024708952241511064189974107Q-1,  0.673550059453815551539866908570375889Q-2, &
0.167191921974188773171133305525295945Q0,   0.353953006033743966537619131807997707Q0,  &
0.163036288715636535656734012694500148Q0,  -0.141906949311411429641535704761714564Q-1, &
0.177482572254522611843442956460569292Q0,   0.313445114741868346798411144814382203Q0,  &
0.352676757516271864626853155865953406Q0,   0.869637112843634643432659873054998518Q-1 /), (/s,s/))
real, parameter ::   b(s) = (/ &
0.173927422568726928686531974610999704Q0,   0.326072577431273071313468025389000296Q0,  &
0.326072577431273071313468025389000296Q0,   0.173927422568726928686531974610999704Q0  /)

! iterate trial steps
g = 0.0; do k = 1,16
g = matmul(g,a)
do i = 1,s
call evalf(y + g(:,i)*dt, g(:,i))
end do
end do

! update the solution
y = y + matmul(g,b)*dt
end subroutine gl8_background

function hubble(y)
real y(2), hubble
hubble = y(2)**3/2.d0-sqrt(y(2)**6/4.d0+k(y)*y(2)**2/6.d0+q(y)*y(2)**4/4.d0)
end function hubble

function den(y)
real y(2), den
den = k(y) + 3.d0*q(y)*y(2)**2+6.d0*hubble(y)*y(2)+15.d-1*y(2)**4
end function den

function k(y)
real y(2), k, ts, a, gamma, Theta, F, hubb, w, omega, logf
F = 9.d-5; Theta = 46.d-4; ts = 5.d-1; gamma = -44.d-4; hubb = 3.d-5
a = (-2.d0*gamma)**(1.d0/3.d0)/Theta
logf = log(y(1)/a)
w = exp(F/Theta**2*(logf-Theta*ts)**2)
omega = w*Theta**6*y(1)**3+hubb*(4.d0*F*(logf)**2-4.d0*Theta*F*ts*logf+&
2.d0*Theta**4*y(1)**3*logf-2.d0*Theta**2)
k = - (12.d0*logf**2*hubb**2-3.d0*w*(omega-hubb*Theta**4*y(1)**3*logf))/(w*w*Theta**4*y(1)**2)
end function k

function q(y)
real y(2), q, ts, a, gamma, Theta, F, hubb, w, omega, logf
F = 9.d-5; Theta = 46.d-4; ts = 5.d-1; gamma = -44.d-4; hubb = 3.d-5
a = (-2.d0*gamma)**(1.d0/3.d0)/Theta
logf = log(y(1)/a)
w = exp(F/Theta**2*(logf-Theta*ts)**2)
omega = w*Theta**6*y(1)**3+hubb*(4.d0*F*(logf)**2-4.d0*Theta*F*ts*logf+&
2.d0*Theta**4*y(1)**3*logf-2.d0*Theta**2)
q = (12.d0*logf**2*hubb**2-2.d0*w*(omega+hubb*Theta**4*y(1)**3*logf))/(w*w*Theta**6*y(1)**4)
end function q

function kprime(y)
real y(2), kprime, ts, a, gamma, Theta, F, hubb, w, omega, logf, gaussf
F = 9.d-5; Theta = 46.d-4; ts = 5.d-1; gamma = -44.d-4; hubb = 3.d-5
a = (-2.d0*gamma)**(1.d0/3.d0)/Theta
logf = log(y(1)/a)
gaussf = exp(-2.d0*F*(logf-Theta*ts)**2/Theta**2)

kprime = (3.d0*(Theta**8 + (8.d0*hubb**2*logf*(-Theta**2 + Theta*(-2.d0*F*ts + Theta)*logf +&
        2.d0*F*logf**2))/(exp((2.d0*F*(-(ts*Theta) + logf)**2)/Theta**2)*y(1)**3) + (hubb*(Theta**3*&
        (-8.d0*F*ts + 4.d0*Theta + Theta**3*y(1)**3) + logf*(4.d0*F*Theta**2*(3.d0 + 2.d0*ts*(-(F*ts)+&
        Theta)) + Theta**5*(2.d0*F*ts + Theta)*y(1)**3 - 2.d0*F*logf*(Theta*(-8.d0*F*ts + 4.d0*Theta +&
        Theta**3*y(1)**3) + 4.d0*F*logf))))/(exp((F*(-(ts*Theta) + logf)**2)/Theta**2)*y(1)**3)))/Theta**6

end function kprime

function qprime(y)
real y(2), qprime, ts, a, gamma, Theta, F, hubb, w, omega, logf, gaussf
F = 9.d-5; Theta = 46.d-4; ts = 5.d-1; gamma = -44.d-4; hubb = 3.d-5
a = (-2.d0*gamma)**(1.d0/3.d0)/Theta
logf = log(y(1)/a)
gaussf = exp(-2.d0*F*(logf-Theta*ts)**2/Theta**2)

qprime = (2.d0*(y(1)**3 + (12.d0*hubb**2*logf*(Theta**2 - 2.d0*logf*(Theta*(-(F*ts) + Theta) + F*logf)))/&
         (exp((2.d0*F*(-(ts*Theta) + logf)**2)/Theta**2)*Theta**8) + (hubb*(Theta**3*(8.d0*F*ts - 8.d0*Theta&
         - 3*Theta**3*y(1)**3) + logf*(4.d0*F*Theta**2*(-3.d0 + 2.d0*F*ts**2 - 4.d0*ts*Theta) +3.d0*Theta**5*(&
        -2*F*ts + Theta)*y(1)**3 + 2.d0*F*logf*(Theta*(-8.d0*F*ts + 8.d0*Theta +3.d0*Theta**3*y(1)**3) + 4.d0*F*logf))))/&
        (exp((F*(-(ts*Theta) + logf)**2)/Theta**2)*Theta**8)))/y(1)**5
end function qprime


end
