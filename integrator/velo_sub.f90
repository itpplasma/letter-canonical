module velo_sub
    use, intrinsic :: iso_fortran_env, only: dp => real64


    abstract interface
        subroutine magfie_p(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
            import :: dp
            real(dp), intent(in) :: x(3)
            real(dp), intent(out) :: bmod, sqrtg
            real(dp), dimension(3), intent(out) :: bder, hcovar, hctrvr, hcurl
        end subroutine magfie_p
    end interface

    abstract interface
        subroutine elefie_p(x, derphi)
            import :: dp
            real(dp), intent(in) :: x(3)
            real(dp), intent(out) :: derphi(3)
        end subroutine elefie_p
    end interface

    contains

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine velo(tau,z,vz,rmu,ro0,magfie,elefie)
    !  Computes the components of the 5-D velocity          -  vz
    !  for given set of phase space coordinates             -  z.
    !
    !  Warning !!!  The dimensionless time is used (see below)
    !
    !  Order of the coordinates is the following:
    !  z(i) = x(i)  for i=1,2,3     - spatial coordinates with real
    !                                 dimension of the general covariant
    !                                 space coordinate system
    !  z(4) = p                     - momentum  module normalized to
    !                                 thermal momentum and sqrt(2);
    !  z(5) = alambd                - cosine of the pitch-angle
    !
    !  Input parameters:
    !            formal:  tau    -   dimensionless time: tau=sqrt(2*T/m)*t
    !                     z      -   see above
    !  Output parameters:
    !            formal:  vz     -   see above
    !
    !
    !  Extra parameters:  rmu    -   inverse relativistic temperature
    !                     ro0    -   Larmor radius for the reference
    !                                magnetic field and temperature:
    !                                ro0=sqrt(2*T/m)/(e*B_ref/m*c)
    !                     magfie -   magnetic field routine
    !                     elefie -   electric field routine

        integer :: i

        real(dp), intent(in)  :: tau, z
        real(dp), intent(out) :: vz

        real(dp), intent(in)  :: rmu, ro0
        procedure(magfie_p) :: magfie
        procedure(elefie_p), optional :: elefie

        real(dp) x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
        real(dp) derphi
        real(dp) p,alambd,p2,ovmu,gamma2,gamma,ppar,vpa,coala
        real(dp) rmumag,rovsqg,rosqgb,rovbm
        real(dp) a_phi,a_b,a_c,hstar
        real(dp) s_hc,hpstar,phidot,blodot,bra
        real(dp) pardeb

        dimension z(5), vz(5)
        dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
        dimension derphi(3)
        dimension a_phi(3),a_b(3),a_c(3),hstar(3)

        do 1 i=1,3
            x(i)=z(i)
    1    continue

    ! in magfie: x(i)   - set of 3 curvilinear space coordinates (input)
    !            bmod   - dimensionless magnetic field module: bmod=B/B_ref
    !            sqrtg  - Jacobian of space coordinates (square root of
    !                     metric tensor
    !            bder   - derivatives of logarithm of bmod over space coords
    !                     (covariant vector)
    !            hcovar - covariant components of the unit vector along
    !                     the magnetic field
    !            hctrvr - contravariant components of the unit vector along
    !                     the magnetic field
    !            hcurl  - contravariant components of the curl of this vector
    !
        call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)

    ! TODO: error handling magfie
    !      if(ierrfield.ne.0) then
    !        vz=0.d0
    !        return
    !      endif
    ! in elefie: x(i)   - space coords (input, see above)
    !            derphi - derivatives of the dimensionless electric potential
    !                     phihat=e*phi/T over space coords (covar. vector)

        if (present(elefie)) then
            call elefie(x,derphi)
        else
            derphi = 0d0
        end if

        p=z(4)
        alambd=z(5)

        p2=p*p
        ovmu=2.d0/rmu
        gamma2=p2*ovmu+1.d0
        gamma=dsqrt(gamma2)
        ppar=p*alambd
    ! vpa - dimensionless parallel velocity: vpa=v_parallel/sqrt(2*T/m)
        vpa=ppar/gamma
        coala=(1.d0-alambd**2)
    ! rmumag - magnetic moment
        rmumag=.5d0*p2*coala/bmod

        rovsqg=ro0/sqrtg
        rosqgb=.5d0*rovsqg/bmod
        rovbm=ro0/bmod

        a_phi(1)=(hcovar(2)*derphi(3)-hcovar(3)*derphi(2))*rosqgb
        a_b(1)=(hcovar(2)*bder(3)-hcovar(3)*bder(2))*rovsqg
        a_phi(2)=(hcovar(3)*derphi(1)-hcovar(1)*derphi(3))*rosqgb
        a_b(2)=(hcovar(3)*bder(1)-hcovar(1)*bder(3))*rovsqg
        a_phi(3)=(hcovar(1)*derphi(2)-hcovar(2)*derphi(1))*rosqgb
        a_b(3)=(hcovar(1)*bder(2)-hcovar(2)*bder(1))*rovsqg

        s_hc=0.d0
        do i=1,3
            a_c(i)=hcurl(i)*rovbm
            s_hc=s_hc+a_c(i)*hcovar(i)
            hstar(i)=hctrvr(i)+ppar*a_c(i)
        enddo
        hpstar=1.d0+ppar*s_hc

    ! velocities in the coordinate space
    !
    ! phidot - derivative of the dmls el. potential over dmls time
    ! blodot - derivative of the logarith of the mag. field module over dmls time
        phidot=0.d0
        blodot=0.d0
        do i=1,3
            bra=vpa*hstar(i)+a_phi(i)+a_b(i)*rmumag/gamma
            vz(i)=bra/hpstar
            phidot=phidot+vz(i)*derphi(i)
            blodot=blodot+vz(i)*bder(i)
        enddo

    ! velocities in the phase space

        vz(4)=-0.5d0*gamma*phidot/p
        vz(5)=-(0.5d0*coala/hpstar)*(sum(hstar*derphi)/p                 &
                + p*sum(hstar*bder)/gamma+alambd*sum(a_phi*bder))

    end subroutine velo

end module velo_sub
