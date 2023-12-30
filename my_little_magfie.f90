module my_little_magfie
    implicit none

contains

    subroutine init
        use field_mod, only : icall, ipert, ampl
        use field_c_mod, only : icall_c
        use libneo_kinds, only : real_kind

        if (icall < 1) then
            call read_field_input
            icall = 1
            print *, 'Perturbation field', ipert, 'Ampl', ampl
        end if

        if (ipert > 0 .and. icall_c < 1) then
            call read_field_input_c
            icall_c = 1
        end if
    end subroutine init

    ! Compute physical components of B in cylindrical coordinates x = (r,phi,z).
    subroutine eval_field_B(x, B)
        real(8), intent(in) :: x(3)
        real(8), intent(inout) :: B(3)
        real(8) :: dummy

        call my_field(x(1),x(2),x(3),B(1),B(2),B(3),dummy,dummy,dummy, &
            dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy)
    end subroutine eval_field_B

    ! Compute physical components of B and A in cylindrical coordinates x = (r,phi,z).
    ! This routine derives from GORILLA vector_potential_rphiz (there it's covariant).
    subroutine eval_field_B_and_A(x, B, A)

        real(8), intent(in) :: x(3)
        real(8), intent(inout) :: B(3), A(3)
        real(8) :: dummy

        call my_field(x(1),x(2),x(3),B(1),B(2),B(3),dummy,dummy,dummy, &
                   dummy,dummy,dummy,dummy,dummy,dummy,A(1),A(2),A(3))
    end subroutine eval_field_B_and_A

    subroutine my_field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
                    ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ,Ar,Ap,Az)

    use field_mod, only : ipert, iequil, ampl
    use field_eq_mod, only : rtf,btf,psif
    use inthecore_mod, only : incore
    use libneo_kinds, only : real_kind

    real(kind=real_kind), intent(in) :: r, z, p
    real(kind=real_kind), intent(inout) :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ         &
                        ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ,Ar,Ap,Az
    real(kind=real_kind) :: rm,zm,Brc,Bpc,Bzc,dBrdRc,dBrdpc,dBrdZc &
                        ,dBpdRc,dBpdpc,dBpdZc,dBzdRc,dBzdpc,dBzdZc,Arc,Apc,Azc

    call stretch_coords(r,z,rm,zm)

    if(iequil.eq.0) then
        call set_zero(Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,   &
                    dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
        Ar=0.d0
        Ap=0.d0
        Az=0.d0
    else
        call field_eq(rm,p,zm,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
                ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

            Ar=0.d0
            Ap=psif/r
            Az=-rtf*btf*log(r)
    end if

    if(ipert.gt.0) then
        ! vacuum perturbation coil field
        incore=-1

        call my_field_divfree(rm,p,zm,Brc,Bpc,Bzc,dBrdRc,dBrdpc,dBrdZc,   &
                    dBpdRc,dBpdpc,dBpdZc,dBzdRc,dBzdpc,dBzdZc,Arc,Apc,Azc)

        call add_scaled(  &
            Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ,  &
            Brc,Bpc,Bzc,dBrdRc,dBrdpc,dBrdZc,dBpdRc,dBpdpc,dBpdZc,dBzdRc,dBzdpc,  &
            dBzdZc,ampl)

        Ar = Ar + ampl*Arc
        Ap = Ap + ampl*Apc
        Az = Az + ampl*Azc

    end if

    end subroutine my_field


    subroutine my_field_divfree(r,phi,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ    &
                            ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ,Ar,Ap,Az)

    use bdivfree_mod, only : nr,nz,ntor,icp,ipoint,hr,hz,pfac,&
        & rpoi,zpoi,apav
    use libneo_kinds, only : complex_kind, real_kind

    implicit none

    double precision, intent(in) :: r, phi, z
    double precision, intent(inout) :: Br, Bp, Bz, dBrdR, dBrdp, dBrdZ,   &
                                     dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ,   &
                                     Ar, Ap, Az

    integer :: n,ierr
    real(kind=real_kind) :: f,fr,fz,frr,frz,fzz
    real(kind=real_kind) :: delBr,delBz,delBp
    real(kind=real_kind) :: deldBrdR,deldBrdp,deldBrdZ
    real(kind=real_kind) :: deldBpdR,deldBpdp,deldBpdZ
    real(kind=real_kind) :: deldBzdR,deldBzdp,deldBzdZ
    complex(kind=complex_kind) :: expon,anr,anz,anr_r,anr_z,anz_r,anz_z
    complex(kind=complex_kind) :: anr_rr,anr_rz,anr_zz,anz_rr,anz_rz,anz_zz
    complex(kind=complex_kind) :: pfac_imaginary


    Bp = 0.0d0
    ! TODO: Solve this with a vector potential being Az = -int Bp dr with gauge Ar=0
    !
    ! call spline(nr,nz,rpoi,zpoi,hr,hz,icp,rbpav_coef,ipoint,r,z,         &
    !             f,fr,fz,frr,frz,fzz,ierr)
    ! Bp = f/r
    ! dBpdR = fr/r - Bp/r
    ! dBpdZ = fz/r
    ! dBpdp = 0.d0

    call spline(nr,nz,rpoi,zpoi,hr,hz,icp,apav,ipoint,r,z,         &
                f,fr,fz,frr,frz,fzz,ierr)

    Br=-fz/r
    Bz=fr/r
    dBrdR=fz/r**2-frz/r
    dBrdZ=-fzz/r
    dBzdR=-fr/r**2+frr/r
    dBzdZ=frz/r
    dBrdp=0.d0
    dBzdp=0.d0

    Ar = 0.d0
    Ap = 0.d0
    Az = 0.d0

    do n=1,ntor
        call spline_vector_potential_n(n, r, z, anr,anz,anr_r,anr_z,anz_r,anz_z, &
        anr_rr,anr_rz,anr_zz,anz_rr,anz_rz,anz_zz)

        pfac_imaginary = cmplx(0.0, pfac, kind=complex_kind)

        expon=exp(cmplx(0.d0,n*pfac*phi, kind=complex_kind))
        delBr=2.d0*real(pfac_imaginary*n*anz*expon/r, kind=real_kind)
        delBz=-2.d0*real(pfac_imaginary*n*anr*expon/r, kind=real_kind)
        delBp=2.d0*real((anr_z-anz_r)*expon, kind=real_kind)
        deldBrdR=-delbr/r+2.d0*real(pfac_imaginary*n*anz_r*expon/r, kind=real_kind)
        deldBrdZ=2.d0*real(pfac_imaginary*n*anz_z*expon/r, kind=real_kind)
        deldBrdp=-2.d0*real((pfac*n)**2*anz*expon/r, kind=real_kind)
        deldBzdR=-delbz/r-2.d0*real(pfac_imaginary*n*anr_r*expon/r, kind=real_kind)
        deldBzdZ=-2.d0*real(pfac_imaginary*n*anr_z*expon/r, kind=real_kind)
        deldBzdp=2.d0*real((pfac*n)**2*anr*expon/r, kind=real_kind)
        deldBpdR=2.d0*real((anr_rz-anz_rr)*expon, kind=real_kind)
        deldBpdZ=2.d0*real((anr_zz-anz_rz)*expon, kind=real_kind)
        deldBpdp=2.d0*real(pfac_imaginary*n*(anr_z-anz_r)*expon, kind=real_kind)

        br=br+delbr
        bz=bz+delbz
        bp=bp+delbp
        dBrdR=dBrdR+deldBrdR
        dBrdZ=dBrdZ+deldBrdZ
        dBrdp=dBrdp+deldBrdp
        dBzdR=dBzdR+deldBzdR
        dBzdZ=dBzdZ+deldBzdZ
        dBzdp=dBzdp+deldBzdp
        dBpdR=dBpdR+deldBpdR
        dBpdZ=dBpdZ+deldBpdZ
        dBpdp=dBpdp+deldBpdp

        Ar = Ar + 2.d0*real(anr*expon, kind=real_kind)
        Az = Az + 2.d0*real(anz*expon, kind=real_kind)
    end do

    end subroutine my_field_divfree


    subroutine read_field_input_c
        use libneo_kinds, only : real_kind
        use field_c_mod, only : ntor,nr,np,nz,icftype,rmin,pmin,zmin,rmax,pmax,zmax

        real(kind=real_kind), dimension(:,:,:), allocatable :: Br,Bp,Bz

        print *,'coils: file type = ',icftype
        if(icftype.eq.4) then
            call read_sizes(nr,np,nz)
        else
            print *,'field_c: unknown coil file type'
            stop
        end if
        allocate(Br(nr,np,nz),Bp(nr,np,nz),Bz(nr,np,nz))

        call read_field4(nr,np,nz,rmin,rmax,pmin,pmax,zmin,zmax,Br,Bp,Bz)

        print *,'coils: nr,np,nz = ',nr,np,nz
        print *,'coils: rmin,rmax = ',rmin,rmax
        print *,'coils: zmin,zmax = ',zmin,zmax
        print *,'coils: pmin,pmax = ',pmin,pmax

        call vector_potentials(nr,np,nz,ntor,rmin,rmax,pmin,pmax,zmin,zmax,br,bp,bz)

        deallocate(Br,Bp,Bz)
    end subroutine read_field_input_c


end module my_little_magfie
