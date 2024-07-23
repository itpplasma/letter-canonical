module magfie_tok
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use magfie, only: field_t

    implicit none

    complex(8), parameter :: imun = dcmplx(0.0d0, 1.0d0)
    character(len=1024) :: input_file = 'field_divB0.inp'

    type, extends(field_t) :: tok_field_t
        type(FourierPerturbation), allocatable :: perturbations(:)
        contains
            procedure :: init_magfie
            procedure :: compute_bfield
            procedure :: compute_abfield

            procedure :: add_perturbation
    end type tok_field_t

    type :: FourierPerturbation
        integer :: mpol
        integer :: ntor
        complex(8) :: Amn(3)
    end type FourierPerturbation

    contains

    subroutine init_magfie(self)
        use field_mod, only : icall, ipert, ampl
        use field_eq_mod, only : icall_eq
        use field_c_mod, only : icall_c
        use libneo_kinds, only : real_kind
        use field_sub, only : read_field_input

        class(tok_field_t), intent(in) :: self

        if (icall < 1) then
            call read_field_input(input_file)
            print *, 'Perturbation field', ipert, 'Ampl', ampl
            call init_stretch_coords
            icall = 1
            call init_field_eq
            icall_eq = 1
        end if

        if (ipert > 0 .and. icall_c < 1) then
            call read_field_input_c
            icall_c = 1
        end if
    end subroutine init_magfie


    subroutine compute_abfield(self, R, phi, Z, AR, Aphi, AZ, BR, Bphi, BZ)
        class(tok_field_t), intent(in) :: self
        real(dp), intent(in) :: R, phi, Z
        real(dp), intent(out) :: AR, Aphi, AZ, BR, Bphi, BZ

        real(dp) :: dummy

        call my_field(R, phi, Z, BR, Bphi, BZ,dummy,dummy,dummy, &
                   dummy,dummy,dummy,dummy,dummy,dummy,AR, Aphi, AZ)

        call perturb_afield(self, R, phi, Z, AR, Aphi, AZ)
        call perturb_bfield(self, R, phi, Z, BR, Bphi, BZ)
    end subroutine compute_abfield


    subroutine compute_bfield(self, R, phi, Z, BR, Bphi, BZ)
        class(tok_field_t), intent(in) :: self
        real(dp), intent(in) :: R, phi, Z
        real(dp), intent(out) :: BR, Bphi, BZ

        real(dp) :: dummy

        call my_field(R, phi, Z, BR, Bphi, BZ,dummy,dummy,dummy, &
            dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy)

        call perturb_bfield(self, R, phi, Z, BR, Bphi, BZ)
    end subroutine compute_bfield


    subroutine perturb_afield(self, R, phi, Z, AR, Aphi, AZ)
        class(tok_field_t), intent(in) :: self
        real(dp), intent(in) :: R, phi, Z
        real(dp), intent(out) :: AR, Aphi, AZ

        real(dp) :: theta
        complex(8) :: expfac, ARmn, Aphimn, AZmn
        integer :: i

        if (.not. allocated(self%perturbations)) return

        ARmn = 0.d0
        APhimn = 0.d0
        AZmn = 0.d0

        theta = atan2(Z, R)

        do i = 1, size(self%perturbations)
            associate(p => self%perturbations(i))
                expfac = exp(dcmplx(0.d0, p%mpol*theta + p%ntor*phi))
                ARmn = ARmn + p%Amn(1)*expfac
                Aphimn = Aphimn + p%Amn(2)*expfac
                AZmn = AZmn + p%Amn(3)*expfac
            end associate
        end do

        AR = AR + real(ARmn, kind=8)
        Aphi = Aphi + real(Aphimn, kind=8)
        AZ = AZ + real(AZmn, kind=8)
    end subroutine perturb_afield


    subroutine perturb_bfield(self, R, phi, Z, BR, Bphi, BZ)
        class(tok_field_t), intent(in) :: self
        real(dp), intent(in) :: R, phi, Z
        real(dp), intent(inout) :: BR, Bphi, BZ

        real(dp) :: theta
        complex(8) :: expfac, BRmn, Bphimn, BZmn
        integer :: i

        if (.not. allocated(self%perturbations)) return

        BRmn = 0.d0
        BPhimn = 0.d0
        BZmn = 0.d0

        theta = atan2(Z, R)


        do i = 1, size(self%perturbations)
            associate(p => self%perturbations(i))
                expfac = exp(dcmplx(0.d0, p%mpol*theta + p%ntor*phi))
                BRmn = BRmn + imun *expfac * &
                    (p%Amn(3)*p%ntor/R - p%mpol*R/(R**2 + Z**2))

                Bphimn = Bphimn + imun/(R**2 + Z**2)*expfac * &
                    p%mpol*(p%Amn(1)*R + p%Amn(3)*Z)

                BZmn = BZmn + expfac/R * &
                    (p%Amn(2) - imun*p%ntor*p%Amn(1)) &
                    - imun*Z/(R**2 + Z**2)*expfac * p%mpol*p%Amn(2)
            end associate
        end do

        BR = BR + real(BRmn, kind=8)
        Bphi = Bphi + real(Bphimn, kind=8)
        BZ = BZ + real(BZmn, kind=8)
    end subroutine perturb_bfield


    subroutine add_perturbation(self, mpol, ntor, Amn)
        class(tok_field_t), intent(inout) :: self
        integer, intent(in) :: mpol, ntor
        complex(8), intent(in) :: Amn(3)

        type(FourierPerturbation) :: new_perturbation

        new_perturbation%mpol = mpol
        new_perturbation%ntor = ntor
        new_perturbation%Amn = Amn

        if (.not. allocated(self%perturbations)) then
            allocate(self%perturbations(1))
            self%perturbations(1) = new_perturbation
        else
            self%perturbations = [self%perturbations, new_perturbation]
        end if

    end subroutine add_perturbation


    subroutine my_field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
                    ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ,Ar,Ap,Az)

    use field_mod, only : ipert, iequil, ampl
    use field_eq_mod, only : rtf,btf,psif
    use inthecore_mod, only : incore
    use libneo_kinds, only : real_kind
    use field_sub

    real(kind=real_kind), intent(in) :: r, p, z
    real(kind=real_kind), intent(inout) :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ         &
                        ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
    real(kind=real_kind), intent(inout), optional :: Ar,Ap,Az

    real(kind=real_kind) :: rm,zm,Brc,Bpc,Bzc,dBrdRc,dBrdpc,dBrdZc &
                        ,dBpdRc,dBpdpc,dBpdZc,dBzdRc,dBzdpc,dBzdZc,Arc,Apc,Azc

    call stretch_coords(r,z,rm,zm)

    if(iequil.eq.0) then
        call set_zero(Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,   &
                    dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
        if (present(Ar)) Ar=0.d0
        if (present(Ap)) Ap=0.d0
        if (present(Az)) Az=0.d0
    else
        call field_eq(rm,p,zm,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
                ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

            if (present(Ar)) Ar=0.d0
            if (present(Ap)) Ap=psif/r
            if (present(Az)) Az=-rtf*btf*log(r)
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

        if (present(Ar)) Ar = Ar + ampl*Arc
        if (present(Ap)) Ap = Ap + ampl*Apc
        if (present(Az)) Az = Az + ampl*Azc

    end if

    end subroutine my_field


    subroutine my_field_divfree(r,phi,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ    &
                            ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ,Ar,Ap,Az)

        use bdivfree_mod, only : nr,nz,ntor,icp,ipoint,hr,hz,pfac,&
            & rpoi,zpoi,apav
        use libneo_kinds, only : complex_kind, real_kind

        implicit none

        real(dp), intent(in) :: r, phi, z
        real(dp), intent(inout) :: Br, Bp, Bz, dBrdR, dBrdp, dBrdZ,   &
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


    subroutine init_stretch_coords
        use input_files, only : iunit, convexfile
        use libneo_kinds, only : real_kind
        use math_constants, only : TWOPI

        integer, parameter :: nrzmx=100 ! possible max. of nrz
        integer i, j, nrz ! number of points "convex wall" in input file
        integer, parameter :: nrhotht=360
        integer :: iflag

        ! points "convex wall"
        real(kind=real_kind), dimension(0:1000) :: rad_w, zet_w
        real(kind=real_kind), dimension(:), allocatable :: rho_w, tht_w
        real(kind=real_kind), dimension(nrhotht) :: rho_wall, tht_wall ! polar coords

        real(kind=real_kind) R0,Rw, Zw, htht, a, b, dummy

        nrz = 0
        rad_w = 0.
        zet_w = 0.
        open(iunit,file=trim(convexfile))
        do i=1,nrzmx
            read(iunit,*,END=10)rad_w(i),zet_w(i)
            nrz = nrz + 1
        end do
    10  continue
        close(iunit)

        allocate(rho_w(0:nrz+1), tht_w(0:nrz+1))
        R0 = (maxval(rad_w(1:nrz)) +  minval(rad_w(1:nrz)))*0.5
        do i=1,nrz
            rho_w(i) = sqrt( (rad_w(i)-R0)**2 + zet_w(i)**2 )
            tht_w(i) = atan2(zet_w(i),(rad_w(i)-R0))
            if(tht_w(i) .lt. 0.) tht_w(i) = tht_w(i) + TWOPI
        end do

        ! make sure points are ordered according to tht_w.
        do
            iflag = 0
            do i=1,nrz-1
                if (tht_w(i) > tht_w(i+1)) then
                iflag = 1
                dummy = rad_w(i+1)
                rad_w(i+1) = rad_w(i)
                rad_w(i) = dummy
                dummy = zet_w(i+1)
                zet_w(i+1) = zet_w(i)
                zet_w(i) = dummy
                dummy = rho_w(i+1)
                rho_w(i+1) = rho_w(i)
                rho_w(i) = dummy
                dummy = tht_w(i+1)
                tht_w(i+1) = tht_w(i)
                tht_w(i) = dummy
                end if
            end do
            if (iflag == 0) exit
        end do
        rad_w(0) = rad_w(nrz)
        zet_w(0) = zet_w(nrz)
        tht_w(0) = tht_w(nrz) - TWOPI
        rho_w(0) = rho_w(nrz)
        rad_w(nrz+1) = rad_w(1)
        zet_w(nrz+1) = zet_w(1)
        tht_w(nrz+1) = tht_w(1) + TWOPI
        rho_w(nrz+1) = rho_w(1)

        htht = TWOPI/(nrhotht-1)
        Rw = 0.d0
        Zw = 0.d0
        do i=2,nrhotht
            tht_wall(i) = htht*(i-1)
            do j=0,nrz
                if(tht_wall(i).ge.tht_w(j) .and. tht_wall(i).le.tht_w(j+1)) then
                if( abs((rad_w(j+1) - rad_w(j))/rad_w(j)) .gt. 1.e-3) then
                    a = (zet_w(j+1) - zet_w(j))/(rad_w(j+1) - rad_w(j))
                    b = zet_w(j) - a*(rad_w(j) - R0)
                    Rw = b/(tan(tht_wall(i)) - a) + R0
                    Zw = a*(Rw - R0) + b
                else
                    a = (rad_w(j+1) - rad_w(j))/(zet_w(j+1) - zet_w(j))
                    b = rad_w(j)-R0 - a*zet_w(j)
                    Zw = b/(1./tan(tht_wall(i)) - a)
                    Rw = a*Zw + b + R0
                end if
                end if
            end do
            rho_wall(i) = sqrt((Rw-R0)**2 + Zw**2)
        end do
        tht_wall(1) = 0.
        rho_wall(1) = rho_wall(nrhotht)

    end subroutine init_stretch_coords


    subroutine init_field_eq
        use input_files, only : ieqfile
        use field_eq_mod, only : use_fpol,skip_read,nrad,nzet,icp,nwindow_r,&
          nwindow_z,psib,btf,rtf,hrad,hzet,psi_axis,psi_sep,&
          psi,psi0,splfpol,splpsi,rad,zet,imi,ima,jmi,jma,ipoint
        use libneo_kinds, only : real_kind
        use field_sub

        integer :: i

        if (.not. skip_read) then
        if (ieqfile == 0) then
            call read_dimeq0(nrad, nzet)
        elseif (ieqfile == 2) then
            call read_dimeq_west(nrad, nzet)
        else
            call read_dimeq1(nrad, nzet)
        end if
        allocate(rad(nrad), zet(nzet))
        allocate(psi0(nrad, nzet), psi(nrad, nzet))
        end if
        if (use_fpol) then
        if (.not. skip_read) then
            allocate(splfpol(0:5, nrad))
            call read_eqfile2(nrad, nzet, psi_axis, psi_sep, btf, rtf, &
                            splfpol(0, :), rad, zet, psi)
        end if
        psib = -psi_axis
        psi_sep = (psi_sep - psi_axis) * 1.d8
        splfpol(0, :) = splfpol(0, :) * 1.d6
        call spline_fpol
        else
        if (.not. skip_read) then
            if (ieqfile == 0) then
            call read_eqfile0(nrad, nzet, psib, btf, rtf, rad, zet, psi)
            elseif (ieqfile == 2) then
            call read_eqfile_west(nrad, nzet, psib, btf, rtf, rad, zet, psi)
            else
            call read_eqfile1(nrad, nzet, psib, btf, rtf, rad, zet, psi)
            end if
        end if
        end if

        ! Filtering:
        do i=1,nzet
        call window_filter(nrad,nwindow_r,psi(:,i),psi0(:,i))
        end do

        do i=1,nrad
        call window_filter(nzet,nwindow_z,psi0(i,:),psi(i,:))
        end do
        ! End filtering

        rad = rad*1.d2 ! cm
        zet = zet*1.d2 ! cm
        rtf = rtf*1.d2 ! cm
        psi = psi*1.d8
        psib= psib*1.d8
        btf = btf*1.d4

        psi=psi+psib

        hrad = rad(2) - rad(1)
        hzet = zet(2) - zet(1)

        ! rectangular domain:
        allocate( imi(nzet),ima(nzet),jmi(nrad),jma(nrad) )
        imi = 1
        ima = nrad
        jmi = 1
        jma = nzet

        !  Computation of the number of data in splpsi
        icp = 0
        do i=1,nzet
        if ( imi(i) .gt. 0 .and. ima(i) .gt. 0 ) then
            icp = icp + ima(i) - imi(i) + 1
        end if
        end do
        write(6,*) 'number of points in the table:  ',icp

        allocate( splpsi(6,6,icp), ipoint(nrad,nzet) )

        call s2dcut(nrad,nzet,hrad,hzet,psi,imi,ima,jmi,jma,icp,splpsi,ipoint)
    end subroutine init_field_eq


    subroutine read_field_input_c
        use libneo_kinds, only : real_kind
        use field_c_mod, only : ntor,nr,np,nz,icftype,rmin,pmin,zmin,rmax,pmax,zmax
        use field_sub

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

end module magfie_tok
