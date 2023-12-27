module my_little_magfie
    implicit none

contains
    ! Compute physical components of B in cylindrical coordinates x = (r,phi,z).
    subroutine eval_field_B(x, B)
        real(8), intent(in) :: x(3)
        real(8), intent(inout) :: B(3)
        real(8) :: dummy

        call my_field(x(1),x(2),x(3),B(1),B(2),B(3),dummy,dummy,dummy, &
            dummy,dummy,dummy,dummy,dummy,dummy)
    end subroutine eval_field_B

    ! Compute physical components of B and A in cylindrical coordinates x = (r,phi,z).
    ! This routine derives from GORILLA vector_potential_rphiz (there it's covariant).
    subroutine eval_field_B_and_A(x, B, A)
        use field_eq_mod, only : rtf,btf,psif

        real(8), intent(in) :: x(3)
        real(8), intent(inout) :: B(3), A(3)
        real(8) :: dummy, Bmod

        call my_field(x(1),x(2),x(3),B(1),B(2),B(3),dummy,dummy,dummy, &
                   dummy,dummy,dummy,dummy,dummy,dummy)

        bmod=sqrt(B(1)**2+B(2)**2+B(3)**2)

        A(1)=0.d0
        A(2)=psif/x(1)
        A(3)=-rtf*btf*log(x(1))

        ! TODO: Add perturbation field
    end subroutine eval_field_B_and_A

    subroutine my_field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
                    ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

    use field_mod, only : icall, ipert, iequil, ampl
    use inthecore_mod, only : incore
    use libneo_kinds, only : real_kind

    real(kind=real_kind), intent(in) :: r, z
    real(kind=real_kind) :: p,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ         &
                        ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
    real(kind=real_kind) :: rm,zm,Brc,Bpc,Bzc,dBrdRc,dBrdpc,dBrdZc &
                        ,dBpdRc,dBpdpc,dBpdZc,dBzdRc,dBzdpc,dBzdZc

    if(icall .eq. 0) then
        call read_field_input
        icall = 1
        print *, 'Perturbation field', ipert, 'Ampl', ampl
    end if

    call stretch_coords(r,z,rm,zm)

    if(iequil.eq.0) then
        call set_zero(Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,   &
                    dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
    else
        call field_eq(rm,p,zm,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
                ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
    end if

    if(ipert.gt.0) then
        ! vacuum perturbation coil field
        incore=-1

        call my_field_c(rm,p,zm,Brc,Bpc,Bzc,dBrdRc,dBrdpc,dBrdZc,   &
                        dBpdRc,dBpdpc,dBpdZc,dBzdRc,dBzdpc,dBzdZc)

        call add_scaled(  &
            Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ,  &
            Brc,Bpc,Bzc,dBrdRc,dBrdpc,dBrdZc,dBpdRc,dBpdpc,dBpdZc,dBzdRc,dBzdpc,  &
            dBzdZc,ampl)

    end if

    end subroutine my_field

    ! ===========================================================================
    subroutine my_field_c(rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ  &
                    ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

    use field_c_mod, only : icall_c
    use libneo_kinds, only : real_kind

    implicit none

    real(kind=real_kind) :: rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ  &
                        ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ

    !-------first call: read data from disk-------------------------------
    if(icall_c .lt. 1) then
        call read_field_input_c
        icall_c = 1
    end if
    !------- end first call ----------------------------------------------

    call field_divfree(rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ    &
                        ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

    end subroutine my_field_c


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
