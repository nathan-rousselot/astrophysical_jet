module mod_usr
  use mod_hd
  implicit none
  double precision :: Rjet, Zjet, rhoj, eta, vj, M, n1, n2

contains

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ Rjet, Zjet

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine usr_params_read

  subroutine usr_init()

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_special_bc    => specialbound_usr
    usr_refine_grid   => specialrefine_grid
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output

    !call set_coordinate_system("cylindrical_2D")
    call hd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    
    hd_gamma=5.0d0/3.0d0
    ! rhoj=hd_gamma
    rhoj=1.0d0
    if(iprob==1)then
        eta=10.d0
    endif
    if(iprob==2)then
        eta=1.d0
    endif
    if(iprob==3)then
        eta=0.01d0
    endif
    vj=5.212261824420211d0
    !vj=4.0d0
    n1 = 8d0
    n2 = 6d0

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)

    ! initialize one grid 

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision :: R(ixG^S),Z(ixG^S)

    ! radial direction R/Rj
    R(ix^S)=x(ix^S,1)/0.1d0
    ! axial direction Z/Zj
    Z(ix^S)=x(ix^S,2)/0.1d0

    where(R(ix^S)<one .and. Z(ix^S)<one)
       w(ix^S,rho_)=rhoj
       w(ix^S,mom(1))=zero
       w(ix^S,mom(2))=rhoj*vj
       w(ix^S,p_)=(rhoj*vj**2.0d0)/hd_gamma !one/(hd_gamma-one)+0.5d0*rhoj*vj**2.0d0  !(rhoj*vj**2.0d0)/hd_gamma
    else where
       w(ix^S,rho_)=rhoj/(2.0d0*eta)  
       w(ix^S,p_)=(rhoj*vj**2.0d0)/hd_gamma
       w(ix^S,mom(1))=zero
       w(ix^S,mom(2))=zero
    end where

    !w(ix^S,mom(1)) = 0.0d0

    !where(dabs(x(ix^S,1))<=0.1d0.and.x(ix^S,2)<=0.1d0)
    !   !w(ix^S,rho_)=rhoj
    !   !w(ix^S,mom(1))=zero
    !   w(ix^S,mom(2))=vj
    !   w(ix^S,p_)= (rhoj*vj**2.0d0)/hd_gamma
    !else where
    !   !w(ix^S,rho_) = rhoj/(eta*2.0d0)
    !   w(ix^S,p_) = (rhoj*vj**2.0d0)/hd_gamma  ! one/(hd_gamma) ! -one)+0.5d0*rhoj*vj**2.0d0
    !   !w(ix^S,mom(1)) = zero
    !   w(ix^S,mom(2)) = zero
    !end where

  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixG^L,ixO^L,iB,w,x)

    ! special boundary types, user defined

    integer, intent(in) :: ixG^L, ixO^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    integer :: ixOInt^L, ix2
    double precision :: R(ixG^S),Z(ixG^S)

    R(ixO^S)=x(ixO^S,1)/0.1d0
    Z(ixO^S)=x(ixO^S,2)/0.1d0

   select case(iB)
     ! implementation of special bottom boundary
     case(3)
      ixOint^L=ixO^L;
      ixOIntmin2=ixOmax2+1
      ixOIntmax2=ixOmax2+1
      call hd_to_primitive(ixG^L,ixOInt^L,w,x)
      do ix2 = ixOmin2,ixOmax2
         w(ixOmin1:ixOmax1,ix2,rho_)  = w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,rho_)
         w(ixOmin1:ixOmax1,ix2,mom(1))= w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,mom(1))
         w(ixOmin1:ixOmax1,ix2,mom(2))= -w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,mom(2))
         w(ixOmin1:ixOmax1,ix2,p_)    = w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,p_)
      enddo
     !case(3)
      ! extrapolate all primitives continuously, and in jet region: fix jet
      !
      ! first switch internal zone above boundary zone to primitive variables
      ! extrapolate primitives, first everywhere on boundary

      ! in jet zone: fix all primitives to the jet values
      where(R(ixO^S)<=1.0d0)!.and.Z(ixO^S)<=1.0d0)
         w(ixO^S,rho_)=w(ixO^S,rho_)+(rhoj-w(ixO^S,rho_))/cosh(((dabs(x(ixO^S,1)-3.2d0))/(1.4d0))**n1)
         w(ixO^S,mom(1))=w(ixO^S,mom(1))+(zero-w(ixO^S,mom(1)))/cosh(dabs(x(ixO^S,1)-3.2d0)**n2)
         w(ixO^S,mom(2))=vj
         w(ixO^S,p_)=w(ixO^S,p_)+((rhoj*vj**2.0d0)/hd_gamma-w(ixO^S,p_))/cosh(dabs(x(ixO^S,1)-3.2d0)**n2)
      !else where
      !   w(ixO^S,rho_)=w(ixO^S,rho_)
      !   w(ixO^S,mom(1))=w(ixO^S,mom(1))
      !   w(ixO^S),mom(2))=-w(ixO^S,mom(2))
      !   w(ixO^S,p_)=w(ixO^S,p_)
      end where

      !where(dabs(x(ixO^S,1)-3.2d0)>0.1d0)
      !   w(ixO^S,rho_)=w(ixO^S,rho_)
      !   w(ixO^S,mom(1))=w(ixO^S,mom(1))
      !   w(ixO^S,mom(2))=-w(ixO^S,mom(2))
      !   w(ixO^S,p_)=w(ixO^S,p_)
      !end where
      ! switch to conservative variables in internal zone
      call hd_to_conserved(ixG^L,ixOInt^L,w,x)
      ! switch to conservative variables in ghost cells
      call hd_to_conserved(ixG^L,ixO^L,w,x)
    case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine specialbound_usr

  subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)

    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

    ! you must set consistent values for integers refine/coarsen:

    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen

    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen
    double precision :: R(ixG^S), Z(ixG^S)

    ! always refine the jet inlet zone
    !if (minval(dabs(x(ix^S,1))) <= 0.3.and.minval(dabs(x(ix^S,2))) <= 0.3) refine=1

    R(ix^S)=x(ix^S,1)/0.1d0
    Z(ix^S)=x(ix^S,2)/0.1d0
    if (any((R(ix^S) <= 3.0d0).and.(Z(ix^S) <= 3.0d0))) refine=1

  end subroutine specialrefine_grid
  
  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    
    w(ixO^S,nw+1)=dlog10(w(ixO^S,rho_))

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames

    varnames='log10rho'
  end subroutine specialvarnames_output

end module mod_usr

