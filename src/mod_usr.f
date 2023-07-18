module mod_usr
  use mod_hd
  implicit none
  double precision :: rhoj, eta, vj

contains

  subroutine usr_init()

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_special_bc    => specialbound_usr
    usr_refine_grid   => specialrefine_grid

    call hd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    
    hd_gamma=1.4d0
    ! hd_gamma=1.d0
    ! rhoj=hd_gamma
    rhoj=1.d0
    if(iprob==1)then
        eta=10.d0
    endif
    if(iprob==2)then
        eta=1.d0
    endif
    if(iprob==3)then
        eta=0.01d0
    endif
    ! vj=100.d0
    vj = 5.212261824420211d0 ! conversion equation : v_j*cm_to_kpc/s_to_t 

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,&
     ixmax1,ixmax2,w,x)

    ! initialize one grid 

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
       ixmax1,ixmax2
    double precision, intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)

    where(dabs(x(ixmin1:ixmax1,ixmin2:ixmax2,1))<0.015d0.and.x(ixmin1:ixmax1,&
       ixmin2:ixmax2,2)<0.015d0)
       w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)=rhoj
       w(ixmin1:ixmax1,ixmin2:ixmax2,mom(1))=zero
       w(ixmin1:ixmax1,ixmin2:ixmax2,mom(2))=rhoj*vj
       w(ixmin1:ixmax1,ixmin2:ixmax2,e_)=one/(hd_gamma-one)+&
          0.5d0*rhoj*vj**2.0d0
    else where
       w(ixmin1:ixmax1,ixmin2:ixmax2,rho_) = rhoj/eta
       w(ixmin1:ixmax1,ixmin2:ixmax2,e_) = one/(hd_gamma-one)
       w(ixmin1:ixmax1,ixmin2:ixmax2,mom(1)) = zero
       w(ixmin1:ixmax1,ixmin2:ixmax2,mom(2)) = zero
    end where

  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,iB,w,x)

    ! special boundary types, user defined

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, iB
    double precision, intent(in) :: qt, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)
    integer :: ixOIntmin1,ixOIntmin2,ixOIntmax1,ixOIntmax2, ix2

   select case(iB)
     ! implementation of special bottom boundary
     case(3)
      ! extrapolate all primitives continuously, and in jet region: fix jet
      !
      ! first switch internal zone above boundary zone to primitive variables
      ixOIntmin1=ixOmin1;ixOIntmin2=ixOmin2;ixOIntmax1=ixOmax1
      ixOIntmax2=ixOmax2;
      ixOIntmin2=ixOmax2+1
      ixOIntmax2=ixOmax2+1
      call hd_to_primitive(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOIntmin1,&
         ixOIntmin2,ixOIntmax1,ixOIntmax2,w,x)
      ! extrapolate primitives, first everywhere on boundary
      do ix2 = ixOmin2,ixOmax2
         w(ixOmin1:ixOmax1,ix2,rho_)  = w(ixOmin1:ixOmax1,ixOmax2+1,rho_)
         w(ixOmin1:ixOmax1,ix2,mom(1))= w(ixOmin1:ixOmax1,ixOmax2+1,mom(1))
         w(ixOmin1:ixOmax1,ix2,mom(2))= w(ixOmin1:ixOmax1,ixOmax2+1,mom(2))
         w(ixOmin1:ixOmax1,ix2,e_)    = w(ixOmin1:ixOmax1,ixOmax2+1,e_)
      enddo
      ! in jet zone: fix all primitives to the jet values
      where(dabs(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1))<0.015d0)
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)=rhoj
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1))=zero
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(2))=vj
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=one
      endwhere
      ! switch to conservative variables in internal zone
      call hd_to_conserved(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOIntmin1,&
         ixOIntmin2,ixOIntmax1,ixOIntmax2,w,x)
      ! switch to conservative variables in ghost cells
      call hd_to_conserved(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,w,x)
    case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine specialbound_usr

  subroutine specialrefine_grid(igrid,level,ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
     ixmin1,ixmin2,ixmax1,ixmax2,qt,w,x,refine,coarsen)

    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

    ! you must set consistent values for integers refine/coarsen:

    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen

    integer, intent(in) :: igrid, level, ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
        ixmin1,ixmin2,ixmax1,ixmax2
    double precision, intent(in) :: qt, w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       1:nw), x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndim)
    integer, intent(inout) :: refine, coarsen

    ! always refine the jet inlet zone
    if (minval(dabs(x(ixmin1:ixmax1,ixmin2:ixmax2,&
       1))) < 0.1.and.minval(dabs(x(ixmin1:ixmax1,ixmin2:ixmax2,&
       2))) < 0.1) refine=1

  end subroutine specialrefine_grid

end module mod_usr
