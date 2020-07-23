!WRF:model_layer:physics
!
!Usama Anber - Brookhaven National Lab.
!uanber@bnl.gov
!
!
!
module module_bl_annpbl
contains
!
!
!-------------------------------------------------------------------------------
!
   subroutine annpbl(pbl_input_mean, pbl_input_scale, pbl_diagnostic_mean,    &
                  pbl_diagnostic_scale,                                       &
                  state_mean, state_scale, sl_real_mean, rt_real_mean,        &
                  w_real_mean, sl_latent_to_real,                             &
                  rt_latent_to_real, w_latent_to_real,                        &
                  rcld_latent_to_real,   &
                  rrain_latent_to_real,                                       &
                  cld_latent_to_real,                                         &
                  pbl_encoder_W, pbl_encoder_b,                               &
                  pbl_hidden_W, pbl_hidden_b,                                 &
                  pbl_tend_decoder_W, pbl_tend_decoder_b, pbl_diag_decoder_W, &
                  pbl_diag_decoder_b,                                         &
                !Arrays Dimensions
                  n_input, n_diag, n_state, n_real, n_latent, n_code,         &
                  n_rcld, n_cld, n_rrain, n_sl, n_rt, n_w,                    &
                !Real WRF variables
                  dt,                                                         & !time step
                  th3d, qv3d, rthblten, rqvblten,rqcblten,w3d,                & !Theta, qv, and tend.
                  lh, hfx, tsk, swdnb, psfc,                                  & ! 2d var.
                  cldfra, swdnflx, qc3d, qr3d,                                & ! cld, sw rad, liquid, rain
                 ! znu,                                                        &
                   
                  ids,ide,jds,jde,kds,kde,                                    &
                  ims,ime,jms,jme,kms,kme,                                    &
                  its,ite, jts,jte, kts,kte,                                  &
                  !z,                                                          &
                  d_top,  d_mid, d_high )


!-----
! use module_radiation_driver
! use module_cumulus_driver
!-----


!-- th3d        3d potential temperature (k)
!-- qv3d        3d water vapor mixing ratio (kg/kg)
!-- qc3d        3d cloud mixing ratio (kg/kg)
!-- qi3d        3d ice mixing ratio (kg/kg)
!               (note: if P_QI<PARAM_FIRST_SCALAR this should be zero filled)
!-- p3d         3d pressure (pa)
!-- p3di        3d pressure (pa) at interface level
!-- pi3d        3d exner function (dimensionless)
!-- rr3d        3d dry air density (kg/m^3)
!-- rublten     u tendency due to
!               pbl parameterization (m/s/s)
!-- rvblten     v tendency due to
!               pbl parameterization (m/s/s)
!-- rthblten    theta tendency due to
!               pbl parameterization (K/s)
!-- rqvblten    qv tendency due to
!               pbl parameterization (kg/kg/s)
!-- rqcblten    qc tendency due to
!               pbl parameterization (kg/kg/s)
!-- rqiblten    qi tendency due to
!               pbl parameterization (kg/kg/s)
!-- cp          heat capacity at constant pressure for dry air (j/kg/k)
!-- g           acceleration due to gravity (m/s^2)
!-- rovcp       r/cp
!-- rd          gas constant for dry air (j/kg/k)
!-- rovg        r/g
!-- dz8w        dz between full levels (m)
!-- lh          latent heat of vaporization (W/m^2)
!-- rv          gas constant for water vapor (j/kg/k)
!-- psfc        pressure at the surface (pa)
!-- swdnflx     downward total shortwave flux 3D (W/m^2)
!-- xland       land mask (1 for land, 2 for water)
!-- hfx         upward heat flux at the surface (w/m^2)
!-- qfx         upward moisture flux at the surface (kg/m^2/s)
!-- swdnb       downward shortwave flux at the surface (W/m^2)
!-- dt          time step (s)
!-- rvovrd      r_v divided by r_d (dimensionless)
!-- ids         start index for i in domain
!-- ide         end index for i in domain
!-- jds         start index for j in domain
!-- jde         end index for j in domain
!-- kds         start index for k in domain
!-- kde         end index for k in domain
!-- ims         start index for i in memory
!-- ime         end index for i in memory
!-- jms         start index for j in memory
!-- jme         end index for j in memory
!-- kms         start index for k in memory
!-- kme         end index for k in memory
!-- its         start index for i in tile
!-- ite         end index for i in tile
!-- jts         start index for j in tile
!-- jte         end index for j in tile
!-- kts         start index for k in tile
!-- kte         end index for k in tile
!----------------------------------------------------------------------

  REAL,       INTENT(IN   ),OPTIONAL    ::     dt

  integer,  intent(in   )   ::      ids,ide, jds,jde, kds,kde,                &
                                     ims,ime, jms,jme, kms,kme,                &
                                     its,ite, jts,jte, kts,kte

  !REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN   ) :: z

  real,     dimension( ims:ime, kms:kme, jms:jme )                           , &
             intent(in   )   ::                                          qv3d, &
                                                                         qc3d, &
                                                                         qr3d, &
                                                                         th3d, &
                                                                         cldfra
                                                          !           rthraten

   real,     dimension( ims:ime, kms:kme, jms:jme )                          , &
              intent(inout) ::                                       rthblten, &
                                                                     rqvblten, &
                                                                     rqcblten


   real,     dimension( ims:ime, kms:kme+1, jms:jme )                       , &
              intent(in)::                                               w3d



   real,     dimension( ims:ime, jms:jme )                                   , &
             intent(in   )   ::                                         swdnb, &
                                                                          hfx, &
                                                                           lh, &
                                                                           tsk,&
                                                                           psfc



  INTEGER, INTENT(IN) ::  d_top, d_mid, d_high
 
  REAL,       DIMENSION( ims:ime, kms:kme ,jms:jme )           , &
               INTENT(INOUT)    ::                           swdnflx

   integer, parameter :: ndiff = 1
   integer ::  i,j,k
   real,     dimension( its:ite, kts:kte*ndiff )  ::                 rqvbl2dt, &
                                                                         qv2d

 

 integer , intent(in) :: n_input  !=33
 integer , intent(in) :: n_diag   !=34
 integer , intent(in) :: n_state  !=17
 integer , intent(in) :: n_real   !=20 should = d_top
 integer , intent(in) :: n_latent !=20
 integer , intent(in) :: n_code   !=128

!
 integer , intent(in) :: n_rcld  !=10
 integer , intent(in) :: n_cld   !=10
 integer , intent(in) :: n_rrain !=2
 integer , intent(in) :: n_sl    !=9
 integer , intent(in) :: n_rt    !=8
 integer , intent(in) :: n_w     !=9

 real, intent(in), dimension(n_input)          :: pbl_input_mean,pbl_input_scale
 real, intent(in), dimension(n_diag)           :: pbl_diagnostic_mean,pbl_diagnostic_scale
 real, intent(in), dimension(n_state)          :: state_mean, state_scale
 real, intent(in), dimension(n_real)           :: sl_real_mean, rt_real_mean,w_real_mean

 real, intent(in), dimension(n_sl, n_latent)   :: sl_latent_to_real
 real, intent(in), dimension(n_rt, n_latent)   :: rt_latent_to_real
 real, intent(in), dimension(n_w, n_latent)    :: w_latent_to_real

 real, intent(in), dimension(n_rcld, n_latent) :: rcld_latent_to_real
 real, intent(in), dimension(n_rrain,n_latent) :: rrain_latent_to_real
 real, intent(in), dimension(n_cld, n_latent)  :: cld_latent_to_real


 real, intent(in), dimension(n_input, n_code)  :: pbl_encoder_W
 real, intent(in), dimension(n_code)           :: pbl_hidden_b
 real, intent(in), dimension(n_code, n_code)   :: pbl_hidden_W
 real, intent(in), dimension(n_code)           :: pbl_encoder_b
 real, intent(in), dimension(n_code, n_state)  :: pbl_tend_decoder_W
 real, intent(in), dimension(n_state)          :: pbl_tend_decoder_b
 real, intent(in), dimension(n_code, n_diag)   :: pbl_diag_decoder_W
 real, intent(in), dimension(n_diag)           :: pbl_diag_decoder_b




!---------------------------------------------------------------------
!
! populate 3D arrays using 2D/column PBL.
!

   qv2d(its:ite,:) = 0.0

!
 do j = jts,jte
!
      do k = kts,kte
        do i = its,ite
          qv2d(i,k) = qv3d(i,k,j)
          qv2d(i,k+kte) = qc3d(i,k,j)
          !if(present(rqiblten)) qv2d(i,k+kte+kte) = qi3d(i,k,j)
        enddo
      enddo
!
      call annpbl2d(pbl_input_mean=pbl_input_mean,pbl_input_scale=pbl_input_scale, &
                    pbl_diagnostic_mean=pbl_diagnostic_mean,           &
                    pbl_diagnostic_scale=pbl_diagnostic_scale,&
                    state_mean=state_mean, state_scale=state_scale,&
                    sl_real_mean=sl_real_mean, rt_real_mean=rt_real_mean,&
                    w_real_mean=w_real_mean,&
                    sl_latent_to_real=sl_latent_to_real,&
                    rt_latent_to_real=rt_latent_to_real,&
                    w_latent_to_real=w_latent_to_real,&
                    rcld_latent_to_real=rcld_latent_to_real,&
                    rrain_latent_to_real=rrain_latent_to_real,&
                    cld_latent_to_real=cld_latent_to_real,&
                    pbl_encoder_W=pbl_encoder_W,&
                    pbl_encoder_b=pbl_encoder_b,&
                    pbl_hidden_W= pbl_hidden_W, pbl_hidden_b=pbl_hidden_b,&
                    pbl_tend_decoder_W=pbl_tend_decoder_W, &
                    pbl_tend_decoder_b=pbl_tend_decoder_b, &
                    pbl_diag_decoder_W=pbl_diag_decoder_W, &
                    pbl_diag_decoder_b=pbl_diag_decoder_b, &
!----
                    n_input=n_input, n_diag=n_diag, n_state=n_state, n_real=n_real,    &
                    n_latent=n_latent, n_code=n_code,n_rcld=n_rcld,n_cld=n_cld,&
                    n_rrain=n_rrain, n_sl=n_sl, n_rt=n_rt, n_w=n_w,&
!----
                    J=j,wx=w3d(ims,kms,j)                            &
                    ,thx=th3d(ims,kms,j)&
                    ,qvx=qv2d(its,kts), qcx = qc3d(its,kts,j),qrx = qr3d(its,kts,j)&
                   ! ,p2d=pdh(its,kts),p2di=pdhi(its,kts)&
                   ! ,pi2d=pi3d(ims,kms,j)&
                    ,ttnp=rthblten(ims,kms,j),qtnp=rqvbl2dt(its,kts)&!,ndiff=ndiff&
                   ! ,dz8w2d=dz8w(ims,kms,j)&
                    ,psfcpa=psfc(ims,j)&
                    ,hfx=hfx(ims,j),lh=lh(ims,j),tsk=tsk(ims,j),swdnb=swdnb(ims,j)&
                    ,cldfrax=cldfra(ims,kms,j),swdnflxx = swdnflx(ims,kms,j)&
                    ,dt=dt& !,z=z&
                    ,ids=ids,ide=ide, jds=jds,jde=jde, kds=kds,kde=kde&
                    ,ims=ims,ime=ime, jms=jms,jme=jme, kms=kms,kme=kme&
                    ,its=its,ite=ite, jts=jts,jte=jte, kts=kts,kte=kte&
                    ,d_top=d_top,d_mid=d_mid, d_high=d_high) 

!
     do k = kts,kte
       do i = its,ite
         rthblten(i,k,j) = rthblten(i,k,j)!/pi3d(i,k,j)
         rqvblten(i,k,j) = rqvbl2dt(i,k)
         rqcblten(i,k,j) = rqvbl2dt(i,k+kte)
         !if(present(rqiblten)) rqiblten(i,k,j) = rqvbl2dt(i,k+kte+kte)
       enddo
     enddo
!
   enddo
!
   end subroutine annpbl

!
!-------------------------------------------------------------------------------
! 2D (column) PBL
!-------------------------------------------------------------------------------

   SUBROUTINE annpbl2d(pbl_input_mean, pbl_input_scale, pbl_diagnostic_mean,  &
                  pbl_diagnostic_scale,                                       &
                  state_mean, state_scale, sl_real_mean, rt_real_mean,        &
                  w_real_mean, sl_latent_to_real,                             &
                  rt_latent_to_real, w_latent_to_real,  &
                  rcld_latent_to_real,   &
                  rrain_latent_to_real,                                       &
                  cld_latent_to_real,   &
                  pbl_encoder_W, pbl_encoder_b,           & 
                  pbl_hidden_W, pbl_hidden_b,                                 &
                  pbl_tend_decoder_W, pbl_tend_decoder_b, pbl_diag_decoder_W, &
                  pbl_diag_decoder_b,                                         &
                !Arrays Dimensions
                  n_input, n_diag, n_state, n_real, n_latent, n_code,         &
                  n_rcld, n_cld, n_rrain, n_sl, n_rt, n_w,                    &
                !Real WRF variables
                  dt,j,                                                       & !time step
                  thx, qvx, ttnp, qtnp, wx,                                   & !Theta, qv, and tend.
                  lh, hfx, tsk, swdnb, psfcpa,                                & ! 2d var.
                  cldfrax, swdnflxx, qcx, qrx,                                & ! cld, sw rad, liquid, rain
                 ! znu,                                                        & 
                  ids,ide,jds,jde,kds,kde,                                    &
                  ims,ime,jms,jme,kms,kme,                                    &
                  its,ite, jts,jte, kts,kte,                                  &
                 ! znu, z,                                                     &
                  d_top,  d_mid, d_high )

!
!--------------------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------------------
!
! Constants             
 real, parameter:: g= 9.81
 real, parameter:: Cpd= 1005.
 real, parameter:: Rd= 287.
 real, parameter:: Lv= 2.5e6

!
 real,     intent(in   )   ::     dt ! time step

 
 integer,  intent(in   )   ::       ids,ide, jds,jde, kds,kde,                 &
                                    ims,ime, jms,jme, kms,kme,                 &
                                    its,ite, jts,jte, kts,kte,                 &
                                    j

 integer :: i, k

!
! Arrays dimensions from namelist
!
 integer , intent(in) :: n_input  !=33
 integer , intent(in) :: n_diag   !=34
 integer , intent(in) :: n_state  !=17
 integer , intent(in) :: n_real   !=20
 integer , intent(in) :: n_latent !=20
 integer , intent(in) :: n_code   !=128

!
 integer , intent(in) :: n_rcld  !=10
 integer , intent(in) :: n_cld   !=10
 integer , intent(in) :: n_rrain !=2
 integer , intent(in) :: n_sl    !=9
 integer , intent(in) :: n_rt    !=8
 integer , intent(in) :: n_w     !=9

!
! ANN Variables: These come from an external .nc file provided in the namelist.
!
 real, intent(in), dimension(n_input)          :: pbl_input_mean, pbl_input_scale
 real, intent(in), dimension(n_diag)           :: pbl_diagnostic_mean, pbl_diagnostic_scale
 real, intent(in), dimension(n_state)          :: state_mean, state_scale
 real, intent(in), dimension(n_real)           :: sl_real_mean, rt_real_mean, w_real_mean

 real, intent(in), dimension(n_sl, n_latent)   :: sl_latent_to_real
 real, intent(in), dimension(n_rt, n_latent)   :: rt_latent_to_real
 real, intent(in), dimension(n_w, n_latent)    :: w_latent_to_real

 real, intent(in), dimension(n_rcld, n_latent) :: rcld_latent_to_real
 real, intent(in), dimension(n_rrain,n_latent) :: rrain_latent_to_real
 real, intent(in), dimension(n_cld, n_latent)  :: cld_latent_to_real


 real, intent(in), dimension(n_input, n_code)  :: pbl_encoder_W
 real, intent(in), dimension(n_code)           :: pbl_hidden_b
 real, intent(in), dimension(n_code, n_code)   :: pbl_hidden_W
 real, intent(in), dimension(n_code)           :: pbl_encoder_b
 real, intent(in), dimension(n_code, n_state)  :: pbl_tend_decoder_W
 real, intent(in), dimension(n_state)          :: pbl_tend_decoder_b
 real, intent(in), dimension(n_code, n_diag)   :: pbl_diag_decoder_W
 real, intent(in), dimension(n_diag)           :: pbl_diag_decoder_b

!
!------ 2D variables (1 Column) ------------------------------------------------
!

 integer, parameter :: ndiff = 1
 
 real,     dimension( ims:ime, kms:kme )                                   , & ! potential temperature 
             intent(in   )   ::                                            thx,&
                                                                           cldfrax, &! cld fraction 2d 
                                                                           swdnflxx! 3D SW rad flux

 real,     dimension( its:ite, kts:kte*ndiff )                             , & !water vapor
             intent(in   )   ::                                            qvx,&
                                                                           qcx,& ! qcloud
                                                                           qrx ! rain
                                                                                        
 real,     dimension( ims:ime, kms:kme+1 )                                , & ! vertical velocity 
             intent(in   )   ::                                            wx


 real,     dimension( ims:ime, kms:kme )                                   , &
             intent(inout)   ::                                          ttnp ! temperature tend.

 real,     dimension( its:ite, kts:kte*ndiff )                             , &
             intent(inout)   ::                                          qtnp ! water vapor tend.

 real,     dimension( ims:ime )                                            , &
             intent(in   )   ::                                          tsk, & 
                                                                         hfx, &
                                                                          lh, &
                                                                        swdnb,& !sfc sw down
                                                                       psfcpa !sfc pressure

 !real,     dimension( kms:kme )                                            ,&
 !            intent(in   )   ::                                           znu !vertical levels. 



 INTEGER, INTENT(IN) ::  d_top ! top of the BPL say 3 km to be specified in the namelist
 INTEGER, INTENT(IN) ::  d_mid, d_high ! index of mid, and high clouds to be specified in the namelist


! local  
 real, dimension(ims:ime,kms:n_real)  :: sl_real, rt_real, w_real
! real, dimension(ims:ime,kms:n_latent):: sl_latent, rt_latent, w_latent
 real, dimension(:,:), allocatable :: sl_latent, rt_latent, w_latent
 real, dimension(ims:ime)  :: cldlow, cldmid, cldhigh, swdn_tod 
 real, dimension(ims:ime )            :: sl_domain_top, rt_domain_top, rrain_domain_top

 real, dimension(ims:ime,kms:n_real)  :: sl_avg, rt_avg, w_avg
 real, dimension(kms:n_real) :: sl_diff, rt_diff, w_diff
 real, dimension(ims:ime, kms:n_real) ::  pbl_input, tend_pbl, sl_pbl_tend_latent
 real, dimension(ims:ime, kms:n_real) ::  rt_pbl_tend_latent
 real, dimension(:,:), allocatable :: X,XX, Y, SL, RT 

 real, dimension(n_code) :: X_ARR ! dummy variable 
 
 real, dimension(ims:ime,n_input) :: pbl_input_raw

 real  ::  ttend, qtend

 CHARACTER*256                :: message
 
! OUTPUT goes in tendencies.
 real, dimension (ims:ime, n_real) :: sl_pbl_tend_real, rt_pbl_tend_real
!

 real    ::  dt2,rdt
 real    ::  dtstep

!
!----------------------------
!

   dtstep = dt
   dt2 = 2.*dtstep
   rdt = 1./dt2


!            
! calculating liquid water static energy sl_real and rt_real and w_real!             
!
! use potential temperature instead of sl

 do k= kts, kte
   do i= its,ite
          sl_real(i,k) = thx(i,k)+300
          rt_real(i,k) = qvx(i,k) + qcx(i,k)
          w_real(i,k)  = 0.5*(wx(i,k) + wx(i,k+1))
    enddo
 enddo


!
! variables at low, mid, high, and domain top levels !
!


 do i= its, ite 
     rrain_domain_top(i) = qrx(i, d_top)               !
     sl_domain_top(i)    = sl_real(i, d_top)           !
     rt_domain_top(i)    = rt_real(i, d_top)           !
     rrain_domain_top(i) = qrx(i,d_top)                !
     cldmid(i)           = cldfrax(i,d_mid )                !
     cldhigh(i)          = cldfrax(i,d_high )                !
     swdn_tod(i)         = swdnflxx(i, d_top)          !
 enddo


!
!------- Interpolate variables onto weights heights ------------------------------
!
 do k=kts,d_top ! d_top=20 at about 3 km
   do i = its, ite
    sl_avg(i,k) = sl_real(i,k) !interp_0(sl_real(i,k), z_force, zzz_avg(k), n_real)
    rt_avg(i,k) = rt_real(i,k) !interp_0(rt_real(i,k), z_force, zzz_avg(k), n_real) 
    w_avg(i,k)  = w_real(i,k)  !interp_0(w_real(i,k) , z_force, zzz_avg(k), n_real)  
   enddo
 enddo

!
!----- Transform variables from physical space to to latent space ---------------
!

 do i= its, ite
       sl_diff = sl_avg(i,:) - sl_real_mean(:)
       rt_diff = rt_avg(i,:) - rt_real_mean(:)
       w_diff  = w_avg(i,:)  - w_real_mean(:)
       !sl_latent(i,:) = matmul(sl_avg(i,:) - sl_real_mean(:) , transpose(sl_latent_to_real) )
       sl_latent = matmul(reshape(sl_diff,(/1,shape(sl_diff)/)) ,transpose(reshape(sl_latent_to_real,[n_sl, n_latent]) ))
                   
! note sl_avg and sl_real_mean are NOT the same size unless interpolated.
       rt_latent = matmul(reshape(rt_diff,(/1,shape(rt_diff)/))  , transpose(reshape(rt_latent_to_real,(/n_rt, n_latent/) ) ))
       w_latent  = matmul(reshape(w_diff,(/1,shape(w_diff)/))  , transpose(reshape(w_latent_to_real, (/n_w, n_latent/) ) ))
!
!-------- Concatenate variables into one input array ---------------------------
!
       pbl_input_raw(i,:) =[sl_latent(1,:), rt_latent(1,:), w_latent(1,:),        &
                            sl_domain_top(i), rt_domain_top(i),                   & 
                            lh(i),hfx(i) , tsk(i), cldmid(i), cldhigh(i),         &
                            swdnb(i), swdn_tod(i), psfcpa(i), rrain_domain_top(i)]

! Here make sure number of elements in input array matches that of modes n_input
! (in the namelist):

       if (size(pbl_input_raw,2).NE. n_input) then
          CALL wrf_message("PBL input array does not match n_input")
          WRITE( message , * ) 'shape of pbl_input_raw = ', shape(pbl_input_raw),&
          ', while the number of inputs (n_input) = ', n_input
          CALL wrf_error_fatal( message )
       endif  
!
! Perform normalization

       pbl_input(i,:) = (pbl_input_raw(i,:) - pbl_input_mean) / pbl_input_scale ! dims. don't match!

! Forward propagation. The following is: ReL (W2 * ReL(W3 * x_in + b3) + b2)

       X = reshape(pbl_input(i,:), [1,n_input] ) ! Make inputs a matrix of dim (1, shape(input)) 

       XX = matmul(X, reshape(pbl_encoder_W, [n_input,n_code])) + reshape(pbl_encoder_b,[1,n_code])
       
       call relu(XX,size(XX), X_ARR )  ! X_ARR is an array (dimension(n_code))
        
       X = reshape(X_ARR,[1,size(XX)]) ! size = n_code 

       XX = matmul(X, reshape(pbl_hidden_W, [n_code,n_code])) + reshape(pbl_hidden_b,[1,n_code])     
       
       call relu(XX, size(XX),X_ARR )  ! X_ARR is an array
       
       X = reshape(X_ARR,[1,size(XX)]) ! size = n_code

       Y = matmul(X, reshape(pbl_tend_decoder_W, [n_code,n_state])) + reshape(pbl_tend_decoder_b,[1,n_state])
       
       tend_pbl(i,:) = Y(1,:) * state_scale

       sl_pbl_tend_latent(i,:) = tend_pbl(i, 1:n_sl) ! indexing in Fortrans starts with 1 not 0

       rt_pbl_tend_latent(i,:) = tend_pbl(i, n_sl+1 : size(tend_pbl))

!
!---- Transform variables from latent space back to physical space -------------
!
       SL = matmul(reshape(sl_pbl_tend_latent(i,:), [1,n_real]) , reshape(sl_latent_to_real, [n_sl,n_latent]))
       RT = matmul(reshape(rt_pbl_tend_latent(i,:), [1,n_real]) , reshape(rt_latent_to_real, [n_rt,n_latent]))

       sl_pbl_tend_real(i,:) = SL(1,:) 
       rt_pbl_tend_real(i,:) = RT(1,:)

 enddo  ! enddo for i

!
!-- initialize tendencies ----------
!
 ttnp(its:ite,:) = 0.
 qtnp(its:ite,:) = 0.

 do k = kts, kte
    do i = its,ite
       if (k .le. d_top) then
          ttend = sl_real(i,k)*rdt + sl_pbl_tend_real(i,k)
          qtend = qvx(i,k)*rdt + rt_pbl_tend_real(i,k)
       else
          ttend = sl_real(i,k)*rdt
          qtend = qvx(i,k)*rdt
       endif
       ttnp(i,k) = ttnp(i,k)+ttend
       qtnp(i,k) = qtnp(i,k)+qtend 
    enddo
 enddo



 deallocate(sl_latent)
 deallocate(rt_latent)
 deallocate(w_latent)
 deallocate(X)
 deallocate(XX)
 !deallocate(RR)
 deallocate(Y)
 deallocate(SL)
 deallocate(RT)

  END SUBROUTINE annpbl2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! defining linear activation function                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 subroutine relu(xx,n, res) !result(res)
    !! Rectified Linear Unit (RELU) activation function.
    IMPLICIT NONE
    integer,intent(in) :: n
    real , dimension(n), intent(in) :: xx
    real, dimension(n), intent(inout) :: res
    res = max(0., xx)
 end subroutine !relu


!pure function relu(x) result(res)
!    !! Rectified Linear Unit (RELU) activation function.
!    real , intent(in) :: x(:)
!    real :: res(size(x))
!    res = max(0., x)
!  end function relu




!
!--------- Constructor----------------------------------------------------------
!
   subroutine annpblinit(rthblten,rqvblten,                       &
                      rqcblten,rqiblten,p_qi,p_first_scalar,                   &
                      restart, allowed_to_read,                                &
                      ids, ide, jds, jde, kds, kde,                            &
                      ims, ime, jms, jme, kms, kme,                            &
                      its, ite, jts, jte, kts, kte                 )
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
   logical , intent(in)          :: restart, allowed_to_read
   integer , intent(in)          ::  ids, ide, jds, jde, kds, kde,             &
                                     ims, ime, jms, jme, kms, kme,             &
                                     its, ite, jts, jte, kts, kte
   integer , intent(in)          ::  p_qi,p_first_scalar
   real , dimension( ims:ime , kms:kme , jms:jme ), intent(out) ::             &
                                                                   !   rublten, &
                                                                   !   rvblten, &
                                                                     rthblten, &
                                                                     rqvblten, &
                                                                     rqcblten, &
                                                                     rqiblten
   integer :: i, j, k, itf, jtf, ktf
!
   jtf = min0(jte,jde-1)
   ktf = min0(kte,kde-1)
   itf = min0(ite,ide-1)
!
   if(.not.restart)then
     do j = jts,jtf
       do k = kts,ktf
         do i = its,itf
            !rublten(i,k,j) = 0.
            !rvblten(i,k,j) = 0.
            rthblten(i,k,j) = 0.
            rqvblten(i,k,j) = 0.
            rqcblten(i,k,j) = 0.
         enddo
       enddo
     enddo
   endif
!
   if (p_qi .ge. p_first_scalar .and. .not.restart) then
     do j = jts,jtf
       do k = kts,ktf
         do i = its,itf
           rqiblten(i,k,j) = 0.
         enddo
       enddo
     enddo
   endif
!
   end subroutine annpblinit
!-------------------------------------------------------------------------------
end module module_bl_annpbl
!-------------------------------------------------------------------------------

