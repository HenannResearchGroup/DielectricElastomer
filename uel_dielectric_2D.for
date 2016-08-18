!************************************************************************
! User element (UEL) for coupled large-deformation elasticity and  
!  dielectric behavior in two-dimensions -- plane-strain and axisymmetric
!************************************************************************
! Element details:
!************************************************************************
!
! Solution variables (or nodal variables) are the displacements (DOFs 1-2)
!  and the electric potential (DOF 11).
!
! Material behavior is Gent rubber elasticity and ideal 
!  (constant permittivity) dielectric behavior.
! 
! This subroutine is for a two-dimensional 4-node isoparametric
!  quadrilateral element as shown below with 4pt (full) or 1pt 
!  (reduced) integration, as well as plane-strain or
!  axisymmetric settings.
!
! In order to avoid locking for the fully-integrated element, we
!  use the F-bar method of de Souza Neto (1996).
!
! Surface charge density boundary conditions in the referential 
!  and current configurations are supported in this element.  Based 
!  on our convention, the face on which the charge density
!  is applied is the "label", i.e.
!  - U1,U2,U3,U4 refer to referential charge densities applied
!    to faces 1,2,3,4, respectively,
!  - U11,U12,U13,U14 refer to spatial charge densities applied
!    to faces 1,2,3,4, respectively.
!  Mechanical, traction- and pressure-type boundary conditions 
!  may be applied to the dummy mesh using the Abaqus built-in 
!  commands *Dload or *Dsload.
!     
!
!              A eta (=xi_2)
!  4-node      |
!   quad       | Face 3
!        4-----------3
!        |     |     |
!        |     |     |
! Face 4 |     ------|---> xi (=xi_1)
!        |           |
!        |           |  Face 2
!        1-----------2
!            Face 1
!
! Shawn A. Chester, July 2011
! David L. Henann, August 2011
!                  September 2012
!
!************************************************************************
! Usage:
!************************************************************************
!
! User element statement in the input file:
!  *User Element,Nodes=4,Type=U1,Iproperties=1,Properties=4,Coordinates=2,Variables=1,Unsymm
!  1,2,11
!
! Note: No local state variables are used in this element, thus we may set the above 
!  parameter 'Variables' to any non-zero integer.
!
! In the subroutine UEL, set 'nInt' = number of integration points
!  Options are nInt=4 (full integration) or nInt=1 (reduced integration)
!
! In the input file, set the parameter pe=1 for plane-strain or 
!  pe=0 for axisymmetric
!
!
!     Material Properties Vector
!     --------------------------------------------------------------
!     Gshear  = props(1)  ! Shear modulus
!     Kbulk   = props(2)  ! Bulk modulus
!     Imax    = props(3)  ! Max value of (I1bar-3) (Gent parameter)
!     epsilon = props(4)  ! Permittivity
!     pe      = jprops(1) ! Plane strain=1, axisymmetric=0
!
!************************************************************************

      subroutine UEL(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,
     +     props,nprops,coords,mcrd,nnode,uall,duall,vel,accn,jtype,
     +     time,dtime,kstep,kinc,jelem,params,ndload,jdltyp,adlmag,
     +     predef,npredf,lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,
     +     njprop,period)

      implicit none
      !
      ! variables defined in uel, passed back to Abaqus
      !
      real*8 rhs(mlvarx,*),amatrx(ndofel,ndofel),svars(*),energy(8),
     +  pnewdt
      !
      ! variables passed into UEL
      !
      integer ndofel,nrhs,nsvars,nprops,mcrd,nnode,jtype,kstep,kinc,
     +  jelem,ndload,jdltyp(mdload,*),npredf,lflags(*),mlvarx,mdload,
     +  jprops(*),njprop
      !
      real*8 props(*),coords(mcrd,nnode),uall(ndofel),duall(mlvarx,*),
     +  vel(ndofel),accn(ndofel),time(2),dtime,params(*),
     +  adlmag(mdload,*),predef(2,npredf,nnode),ddlmag(mdload,*),period
      !
      ! variables defined and used in the UEL
      !
      integer i,j,k,A11,B11,A12,B12,nInt,nIntPt,intpt,nDim,stat,pe,face
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      parameter(nDim=2) ! number of spatial dimensions, do not change
      parameter(nInt=4) ! number of integration points
      !
      ! nInt=4: fully-integrated, nInt=1: reduced integration
      ! 
      ! When reduced integration is selected, be sure to
      !  specify an hourglass stiffness of about 0.005G
      !  for the dummy mesh in the input file.  Also, change
      !  the dummy mesh element to reduced.
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      real*8 u(nNode,2),du(nNode,ndofel),phi(nNode),coordsC(mcrd,nNode),
     +  Ru(2*nNode,1),Rphi(nNode,1),Kuu(2*nNode,2*nNode),
     +  Kuphi(2*nNode,nNode),Kphiu(nNode,2*nNode),Kphiphi(nNode,nNode),
     +  Iden(3,3),xi(nInt,2),w(nInt),sh0(nNode),sh(nNode),dsh0(nNode,2),
     +  dshC0(nNode,2),dsh(nNode,2),dshC(nNode,2),dshxi(nNode,2),
     +  detMapJ0,detMapJ0C,detMapJ,detMapJC,Fc_tau(3,3),detFc_tau,
     +  F_tau(3,3),detF_tau,ER(3,1),T_tau(3,3),D_tau(3,1),
     +  SpUUMod(3,3,3,3),SpPhiPhiMod(3,3),SpPhiUMod(3,3,3),
     +  SpUPhiMod(3,3,3),Smat(3,1),Bmat(3,2*nNode),Gmat(4,2*nNode),
     +  G0mat(4,2*nNode),Amat(4,4),Qmat(4,4),AmatPhiU(2,4),
     +  AmatUPhi(4,2),Gmatphi(2,nNode),body(3),BodyForceRes(2*nNode,1),
     +  bodyCharge,BodyChargeRes(nNode,1),Smatphi(2,1),Amatphi(2,2),
     +  SmatAx(4,1),BmatAx(4,2*nNode),GmatAx(5,2*nNode),
     +  G0MatAx(5,2*nNode),AmatAx(5,5),QmatAx(5,5),AmatPhiUAx(2,5),
     +  AmatUPhiAx(5,2),BodyForceResAx(2*nNode,1),
     +  BodyChargeResAx(nNode,1),AR,AR0,ARc,flux,Le
      !
      real*8 zero,one,two,half,Pi,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,Pi=3.141592653d0,
     +     three=3.d0,third=1.d0/3.d0)


      ! Check the procedure type, this should be a coupled
      !  temperature displacement, which is either 72 or 73
      !
      if((lflags(1).eq.72).or.(lflags(1).eq.73)) then
         !
         ! correct procedure specified
         !
      else
         write(*,*) 'Abaqus does not have the right procedure'
         write(*,*) 'go back and check the procedure type'
         write(*,*) 'lflags(1)=',lflags(1)
         call xit
      endif


      ! Make sure Abaqus knows you are doing a large
      !  deformation problem
      !
      if(lflags(2).eq.0) then
         !
         ! lflags(2)=0 -> small disp.
         ! lflags(2)=1 -> large disp.
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a small displacement analysis'
         write(*,*) 'go in and set nlgeom=yes'
         call xit
      endif


      ! Check to see if you are doing a general
      !  step or a linear purturbation step
      !
      if(lflags(4).eq.1) then
         !
         ! lflags(4)=0 -> general step
         ! lflags(4)=1 -> linear purturbation step
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a linear purturbation step'
         call xit         
      endif


      ! Do nothing if a ``dummy'' step
      !
      if(dtime.eq.zero) return


      ! Get flag for plane strain or axisymmetric
      !
      pe = jprops(1)


      ! Identity tensor
      !
      call onem(Iden)


      ! Initialize the residual and tangent matrices to zero.
      !
      Ru  = zero
      Rphi  = zero
      Kuu = zero
      Kphiphi = zero
      Kuphi = zero
      Kphiu = zero
      Energy = zero


      ! Obtain nodal displacements and potentials
      !
      k = 0
      do i=1,nNode
         do j=1,nDim
            k = k + 1
            u(i,j) = Uall(k)
            du(i,j) = DUall(k,1)
         enddo
         k = k + 1
         phi(i) = Uall(k)
      enddo


      ! Obtain current nodal coordinates
      !
      do i=1,nNode
         do j=1,nDim
            coordsC(j,i) = coords(j,i) + u(i,j)
         enddo
      enddo


      ! Impose the heuristic restriction that the displacement 
      !  increment is less than ten times the element diagonal,
      !  which is taken as an approximate element size.
      !
      Le = dsqrt(((coordsC(1,1)-coordsC(1,3))**two) + 
     +                    ((coordsC(2,1)-coordsC(2,3))**two))
      !
      do i=1,nNode
         do j=1,nDim
            if(dabs(du(i,j)).gt.10.d0*Le) then
               pnewdt = 0.5
               return
            endif
         enddo
      enddo


      !----------------------------------------------------------------
      ! 
      ! Take this opportunity to perform calculations at the element
      !  centroid.  Get the deformation gradient for use in the 'F-bar' method.
      !
      ! Reference for the F-bar method:
      !  de Souza Neto, E.A., Peric, D., Dutko, M., Owen, D.R.J., 1996.
      !  Design of simple low order finite elements for large strain
      !  analysis of nearly incompressible solids. International Journal
      !  of Solids and Structures 33, 3277-3296.
      !
      ! Obtain shape functions and their local gradients at the element
      !  centroid, that means xi=eta=zeta=0.0, and nIntPt=1
      !
      if(nNode.eq.4) then
         call calcShape2DLinear(1,zero,1,sh0,dshxi)
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif


      ! Map shape functions from local to global reference coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif


      ! Map shape functions from local to global current coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif


      AR0  = one ! for plane-strain you could put the depth here
      ARc  = one !
      if(pe.ne.1) then ! this is axisymmetric
         AR0 = zero
         ARc = zero
         do i=1,nNode
            ! radial coord in ref config at centroid
            AR0 = AR0 + sh0(i)*coords(1,i)
            ! radial coord in current config at centroid
            ARc = ARc + sh0(i)*(coords(1,i) + u(i,1))
         enddo
      endif


      ! Calculate the deformation gradient at the element centroid
      !  at the end of the increment for use in the 'F-bar' method.
      !  The subscript tau denotes the time at the end of the increment.
      !
      Fc_tau = Iden
      do i=1,nDim
         do j=1,nDim
            do k=1,nNode
               Fc_tau(i,j) = Fc_tau(i,j) + dsh0(k,j)*u(k,i)
            enddo
         enddo
      enddo
      if(pe.ne.1) then
         !
         ! modify for axisymmetric
         !
         Fc_tau(3,3) = ARc/AR0 ! 'hoop' stretch
         !
         ! axisymmetric implementation of detF
         !
         call mdet(Fc_tau,detFc_tau)
      else
         !
         ! modify for plane-strain
         !
         Fc_tau(3,3) = one
         !
         ! 2D plane-strain implementation detF
         !
         detFc_tau = Fc_tau(1,1)*Fc_tau(2,2) - Fc_tau(1,2)*Fc_tau(2,1)
      endif
      !
      ! With the deformation gradient known at the element centroid
      !  we are now able to implement the `F-bar' method later
      !
      !----------------------------------------------------------------




      !----------------------------------------------------------------
      ! Begin the loop over integration points
      !
      ! Obtain integration point local coordinates and weights
      !
      if(nNode.eq.4) then
         !
         ! gauss integration for a rectangular element
         !
         if(nInt.eq.4) then
            call xint2D4pt(xi,w,nIntPt) ! 4-pt integration, nInt=4 above
         elseif(nInt.eq.1) then
            call xint2D1pt(xi,w,nIntPt) ! 1-pt integration, nInt=1 above
         else
            write(*,*) 'Invalid number of int points, nInt=',nInt
            call xit
         endif
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif



      ! Loop over integration points
      !
      do intpt=1,nIntPt


         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.4) then
            call calcShape2DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.4'
            call xit
         endif
         

         ! Map shape functions from local to global reference coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            call xit 
         endif


         ! Map shape functions from local to global current coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coordsC,dshC,detMapJC,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coordsC,dshC,detMapJC,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            call xit 
         endif


         AR0  = one ! for plane-strain you could put the depth here
         AR   = one !
         if(pe.ne.1) then ! this is axisymmetric
            AR0  = zero
            AR   = zero
            do i=1,nNode
               ! radial coord in reference config
               AR0 = AR0 + sh(i)*coords(1,i)
               ! radial coord in current config
               AR  = AR  + sh(i)*(coords(1,i) + u(i,1))
            enddo
            AR0  = two*Pi*AR0
            AR   = two*Pi*AR
         endif


         ! Obtain the referential electric field at this integration point
         !
         ER = zero
         do k=1,nNode
            do i=1,nDim
               ER(i,1) = ER(i,1) + phi(k)*dsh(k,i)
            enddo
         enddo


         ! Obtain the deformation gradient at this integration point.
         ! The subscript tau denotes the time at the end of the 
         !  increment.
         !
         F_tau = Iden
         do i=1,nDim
            do j=1,nDim
               do k=1,nNode
                  F_tau(i,j) = F_tau(i,j) + dsh(k,j)*u(k,i)
               enddo
            enddo
         enddo
         !
         ! modify F(3,3) for plane-strain or axisymmetric
         !
         if(pe.ne.1) then
            !
            ! this is axisymmetric
            !
            F_tau(3,3) = AR/AR0
         else
            !
            ! this is plane strain
            !
            F_tau(3,3) = one
         endif
         !
         ! Modify the deformation gradient for the `F-bar' method
         !  only when using the 4 node fully integrated linear
         !  element, do not use the `F-bar' method for any other element
         !
         if((nNode.eq.4).and.(nInt.eq.4)) then
            if(pe.eq.1) then
               !
               !  2D plane-strain implementation
               !
               detF_tau = F_tau(1,1)*F_tau(2,2) - F_tau(1,2)*F_tau(2,1)
               do i=1,nDim
                  do j=1,nDim
                     F_tau(i,j) =((detFc_tau/detF_tau)**half)*F_tau(i,j)
                  enddo
               enddo
            else
               !
               ! 2D axisymmetric implementation
               !
               call mdet(F_tau,detF_tau)
               F_tau = ((detFc_tau/detF_tau)**third)*F_tau
            endif
         endif


         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         !
         ! Perform the constitutive update at this integ. point
         !
         call Gent(props,nprops,F_tau,ER,T_tau,D_tau,
     +        SpUUMod,SpUPhiMod,SpPhiUMod,SpPhiPhiMod,stat)
         !
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
         !
         !
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


         ! Compute/update the displacement residual vector
         !
         if(pe.eq.1) then
            !
            ! this is plane strain
            !
            Smat(1,1) = T_tau(1,1)
            Smat(2,1) = T_tau(2,2)
            Smat(3,1) = T_tau(1,2)
            !
            Bmat = zero
            do k=1,nNode
               Bmat(1,1+nDim*(k-1)) = dshC(k,1)
               Bmat(2,2+nDim*(k-1)) = dshC(k,2)
               Bmat(3,1+nDim*(k-1)) = dshC(k,2)
               Bmat(3,2+nDim*(k-1)) = dshC(k,1)
            enddo
            !
            body = zero ! The body force vector may be specified here
            !
            BodyForceRes = zero
            do k=1,nNode
               BodyForceRes(1+nDim*(k-1),1) = sh(k)*body(1)
               BodyForceRes(2+nDim*(k-1),1) = sh(k)*body(2)
            enddo
            !
            Ru = Ru 
     +           - matmul(transpose(Bmat),Smat)*detmapJC*w(intpt)*AR
     +           + BodyForceRes*detmapJC*w(intpt)*AR
         else
            !
            ! this is axisymetric
            !
            SmatAx(1,1) = T_tau(1,1)
            SmatAx(2,1) = T_tau(2,2)
            SmatAx(3,1) = T_tau(1,2)
            SmatAx(4,1) = T_tau(3,3)
            !
            BmatAx = zero
            do k=1,nNode
               BmatAx(1,1+nDim*(k-1)) = dshC(k,1)
               BmatAx(2,2+nDim*(k-1)) = dshC(k,2)
               BmatAx(3,1+nDim*(k-1)) = dshC(k,2)
               BmatAx(3,2+nDim*(k-1)) = dshC(k,1)
               BmatAx(4,1+nDim*(k-1)) = sh(k)/(AR/(two*Pi))
            enddo
            !
            body = zero ! The body force vector may be specified here
            !
            BodyForceResAx = zero
            do k=1,nNode
               BodyForceResAx(1+nDim*(k-1),1) = sh(k)*body(1)
               BodyForceResAx(2+nDim*(k-1),1) = sh(k)*body(2)
            enddo
            !
            Ru = Ru
     +          - matmul(transpose(BmatAx),SmatAx)*detmapJC*w(intpt)*AR
     +          + BodyForceResAx*detmapJC*w(intpt)*AR
         endif


         ! Compute/update the electric potential residual vector
         !
         Smatphi(1,1) = D_tau(1,1)
         Smatphi(2,1) = D_tau(2,1)
         !
         Gmatphi = zero
         do i = 1,nDim
            do k = 1,nNode
               Gmatphi(i,k) = dshC(k,i)
            end do
         end do
         !
         bodyCharge = zero ! the free charge density may be specified here
         !
         BodyChargeRes = zero
         do k=1,nNode
            BodyChargeRes(k,1) = sh(k)*bodyCharge
         end do
         !
         Rphi = Rphi + matmul(transpose(Gmatphi),Smatphi)
     +                 *detmapJC*w(intpt)*AR
     +               + BodyChargeRes*detmapJC*w(intpt)*AR


         ! Compute/update the displacement tangent matrix
         !
         if(pe.eq.1) then
            !
            ! this is plane strain
            !
            Gmat = zero
            do k=1,nNode
               Gmat(1,1+nDim*(k-1)) = dshC(k,1)
               Gmat(2,2+nDim*(k-1)) = dshC(k,1)
               Gmat(3,1+nDim*(k-1)) = dshC(k,2)
               Gmat(4,2+nDim*(k-1)) = dshC(k,2)
            enddo

            G0mat = zero
            do k=1,nNode
               G0mat(1,1+nDim*(k-1)) = dshC0(k,1)
               G0mat(2,2+nDim*(k-1)) = dshC0(k,1)
               G0mat(3,1+nDim*(k-1)) = dshC0(k,2)
               G0mat(4,2+nDim*(k-1)) = dshC0(k,2)
            enddo
            !
            Amat = zero
            Amat(1,1) = SpUUMod(1,1,1,1)
            Amat(1,2) = SpUUMod(1,1,2,1)
            Amat(1,3) = SpUUMod(1,1,1,2)
            Amat(1,4) = SpUUMod(1,1,2,2)
            Amat(2,1) = SpUUMod(2,1,1,1)
            Amat(2,2) = SpUUMod(2,1,2,1)
            Amat(2,3) = SpUUMod(2,1,1,2)
            Amat(2,4) = SpUUMod(2,1,2,2)
            Amat(3,1) = SpUUMod(1,2,1,1)
            Amat(3,2) = SpUUMod(1,2,2,1)
            Amat(3,3) = SpUUMod(1,2,1,2)
            Amat(3,4) = SpUUMod(1,2,2,2)
            Amat(4,1) = SpUUMod(2,2,1,1)
            Amat(4,2) = SpUUMod(2,2,2,1)
            Amat(4,3) = SpUUMod(2,2,1,2)
            Amat(4,4) = SpUUMod(2,2,2,2)
            !
            Qmat = zero
            Qmat(1,1) = half*(Amat(1,1)+Amat(1,4)) - half*T_tau(1,1)
            Qmat(2,1) = half*(Amat(2,1)+Amat(2,4)) - half*T_tau(1,2)
            Qmat(3,1) = half*(Amat(3,1)+Amat(3,4)) - half*T_tau(1,2)
            Qmat(4,1) = half*(Amat(4,1)+Amat(4,4)) - half*T_tau(2,2)
            Qmat(1,4) = half*(Amat(1,1)+Amat(1,4)) - half*T_tau(1,1)
            Qmat(2,4) = half*(Amat(2,1)+Amat(2,4)) - half*T_tau(1,2)
            Qmat(3,4) = half*(Amat(3,1)+Amat(3,4)) - half*T_tau(1,2)
            Qmat(4,4) = half*(Amat(4,1)+Amat(4,4)) - half*T_tau(2,2)
            !
            if((nNode.eq.4).and.(nInt.eq.4)) then
               !
               ! This is the tangent using the F-bar method with the
               !  4 node fully integrated linear element
               !
               Kuu = Kuu
     +              + matmul(matmul(transpose(Gmat),Amat),Gmat)
     +              *detMapJC*w(intpt)*AR
     +              +matmul(transpose(Gmat),matmul(Qmat,(G0mat - Gmat)))
     +              *detMapJC*w(intpt)*AR
            else
               !
               ! This is the tangent not using the F-bar method with all
               !  other elements
               !
               Kuu = Kuu
     +              + matmul(matmul(transpose(Gmat),Amat),Gmat)
     +              *detMapJC*w(intpt)*AR
            endif
            !
         else
            !
            ! this is axisymmetric
            !
            GmatAx = zero
            do k=1,nNode
               GmatAx(1,1+nDim*(k-1)) = dshC(k,1)
               GmatAx(2,2+nDim*(k-1)) = dshC(k,1)
               GmatAx(3,1+nDim*(k-1)) = dshC(k,2)
               GmatAx(4,2+nDim*(k-1)) = dshC(k,2)
               GmatAx(5,1+nDim*(k-1)) = sh(k)/(AR/(two*Pi))
            enddo
            !
            G0matAx = zero
            do k=1,nNode
               G0matAx(1,1+nDim*(k-1)) = dshC0(k,1)
               G0matAx(2,2+nDim*(k-1)) = dshC0(k,1)
               G0matAx(3,1+nDim*(k-1)) = dshC0(k,2)
               G0matAx(4,2+nDim*(k-1)) = dshC0(k,2)
               G0matAx(5,1+nDim*(k-1)) = sh0(k)/ARc
            enddo
            !
            AmatAx = zero
            AmatAx(1,1) = SpUUMod(1,1,1,1)
            AmatAx(1,2) = SpUUMod(1,1,2,1)
            AmatAx(1,3) = SpUUMod(1,1,1,2)
            AmatAx(1,4) = SpUUMod(1,1,2,2)
            AmatAx(1,5) = SpUUMod(1,1,3,3)
            AmatAx(2,1) = SpUUMod(2,1,1,1)
            AmatAx(2,2) = SpUUMod(2,1,2,1)
            AmatAx(2,3) = SpUUMod(2,1,1,2)
            AmatAx(2,4) = SpUUMod(2,1,2,2)
            AmatAx(2,5) = SpUUMod(2,1,3,3)
            AmatAx(3,1) = SpUUMod(1,2,1,1)
            AmatAx(3,2) = SpUUMod(1,2,2,1)
            AmatAx(3,3) = SpUUMod(1,2,1,2)
            AmatAx(3,4) = SpUUMod(1,2,2,2)
            AmatAx(3,5) = SpUUMod(1,2,3,3)
            AmatAx(4,1) = SpUUMod(2,2,1,1)
            AmatAx(4,2) = SpUUMod(2,2,2,1)
            AmatAx(4,3) = SpUUMod(2,2,1,2)
            AmatAx(4,4) = SpUUMod(2,2,2,2)
            AmatAx(4,5) = SpUUMod(2,2,3,3)
            AmatAx(5,1) = SpUUMod(3,3,1,1)
            AmatAx(5,2) = SpUUMod(3,3,2,1)
            AmatAx(5,3) = SpUUMod(3,3,1,2)
            AmatAx(5,4) = SpUUMod(3,3,2,2)
            AmatAx(5,5) = SpUUMod(3,3,3,3)
            !
            QmatAx = zero
            QmatAx(1,1) = third*(AmatAx(1,1)+AmatAx(1,4)+AmatAx(1,5)) 
     +           - (two/three)*T_tau(1,1)
            QmatAx(2,1) = third*(AmatAx(2,1)+AmatAx(2,4)+AmatAx(2,5))
     +           - (two/three)*T_tau(1,2)
            QmatAx(3,1) = third*(AmatAx(3,1)+AmatAx(3,4)+AmatAx(3,5))
     +           - (two/three)*T_tau(1,2)
            QmatAx(4,1) = third*(AmatAx(4,1)+AmatAx(4,4)+AmatAx(4,5))
     +           - (two/three)*T_tau(2,2)
            QmatAx(5,1) = third*(AmatAx(5,1)+AmatAx(5,4)+AmatAx(5,5))
     +           - (two/three)*T_tau(3,3)
            QmatAx(1,4) = QmatAx(1,1)
            QmatAx(2,4) = QmatAx(2,1)
            QmatAx(3,4) = QmatAx(3,1)
            QmatAx(4,4) = QmatAx(4,1)
            QmatAx(5,4) = QmatAx(5,1)
            QmatAx(1,5) = QmatAx(1,1)
            QmatAx(2,5) = QmatAx(2,1)
            QmatAx(3,5) = QmatAx(3,1)
            QmatAx(4,5) = QmatAx(4,1)
            QmatAx(5,5) = QmatAx(5,1)
            !
            if((nNode.eq.4).and.(nInt.eq.4)) then
               !
               ! This is the tangent using the F-bar method with the
               !  4 node fully integrated linear element
               !
               Kuu = Kuu
     +              + matmul(matmul(transpose(GmatAx),AmatAx),GmatAx)
     +              *detMapJC*w(intpt)*AR
     +              + matmul(transpose(GmatAx),matmul(QmatAx,
     +              (G0matAx-GmatAx)))*detMapJC*w(intpt)*AR
            else
               !
               ! This is the tangent NOT using the F-bar method with all
               !  other elements
               !
               Kuu = Kuu
     +              + matmul(matmul(transpose(GmatAx),AmatAx),GmatAx)
     +              *detMapJC*w(intpt)*AR
            endif
            !
         endif


         ! Compute/update the electric potential tangent matrix
         !
         Amatphi(1,1) = SpPhiPhiMod(1,1)
         Amatphi(1,2) = SpPhiPhiMod(1,2)
         Amatphi(2,1) = SpPhiPhiMod(2,1)
         Amatphi(2,2) = SpPhiPhiMod(2,2)
         !
         Kphiphi = Kphiphi - 
     +             (matmul(matmul(transpose(Gmatphi),Amatphi),Gmatphi)
     +                    *detMapJC*w(intpt))*AR


         ! Compute/update the electric potential/displacement tangent matrix
         !   Strictly speaking, use of 'F-bar' will affect this tangent;
         !   however, the effect is expected to be small and is neglected here.
         !
         if(pe.eq.1) then
            !
            ! this is plane strain
            !
            AmatPhiU = zero
            AmatPhiU(1,1) = SpPhiUMod(1,1,1)
            AmatPhiU(1,2) = SpPhiUMod(1,2,1)
            AmatPhiU(1,3) = SpPhiUMod(1,1,2)
            AmatPhiU(1,4) = SpPhiUMod(1,2,2)
            AmatPhiU(2,1) = SpPhiUMod(2,1,1)
            AmatPhiU(2,2) = SpPhiUMod(2,2,1)
            AmatPhiU(2,3) = SpPhiUMod(2,1,2)
            AmatPhiU(2,4) = SpPhiUMod(2,2,2)
            !
            Kphiu = Kphiu - matmul(matmul(transpose(Gmatphi),AmatPhiU),
     +                          Gmat)*detMapJC*w(intpt)*AR
            !
         else
            !
            ! this is axisymmetric
            !
            AmatPhiUAx = zero
            AmatPhiUAx(1,1) = SpPhiUMod(1,1,1)
            AmatPhiUAx(1,2) = SpPhiUMod(1,2,1)
            AmatPhiUAx(1,3) = SpPhiUMod(1,1,2)
            AmatPhiUAx(1,4) = SpPhiUMod(1,2,2)
            AmatPhiUAx(1,5) = SpPhiUMod(1,3,3)
            AmatPhiUAx(2,1) = SpPhiUMod(2,1,1)
            AmatPhiUAx(2,2) = SpPhiUMod(2,2,1)
            AmatPhiUAx(2,3) = SpPhiUMod(2,1,2)
            AmatPhiUAx(2,4) = SpPhiUMod(2,2,2)
            AmatPhiUAx(2,5) = SpPhiUMod(2,3,3)
            !
            Kphiu = Kphiu - matmul(matmul(transpose(Gmatphi),
     +                        AmatPhiUAx),GmatAx)*detMapJC*w(intpt)*AR
            !
         end if


         ! Compute/update the displacement/electric potential tangent matrix
         !
         if(pe.eq.1) then
            !
            ! this is plane strain
            !
            AmatUPhi = zero
            AmatUPhi(1,1) = SpUPhiMod(1,1,1)
            AmatUPhi(2,1) = SpUPhiMod(2,1,1)
            AmatUPhi(3,1) = SpUPhiMod(1,2,1)
            AmatUPhi(4,1) = SpUPhiMod(2,2,1)
            AmatUPhi(1,2) = SpUPhiMod(1,1,2)
            AmatUPhi(2,2) = SpUPhiMod(2,1,2)
            AmatUPhi(3,2) = SpUPhiMod(1,2,2)
            AmatUPhi(4,2) = SpUPhiMod(2,2,2)
            !
            Kuphi = Kuphi + matmul(matmul(transpose(Gmat),AmatUPhi),
     +                          Gmatphi)*detMapJC*w(intpt)*AR
            !
         else
            !
            ! this is axisymmetric
            !
            AmatUPhiAx = zero
            AmatUPhiAx(1,1) = SpUPhiMod(1,1,1)
            AmatUPhiAx(2,1) = SpUPhiMod(2,1,1)
            AmatUPhiAx(3,1) = SpUPhiMod(1,2,1)
            AmatUPhiAx(4,1) = SpUPhiMod(2,2,1)
            AmatUPhiAx(5,1) = SpUPhiMod(3,3,1)
            AmatUPhiAx(1,2) = SpUPhiMod(1,1,2)
            AmatUPhiAx(2,2) = SpUPhiMod(2,1,2)
            AmatUPhiAx(3,2) = SpUPhiMod(1,2,2)
            AmatUPhiAx(4,2) = SpUPhiMod(2,2,2)
            AmatUPhiAx(5,2) = SpUPhiMod(3,3,2)
            !
            Kuphi = Kuphi + matmul(matmul(transpose(GmatAx),
     +                       AmatUPhiAx),Gmatphi)*detMapJC*w(intpt)*AR
            !
         end if

      enddo
      !
      ! End the loop over integration points
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! Start loop over surface charge terms; surface charge densities
      !  in the referential and current configurations are supported.
      !
      ! Based on our convention, the face on which the charge density
      !  is applied is the "label", i.e.
      !  - U1,U2,U3,U4 refer to referential charge densities applied
      !    to faces 1,2,3,4, respectively,
      !  - U11,U12,U13,U14 refer to spatial charge densities applied
      !    to faces 1,2,3,4, respectively.
      !
      if (pe.eq.1) then
        !
        ! This is plane strain
        !
        if(ndload.gt.0) then
          !
          ! loop over faces and make proper modifications to
          !  residuals and tangents if needed
          !
          do i=1,ndload
            !
            face = jdltyp(i,1) ! label
            flux = adlmag(i,1) ! charge density
            !
            if(face.eq.1) then
               !
               ! charge density (reference) on face 1 of the element
               !
               ! Modify the displacement residual, loop over nodes
               !
               Le = dsqrt(((coords(1,1)-coords(1,2))**two) + 
     +                    ((coords(2,1)-coords(2,2))**two))
               Rphi(1,1) = Rphi(1,1) + half*flux*Le
               Rphi(2,1) = Rphi(2,1) + half*flux*Le
               !
               ! No modification to the tangent matrix
               !
            elseif(face.eq.2) then
               !
               ! charge density (reference) on face 2 of the element
               !
               ! Modify the displacement residual, loop over nodes
               !
               Le = dsqrt(((coords(1,2)-coords(1,3))**two) + 
     +                    ((coords(2,2)-coords(2,3))**two))
               Rphi(2,1) = Rphi(2,1) + half*flux*Le
               Rphi(3,1) = Rphi(3,1) + half*flux*Le
               !
               ! No modification to the tangent matrix
               !
            elseif(face.eq.3) then
               !
               ! charge density (reference) on face 3 of the element
               !
               ! Modify the displacement residual, loop over nodes
               !
               Le = dsqrt(((coords(1,3)-coords(1,4))**two) + 
     +                    ((coords(2,3)-coords(2,4))**two))
               Rphi(3,1) = Rphi(3,1) + half*flux*Le
               Rphi(4,1) = Rphi(4,1) + half*flux*Le
               !
               ! No modification to the tangent matrix
               !
            elseif(face.eq.4) then
               !
               ! charge density (reference) on face 4 of the element
               !
               ! Modify the displacement residual, loop over nodes
               !
               Le = dsqrt(((coords(1,4)-coords(1,1))**two) + 
     +                    ((coords(2,4)-coords(2,1))**two))
               Rphi(4,1) = Rphi(4,1) + half*flux*Le
               Rphi(1,1) = Rphi(1,1) + half*flux*Le
               !
               ! No modification to the tangent matrix
               !
            elseif(face.eq.11) then
               !
               ! charge density (spatial) on face 1 of the element
               !
               ! Modify the displacement residual, loop over nodes
               !
               Le = dsqrt(((coordsC(1,1)-coordsC(1,2))**two) + 
     +                    ((coordsC(2,1)-coordsC(2,2))**two))
               Rphi(1,1) = Rphi(1,1) + half*flux*Le
               Rphi(2,1) = Rphi(2,1) + half*flux*Le
               !
               ! Modify the tangent matrix
               !
               Kphiu(1,1) = Kphiu(1,1) - 
     +               (coordsC(1,1)-coordsC(1,2))*flux/(two*Le)
               Kphiu(1,2) = Kphiu(1,2) -
     +               (coordsC(2,1)-coordsC(2,2))*flux/(two*Le)
               Kphiu(1,3) = Kphiu(1,3) -
     +               (coordsC(1,2)-coordsC(1,1))*flux/(two*Le)
               Kphiu(1,4) = Kphiu(1,4) -
     +               (coordsC(2,2)-coordsC(2,1))*flux/(two*Le)
               Kphiu(2,1) = Kphiu(2,1) -
     +               (coordsC(1,1)-coordsC(1,2))*flux/(two*Le)
               Kphiu(2,2) = Kphiu(2,2) -
     +               (coordsC(2,1)-coordsC(2,2))*flux/(two*Le)
               Kphiu(2,3) = Kphiu(2,3) -
     +               (coordsC(1,2)-coordsC(1,1))*flux/(two*Le)
               Kphiu(2,4) = Kphiu(2,4) -
     +               (coordsC(2,2)-coordsC(2,1))*flux/(two*Le)
               ! 
            elseif(face.eq.12) then
               !
               ! charge density (spatial) on face 2 of the element
               !
               ! Modify the displacement residual, loop over nodes
               !
               Le = dsqrt(((coordsC(1,2)-coordsC(1,3))**two) + 
     +                    ((coordsC(2,2)-coordsC(2,3))**two))
               Rphi(2,1) = Rphi(2,1) + half*flux*Le
               Rphi(3,1) = Rphi(3,1) + half*flux*Le
               !
               ! Modify the tangent matrix
               !
               Kphiu(2,3) = Kphiu(2,3) - 
     +               (coordsC(1,2)-coordsC(1,3))*flux/(two*Le)
               Kphiu(2,4) = Kphiu(2,4) -
     +               (coordsC(2,2)-coordsC(2,3))*flux/(two*Le)
               Kphiu(2,5) = Kphiu(2,5) -
     +               (coordsC(1,3)-coordsC(1,2))*flux/(two*Le)
               Kphiu(2,6) = Kphiu(2,6) - 
     +               (coordsC(2,3)-coordsC(2,2))*flux/(two*Le)
               Kphiu(3,3) = Kphiu(3,3) -
     +               (coordsC(1,2)-coordsC(1,3))*flux/(two*Le)
               Kphiu(3,4) = Kphiu(3,4) - 
     +               (coordsC(2,2)-coordsC(2,3))*flux/(two*Le)
               Kphiu(3,5) = Kphiu(3,5) -
     +               (coordsC(1,3)-coordsC(1,2))*flux/(two*Le)
               Kphiu(3,6) = Kphiu(3,6) -
     +               (coordsC(2,3)-coordsC(2,2))*flux/(two*Le)
               !
            elseif(face.eq.13) then
               !
               ! charge density (spatial) on face 3 of the element
               !
               ! Modify the displacement residual, loop over nodes
               !
               Le = dsqrt(((coordsC(1,3)-coordsC(1,4))**two) + 
     +                    ((coordsC(2,3)-coordsC(2,4))**two))
               Rphi(3,1) = Rphi(3,1) + half*flux*Le
               Rphi(4,1) = Rphi(4,1) + half*flux*Le
               !
               ! Modify the tangent matrix
               !
               Kphiu(3,5) = Kphiu(3,5) -
     +               (coordsC(1,3)-coordsC(1,4))*flux/(two*Le)
               Kphiu(3,6) = Kphiu(3,6) - 
     +               (coordsC(2,3)-coordsC(2,4))*flux/(two*Le)
               Kphiu(3,7) = Kphiu(3,7) -
     +               (coordsC(1,4)-coordsC(1,3))*flux/(two*Le)
               Kphiu(3,8) = Kphiu(3,8) -
     +               (coordsC(2,4)-coordsC(2,3))*flux/(two*Le)
               Kphiu(4,5) = Kphiu(4,5) -
     +               (coordsC(1,3)-coordsC(1,4))*flux/(two*Le)
               Kphiu(4,6) = Kphiu(4,6) -
     +               (coordsC(2,3)-coordsC(2,4))*flux/(two*Le)
               Kphiu(4,7) = Kphiu(4,7) -
     +               (coordsC(1,4)-coordsC(1,3))*flux/(two*Le)
               Kphiu(4,8) = Kphiu(4,8) - 
     +               (coordsC(2,4)-coordsC(2,3))*flux/(two*Le)
               !
            elseif(face.eq.14) then
               !
               ! charge density (spatial) on face 4 of the element
               !
               ! Modify the displacement residual, loop over nodes
               !
               Le = dsqrt(((coordsC(1,4)-coordsC(1,1))**two) + 
     +                    ((coordsC(2,4)-coordsC(2,1))**two))
               Rphi(4,1) = Rphi(4,1) + half*flux*Le
               Rphi(1,1) = Rphi(1,1) + half*flux*Le
               !
               ! Modify the tangent matrix
               !
               Kphiu(4,7) = Kphiu(4,7) - 
     +               (coordsC(1,4)-coordsC(1,1))*flux/(two*Le)
               Kphiu(4,8) = Kphiu(4,8) - 
     +               (coordsC(2,4)-coordsC(2,1))*flux/(two*Le)
               Kphiu(4,1) = Kphiu(4,1) - 
     +               (coordsC(1,1)-coordsC(1,4))*flux/(two*Le)
               Kphiu(4,2) = Kphiu(4,2) - 
     +               (coordsC(2,1)-coordsC(2,4))*flux/(two*Le)
               Kphiu(1,7) = Kphiu(1,7) - 
     +               (coordsC(1,4)-coordsC(1,1))*flux/(two*Le)
               Kphiu(1,8) = Kphiu(1,8) - 
     +               (coordsC(2,4)-coordsC(2,1))*flux/(two*Le)
               Kphiu(1,1) = Kphiu(1,1) - 
     +               (coordsC(1,1)-coordsC(1,4))*flux/(two*Le)
               Kphiu(1,2) = Kphiu(1,2) - 
     +               (coordsC(2,1)-coordsC(2,4))*flux/(two*Le)
               !
            else
               write(*,*) 'Incorrect dload type',face
               call xit
            endif
            !
          end do
          !
        endif
        !
      else
        !
        ! This is axisymmetric
        !
        if(ndload.gt.0) then
          !
          ! loop over faces and make proper modifications to
          !  residuals and tangents if needed
          !
          do i=1,ndload
            !
            face = jdltyp(i,1) ! label
            flux = adlmag(i,1) ! charge density
            !
	    if(face.eq.1) then
               !
               ! charge density (reference) on face 1 of the element
               !
               ! Modify the displacement residual, loop over nodes
               !
               Le = dsqrt(((coords(1,1)-coords(1,2))**two) + 
     +                    ((coords(2,1)-coords(2,2))**two))
               Rphi(1,1) = Rphi(1,1) + pi*flux*Le*
     +                   third*(two*coords(1,1)+coords(1,2))
               Rphi(2,1) = Rphi(2,1) + pi*flux*Le*
     +                   third*(coords(1,1)+two*coords(1,2))
               !
               ! No modification to the tangent matrix
               !
            elseif(face.eq.2) then
               !
               ! charge density (reference) on face 2 of the element
               !
               ! Modify the displacement residual, loop over nodes
               !
               Le = dsqrt(((coords(1,2)-coords(1,3))**two) + 
     +                    ((coords(2,2)-coords(2,3))**two))
               Rphi(2,1) = Rphi(2,1) + pi*flux*Le*
     +                   third*(two*coords(1,2)+coords(1,3))
               Rphi(3,1) = Rphi(3,1) + pi*flux*Le*
     +                   third*(coords(1,2)+two*coords(1,3))
               !
               ! No modification to the tangent matrix
               !
            elseif(face.eq.3) then
               !
               ! charge density (reference) on face 3 of the element
               !
               ! Modify the displacement residual, loop over nodes
               !
               Le = dsqrt(((coords(1,3)-coords(1,4))**two) + 
     +                    ((coords(2,3)-coords(2,4))**two))
               Rphi(3,1) = Rphi(3,1) + pi*flux*Le*
     +                   third*(two*coords(1,3)+coords(1,4))
               Rphi(4,1) = Rphi(4,1) + pi*flux*Le*
     +                   third*(coords(1,3)+two*coords(1,4))
               !
               ! No modification to the tangent matrix
               !
            elseif(face.eq.4) then
               !
               ! charge density (reference) on face 4 of the element
               !
               ! Modify the displacement residual, loop over nodes
               !
               Le = dsqrt(((coords(1,4)-coords(1,1))**two) + 
     +                    ((coords(2,4)-coords(2,1))**two))
               Rphi(4,1) = Rphi(4,1) + pi*flux*Le*
     +                   third*(two*coords(1,4)+coords(1,1))
               Rphi(1,1) = Rphi(1,1) + pi*flux*Le*
     +                   third*(coords(1,4)+two*coords(1,1))
               !
               ! No modification to the tangent matrix
               !
            elseif(face.eq.11) then
               !
               ! charge density (spatial) on face 1 of the element
               !
               ! Modify the displacement residual, loop over nodes
               !
               Le = dsqrt(((coordsC(1,1)-coordsC(1,2))**two) + 
     +                    ((coordsC(2,1)-coordsC(2,2))**two))
               Rphi(1,1) = Rphi(1,1) + pi*flux*Le*
     +                   third*(two*coordsC(1,1)+coordsC(1,2))
               Rphi(2,1) = Rphi(2,1) + pi*flux*Le*
     +                   third*(coordsC(1,1)+two*coordsC(1,2))
               !
               ! Modify the tangent matrix
               !
               Kphiu(1,1) = Kphiu(1,1) - (pi*flux/Le)*
     +               (coordsC(1,1)-coordsC(1,2))*
     +                third*(two*coordsC(1,1)+coordsC(1,2)) -
     +                pi*flux*Le*third*two
               Kphiu(1,2) = Kphiu(1,2) - (pi*flux/Le)*
     +               (coordsC(2,1)-coordsC(2,2))*
     +                third*(two*coordsC(1,1)+coordsC(1,2))
               Kphiu(1,3) = Kphiu(1,3) - (pi*flux/Le)*
     +               (coordsC(1,2)-coordsC(1,1))*
     +                third*(two*coordsC(1,1)+coordsC(1,2)) -
     +                pi*flux*Le*third*one
               Kphiu(1,4) = Kphiu(1,4) - (pi*flux/Le)*
     +               (coordsC(2,2)-coordsC(2,1))*
     +                third*(two*coordsC(1,1)+coordsC(1,2))
               Kphiu(2,1) = Kphiu(2,1) - (pi*flux/Le)*
     +               (coordsC(1,1)-coordsC(1,2))*
     +                third*(coordsC(1,1)+two*coordsC(1,2)) -
     +                pi*flux*Le*third*one
               Kphiu(2,2) = Kphiu(2,2) - (pi*flux/Le)*
     +               (coordsC(2,1)-coordsC(2,2))*
     +                third*(coordsC(1,1)+two*coordsC(1,2))
               Kphiu(2,3) = Kphiu(2,3) - (pi*flux/Le)*
     +               (coordsC(1,2)-coordsC(1,1))*
     +                third*(coordsC(1,1)+two*coordsC(1,2)) -
     +                pi*flux*Le*third*two
               Kphiu(2,4) = Kphiu(2,4) - (pi*flux/Le)*
     +               (coordsC(2,2)-coordsC(2,1))*
     +                third*(coordsC(1,1)+two*coordsC(1,2))
               !
            elseif(face.eq.12) then
               !
               ! charge density (spatial) on face 2 of the element
               !
               ! Modify the displacement residual, loop over nodes
               !
               Le = dsqrt(((coordsC(1,2)-coordsC(1,3))**two) + 
     +                    ((coordsC(2,2)-coordsC(2,3))**two))
               Rphi(2,1) = Rphi(2,1) + pi*flux*Le*
     +                   third*(two*coordsC(1,2)+coordsC(1,3))
               Rphi(3,1) = Rphi(3,1) + pi*flux*Le*
     +                   third*(coordsC(1,2)+two*coordsC(1,3))
               !
               ! Modify the tangent matrix
               !
               Kphiu(2,3) = Kphiu(2,3) - (pi*flux/Le)*
     +               (coordsC(1,2)-coordsC(1,3))*
     +                third*(two*coordsC(1,2)+coordsC(1,3)) -
     +                pi*flux*Le*third*two
               Kphiu(2,4) = Kphiu(2,4) - (pi*flux/Le)*
     +               (coordsC(2,2)-coordsC(2,3))*
     +                third*(two*coordsC(1,2)+coordsC(1,3))
               Kphiu(2,5) = Kphiu(2,5) - (pi*flux/Le)*
     +               (coordsC(1,3)-coordsC(1,2))*
     +                third*(two*coordsC(1,2)+coordsC(1,3)) -
     +                pi*flux*Le*third*one
               Kphiu(2,6) = Kphiu(2,6) - (pi*flux/Le)*
     +               (coordsC(2,3)-coordsC(2,2))*
     +                third*(two*coordsC(1,2)+coordsC(1,3))
               Kphiu(3,3) = Kphiu(3,3) - (pi*flux/Le)*
     +               (coordsC(1,2)-coordsC(1,3))*
     +                third*(coordsC(1,2)+two*coordsC(1,3)) -
     +                pi*flux*Le*third*one
               Kphiu(3,4) = Kphiu(3,4) - (pi*flux/Le)*
     +               (coordsC(2,2)-coordsC(2,3))*
     +                third*(coordsC(1,2)+two*coordsC(1,3))
               Kphiu(3,5) = Kphiu(3,5) - (pi*flux/Le)*
     +               (coordsC(1,3)-coordsC(1,2))*
     +                third*(coordsC(1,2)+two*coordsC(1,3)) -
     +                pi*flux*Le*third*two
               Kphiu(3,6) = Kphiu(3,6) - (pi*flux/Le)*
     +               (coordsC(2,3)-coordsC(2,2))*
     +                third*(coordsC(1,2)+two*coordsC(1,3))
               !
            elseif(face.eq.13) then
               !
               ! charge density (spatial) on face 3 of the element
               !
               ! Modify the displacement residual, loop over nodes
               !
               Le = dsqrt(((coordsC(1,3)-coordsC(1,4))**two) + 
     +                    ((coordsC(2,3)-coordsC(2,4))**two))
               Rphi(3,1) = Rphi(3,1) + pi*flux*Le*
     +                   third*(two*coordsC(1,3)+coordsC(1,4))
               Rphi(4,1) = Rphi(4,1) + pi*flux*Le*
     +                   third*(coordsC(1,3)+two*coordsC(1,4))
               !
               ! Modify the tangent matrix
               !
               Kphiu(3,5) = Kphiu(3,5) - (pi*flux/Le)*
     +               (coordsC(1,3)-coordsC(1,4))*
     +                third*(two*coordsC(1,3)+coordsC(1,4)) -
     +                pi*flux*Le*third*two
               Kphiu(3,6) = Kphiu(3,6) - (pi*flux/Le)*
     +               (coordsC(2,3)-coordsC(2,4))*
     +                third*(two*coordsC(1,3)+coordsC(1,4))
               Kphiu(3,7) = Kphiu(3,7) - (pi*flux/Le)*
     +               (coordsC(1,4)-coordsC(1,3))*
     +                third*(two*coordsC(1,3)+coordsC(1,4)) -
     +                pi*flux*Le*third*one
               Kphiu(3,8) = Kphiu(3,8) - (pi*flux/Le)*
     +               (coordsC(2,4)-coordsC(2,3))*
     +                third*(two*coordsC(1,3)+coordsC(1,4))
               Kphiu(4,5) = Kphiu(4,5) - (pi*flux/Le)*
     +               (coordsC(1,3)-coordsC(1,4))*
     +                third*(coordsC(1,3)+two*coordsC(1,4)) -
     +                pi*flux*Le*third*one
               Kphiu(4,6) = Kphiu(4,6) - (pi*flux/Le)*
     +               (coordsC(2,3)-coordsC(2,4))*
     +                third*(coordsC(1,3)+two*coordsC(1,4))
               Kphiu(4,7) = Kphiu(4,7) - (pi*flux/Le)*
     +               (coordsC(1,4)-coordsC(1,3))*
     +                third*(coordsC(1,3)+two*coordsC(1,4)) -
     +                pi*flux*Le*third*two
               Kphiu(4,8) = Kphiu(4,8) - (pi*flux/Le)*
     +               (coordsC(2,4)-coordsC(2,3))*
     +                third*(coordsC(1,3)+two*coordsC(1,4))
               !
            elseif(face.eq.14) then
               !
               ! charge density (spatial) on face 4 of the element
               !
               ! Modify the displacement residual, loop over nodes
               !
               Le = dsqrt(((coordsC(1,4)-coordsC(1,1))**two) + 
     +                    ((coordsC(2,4)-coordsC(2,1))**two))
               Rphi(4,1) = Rphi(4,1) + pi*flux*Le*
     +                   third*(two*coordsC(1,4)+coordsC(1,1))
               Rphi(1,1) = Rphi(1,1) + pi*flux*Le*
     +                   third*(coordsC(1,4)+two*coordsC(1,1))
               !
               ! Modify the tangent matrix
               !
               Kphiu(4,7) = Kphiu(4,7) - (pi*flux/Le)*
     +               (coordsC(1,4)-coordsC(1,1))*
     +                third*(two*coordsC(1,4)+coordsC(1,1)) -
     +                pi*flux*Le*third*two
               Kphiu(4,8) = Kphiu(4,8) - (pi*flux/Le)*
     +               (coordsC(2,4)-coordsC(2,1))*
     +                third*(two*coordsC(1,4)+coordsC(1,1))
               Kphiu(4,1) = Kphiu(4,1) - (pi*flux/Le)*
     +               (coordsC(1,1)-coordsC(1,4))*
     +                third*(two*coordsC(1,4)+coordsC(1,1)) -
     +                pi*flux*Le*third*one
               Kphiu(4,2) = Kphiu(4,2) - (pi*flux/Le)*
     +               (coordsC(2,1)-coordsC(2,4))*
     +                third*(two*coordsC(1,4)+coordsC(1,1))
               Kphiu(1,7) = Kphiu(1,7) - (pi*flux/Le)*
     +               (coordsC(1,4)-coordsC(1,1))*
     +                third*(coordsC(1,4)+two*coordsC(1,1)) -
     +                pi*flux*Le*third*one
               Kphiu(1,8) = Kphiu(1,8) - (pi*flux/Le)*
     +               (coordsC(2,4)-coordsC(2,1))*
     +                third*(coordsC(1,4)+two*coordsC(1,1))
               Kphiu(1,1) = Kphiu(1,1) - (pi*flux/Le)*
     +               (coordsC(1,1)-coordsC(1,4))*
     +                third*(coordsC(1,4)+two*coordsC(1,1)) -
     +                pi*flux*Le*third*two
               Kphiu(1,2) = Kphiu(1,2) - (pi*flux/Le)*
     +               (coordsC(2,1)-coordsC(2,4))*
     +                third*(coordsC(1,4)+two*coordsC(1,1))
               !
            else
               write(*,*) 'Incorrect dload type',face
               call xit
            endif
            !
          end do
          !
        endif
        !
      endif
      !
      ! End loop over flux and traction terms
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the Stiffness matrix.  This
      !  is essentially giving Abaqus the residual and the tangent matrix.
      !
      ! Return Abaqus the right hand side vector
      !
      do i=1,nNode
         A11 = (nDim+1)*(i-1)+1
         A12 = nDim*(i-1)+1
         !
         ! displacement
         !
         rhs(A11,1) = Ru(A12,1)
         rhs(A11+1,1) = Ru(A12+1,1)
         !
         ! electric potential
         !
         rhs(A11+2,1) = Rphi(i,1)
      enddo
      !
      ! Return Abaqus the tangent matrix
      !
      amatrx = zero
      do i=1,nNode
         do j=1,nNode
            A11 = (nDim+1)*(i-1)+1
            A12 = nDim*(i-1)+1
            B11 = (nDim+1)*(j-1)+1
            B12 = nDim*(j-1)+1
            !
            ! displacement
            !
            amatrx(A11,B11) = Kuu(A12,B12)
            amatrx(A11,B11+1) = Kuu(A12,B12+1)
            amatrx(A11+1,B11) = Kuu(A12+1,B12)
            amatrx(A11+1,B11+1) = Kuu(A12+1,B12+1)
            !
            ! electric potential
            !
            amatrx(A11+2,B11+2) = Kphiphi(i,j)
            !
            ! displacement/electric potential
            !
            amatrx(A11,B11+2) = Kuphi(A12,j)
            amatrx(A11+1,B11+2) = Kuphi(A12+1,j)
            !
            ! electric potential/displacement
            !
            amatrx(A11+2,B11) = Kphiu(i,B12)
            amatrx(A11+2,B11+1) = Kphiu(i,B12+1)
            !
         enddo
      enddo
      !
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------


      return
      end subroutine uel

!************************************************************************
!     Material subroutine
!************************************************************************

      subroutine Gent(props,nprops,F_tau,ER,T_tau,D_tau,
     +        SpUUMod,SpUPhiMod,SpPhiUMod,SpPhiPhiMod,stat)

      implicit none
      !
      integer i,j,k,l,m,n,nprops,stat
      !
      real*8 props(nprops),F_tau(3,3),ER(3,1),T_tau(3,3),D_tau(3,1),
     +  SpUUMod(3,3,3,3),SpUPhiMod(3,3,3),SpPhiUMod(3,3,3),
     +  SpPhiPhiMod(3,3),Iden(3,3),Gshear,Kbulk,Imax,permit,detF,
     +  Finv(3,3),FT(3,3),FinvT(3,3),C_tau(3,3),Cinv(3,3),detC,trC,
     +  I1bar,fac,GShearGent,E(3,1),vec(3,1),I6,TR_tau(3,3),dGdF(3,3),
     +  dTRdF(3,3,3,3)
      !
      real*8 zero,one,two,three,fourth,third,half
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,fourth=1.d0/4.d0,
     +     third=1.d0/3.d0,half=1.d0/2.d0)
 

      ! Identity matrix
      !
      call onem(Iden)
 

      ! Obtain relevant material properties
      !
      Gshear = props(1) ! Shear modulus
      Kbulk  = props(2) ! Bulk modulus
      Imax   = props(3) ! Max I1bar
      permit = props(4) ! Permittivity


      ! Compute the determinant of the deformation gradient
      !
      call mdet(F_tau,detF)


      ! Compute the determinant, the inverse, the transpose, 
      !  and the inverse transpose of the deformation gradient
      !
      call matInv3D(F_tau,Finv,detF,stat)
      FT = transpose(F_tau)
      FinvT = transpose(Finv)


      ! Compute the right Cauchy-Green tensor, its inverse, 
      !  its determinant, and its trace
      !
      C_tau = matmul(transpose(F_tau),F_tau)
      call matInv3D(C_tau,Cinv,detC,stat)
      trC = C_tau(1,1) + C_tau(2,2) + C_tau(3,3)
 

      ! Compute the trace of the distortional right Cauchy-Green tensor
      !
      I1bar = (detF**(-two/three))*trC
 

      ! Compute the ``I-3'' factor appearing in the Gent model
      !
      fac = (I1bar - three)/Imax
      if(fac.gt.0.95d0) fac = 0.95d0
      fac = one/(one - fac)
 

      ! Compute the ``shear'' modulus. Note: this is not really the shear
      !  modulus, but it will help when computing the material tangent later
      !
      GShearGent = Gshear*fac
 

      ! Compute the derivative of the ``shear modulus'' with respect
      !  to the deformation gradient for use in the material tangent
      !
      dGdF = two*(Gshear/Imax)*(detF**(-two/three))*
     +     fac*fac*(F_tau - third*trC*FinvT)
      
      
      ! Compute the current electric field
      !
      E = matmul(FinvT,ER)
      
      
      ! Compute quantities related to electric field
      !
      vec = matmul(Cinv,ER)
      I6 = E(1,1)*E(1,1) + E(2,1)*E(2,1) + E(3,1)*E(3,1)
 

      ! Compute the 1st Piola stress
      !
      TR_tau = (detF**(-two/three))*GShearGent*(F_tau-third*trC*FinvT)
     +     + Kbulk*detF*(detF - one)*FinvT
     +     + permit*detF*(matmul(E,transpose(matmul(Cinv,ER))) 
     +                    - half*I6*FinvT)
 

      ! Compute the Cauchy stress
      !
      T_tau = (one/detF)*matmul(TR_tau,transpose(F_tau))


      ! Compute the current electric displacement
      !
      D_tau = permit*E


      ! Calculate the material tangent modulus
      !
      dTRdF = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  dTRdF(i,j,k,l) = dTRdF(i,j,k,l)
     +                 + (detF**(-two/three))*dGdF(k,l)*
     +                 (
     +                 F_tau(i,j) - third*trC*Finv(j,i)
     +                 )
     +                 + (detF**(-two/three))*GshearGent*
     +                 (
     +                 (-two/three)*F_tau(i,j)*Finv(l,k)
     +                 + (two/9.d0)*trC*Finv(j,i)*Finv(l,k)
     +                 + Iden(i,k)*Iden(j,l)
     +                 + third*trC*Finv(l,i)*Finv(j,k)
     +                 - (two/three)*Finv(j,i)*F_tau(k,l)
     +                 )
     +                 + detF*Kbulk*
     +                 (
     +                 (detF-one)*Finv(j,i)*Finv(l,k)
     +                 + detF*Finv(j,i)*Finv(l,k)
     +                 - (detF-one)*Finv(l,i)*Finv(j,k)
     +                 )
     +                 + detF*permit*
     +                 (
     +                 - Finv(l,i)*vec(j,1)*E(k,1) 
     +                 - E(i,1)*Finv(j,k)*vec(l,1)
     +                 - E(i,1)*E(k,1)*Cinv(j,l)
     +                 + Finv(j,i)*E(k,1)*vec(l,1)
     +                 + E(i,1)*vec(j,1)*Finv(l,k)
     +                 + half*I6*Finv(l,i)*Finv(j,k)
     +                 - half*I6*Finv(j,i)*Finv(l,k)
     +                 )
               enddo
            enddo
         enddo
      enddo
      
      
      ! Calculate the spatial tangent modulus
      !
      SpUUMod = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  do m=1,3
                     do n=1,3
                        SpUUMod(i,j,k,l) = SpUUMod(i,j,k,l) +
     +                      (dTRdF(i,m,k,n)*F_tau(j,m)*F_tau(l,n))/detF
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo


      ! Calculate the spatial stress/electric potential modulus modulus
      !
      SpUPhiMod = zero
      do i=1,3
         do j=1,3
            do l=1,3
               SpUPhiMod(i,j,l) = SpUPhiMod(i,j,l) +
     +               permit*(iden(i,l)*E(j,1) + 
     +                       E(i,1)*iden(j,l) - 
     +                       iden(i,j)*E(l,1))
            end do
         end do
      end do


      ! Calculate the spatial electric displacement/electric potential modulus
      !
      SpPhiPhiMod = permit*iden


      ! Calculate the spatial electric displacement/strain modulus
      !
      SpPhiUMod = zero
      do j=1,3
         do k=1,3
            do l=1,3
               SpPhiUMod(j,k,l) = SpPhiUMod(j,k,l) - 
     +               permit*(iden(j,k)*E(l,1) +
     +                       E(k,1)*iden(j,l) -
     +                       E(j,1)*iden(l,k))
            end do
         end do
      end do


      return
      end subroutine Gent 
      
!****************************************************************************
!     Element subroutines
!****************************************************************************

      subroutine xint2D1pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 1 gauss point for integration
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(1,2), w(1)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 1


      ! Gauss weights
      !
      w = 4.d0
      

      ! Gauss pt location in master element
      !
      xi(1,1) = 0.d0
      xi(1,2) = 0.d0


      return
      end subroutine xint2D1pt
      
!************************************************************************

      subroutine xint2D4pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 4 gauss points for integration
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(4,2), w(4)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 4


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)


      return
      end subroutine xint2D4pt
      
!************************************************************************

      subroutine calcShape2DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      !
      implicit none
      !
      integer intpt,nDim,nIntPt
      !
      real*8 xi_int(nIntPt,2),sh(4),dshxi(4,2),xi,eta
      !
      real*8 zero,one,fourth
      parameter(zero=0.d0,one=1.d0,fourth=1.d0/4.d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      
      
      ! The shape functions
      !
      sh(1) = fourth*(one - xi)*(one - eta)
      sh(2) = fourth*(one + xi)*(one - eta)
      sh(3) = fourth*(one + xi)*(one + eta)
      sh(4) = fourth*(one - xi)*(one + eta)
      
      
      ! The first derivatives
      !
      dshxi(1,1) = -fourth*(one - eta)
      dshxi(1,2) = -fourth*(one - xi)
      dshxi(2,1) = fourth*(one - eta)
      dshxi(2,2) = -fourth*(one + xi)
      dshxi(3,1) = fourth*(one + eta)
      dshxi(3,2) = fourth*(one + xi)
      dshxi(4,1) = -fourth*(one + eta)
      dshxi(4,2) = fourth*(one - xi)
      

      return
      end subroutine calcShape2DLinear

!************************************************************************

      subroutine mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(3,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))
      

      return
      end subroutine mapShape2D

!*************************************************************************

      subroutine mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      ! This subroutine is exactly the same as the regular mapShape2D
      !  with the exception that coords(2,nNode) here and coords(3,nNode)
      !  in the regular.  I have noticed that a "heat transfer" and 
      !  "static" step uses MCRD=2, but for "coupled-temperature-displacement"
      !  you will get MCRD=3, even for a plane analysis.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(2,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))
      

      return
      end subroutine mapShape2Da

!****************************************************************************
!     Utility subroutines
!****************************************************************************

      subroutine matInv3D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(3,3),A_inv(3,3),det_A,det_A_inv


      istat = 1
      
      det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
      
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv3D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
          
      det_A_inv = 1.d0/det_A
        
      A_inv(1,1) = det_A_inv*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_inv(1,2) = det_A_inv*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_inv(1,3) = det_A_inv*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_inv(2,1) = det_A_inv*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_inv(2,2) = det_A_inv*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_inv(2,3) = det_A_inv*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_inv(3,1) = det_A_inv*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_inv(3,2) = det_A_inv*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_inv(3,3) = det_A_inv*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
      

      return
      end subroutine matInv3D

!****************************************************************************

      subroutine matInv2D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse, and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(2,2),A_inv(2,2),det_A,det_A_inv

      
      istat = 1
      
      det_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)
        
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv2D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
            
      det_A_inv = 1.d0/det_A
          
      A_inv(1,1) =  det_A_inv*A(2,2)
      A_inv(1,2) = -det_A_inv*A(1,2)
      A_inv(2,1) = -det_A_inv*A(2,1)
      A_inv(2,2) =  det_A_inv*A(1,1)


      return
      end subroutine matInv2D

!****************************************************************************

      subroutine mdet(A,det)
      !
      ! This subroutine calculates the determinant
      ! of a 3 by 3 matrix [A]
      !
      implicit none
      !
      real*8  A(3,3),det


      det = A(1,1)*A(2,2)*A(3,3) 
     +	  + A(1,2)*A(2,3)*A(3,1)
     +	  + A(1,3)*A(2,1)*A(3,2)
     +	  - A(3,1)*A(2,2)*A(1,3)
     +	  - A(3,2)*A(2,3)*A(1,1)
     +	  - A(3,3)*A(2,1)*A(1,2)


      return
      end subroutine mdet
	
!****************************************************************************

      subroutine onem(A)
      !
      ! This subroutine stores the identity matrix in the
      ! 3 by 3 matrix [A]
      !
      implicit none
      !
      integer i,j
      !
      real*8 A(3,3)


      do i=1,3
         do J=1,3
	    if (i .eq. j) then
              A(i,j) = 1.0
            else
              A(i,j) = 0.0
            end if
         end do
      end do


      return
      end subroutine onem

!****************************************************************************
