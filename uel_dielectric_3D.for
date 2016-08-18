!************************************************************************
! User element (UEL) for coupled large-deformation elasticity and  
!  dielectric behavior in three-dimensions
!************************************************************************
! Element details:
!************************************************************************
!
! Solution variables (or nodal variables) are the displacements (DOFs 1-3)
!  and the electric potential (DOF 11).
!
! Material behavior is Gent rubber elasticity and ideal 
!  (constant permittivity) dielectric behavior.
! 
! This subroutine is for a three-dimensional 8-node isoparametric
!  brick element as shown below with 8pt (full) or 1pt (reduced) 
!  integration.
!
! In order to avoid locking for the fully-integrated element, we
!  use the F-bar method of de Souza Neto (1996).
!
! Surface charge density boundary conditions are NOT supported in 
!  this element at this time.  Mechanical, traction- and pressure-type 
!  boundary conditions may be applied to the dummy mesh using the 
!  Abaqus built-in commands *Dload or *Dsload.
!
!
!  8-node     8-----------7
!  brick     /|          /|       zeta
!           / |         / |       
!          5-----------6  |       |     eta
!          |  |        |  |       |   /
!          |  |        |  |       |  /
!          |  4--------|--3       | /
!          | /         | /        |/
!          |/          |/         O--------- xi
!          1-----------2        origin at cube center
!
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
!  *User Element,Nodes=8,Type=U1,Iproperties=0,Properties=4,Coordinates=3,Variables=1,Unsymm
!  1,2,3,11
!
! Note: No local state variables are used in this element, thus we may set the above 
!  parameter 'Variables' to any non-zero integer.
!
! In the subroutine UEL, set 'nInt' = number of integration points
!  Options are nInt=8 (full integration) or nInt=1 (reduced integration)
!
!
!     Material Properties Vector
!     --------------------------------------------------------------
!     Gshear  = props(1)  ! Shear modulus
!     Kbulk   = props(2)  ! Bulk modulus
!     Imax    = props(3)  ! Max value of (I1bar-3) (Gent parameter)
!     epsilon = props(4)  ! Permittivity
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
      integer i,j,k,A11,B11,A12,B12,nInt,nIntPt,intpt,nDim,stat
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      parameter(nDim=3) ! number of spatial dimensions, do not change
      parameter(nInt=8) ! number of integration points
      !
      ! nInt=8: fully-integrated, nInt=1: reduced integration
      ! 
      ! When reduced integration is selected, be sure to
      !  specify an hourglass stiffness of about 0.005G
      !  for the dummy mesh in the input file.  Also, change
      !  the dummy mesh element to reduced.
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      real*8 u(nNode,3),du(nNode,ndofel),phi(nNode),coordsC(mcrd,nNode),
     +  Ru(3*nNode,1),Rphi(nNode,1),Kuu(3*nNode,3*nNode),
     +  Kuphi(3*nNode,nNode),Kphiu(nNode,3*nNode),Kphiphi(nNode,nNode),
     +  Iden(3,3),xi(nInt,3),w(nInt),sh0(nNode),sh(nNode),dsh0(nNode,3),
     +  dshC0(nNode,3),dsh(nNode,3),dshC(nNode,3),dshxi(nNode,3),
     +  detMapJ0,detMapJ0C,detMapJ,detMapJC,Fc_tau(3,3),detFc_tau,
     +  F_tau(3,3),detF_tau,ER(3,1),T_tau(3,3),D_tau(3,1),
     +  SpUUMod(3,3,3,3),SpUPhiMod(3,3,3),SpPhiUMod(3,3,3),
     +  SpPhiPhiMod(3,3),Smat(6,1),Bmat(6,3*nNode),Gmat(9,3*nNode),
     +  G0mat(9,3*nNode),Amat(9,9),Qmat(9,9),AmatPhiU(3,9),
     +  AmatUPhi(9,3),Gmatphi(3,nNode),body(3),BodyForceRes(3*nNode,1),
     +  bodyCharge,BodyChargeRes(nNode,1),Le
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


      ! Identity tensor
      !
      call onem(Iden)


      ! Initialize the residual and tangent matrices to zero.
      !
      Ru  = zero
      Rphi  = zero
      Kuu = zero
      Kuphi = zero
      Kphiphi = zero
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
      Le = dsqrt(((coordsC(1,1)-coordsC(1,7))**two) + 
     +           ((coordsC(2,1)-coordsC(2,7))**two) +
     +           ((coordsC(3,1)-coordsC(3,7))**two))
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
      !  centroid.  Get the deformation gradient for use in the `F-bar' method.
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
      if(nNode.eq.8) then
         call calcShape3DLinear(1,zero,1,sh0,dshxi)
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.8'
         call xit
      endif


      ! Map shape functions from local to global reference coordinate system
      !
      call mapShape3D(nNode,dshxi,coords,dsh0,detMapJ0,stat)
      if(stat.eq.0) then
         pnewdt = 0.5
         return
      endif


      ! Map shape functions from local to global current coordinate system
      !
      call mapShape3D(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
      if(stat.eq.0) then
         pnewdt = 0.5
         return
      endif


      ! Calculate the deformation gradient at the element centroid
      !  at the end of the increment for use in the `F-bar' method
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
      call mdet(Fc_tau,detFc_tau)
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
      if(nNode.eq.8) then
         if(nInt.eq.1) then
            call xint3D1pt(xi,w,nIntPt) ! 1-pt integration, nInt=1 above
         elseif(nInt.eq.8) then
            call xint3D8pt(xi,w,nIntPt) ! 8-pt integration, nInt=8 above
         else
            write(*,*) 'Invalid number of int points, nInt=',nInt
            call xit
         endif
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.8'
         call xit
      endif


      ! Loop over integration points
      !
      do intpt=1,nIntPt


         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.8) then
            call calcShape3DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.8'
            call xit
         endif


         ! Map shape functions from local to global reference coordinate system
         !
         call mapShape3D(nNode,dshxi,coords,dsh,detMapJ,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif


         ! Map shape functions from local to global current coordinate system
         !
         call mapShape3D(nNode,dshxi,coordsC,dshC,detMapJC,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
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
         !  The subscript tau denotes the time at the end of the increment.
         !
         F_tau = Iden
         do i=1,nDim
            do j=1,nDim
               do k=1,nNode
                  F_tau(i,j) = F_tau(i,j) + dsh(k,j)*u(k,i)
               enddo
            enddo
         enddo
         
         
         ! Modify the deformation gradient for the `F-bar' method
         !  only when using the 8 node fully integrated linear
         !  element, do not use the `F-bar' method for reduced element
         !
         if((nNode.eq.8).and.(nInt.eq.8)) then
            call mdet(F_tau,detF_tau)
            F_tau = ((detFc_tau/detF_tau)**third)*F_tau
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
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


         ! Compute/update the displacement residual vector
         !
         Smat(1,1) = T_tau(1,1)
         Smat(2,1) = T_tau(2,2)
         Smat(3,1) = T_tau(3,3)
         Smat(4,1) = T_tau(1,2)
         Smat(5,1) = T_tau(2,3)
         Smat(6,1) = T_tau(1,3)
         !
         Bmat = zero
         do k=1,nNode
            Bmat(1,1+nDim*(k-1)) = dshC(k,1)
            Bmat(2,2+nDim*(k-1)) = dshC(k,2)
            Bmat(3,3+nDim*(k-1)) = dshC(k,3)
            Bmat(4,1+nDim*(k-1)) = dshC(k,2)
            Bmat(4,2+nDim*(k-1)) = dshC(k,1)
            Bmat(5,2+nDim*(k-1)) = dshC(k,3)
            Bmat(5,3+nDim*(k-1)) = dshC(k,2)
            Bmat(6,1+nDim*(k-1)) = dshC(k,3)
            Bmat(6,3+nDim*(k-1)) = dshC(k,1)
         enddo
         !
         body = zero ! The body force vector may be specified here
         !
         BodyForceRes = zero
         do k=1,nNode
            BodyForceRes(1+nDim*(k-1),1) = sh(k)*body(1)
            BodyForceRes(2+nDim*(k-1),1) = sh(k)*body(2)
            BodyForceRes(3+nDim*(k-1),1) = sh(k)*body(3)
         enddo
         !
         Ru = Ru + matmul(transpose(Bmat),Smat)*detmapJC*w(intpt)
     +        - BodyForceRes*detmapJC*w(intpt)


         ! Compute/update the electric potential residual vector
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
         Rphi = Rphi + matmul(transpose(Gmatphi),D_tau)
     +                 *detmapJC*w(intpt) 
     +               + BodyChargeRes*detmapJC*w(intpt)


         ! Compute/update the displacement tangent matrix
         !
         Gmat = zero
         do k=1,nNode
            Gmat(1,1+nDim*(k-1)) = dshC(k,1)
            Gmat(2,2+nDim*(k-1)) = dshC(k,1)
            Gmat(3,3+nDim*(k-1)) = dshC(k,1)
            Gmat(4,1+nDim*(k-1)) = dshC(k,2)
            Gmat(5,2+nDim*(k-1)) = dshC(k,2)
            Gmat(6,3+nDim*(k-1)) = dshC(k,2)
            Gmat(7,1+nDim*(k-1)) = dshC(k,3)
            Gmat(8,2+nDim*(k-1)) = dshC(k,3)
            Gmat(9,3+nDim*(k-1)) = dshC(k,3)
         enddo
         !
         G0mat = zero
         do k=1,nNode
            G0mat(1,1+nDim*(k-1)) = dshC0(k,1)
            G0mat(2,2+nDim*(k-1)) = dshC0(k,1)
            G0mat(3,3+nDim*(k-1)) = dshC0(k,1)
            G0mat(4,1+nDim*(k-1)) = dshC0(k,2)
            G0mat(5,2+nDim*(k-1)) = dshC0(k,2)
            G0mat(6,3+nDim*(k-1)) = dshC0(k,2)
            G0mat(7,1+nDim*(k-1)) = dshC0(k,3)
            G0mat(8,2+nDim*(k-1)) = dshC0(k,3)
            G0mat(9,3+nDim*(k-1)) = dshC0(k,3)
         enddo
         !
         Amat = zero
         Amat(1,1) = SpUUMod(1,1,1,1)
         Amat(1,2) = SpUUMod(1,1,2,1)
         Amat(1,3) = SpUUMod(1,1,3,1)
         Amat(1,4) = SpUUMod(1,1,1,2)
         Amat(1,5) = SpUUMod(1,1,2,2)
         Amat(1,6) = SpUUMod(1,1,3,2)
         Amat(1,7) = SpUUMod(1,1,1,3)
         Amat(1,8) = SpUUMod(1,1,2,3)
         Amat(1,9) = SpUUMod(1,1,3,3)
         Amat(2,1) = SpUUMod(2,1,1,1)
         Amat(2,2) = SpUUMod(2,1,2,1)
         Amat(2,3) = SpUUMod(2,1,3,1)
         Amat(2,4) = SpUUMod(2,1,1,2)
         Amat(2,5) = SpUUMod(2,1,2,2)
         Amat(2,6) = SpUUMod(2,1,3,2)
         Amat(2,7) = SpUUMod(2,1,1,3)
         Amat(2,8) = SpUUMod(2,1,2,3)
         Amat(2,9) = SpUUMod(2,1,3,3)
         Amat(3,1) = SpUUMod(3,1,1,1)
         Amat(3,2) = SpUUMod(3,1,2,1)
         Amat(3,3) = SpUUMod(3,1,3,1)
         Amat(3,4) = SpUUMod(3,1,1,2)
         Amat(3,5) = SpUUMod(3,1,2,2)
         Amat(3,6) = SpUUMod(3,1,3,2)
         Amat(3,7) = SpUUMod(3,1,1,3)
         Amat(3,8) = SpUUMod(3,1,2,3)
         Amat(3,9) = SpUUMod(3,1,3,3)
         Amat(4,1) = SpUUMod(1,2,1,1)
         Amat(4,2) = SpUUMod(1,2,2,1)
         Amat(4,3) = SpUUMod(1,2,3,1)
         Amat(4,4) = SpUUMod(1,2,1,2)
         Amat(4,5) = SpUUMod(1,2,2,2)
         Amat(4,6) = SpUUMod(1,2,3,2)
         Amat(4,7) = SpUUMod(1,2,1,3)
         Amat(4,8) = SpUUMod(1,2,2,3)
         Amat(4,9) = SpUUMod(1,2,3,3)
         Amat(5,1) = SpUUMod(2,2,1,1)
         Amat(5,2) = SpUUMod(2,2,2,1)
         Amat(5,3) = SpUUMod(2,2,3,1)
         Amat(5,4) = SpUUMod(2,2,1,2)
         Amat(5,5) = SpUUMod(2,2,2,2)
         Amat(5,6) = SpUUMod(2,2,3,2)
         Amat(5,7) = SpUUMod(2,2,1,3)
         Amat(5,8) = SpUUMod(2,2,2,3)
         Amat(5,9) = SpUUMod(2,2,3,3)
         Amat(6,1) = SpUUMod(3,2,1,1)
         Amat(6,2) = SpUUMod(3,2,2,1)
         Amat(6,3) = SpUUMod(3,2,3,1)
         Amat(6,4) = SpUUMod(3,2,1,2)
         Amat(6,5) = SpUUMod(3,2,2,2)
         Amat(6,6) = SpUUMod(3,2,3,2)
         Amat(6,7) = SpUUMod(3,2,1,3)
         Amat(6,8) = SpUUMod(3,2,2,3)
         Amat(6,9) = SpUUMod(3,2,3,3)
         Amat(7,1) = SpUUMod(1,3,1,1)
         Amat(7,2) = SpUUMod(1,3,2,1)
         Amat(7,3) = SpUUMod(1,3,3,1)
         Amat(7,4) = SpUUMod(1,3,1,2)
         Amat(7,5) = SpUUMod(1,3,2,2)
         Amat(7,6) = SpUUMod(1,3,3,2)
         Amat(7,7) = SpUUMod(1,3,1,3)
         Amat(7,8) = SpUUMod(1,3,2,3)
         Amat(7,9) = SpUUMod(1,3,3,3)
         Amat(8,1) = SpUUMod(2,3,1,1)
         Amat(8,2) = SpUUMod(2,3,2,1)
         Amat(8,3) = SpUUMod(2,3,3,1)
         Amat(8,4) = SpUUMod(2,3,1,2)
         Amat(8,5) = SpUUMod(2,3,2,2)
         Amat(8,6) = SpUUMod(2,3,3,2)
         Amat(8,7) = SpUUMod(2,3,1,3)
         Amat(8,8) = SpUUMod(2,3,2,3)
         Amat(8,9) = SpUUMod(2,3,3,3)
         Amat(9,1) = SpUUMod(3,3,1,1)
         Amat(9,2) = SpUUMod(3,3,2,1)
         Amat(9,3) = SpUUMod(3,3,3,1)
         Amat(9,4) = SpUUMod(3,3,1,2)
         Amat(9,5) = SpUUMod(3,3,2,2)
         Amat(9,6) = SpUUMod(3,3,3,2)
         Amat(9,7) = SpUUMod(3,3,1,3)
         Amat(9,8) = SpUUMod(3,3,2,3)
         Amat(9,9) = SpUUMod(3,3,3,3)
         !
         Qmat = zero
         Qmat(1,1) = third*(Amat(1,1)+Amat(1,5)+Amat(1,9)) 
     +        - (two/three)*T_tau(1,1)
         Qmat(2,1) = third*(Amat(2,1)+Amat(2,5)+Amat(2,9))
     +        - (two/three)*T_tau(2,1)
         Qmat(3,1) = third*(Amat(3,1)+Amat(3,5)+Amat(3,9))
     +        - (two/three)*T_tau(3,1)
         Qmat(4,1) = third*(Amat(4,1)+Amat(4,5)+Amat(4,9))
     +        - (two/three)*T_tau(1,2)
         Qmat(5,1) = third*(Amat(5,1)+Amat(5,5)+Amat(5,9))
     +        - (two/three)*T_tau(2,2)
         Qmat(6,1) = third*(Amat(6,1)+Amat(6,5)+Amat(6,9))
     +        - (two/three)*T_tau(3,2)
         Qmat(7,1) = third*(Amat(7,1)+Amat(7,5)+Amat(7,9))
     +        - (two/three)*T_tau(1,3)
         Qmat(8,1) = third*(Amat(8,1)+Amat(8,5)+Amat(8,9))
     +        - (two/three)*T_tau(2,3)
         Qmat(9,1) = third*(Amat(9,1)+Amat(9,5)+Amat(9,9))
     +        - (two/three)*T_tau(3,3)
         Qmat(1,5) = Qmat(1,1)
         Qmat(2,5) = Qmat(2,1)
         Qmat(3,5) = Qmat(3,1)
         Qmat(4,5) = Qmat(4,1)
         Qmat(5,5) = Qmat(5,1)
         Qmat(6,5) = Qmat(6,1)
         Qmat(7,5) = Qmat(7,1)
         Qmat(8,5) = Qmat(8,1)
         Qmat(9,5) = Qmat(9,1)
         Qmat(1,9) = Qmat(1,1)
         Qmat(2,9) = Qmat(2,1)
         Qmat(3,9) = Qmat(3,1)
         Qmat(4,9) = Qmat(4,1)
         Qmat(5,9) = Qmat(5,1)
         Qmat(6,9) = Qmat(6,1)
         Qmat(7,9) = Qmat(7,1)
         Qmat(8,9) = Qmat(8,1)
         Qmat(9,9) = Qmat(9,1)
         !
         if((nNode.eq.8).and.(nInt.eq.8)) then
            !
            ! This is the tangent using the F-bar method with the
            !  8 node fully-integrated element
            !
            Kuu = Kuu
     +           - matmul(matmul(transpose(Gmat),Amat),Gmat)
     +           *detMapJC*w(intpt)
     +           - matmul(transpose(Gmat),matmul(Qmat,
     +           (G0mat-Gmat)))*detMapJC*w(intpt)
         else
            !
            ! This is the tangent NOT using the F-bar method with all
            !  other elements
            !
            Kuu = Kuu -
     +           matmul(matmul(transpose(Gmat),Amat),Gmat)
     +           *detMapJC*w(intpt)
         endif


         ! Compute/update the electric potential tangent matrix
         !
         Kphiphi = Kphiphi - 
     +             matmul(matmul(transpose(Gmatphi),SpPhiPhiMod),
     +                    Gmatphi)*detMapJC*w(intpt)


         ! Compute/update the electric potential/displacement tangent matrix
         !   Strictly speaking, use of 'F-bar' will affect this tangent;
         !   however, the effect is expected to be small and is neglected here.
         !
         AmatPhiU = zero
         AmatPhiU(1,1) = SpPhiUMod(1,1,1)
         AmatPhiU(1,2) = SpPhiUMod(1,2,1)
         AmatPhiU(1,3) = SpPhiUMod(1,3,1)
         AmatPhiU(1,4) = SpPhiUMod(1,1,2)
         AmatPhiU(1,5) = SpPhiUMod(1,2,2)
         AmatPhiU(1,6) = SpPhiUMod(1,3,2)
         AmatPhiU(1,7) = SpPhiUMod(1,1,3)
         AmatPhiU(1,8) = SpPhiUMod(1,2,3)
         AmatPhiU(1,9) = SpPhiUMod(1,3,3)
         AmatPhiU(2,1) = SpPhiUMod(2,1,1)
         AmatPhiU(2,2) = SpPhiUMod(2,2,1)
         AmatPhiU(2,3) = SpPhiUMod(2,3,1)
         AmatPhiU(2,4) = SpPhiUMod(2,1,2)
         AmatPhiU(2,5) = SpPhiUMod(2,2,2)
         AmatPhiU(2,6) = SpPhiUMod(2,3,2)
         AmatPhiU(2,7) = SpPhiUMod(2,1,3)
         AmatPhiU(2,8) = SpPhiUMod(2,2,3)
         AmatPhiU(2,9) = SpPhiUMod(2,3,3)
         AmatPhiU(3,1) = SpPhiUMod(3,1,1)
         AmatPhiU(3,2) = SpPhiUMod(3,2,1)
         AmatPhiU(3,3) = SpPhiUMod(3,3,1)
         AmatPhiU(3,4) = SpPhiUMod(3,1,2)
         AmatPhiU(3,5) = SpPhiUMod(3,2,2)
         AmatPhiU(3,6) = SpPhiUMod(3,3,2)
         AmatPhiU(3,7) = SpPhiUMod(3,1,3)
         AmatPhiU(3,8) = SpPhiUMod(3,2,3)
         AmatPhiU(3,9) = SpPhiUMod(3,3,3)
         !
         Kphiu = Kphiu - matmul(matmul(transpose(Gmatphi),AmatPhiU),
     +                          Gmat)*detMapJC*w(intpt)


         ! Compute/update the displacement/electric potential tangent matrix
         !
         AmatUPhi = zero
         AmatUPhi(1,1) = SpUPhiMod(1,1,1)
         AmatUPhi(2,1) = SpUPhiMod(2,1,1)
         AmatUPhi(3,1) = SpUPhiMod(3,1,1)
         AmatUPhi(4,1) = SpUPhiMod(1,2,1)
         AmatUPhi(5,1) = SpUPhiMod(2,2,1)
         AmatUPhi(6,1) = SpUPhiMod(3,2,1)
         AmatUPhi(7,1) = SpUPhiMod(1,3,1)
         AmatUPhi(8,1) = SpUPhiMod(2,3,1)
         AmatUPhi(9,1) = SpUPhiMod(3,3,1)
         AmatUPhi(1,2) = SpUPhiMod(1,1,2)
         AmatUPhi(2,2) = SpUPhiMod(2,1,2)
         AmatUPhi(3,2) = SpUPhiMod(3,1,2)
         AmatUPhi(4,2) = SpUPhiMod(1,2,2)
         AmatUPhi(5,2) = SpUPhiMod(2,2,2)
         AmatUPhi(6,2) = SpUPhiMod(3,2,2)
         AmatUPhi(7,2) = SpUPhiMod(1,3,2)
         AmatUPhi(8,2) = SpUPhiMod(2,3,2)
         AmatUPhi(9,2) = SpUPhiMod(3,3,2)
         AmatUPhi(1,3) = SpUPhiMod(1,1,3)
         AmatUPhi(2,3) = SpUPhiMod(2,1,3)
         AmatUPhi(3,3) = SpUPhiMod(3,1,3)
         AmatUPhi(4,3) = SpUPhiMod(1,2,3)
         AmatUPhi(5,3) = SpUPhiMod(2,2,3)
         AmatUPhi(6,3) = SpUPhiMod(3,2,3)
         AmatUPhi(7,3) = SpUPhiMod(1,3,3)
         AmatUPhi(8,3) = SpUPhiMod(2,3,3)
         AmatUPhi(9,3) = SpUPhiMod(3,3,3)
         !
         Kuphi = Kuphi - matmul(matmul(transpose(Gmat),AmatUPhi),
     +                          Gmatphi)*detMapJC*w(intpt)

      enddo
      !
      ! End the loop over integration points
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! Traction terms and surface charge terms have not been
      !  implemented in the three-dimensional code at this time.
      ! Mechanical, traction- and pressure-type boundary conditions 
      !  may be applied to the dummy mesh using the Abaqus built-in 
      !  commands *Dload or *Dsload.
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
         rhs(A11+2,1) = Ru(A12+2,1)
         !
         ! electric potential
         !
         rhs(A11+3,1) = Rphi(i,1)
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
            amatrx(A11,B11)     = Kuu(A12,B12)
            amatrx(A11,B11+1)   = Kuu(A12,B12+1)
            amatrx(A11,B11+2)   = Kuu(A12,B12+2)
            amatrx(A11+1,B11)   = Kuu(A12+1,B12)
            amatrx(A11+1,B11+1) = Kuu(A12+1,B12+1)
            amatrx(A11+1,B11+2) = Kuu(A12+1,B12+2)
            amatrx(A11+2,B11)   = Kuu(A12+2,B12)
            amatrx(A11+2,B11+1) = Kuu(A12+2,B12+1)
            amatrx(A11+2,B11+2) = Kuu(A12+2,B12+2)
            !
            ! electric potential
            !
            amatrx(A11+3,B11+3) = Kphiphi(i,j)
            !
            ! displacement/electric potential
            !
            amatrx(A11,B11+3) = Kuphi(A12,j)
            amatrx(A11+1,B11+3) = Kuphi(A12+1,j)
            amatrx(A11+2,B11+3) = Kuphi(A12+2,j)
            !
            ! electric potential/displacement
            !
            amatrx(A11+3,B11) = Kphiu(i,B12)
            amatrx(A11+3,B11+1) = Kphiu(i,B12+1)
            amatrx(A11+3,B11+2) = Kphiu(i,B12+2)
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

      subroutine xint3D1pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using a 1 gauss point for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !      
      integer nIntPt,nDim
      !
      real*8 xi(1,3),w(1)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 1


      ! Gauss weights
      !
      w(1) = 8.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = 0.d0
      xi(1,2) = 0.d0
      xi(1,3) = 0.d0


      return
      end subroutine xint3D1pt
      
!************************************************************************

      subroutine xint3D8pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 8 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(8,3),w(8)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 8


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      w(5) = 1.d0
      w(6) = 1.d0
      w(7) = 1.d0
      w(8) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(1,3) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(2,3) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(3,3) = -dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)
      xi(4,3) = -dsqrt(1.d0/3.d0)
      xi(5,1) = -dsqrt(1.d0/3.d0)
      xi(5,2) = -dsqrt(1.d0/3.d0)
      xi(5,3) = dsqrt(1.d0/3.d0)
      xi(6,1) = dsqrt(1.d0/3.d0)
      xi(6,2) = -dsqrt(1.d0/3.d0)
      xi(6,3) = dsqrt(1.d0/3.d0)
      xi(7,1) = -dsqrt(1.d0/3.d0)
      xi(7,2) = dsqrt(1.d0/3.d0)
      xi(7,3) = dsqrt(1.d0/3.d0)
      xi(8,1) = dsqrt(1.d0/3.d0)
      xi(8,2) = dsqrt(1.d0/3.d0)
      xi(8,3) = dsqrt(1.d0/3.d0)


      return
      end subroutine xint3D8pt

!************************************************************************

      subroutine calcShape3DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      !
      implicit none
      !
      integer intpt,nDim,nIntPt,i,j
      !
      real*8 xi_int(nIntPt,3),sh(8),dshxi(8,3),xi,eta,zeta
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      zeta = xi_int(intpt,3)
      
      
      ! The shape functions
      !
      sh(1) = eighth*(one - xi)*(one - eta)*(one - zeta)
      sh(2) = eighth*(one + xi)*(one - eta)*(one - zeta)
      sh(3) = eighth*(one + xi)*(one + eta)*(one - zeta)
      sh(4) = eighth*(one - xi)*(one + eta)*(one - zeta)
      sh(5) = eighth*(one - xi)*(one - eta)*(one + zeta)
      sh(6) = eighth*(one + xi)*(one - eta)*(one + zeta)
      sh(7) = eighth*(one + xi)*(one + eta)*(one + zeta)
      sh(8) = eighth*(one - xi)*(one + eta)*(one + zeta)
      
      
      ! The first derivatives
      !
      dshxi(1,1) = -eighth*(one - eta)*(one - zeta)
      dshxi(1,2) = -eighth*(one - xi)*(one - zeta)
      dshxi(1,3) = -eighth*(one - xi)*(one - eta)
      dshxi(2,1) = eighth*(one - eta)*(one - zeta)
      dshxi(2,2) = -eighth*(one + xi)*(one - zeta)
      dshxi(2,3) = -eighth*(one + xi)*(one - eta)
      dshxi(3,1) = eighth*(one + eta)*(one - zeta)
      dshxi(3,2) = eighth*(one + xi)*(one - zeta)
      dshxi(3,3) = -eighth*(one + xi)*(one + eta)
      dshxi(4,1) = -eighth*(one + eta)*(one - zeta)
      dshxi(4,2) = eighth*(one - xi)*(one - zeta)
      dshxi(4,3) = -eighth*(one - xi)*(one + eta)
      dshxi(5,1) = -eighth*(one - eta)*(one + zeta)
      dshxi(5,2) = -eighth*(one - xi)*(one + zeta)
      dshxi(5,3) = eighth*(one - xi)*(one - eta)
      dshxi(6,1) = eighth*(one - eta)*(one + zeta)
      dshxi(6,2) = -eighth*(one + xi)*(one + zeta)
      dshxi(6,3) = eighth*(one + xi)*(one - eta)
      dshxi(7,1) = eighth*(one + eta)*(one + zeta)
      dshxi(7,2) = eighth*(one + xi)*(one + zeta)
      dshxi(7,3) = eighth*(one + xi)*(one + eta)
      dshxi(8,1) = -eighth*(one + eta)*(one + zeta)
      dshxi(8,2) = eighth*(one - xi)*(one + zeta)
      dshxi(8,3) = eighth*(one - xi)*(one + eta)
      
      
      return
      end subroutine calcShape3DLinear

!*************************************************************************

      subroutine mapShape3D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,3),dsh(nNode,3),coords(3,nNode),mapJ(3,3),
     +  mapJ_inv(3,3),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)
      

      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,3
        do j=1,3
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv3D(mapJ,mapJ_inv,detMapJ,stat)


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))


      return
      end subroutine mapShape3D

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
        write(*,*) 'WARNING: SUBROUTINE matInv3:'
        write(*,*) 'WARNING: DET of MAT=',DET_A
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
