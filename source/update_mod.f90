! Copyright (C) 2007 Barbara Ercolano
!
! Version 3.00
module update_mod
    use constants_mod
    use common_mod
    use composition_mod
    use emission_mod
    use grid_mod
    use interpolation_mod
    use xSec_mod

    contains

      subroutine updateCell(grid, xP, yP, zP)
        implicit none

        real, parameter :: Y0 = 0.5, Y1 = 0.2          ! see Baldwin et al. 1991
        real                           :: comptonHeat, comptonCool, comptonRecoilHeat, &
             & comptonRecoilIonH

        real, dimension(nElements, nstages,4) :: &
             & chex                                    ! ch exchange coeff in cm^3/s

        real, dimension(2,nelements,nstages)       :: collIon      ! collisional ionisation
                                                       ! 1=contribution to ion balance
                                                       ! 2=contribution to thermal balance
        real, dimension(2,10,nelements, nstages)       :: photoIon     ! photoionisation
                                                       ! 1=contribution to ion balance
                                                       ! 2=contribution to thermal balance
        real                           :: deltaXHI        ! delta ionDen of H0
        real                           :: ionJ         !
        real                           :: phXSec      ! ph xSec of ion M at freqency bin j
        real                           :: thResidual   ! thermal balance residual
        real                           :: thResidualHigh! ther bal residual-high lim
        real                           :: thResidualLow! ther bal residual-low lim
        real                           :: Thigh        ! T high limit [K]
        real                           :: Tlow         ! T low limit  [K]
        real                           :: XOldHI       ! old ionDen of H0 at this cell
        real :: elemAbunUsed(nelements)

        ! dust-gas interaction heating and cooling process
        real                           :: grainEmi, grainRec ! grain emissions and recom
        real, pointer                  :: gasDustColl_d(:,:) ! cooling/heating of the
                                                             ! dust through collisions with grains
        real                           :: gasDustColl_g=0. ! cooling/heating of the
                                                       ! gas through collisions with grains
        real,pointer                   :: photoelHeat_d(:,:) ! cooling of dust by photoelectric
                                                             ! emission
        real                           :: photoelHeat_g=0. ! heating of gas by dust photoelctric
                                                       !emission
        real,pointer                   :: grainPot(:,:)  ! [Ryd]

        real, parameter                :: hcRyd_k = &  ! constant: h*cRyd/k (Ryd at inf used) [K]
             & 157893.94
        real, parameter                :: hcRyd = &    ! constant: h*c*Ryd (Ryd at inf used) [erg]
             & 2.1799153e-11
        real, parameter                :: thLimit = 0.02! convergence limit T-iteration



        real, dimension(nElements, nstages) &
             & :: alphaTot     ! total recombination coeffs

        integer, intent(in)            :: xP, yP, zP   ! cell indexes on the Cartesian axes




        integer ,pointer               :: grainPotP(:,:)
        integer                        :: cellP        ! points to this cell
        integer                        :: elem         ! element counter
        integer                        :: err          ! allocation error status
        integer                        :: ion          ! ionization stage counter
        integer                        :: i,j          ! counter
        integer                        :: nIterateGC   ! # of GC iterations
        integer                        :: nIterateT    ! # of T iterations
        integer                        :: nIterateX    ! # of X iterations
        integer                        :: nShell       ! shell counter
!        integer                        :: domCool      ! dominant collants : CEL = 1; RL = 2

        integer, parameter             ::  nTbins=300  ! number of enthalpy bins for T spike
        integer, parameter             :: maxIterateGC&! limit to number of grain charge it
             & = 100
        integer, parameter             :: maxIterateX& ! limit to number of X-iterat ions
             & = 15
        integer, parameter             :: maxIterateT& ! limit to number of T-iterations
             & = 20

        type(grid_type), intent(inout) :: grid         ! the grid

        logical                        :: lgHit        ! has this cell been hit by a photon?
        logical                        :: lgVerbose=.false.

        ! check whether this cell is outside the nebula
        if (grid%active(xP, yP, zP)<=0) return

        cellP = grid%active(xP, yP, zP)

        ! initialise lgBlack
        grid%lgBlack(cellP) = 0

        ! initialise lgHit
        lgHit = .false.



        ! find out if this cell has been hit by at least one photon
        if (.not.lgDebug) then
           do i = 1, nbins
              if ( grid%Jste(cellP,i) > 0.) then
                 lgHit = .true.
                 exit
              end if
           end do
        else
           do i = 1, nbins
              if ( (grid%Jste(cellP,i) > 0.) .or. &
                   & (grid%Jdif(cellP,i) > 0.)) then
                 lgHit = .true.
                 exit
              end if
           end do
        end if


        ionJ = 0.
        do i = lymanP, nbins
           ionJ = ionJ+grid%Jste(cellP,i)
        end do


        ! do not update this cell if there were no hits
        if (.not.lgHit) then

           grid%lgBlack(cellP) = 1

           if (lgTalk) print*, "! updateCell [talk]: no photon hits, returning...", xP,yP,zP

           if (lgDust) TdustTemp(:,:,cellP)       = grid%Tdust(:,:,cellP)

           if (lgGas) then
              ! the grid values stay the same
              TeTemp(cellP)         = grid%Te(cellP)
              NeTemp(cellP)         = grid%Ne(cellP)
              ionDenTemp(cellP,:,:) = grid%ionDen(cellP,:,:)
           end if

           grid%noHit = grid%noHit+1.

           return
        end if

        if (lgGas) then

           ! initialize local Te and Ne

           TeUsed     = grid%Te(cellP)
           NeUsed     = grid%Ne(cellP)
           HdenUsed   = grid%Hden(cellP)
           ionDenUsed = grid%ionDen(cellP, :, :)
           elemAbunUsed = grid%elemAbun(grid%abFileIndex(xp,yp,zp),:)

           ! save present value of H0 abundance in XOldHI
           XOldHI = grid%ionDen(cellP,elementXref(1),1)

           ! zero out T-iteration components
           nIterateT      = 0
           thResidual     = 0.
           thResidualHigh = 0.
           thResidualLow  = 0.
           Thigh          = 0.
           Tlow           = 0.

           err=0

           if (lgDust .and. lgGas .and. lgPhotoelectric) then
              allocate (grainPot(1:nSPecies, 1:nsizes))
              if (err /= 0) then
                 print*, "! updateCell:cannot allocate grid memory,grainpot"
                 stop
              end if
              grainPot=0.
              allocate (grainPotP(1:nSPecies, 1:nsizes))
              if (err /= 0) then
                 print*, "! updateCell:cannot allocate grid memory,grainpotp"
                 stop
              end if
              grainPotP=0
              allocate (photoelHeat_d(1:nSPecies, 1:nsizes))
              if (err /= 0) then
                 print*, "! updateCell:cannot allocate grid memory,photoelHeat_d"
                 stop
              end if
              photoelHeat_d=0.
              allocate ( gasDustColl_d(1:nSPecies, 1:nsizes))
              if (err /= 0) then
                 print*, "! updateCell:cannot allocate grid memory"
                 stop
              end if
              gasDustColl_d=0.
           end if


         ! call the recursive iteration procedure for the gas
         call iterateT()

         ! call the non recursive itaration procedure for the dust
         if (lgDust) call getDustT()

         ! this was added to help implementing MPI comunication
         if (lgDust) TdustTemp(:,:,cellP)          = grid%Tdust(:,:,cellP)

         if (lgNeInput .and. lgTalk) then
            print*, '! updateCell: [talk] cell:', xP,yP,zP, " NeInput: ", &
                 &grid%NeInput(cellP), " NeUsed: ", &
                 &grid%Ne(cellP),  " N_gas: ", &
                 &grid%Hden(cellP)
         end if


      else

         if (.not.lgDust) then
            print*, '! updateCell: no gas or dust present. the grid is empty.'
            stop
         end if

         XOldHI = grid%Tdust(0,0,cellP)

         if (lgRadPress) call getRadAcc()


         call getDustT()

         ! determine if the model has converged at this cell
         ! NOT  variable names refer to gas phase for reasons of laziness
         deltaXHI = (grid%Tdust(0,0,cellP) - XOldHI) / XOldHI
         if ( abs(deltaXHI) <= XHILimit ) then
            grid%lgConverged(cellP) = 1
         else
            grid%lgConverged(cellP) = 0
         end if

         if (lgTalk) &
              & print*, "updateCell: [talk] cell", xP,yP,zP, "; converged?",&
              & grid%lgConverged(cellP), "; mean dust T: ", &
              & grid%Tdust(0,0,cellP), "; &
              & mean dust T old: ", XOldHI,"; dT(dust): ", deltaXHI

         ! this was added to help implementing MPI comunication
         TdustTemp(:,:,cellP)          = grid%Tdust(:,:,cellP)

      end if

    contains


        recursive subroutine iterateT()
            implicit none

            real                          :: coolInt         ! tot cooling integral [erg/s/Hden]
            real                          :: heatInt         ! tot heating integral [erg/s/Hden]

            real, parameter               :: Xmax = 0.9999   ! max rel abundance

            integer                       :: isp, ai         ! counters

            logical                       :: lgGCBConv       ! grain charge converged?
            logical                       :: lgIBConv        ! converged?

            ! step up T-iteration
            nIterateT = nIterateT + 1

            ! calculate the recombination coefficients for this cell
            call calcalpha()

            ! photoionisation banter
            ! calculate photo rates
            call photoionisation()

            ! calculate collisional rates -NOTE: NEED TO MULTIPLY COLLION BY Ne!!!
            call collisionalIonisation()

            ! initialize X iteration counter
            nIterateX = 1

            if (lgCompton) call compton()

            ! calculate the ion abundances for the remaning ions
            call ionBalance(lgIBConv)

            ! calculate dust-gas interaction heating/cooling terms
            if (lgDust .and. lgGas .and. lgPhotoelectric) then
               do isp = 1, nSpecies
                  do ai = 1, nsizes
                     nIterateGC     = 0
                     grainEmi       = 0.
                     grainRec       = 0.

                     call setGrainPotential(isp,ai,lgGCBConv)

                     call locate(nuArray, grainPot(isp,ai), grainPotP(isp,ai))
                     if (grainPotP(isp,ai) == 0) grainPotP = 1

                  end do
               end do
               call setPhotoelHeatCool()


               call setDustGasCollHeatCool()

            end if

            ! solve thermal balance
            call thermBalance(heatInt, coolInt)

            ! Now calculate new Te

            ! calculate thResidual = heatInt - coolInt
            thResidual = heatInt - coolInt

            if ( abs(thResidual) >= thLimit*(abs(coolint)+abs(heatInt))/2.) then ! start convergence condition

                if ( nIterateT < maxIterateT ) then ! start nIterateT condition

                    if ( thResidual < 0) then
                       ! too cool

                        thResidualLow = thResidual
                        Tlow          = TeUsed

                        if ( Thigh /= 0. ) then

                            TeUsed = Tlow-(Thigh-Tlow)*thResidualLow/&
                                 & (thResidualHigh-thResidualLow)

                        else

                            TeUsed = TeUsed/1.2

                        end if

                    else

                        thResidualHigh = thResidual
                        Thigh          = TeUsed

                        if ( Tlow /= 0. ) then

                            TeUsed = Tlow-(Thigh-Tlow)*thResidualLow/&
                                 & (thResidualHigh-thResidualLow)

                        else

                            TeUsed = TeUsed*1.2

                        end if

                    end if


                    if (nIterateT >= 6) then
                        if (abs(Thigh-Tlow) <= (0.002*(Thigh+Tlow)) ) then

                            if ( thResidual < 0 ) then

                                TeUsed = Thigh - 100.
                                 Thigh  = 0.

                            else

                                TeUsed = Tlow+100.
                                Tlow   = 0.

                            end if
                        end if
                    end if

                    if (TeUsed <= 0.) TeUsed = 1.

                    ! next T-iteration
                    call iterateT()
                    return

                else ! if nIterateT > maxIterateT

                    ! after maxIterateT number of iterations the temperature
                    ! is determined as follows

                    if ( Thigh == 0. ) then

                        TeUsed = Tlow

                    else

                        if (Tlow == 0.) then

                            TeUsed = Thigh

                        else

                            TeUsed = (Thigh+Tlow)*0.5

                        end if

                    end if


                    if (lgVerbose) print*, "! iterateT: [warning] no convergence after ", &
                         & nIterateT, " steps. (cell, T)", cellP,xP,yP,zP,TeUsed

                    grid%noTeBal = grid%noTeBal+1.

!CHANGED
!                    grid%lgBlack(cellP) = 1

                end if

            else ! if thResidual < thLimit

                ! the T-iteration has converged

                if (lgTalk) print*, "! iterateT: [talk] convergence achieved after ",&
                     & nIterateT, " steps. (cell, T)", xP,yP,zP,TeUsed, grid%abFileIndex(xP,yP,zP)

            end if ! end convergence condition

            if (.not.lgIBConv) then
               grid%noIonBal = grid%noIonBal+1

!CHANGED
!               grid%lgBlack(cellP) = 1
            end if

            grid%Te(cellP)         = TeUsed

            ! determine if the model has converged at this cell
            ! converge on temperature
            deltaXHI = (grid%ionDen(cellP,elementXref(1),1) - XOldHI) / XOldHI
!            deltaXHI = (grid%Te(cellP) - XOldHI) / XOldHI
            if ( abs(deltaXHI) <= XHILimit ) then
               grid%lgConverged(cellP) = 1
            else
               grid%lgConverged(cellP) = 0
            end if

 !           if (lgTalk) print*, "iterateT: [talk] cell", xP,yP,zP, "; converged?",&
 !                & grid%lgConverged(cellP), "; Te: ", &
 !                & grid%Te(cellP), &
 !                & "; Te old: ", XOldHI,"; dX(H0): ", deltaXHI

            ! this was added to help implementing MPI comunication
            TeTemp(cellP)          = TeUsed

        end subroutine iterateT

        recursive subroutine setGrainPotential(iSp, ai, lgGCBConv)
          implicit none

          real,save            :: delta,delta1
          real,save            :: grainPotOld ! local copy of grai pot
          real                 :: grainEmi, grainRec ! grain emissions and recom
          real,save            :: grainEmiOld, grainRecOld ! grain emissions and recom
          real,parameter       :: errorLim = 0.005, dm = 0.05 ! loop convergence
          real,parameter       :: safeLim = 100 ! loop safety limit
          real                 :: threshold,fac
          real,save :: dVg,slope

          integer, intent(in)  :: iSp, ai

          logical, intent(inout) :: lgGCBConv   ! converged?

          nIterateGC = nIterateGC+1

          if (nIterateGC==1) then
             dVg=0.05
             grainPot(isp,ai) = 0.
             grainPotOld = grainPot(isp,ai)
             lgGCBConv = .true.
             grainEmiOld = getGrainEmission(grainPot(isp,ai), isp, ai)
             grainRecOld = getGrainRecombination(grainPot(isp,ai))
             grainPot(isp,ai) = grainPotOld+0.05
          end if

          threshold = max(grainVn(isp)+grainPot(isp,ai),grainVn(isp))
          grainEmi = getGrainEmission(grainPot(isp,ai), isp, ai)
          grainRec = getGrainRecombination(grainPot(isp,ai))
          delta = grainEmi-grainRec

          delta1 = abs(delta/(0.5*max(1.e-35, grainEmi+grainRec)))


          ! check for convergence
          if (delta1<errorLim) then
             return
          else

             if (grainPot(iSp,ai) /= grainPotOld) then
                fac = (grainEmi-grainEmiOld)-(grainRec-grainRecOld)
                if (fac/=0.) slope = fac/(grainPot(iSp,ai)-grainPotOld)
             end if

             grainPotOld=grainPot(iSp,ai)
             grainRecOld = grainRec
             grainEmiOld = grainEmi

             delta1 = -delta/slope
             delta = abs(delta1)
             delta = min(delta, dm*threshold)
             delta = sign(delta, delta1)

             grainPot(iSp,ai) = grainPot(iSp,ai)+delta

             if (nIterateGC< maxIterateGC) then
                call setGrainPotential(isp,ai,lgGCBConv)
                return
             else
                print*, '! setGrainPotential: no convergence', cellP,grainPot(isp,ai),grainEmi,grainRec
                lgGCBConv=.false.
                return
             end if

          end if

        end subroutine setGrainPotential

        ! calculate the grain recombination
        ! using eqn 18 etc of Baldwin et al. (1991)
        ! cellFactor is dependant on the physical conditions of the gas,
        function getGrainEmission(Vg, isp,ai)
          implicit none


          real :: getGrainEmission
          real, intent(in) :: Vg ! grain potential
          real :: Yn, Yhat

          real :: Qa ! grain absorption efficiency
          real :: thres, photFlux

          integer, intent(in) :: isp, ai
          integer :: ifreq, ip

          getGrainEmission=0.


          ! get the threshold
          thres = max(grainVn(isp)+Vg,grainVn(isp))

          call locate( nuArray, thres, ip)
          ip = ip+1


          do ifreq = ip, nbins


             Yn = min(Y0*(1.-grainVn(isp)/nuArray(ifreq)), Y1)

             Yhat = Yn*min(1., max(0.,1.-Vg/(nuArray(ifreq)-grainVn(isp))))

             if (.not. lgDebug) then
!                photFlux = (grid%JPEots(cellP,ifreq) + grid%Jste(cellP,ifreq))/(hcRyd*nuArray(ifreq))
                photFlux =  grid%Jste(cellP,ifreq)/(hcRyd*nuArray(ifreq))
             else
!                photFlux = (grid%JPEots(cellP,ifreq) + grid%Jste(cellP,ifreq)+grid%Jdif(cellP,ifreq))/&
!                     & (hcRyd*nuArray(ifreq))
             end if

             photFlux = photFlux*fourPi

             Qa = XSecArray(dustAbsXsecP(isp,ai)+ifreq-1)/(1.e-8*grainRadius(ai)**2.)

             getGrainEmission = getGrainEmission+Yhat*photFlux*Qa

          end do

        end function getGrainEmission


        ! calculate the grain recombination
        ! using eqn 23 etcof Baldwin et al. (1991)
        ! cellFactor is dependant on the physical conditions of the gas,
        function getGrainRecombination(Vg)
          implicit none

          real :: getGrainRecombination
          real, intent(in) :: Vg ! grain potential [ryd]

          real :: eta   ! Coulomb correction
          real :: cpDen ! colliding particle number density [cm^-3]
          real :: eightkT_pi ! 8*k * Te/Pi [erg]
          real :: mcp   ! mass of colliding particle in [g]
          real :: kT    ! k*Te [ryd]
          real :: S     ! sticking coefficient
          real :: vmean ! colliding particle mean velocity
          real :: Z     ! colliding particle charge

          integer :: istage

          getGrainRecombination = 0.
          eta = 0.
          kT =6.336e-6*TeUsed
          eightkT_pi = (1.1045e-15)*TeUsed/Pi

          ! add e- collisions contributions

          vmean = sqrt(eightkT_pi/me)

          S  =1. ! electron sticking probability


          ! eq 24 of Baldwin et al 91
          eta = -Vg/kT

          if (eta <= 0.) then
             eta = 1.-eta
          else if (eta >0.) then
!print*, 'eta ', eta
             eta = exp(-eta)
          else
             print*, "! getGrainRecombination: insane eta for e-", eta
          end if

          getGrainRecombination = getGrainRecombination + &
               & NeUsed*vmean*S*eta

          ! add contribution from all other neutral and ionic species
          do elem = 1, nElements
             if (lgElementOn(elem)) then
                do istage = 2, min(elem+1,nstages)
                   ! get cpDen
                   cpDen = grid%ionDen(cellP,elementXref(elem),istage)*&
                        & grid%elemAbun(grid%abFileIndex(xP,yP,zP),elem)*&
                        & grid%Hden(cellP)
                   mcp = (aWeight(elem)*amu)
                   vmean = sqrt(eightkT_pi/mcp)

                   S = 1.
                   Z = real(istage-1)

                   eta = Z*Vg/kT

                   if (eta <= 0.) then
                      eta = 1.-eta
                   else if (eta >0.) then
                      eta = exp(-eta)
                   else
                      print*, "! getGrainRecombination: insane eta", &
                           & eta, elem, istage
                   end if

                   getGrainRecombination = getGrainRecombination - &
                        & cpDen*vmean*S*eta

                end do
             end if
          end do

        end function getGrainRecombination

        ! see Baldwin et al 1991; but beware that we are resolving the size distribution.
        ! so must keep size dependance, so Qa ia really Pi a^2 Qa
        subroutine setPhotoelHeatCool()
          implicit none

          real    :: Qa, photFlux, EY, Yhat,Yn,th

          integer :: ns,na,ifreq,thP

          ! calculate the cooling of dust by photoelectric emission
          !  and heating of gas
          ! Baldwin et al. 1991 eqn 25-27

          photoelHeat_d=0.
          photoelHeat_g=0.
          do ns = 1, nSpecies
             do na = 1, nSizes

                if (grainPotP(ns,na) <= 0) then
                   print*, "! setPhotoelHeatCool irregular grain potential index", &
                        grainPotP(ns,na), grainPot(ns,na)
                   stop
                end if

                th = max(grainVn(ns)+grainPot(ns,na), grainVn(ns))
                call locate(nuArray,th,thP)
                if(thP<=0) thP=1

                do ifreq = thP, nbins

!                do ifreq = grainPotP(ns,na), nbins

                   Yn = min(Y0*(1.-grainVn(ns)/nuArray(ifreq)), Y1)

                   Yhat = Yn*min(1., max(0.,1.-grainPot(ns,na)/(nuArray(ifreq)-grainVn(ns))))

                   if (.not. lgDebug) then
!                      photFlux = (grid%JPEots(cellP,ifreq) + &
!                           & grid%Jste(cellP,ifreq))/(hcRyd*nuArray(ifreq))
                      photFlux = grid%Jste(cellP,ifreq)/(hcRyd*nuArray(ifreq))
                   else
!                      photFlux = (grid%JPEots(cellP,ifreq) + grid%Jste(cellP,ifreq)+&
!                           & grid%Jdif(cellP,ifreq))/(hcRyd*nuArray(ifreq))
                      photFlux = (grid%Jste(cellP,ifreq)+grid%Jdif(cellP,ifreq))/(hcRyd*nuArray(ifreq))

                   end if

                   photFlux = photFLux

                   Qa = XSecArray(dustAbsXsecP(ns,na)+ifreq-1)


                   EY = Yn*0.5* min(nuArray(ifreq)-grainVn(ns),&
                        & max(0., ((nuArray(ifreq)-grainVn(ns))**2.-grainPot(ns,na)**2.)/&
                        & (nuArray(ifreq)-grainVn(ns))))

!                   photoelHeat_d(ns,na) = photoelHeat_d(ns,na)+Qa*photFlux*EY*hcRyd/Pi
                   photoelHeat_d(ns,na) = photoelHeat_d(ns,na)+Qa*photFlux*EY*hcRyd/Pi

                   photoelHeat_g = photoelHeat_g+Qa*photFlux*(EY-Yhat*grainPot(ns,na))*&
                        & grainAbun(ns)*grainWeight(na)

                end do
             end do
          end do

          photoelHeat_g = photoelHeat_g*grid%Ndust(cellP)*hcRyd/grid%Hden(cellP)

        end subroutine setPhotoelHeatCool

        ! sets the cooling and heating rates of gas and dust due to collisions between the two phases
        ! see Baldwin et al 1991
        ! only collisions with up to the 3 times ionised case are considered here
        ! (process becomes unimportant for higher ionisation cases)
        subroutine setDustGasCollHeatCool()
          implicit none

          real :: Z , kT, eta, psi,xi,  S, IPerg
          real :: eightkT_pi ! 8*k * Te/Pi [erg]
          integer :: ss(30)
          real :: vmean, mcp ! colliding particle mean velocity and mean particle mass
          integer :: nelectrons, istage, ns, na

          ! outer shell array
          ss = (/1,1,2,2,3,3,3,3,3,3,4,4,5,5,5,5,5,5,&
               &6,6,6,6,6,6,6,6,6,6,7,7/)

          gasDustColl_g = 0.
          gasDustColl_d = 0.

          kT =6.336e-6*TeUsed ! ryd
          eightkT_pi = (1.1045e-15)*TeUsed/Pi


          do elem = 1, nElements ! 0 for electrons
             if ( lgElementOn(elem)) then

                   do istage = 1, min(nstages,elem+1)
                      Z = real(istage-1)

                      mcp = aWeight(elem)*amu
                      vmean = sqrt(eightkT_pi/mcp)

                      ! number of e-'s of the ionisation stage above
                      nelectrons = elem-istage+1

                      if (istage>1) then
                         if (elem == 1) then
                            IPerg = nuArray(HlevNuP(1))*ryd2erg
                         else if (elem == 2 .and. istage == 2) then
                            IPerg = nuArray(HeIlevNuP(1))*ryd2erg
                         else if (elem == 2 .and. istage == 3) then
                            IPerg = nuArray(HeIIlevNuP(1))*ryd2erg
                         else
                            IPerg = nuArray(elementP(elem,istage-1,nShells(elem,istage-1),1))*ryd2erg
                         end if
                      end if

                      do ns = 1, nSpecies
                         if (istage == 1) then
                            S = 2*mcp*MsurfAtom(ns)/(mcp+MsurfAtom(ns))**2.
                         else
                            S = 1.
                         end if
                         do na = 1, nSizes
                            psi = Z*grainPot(ns,na)/kT
                            if (psi <= 0. ) then
                               eta = 1.-psi
                               xi = 1. - psi/2.
                            else
!print*, 'psi ', psi
                               eta = exp(-psi)
                               xi = (1.+psi/2.) * eta
                            end if

                            gasDustColl_d(ns,na) = gasDustColl_d(ns,na)+&
                                 & ionDenUsed(elementXref(elem),istage)*&
                                 & grid%elemAbun(grid%abFileIndex(xP,yP,zP),elem)*&
                                 & grid%Hden(cellP)*&
                                 & Pi*grainRadius(na)*grainRadius(na)*1.e-8*&
                                 & S*vmean*(2.*kT*Ryd2erg*xi-eta*&
                                 & (Z*grainPot(ns,na)*Ryd2erg-IPerg+&
                                 & 2.*kBoltzmann*grid%Tdust(ns,na,cellP)))

                            gasDustColl_g = gasDustColl_g + &
                                 & ionDenUsed(elementXref(elem),istage)*&
                                 & grainWeight(na)*grainAbun(ns)*grid%Ndust(cellP)*&
                                 & grid%elemAbun(grid%abFileIndex(xP,yP,zP),elem)*&
                                 & Pi*grainRadius(na)*grainRadius(na)*1.e-8*&
                                 & S*vmean*(2.*kT*Ryd2erg*xi-eta*2.*kBoltzmann*&
                                 & grid%Tdust(ns,na,cellP))
                         end do
                      end do
                   end do
                end if
             end do

             ! add contribution of e- collisions
             vmean = sqrt(eightkT_pi/me)
             do ns = 1, nSpecies
                S = 1.
                Z=-1.
                do na = 1, nSizes

                   psi = Z*grainPot(ns,na)/kT
                   if (psi <= 0. ) then
                      eta = 1.-psi
                      xi = 1. - psi/2.
                   else
!print*, 'psi2 ', psi
                      eta = exp(-psi)
                      xi = (1.+psi/2.) * eta
                   end if

                   gasDustColl_d(ns,na) = gasDustColl_d(ns,na)+&
                        & NeUsed* &
                        & Pi*grainRadius(na)*grainRadius(na)*1.e-8*&
                        & S*vmean*(2.*kT*Ryd2erg*xi-eta*&
                        & (Z*grainPot(ns,na)*ryd2erg))

                   gasDustColl_g = gasDustColl_g + &
                        & NeUsed* &
                        & grainWeight(na)*grainAbun(ns)*grid%Ndust(cellP)*&
                        & Pi*grainRadius(na)*grainRadius(na)*1.e-8*&
                        & S*vmean*(2.*kT*Ryd2erg*xi)/&
                        & grid%Hden(cellP)
                end do
             end do

             ! factor of fourpi to make up for the lack of fourpi
             ! in the balance eqns for J
             gasDustColl_d(:,:)= gasDustColl_d(:,:)/fourpi
             gasDustColl_g= gasDustColl_g/fourpi

           end subroutine setDustGasCollHeatCool


        subroutine thermBalance(heatInt, coolInt)
            implicit none

            real, intent(out)      :: heatInt, &    ! total heating and
                 & coolInt       ! cooling integrals
!            integer, intent(out)   :: dc ! dominant coolant ; CEL=1;RL=2
            ! local variables

            integer                :: i,j,k         ! counters
            integer                :: elem, ion     ! counters

            real                   :: betaFF        ! energy loss coeff due to ff rad
            real                   :: betaRec       ! energy loss coeff due to recomb
            real                   :: ch12, ch13, & ! collision excitation of
                 & ex12, ex13,&                     ! hydrogen data
                 & th12, th13                       ! (Mathis, Ly alpha, beta)
            real                   :: coolFF        ! cool due to FF radiation [erg/s/Hden]
            real                   :: coolColl      ! cool due to coll excit   [erg/s/Hden]
            real                   :: coolCollH     ! cool due to coll excit   [erg/s/Hden]
            real                   :: coolRec       ! cool due to recombination[erg/s/Hden]
            real                   :: fcool
            real                   :: heatSte       ! tot heat gain due to stellar phot
            real                   :: heatDif       ! tot heat gain due to diffuse phot
            real                   :: log10Te       ! log10(Te)
            real                   :: Np            ! proton density
            real                   :: Te4           ! Te/10000.

            coolFF  = 0.
            coolRec = 0.

            do i = 1, nElements
               if (lgElementOn(i) .and. nstages > i) then
!                  if (grid%ionDen(cellP,elementXref(i),i+1) > 1.e-5) then

                     log10Te = log10(TeUsed/real(i*i))
                     Te4     = TeUsed / (i*i*1.e4)

                     ! find the N(Xi+)
                     Np = grid%ionDen(cellP,elementXref(i),i+1)*grid%elemAbun(grid%abFileIndex(xP,yP,zP),i)

                     ! cooling of gas due to FF radiation from H-like
                     ! fits to Hummer, MNRAS 268(1994) 109, Table 1. or  least square fitting to m=4

                     betaFF = real(i)*(1.0108464E-11 + 9.7930778E-13*log10Te - &
                          & 6.6433144E-13*log10Te*log10Te + 2.4793747E-13*log10Te*log10Te*log10Te -&
                          & 2.3938215E-14*log10Te*log10Te*log10Te*log10Te)

                     if (.not. isnan(Np*NeUsed*betaFF*kBoltzmann*TeUsed/sqrt(TeUsed)))&
                          & coolFF = coolFF+Np*NeUsed*betaFF*kBoltzmann*TeUsed/sqrt(TeUsed)

                     ! cooling of gas due to recombination of H-like
                     ! fits to Hummer, MNRAS 268(1994) 109, Table 1.
                     ! least square fitting to m=4
                     betaRec = real(i)*(9.4255985E-11 -4.04794384E-12*log10Te &
                          & -1.0055237E-11*log10Te*log10Te +  1.99266862E-12*log10Te*log10Te*log10Te&
                          & -1.06681387E-13*log10Te*log10Te*log10Te*log10Te)

                     if (.not. isnan(Np*NeUsed*betaRec*kBoltzmann*TeUsed/sqrt(TeUsed))) &
                          & coolRec = coolRec+Np*NeUsed*betaRec*kBoltzmann*TeUsed/sqrt(TeUsed)

!                  end if
               end if
            end do

            if (lgTraceHeating.and.taskid==0) then
               if (nIterateT > 1) then
                  do i = 1, 14
!                     backspace 57
                  end do
               end if
               write(57,*) 'Cell,Te: ', xp,yp,zp, TeUsed, NeUsed
               write(57,*) 'Ionising radiation: ', ionJ
               write(57,*) 'FF H-like: ', coolFF
               write(57,*) 'Rec H-like: ', coolrec

            end if

            log10Te = log10(TeUsed)
            Te4     = TeUsed / (1.e4)

            do i = 2, nElements
               if (lgElementOn(i) .and. nstages > i) then
!                  if (grid%ionDen(cellP,elementXref(i),i) > 1.e-5) then

                     log10Te = log10(4.*TeUsed/real(i*i))
                     Te4     = 4.*TeUsed / (i*i*1.e4)

                     ! cooling of gas due to FF radiation from He-like
                     ! fits to Hummer and Storey, MNRAS 297(1998) 1073, Table 6. least square fitting to m=4

                     ! find N(X+)
                     Np = grid%ionDen(cellP,elementXref(i),i)*grid%elemAbun(grid%abFileIndex(xP,yP,zP),i)

!                     Np =  grid%ionDen(cellP,elementXref(2),2)*grid%elemAbun(grid%abFileIndex(xp, yP, zP),2)

                     betaFF =  (real(i)/2.)*( 1.070073e-11    -2.5730207e-13*log10Te + &
                          & 2.109134e-13*log10Te*log10Te )
!                     betaFF = 1.070073e-11    -2.5730207e-13*log10Te + &
!                          & 2.109134e-13*log10Te*log10Te


                     if (.not. isnan(Np*NeUsed*betaFF*kBoltzmann*TeUsed/sqrt(TeUsed))) &
                          & coolFF = coolFF + Np*NeUsed*betaFF*kBoltzmann*TeUsed/sqrt(TeUsed)

                     ! cooling of gas due to recombination of He+
                     ! fits to Hummer and Storey, MNRAS 297(1998) 1073, Table 6. least square fitting to m=4
                     betaRec =    (real(i)/2.)*(  1.06926e-10 -2.41756e-11*log10Te+ &
                          &1.14542e-12 *log10Te*log10Te  &
                          &-5.06535e-13 *log10Te*log10Te*log10Te +&
                          & 9.16745e-14 *log10Te*log10Te*log10Te*log10Te)

                     if (.not. isnan(Np*NeUsed*betaRec*kBoltzmann*TeUsed/sqrt(TeUsed))) &
                          & coolRec = coolRec + Np*NeUsed*betaRec*kBoltzmann*TeUsed/sqrt(TeUsed)

!                  end if
               end if
            end do

            if (lgTraceHeating.and.taskid==0) then
               write(57,*) 'FF H-like + He+: ',coolFF
               write(57,*) 'Rec H-like + He+: ',coolRec
            end if


            ! collisional excitation of Hydrogen
            ! Mathis, Ly alpha, beta
            ch12 = 2.47e-8
            ch13 = 1.32e-8
            ex12 = -0.228
            ex13 = -0.460
            th12 = 118338.
            th13 = 140252.

            if (TeUsed > 5000.) then

!print*, 'th12/TeUsed ', th12/TeUsed, th13/TeUsed
                coolColl = (ch12*exp(-th12/TeUsed)*Te4**ex12 + &
                     & ch13*exp(-th13/TeUsed)*Te4**ex13) * &
                     & hcRyd*grid%ionDen(cellP,elementXref(1),1)*NeUsed
!print*, 'th12/TeUsed '
             else

                coolColl = 0.

             end if

            if (lgTraceHeating.and.taskid==0) then
               write(57,*) 'Coll exc H: ',coolColl, NeUsed
               fcool = 0.
            end if


             ! cooling due to collisional ionisation of heavy metals
            coolcollh=0.
             do elem =1, nElements
                if (lgElementOn(elem)) then
                   do ion = 1, min(elem,nstages-1)
                      coolCollH = coolCollH + grid%elemAbun(grid%abFileIndex(xP,yP,zP),elem)*&
                           &grid%ionDen(cellP,elementXref(elem),ion) *collIon(2,elem,ion)*NeUsed
                   end do
                end if
             end do

            if (lgTraceHeating.and.taskid==0) then
               write(57,*) 'Coll ionisation of heavies: ',coolCollH
               fcool = 0.
            end if


             ! collisional excitation of Heavies

             ! get the emissivities of the forb lines
             call forLines()
!print*,'a'
!             coolCEL = 0.

             ! sum all contributions from the heavies to coolColl
             do elem = 3, nElements
                do ion = 1, min(elem+1, nstages)
                   if (.not.lgElementOn(elem)) exit
                   if (lgDataAvailable(elem, ion)) then

                      do j = 1, nForLevels
                         do k = 1, nForLevels

                            if ( forbiddenLines(elem,ion,j,k) > 1.e-35) then
                               coolColl = coolColl + forbiddenLines(elem,ion,j,k)
!                               coolCEL = coolCEL + forbiddenLines(elem,ion,j,k)
                            else
                               forbiddenLines(elem,ion,j,k) = 0.
                            end if

                            if (lgTraceHeating.and.taskid==0) then
                               fcool = fcool + forbiddenLines(elem,ion,j,k)
                            end if
                         end do
                      end do

                   end if
                end do
             end do

            if (lgTraceHeating.and.taskid==0) then
               write(57,*) 'CELs cool: ',fcool
               fcool = 0.
            end if


            ! add cooling by rec lines

!            call RecLinesEmission()

!            hydrolines = hydrolines*1.e-25
!            HeIRecLines = HeIRecLines*1.e-25

            ! hydrogenic rec lines
!            do izp = 1, 30
!               do iup = 3, 15
!                  do ilow =2, min(8, iup-1)

!                     coolColl = coolColl + hydroLines(izp, iup, ilow)
!                     if (lgTraceHeating.and.taskid==0) fcool = fcool+hydroLines(izp, iup, ilow)
!                  end do
!               end do
!            end do
!            do i = 1, 34
!               coolColl = coolColl + HeIRecLines(i)
!               if (lgTraceHeating.and.taskid==0) fcool = fcool+HeIRecLines(i)
!            end do

!            if (lgTraceHeating.and.taskid==0) then
!               write(57,*) 'rec lines: ',fcool
!               fcool = 0.
!            end if

!            if (coolCEL < coolRL) dc = 2

            if (lgCompton) then

               ! compton heating must be multipied by Ne
               comptonHeat = comptonHeat*NeUsed
               comptonCool = comptonCool*NeUsed

               ! cooling by Compton
               if (lgTraceHeating.and.taskid==0) then
                  write(57,*) 'Compton cool: ',comptonCool
               end if

               ! cooling by Compton
               if (lgTraceHeating.and.taskid==0) then
                  write(57,*) 'Compton heat: ', comptonHeat
                  write(57,*) 'Compton recoil heat: ', comptonRecoilHeat
               end if
            end if

             ! heating due to photoionization

             ! re-initialize heatSte and heatDif
             heatSte = 0.
             if (lgDebug) heatDif = 0.
!print*, 'heatste -photoionisation-'
             do elem = 1, nElements  ! begin element loop

                 do ion = 1, min(elem, nStages-1) ! begin ion loop
 !                   heationste = 0.
                    if(.not.lgElementOn(elem)) exit
                    do nShell = 1, nShells(elem, ion)
!print*, elem, ion, nshell
!print*, photoIon(2,nshell,elem,ion), ionDenUsed(elementXref(elem),ion), grid%elemAbun(grid%abFileIndex(xP,yP,zP),elem)
                       heatSte = heatSte+photoIon(2,nshell,elem,ion)*ionDenUsed(elementXref(elem),ion)*&
                                   & grid%elemAbun(grid%abFileIndex(xP,yP,zP),elem)


!print*, elem, ion, photoIon(2,nshell,elem,ion), ionDenUsed(elementXref(elem),ion), grid%elemAbun(grid%abFileIndex(xP,yP,zP),elem)
!                       heatSte    = heatSte + heatIonSte
                    end do ! end shell loop

!                    heatIonSte = heatIonSte*ionDenUsed(elementXref(elem),ion)*&
!                         & grid%elemAbun(grid%abFileIndex(xP,yP,zP),elem)
!                     if (lgDebug) &
!                          & heatIonDif = heatIonDif*ionDenUsed(elementXref(elem),ion)*&
!                          &grid%elemAbun(grid%abFileIndex(xP,yP,zP),elem)
!                    heatSte    = heatSte + heatIonSte

!                     if (lgDebug) &
!                          & heatDIf    = heatDif + heatIonDif

                  end do ! end ion loop

               end do ! end element loop

               ! calculate the total heating and cooling integrals

!               if (lgDebug) then
!                  heatInt = heatSte + heatDIf
!               else
                  heatInt = heatSte
!               end if

               if (lgTraceHeating.and.taskid==0) then
                  write(57,*) 'Dust gas coll cool: ',gasDustColl_g
                  write(57,*) 'Heat photionization: ', heatInt
                  write(57,*) 'Heat photoelectric: ', photoelHeat_g
               end if

               coolInt = coolFF + coolRec + coolColl + coolCollH

               if (lgCompton) then
                  heatInt = heatInt+comptonHeat+comptonRecoilHeat
                  coolInt = coolInt+comptonCool
               end if

               if (lgDust .and. lgPhotoelectric) then
                  coolInt = coolInt+gasDustColl_g
                  heatInt = heatInt+photoelHeat_g
               end if

               if (lgTraceHeating.and.taskid==0) then
                  write(57,*) 'CoolInt: ',coolInt
                  write(57,*) 'HeatInt: ',heatInt
               end if

             end subroutine thermBalance



        ! this subroutine is the driver for the calculation of the emissivity
        ! from the heavy elements forbidden lines.
        subroutine forLines()
          implicit none

          integer        :: elem, ion ! counters

          ! re-initialize forbiddenLines
          forbiddenLines = 0.

          do elem = 3, nElements
             do ion = 1, min(elem+1, nstages)
                if (.not.lgElementOn(elem)) exit

                if (lgDataAvailable(elem, ion)) then

                   if (ion<min(elem+1,nstages).and.ionDenUsed(elementXref(elem), ion)>0.) then
                      call equilibrium(file_name = dataFile(elem, ion), &
                           &ionDenUp = ionDenUsed(elementXref(elem), ion+1)&
                           &/ionDenUsed(elementXref(elem), ion), &
                           & Te = TeUsed, Ne = NeUsed, &
                           & flineEm = forbiddenLines(elem, ion,:,:),rec = .false.)
                   else
                      call equilibrium(file_name = dataFile(elem, ion), ionDenUp = 0., &
                           & Te = TeUsed, Ne = NeUsed, flineEm = forbiddenLines(elem, ion,:,:),rec = .false.)
                   end if

                   forbiddenLines(elem, ion, :, :) = forbiddenLines(elem, ion, :, :)*&
                        & grid%elemAbun(grid%abFileIndex(xP,yP,zP),elem)*&
                        & ionDenUsed(elementXref(elem), ion)

                end if

             end do
          end do

          ! scale the forbidden lines emissivity to give units of [erg/s/Ngas]
          ! comment: the forbidden line emissivity is so far in units of cm^-1/s/Ngas
          !          the energy [erg] of unit wave number [cm^-1] is 1.9865e-16, hence
          !          the right units are obtained by multiplying by 1.9865e-16
          forbiddenLines = forbiddenLines*1.9865e-16


        end subroutine forLines

        subroutine ionBalance2(lgConv)
            implicit none

            real(kind=8)                   :: ionRatio(2,2), ionProd(2,2), denominator(2)
            real(kind=8)                   :: mat(nElements)
            real(kind=8)                   :: collIonH, expFact
            real(kind=8)                   :: sumMat
            real(kind=8)                   :: ratio
            real(kind=8)                   :: in(nElements) ! pop creation
            real(kind=8)                   :: out(nElements)! pop destruction

            real                           :: correction    ! used in lgNeInput = .t.
            real                           :: ionDenOld

            integer                :: nAugerelec
            integer                :: elem, ion     ! element and ion counters
            integer                :: i             ! counter
            integer                :: maxim
            integer                :: nElec         ! of  e's in the ion

            logical, intent(out)   :: lgConv        ! did ion bal converge?
            logical                :: lgConvEl(nelements)
            logical                :: lgNegative

            ! initialize variables
            correction  = 0.
            ionProd = 1.
            ionRatio = 0.
            denominator = 1.

            ! take into account collisional ionization of H
            ! Drake & Ulrich, ApJS42(1980)351
            expFact = 157893.94/TeUsed

            if (expFact > 75.) then
               ! prevents underflow of exponential factor in collIon
               expFact = 75.
            end if

            collIonH = 2.75E-16*TeUsed*sqrt(TeUsed)*&
                 & (157893.94/TeUsed+2.)*exp(-expFact)

            ! do H and He separately
            do elem = 1, 2
               do ion = 1, min(elem, nstages-1)

                  if (.not.lgElementOn(elem)) exit

                  if (elem>1) then
                     collIonH=0.
                     comptonRecoilIonH = 0.
                  end if

                  ionRatio(elem,ion) =  (photoIon(1,1,elem,ion)+collIonH*NeUsed+comptonRecoilIonH)/&
                       & (NeUsed*alphaTot(elem,ion)&
                       & +chex(elem,ion,1)*grid%Hden(grid%active(xP,yp,zP)) * &
                       & grid%ionDen(cellP,elementXref(1),1))
!print*, elem, ion, photoIon(1,1,elem,ion), collionH, NeUsed, comptonRecoilIonH
!print*, NeUsed, alphaTot(elem,ion), chex(elem,ion,1),grid%Hden(grid%active(xP,yp,zP)), grid%ionDen(cellP,elementXref(1),1)
!print*, ' '


               end do
            end do

            ! calculate the products of ionRatio
            do elem = 1, 2
               do ion = 1, min(elem, nstages-1)
                  if (.not.lgElementOn(elem)) exit

                  ! generate the product
                  do i = 1, ion
                     ionProd(elem, ion) = ionProd(elem,ion)*ionRatio(elem, i)
                  end do

               end do
            end do


            ! calculate denominators for final ion abundances
            do elem = 1, 2
               do ion = 1, min(elem, nstages-1)
                  if (.not.lgElementOn(elem)) exit

                  denominator(elem) = denominator(elem) + ionProd(elem, ion)

               end do
            end do


            ! calculate the new abundances for all ions
            do elem = 1, 2


               do j = 1, maxIterateX
                  lgConvEl = .true.

                  do ion = 1, min(elem+1, nstages)

                     if (.not.lgElementOn(elem)) exit
                     ionDenOld = grid%ionDen(cellP,elementXref(elem),ion)

                     if (ion == 1) then
                        grid%ionDen(cellP,elementXref(elem),ion) = 1./denominator(elem)
                     else

                        grid%ionDen(cellP,elementXref(elem),ion) = &
                             & ionProd(elem, ion-1)/denominator(elem)

                     end if

                     ! take back within the limit
                     if (grid%ionDen(cellP,elementXref(elem),ion) > xMax) &
                          & grid%ionDen(cellP,elementXref(elem),ion) = xMax
                     if (grid%ionDen(cellP,elementXref(elem),ion) < 1.e-20) &
                          & grid%ionDen(cellP,elementXref(elem),ion) = 1.e-20

                     ionDenUsed(elementXref(elem), ion) = &
                          & grid%ionDen(cellP,elementXref(elem),ion)
                     ! this was added to help MPI communication
                     ionDenTemp(cellP,elementXref(elem),ion) = &
                          & grid%ionDen(cellP,elementXref(elem),ion)

!print*,ionDenOld,grid%ionDen(cellP, elementXref(elem),ion), abs( (ionDenOld/grid%ionDen(cellP,&
!                                & elementXref(elem),ion))-1.)

                     if (grid%ionDen(cellP,elementXref(elem),ion)>1.e-10) then
                        if (abs( (ionDenOld/grid%ionDen(cellP,&
                             & elementXref(elem),ion))-1.) > 0.2) lgConvEl(elem) = .false.
                     end if


                  end do
                  if (lgConvEl(elem)) exit

               end do
            end do

            ! calculate the X(i+1)/X(i) ratio
            do elem = 3, nElements
               if (lgElementOn(elem)) then

                  do j = 1, maxIterateX

                     lgConvEl = .true.
                     out = 0.
                     in = 0.


                     do ion = 1, min(elem, nstages-1)
                        out(ion) = out(ion)+collIon(1,elem,ion)*NeUsed

!print*, elem, ion, collIon(1,elem,ion),NeUsed
                        do nshell = 1, nshells(elem,ion)
!print*, nshell, photoIon(1,nshell,elem,ion)


                           out(ion) = out(ion)+photoIon(1,nshell,elem,ion)
                           nAugerelec = nauger(elem,ion,nshell)

                           ! loop over electrons that come out of shell
                           ! with multiple electron ejection
                           do nelec=2,nAugerelec

                              ! this is highest possible stage of ionization -
                              ! do not want to ignore ionization that go beyond this
                              maxim = min( ion+nelec-1, min(nstages-1,elem))

                              if(grid%ionDen(cellP,elementXref(elem),maxim) > 1e-30 ) then
                                 ratio = dble(grid%ionDen(cellP,elementXref(elem),ion)) / &
                                      & dble(grid%ionDen(cellP,elementXref(elem),maxim))
                              else
                                 ratio = 1.
                              endif
                              ! yield here is fraction removing ion electrons
                              out(maxim) = out(maxim) +&
                                   & photoIon(1,nshell,elem,ion) * auger(elem,ion,nshell,nelec) * ratio

                           end do
                        end do

                     end do

                     do ion = 1, min(elem, nstages-1)
!print*, ion, alphaTot(elem,ion), chex(elem,ion,1), grid%Hden(grid%active(xP,yp,zP)),&
!                             & grid%ionDen(cellP,elementXref(1),1)
                        in(ion) = dble(NeUsed*alphaTot(elem,ion)&
                             & +chex(elem,ion,1)*grid%Hden(grid%active(xP,yp,zP)) * &
                             & grid%ionDen(cellP,elementXref(1),1))

                     end do

                     ! invert and solve bidiagonal matrix
                     mat(1) = 1.d0
                     sumMat = 1.d0

                     do ion = 2, min(nstages,elem+1)
                        mat(ion) = mat(ion-1) * out(ion-1)/in(ion-1)
                        sumMat = sumMat+mat(ion)
                     end do

                     sumMat = 1.d0/sumMat

                     lgNegative = .false.

                     grid%ionDen(cellP,elementXref(elem),min(nstages,elem+1)) = mat(min(nstages,elem+1))*sumMat
                     if (grid%ionDen(cellP,elementXref(elem),min(nstages,elem+1))<0.) lgNegative=.true.

                     do ion=1, min(elem,nstages-1)
                        if (nIterateMC < 2 ) then
                           if (grid%ionDen(cellP,elementXref(elem),ion)<=0.) &
                                &grid%ionDen(cellP,elementXref(elem),ion) = 1.e-5
                           ionDenOld = grid%ionDen(cellP,elementXref(elem),ion)
!print*, 'c', grid%ionDen(cellP,elementXref(elem),ion), mat(ion)*sumMat,sumMat, mat(ion), elem, ion

                           grid%ionDen(cellP,elementXref(elem),ion) = 10.**( 0.5*( &
                                & log10(grid%ionDen(cellP,elementXref(elem),ion))+ &
                                & log10(mat(ion)*sumMat)))

                           if (grid%ionDen(cellP,elementXref(elem),i)>1.e-15) then
                              if (abs(ionDenOld/grid%ionDen(cellP,elementXref(elem),i)-1) > 0.2 ) &
                                   &lgConv = .false.
                           end if
                        else
                           ionDenOld = grid%ionDen(cellP,elementXref(elem),ion)
                           grid%ionDen(cellP,elementXref(elem),ion) = mat(ion)*sumMat
!print*,ionDenOld,grid%ionDen(cellP, elementXref(elem),ion), abs( (ionDenOld/grid%ionDen(cellP,&
!                           & elementXref(elem),ion))-1.)
                           if (grid%ionDen(cellP,elementXref(elem),i)>1.e-10) then
                              if (abs( (ionDenOld/grid%ionDen(cellP,&
                                   & elementXref(elem),ion))-1.) > 0.2) lgConvEl(elem) = .false.

                           end if
                        end if
                       if (grid%ionDen(cellP,elementXref(elem),ion) < 0.) lgNegative = .true.
                        if (lgNegative) then
                           print*, '! ionBalance: negative populations [elem,ion]', elem, ion, mat(ion), summat
                        end if


                        if (grid%ionDen(cellP,elementXref(elem),ion) > xMax) &
                             & grid%ionDen(cellP,elementXref(elem),ion) = xMax
                        if (grid%ionDen(cellP,elementXref(elem),ion) < 1.e-20) &
                             & grid%ionDen(cellP,elementXref(elem),ion) = 1.e-20


                     end do


                     ionDenUsed(elementXref(elem), :) = &
                          & grid%ionDen(cellP,elementXref(elem),:)
                     ! this was added to help MPI communication
                     ionDenTemp(cellP,elementXref(elem),:) = &
                          & grid%ionDen(cellP,elementXref(elem),:)

                     if (lgConvEl(elem)) exit

                  end do

               end if
            end do

            ! calculate new Ne
            NeUsed = 0.
            do elem = 1, nElements
               do ion = 2, min(elem+1, nstages)
                  if (lgElementOn(elem)) then
                     if( ionDenUsed(elementXref(elem),ion) >= 1.e-10) &
                          & NeUsed = NeUsed + (ion-1)*&
                          &grid%elemAbun(grid%abFileIndex(xP,yP,zP),elem)*&
                          &ionDenUsed(elementXref(elem), ion)
                  end if
               end do
            end do

            NeUsed = NeUsed * grid%Hden(cellP)


            if (NeUsed==0. .and. lgVerbose) print*, '! ionBalance [warning]: cell ', xP,yP,zP, &
                 &'; NeUsed = ',  NeUsed

            if (NeUsed == 0.) then
               NeUsed = 1.
            end if

            if (LgNeInput) then
               correction = NeUsed/grid%NeInput(cellP)
               grid%Hden(cellP) = grid%Hden(cellP)/correction
            end if

            ! this was added to help MPI implementation
            NeTemp(cellP) = NeUsed

            grid%Ne(cellP) = NeUsed

            lgConv=.true.
            do elem = 1, nelements
               if (.not.lgConvEl(elem)) then
                  print*, "! ionBalance: [warning] convergence not reached after ", &
                       &maxIterateX, " steps [element,cell] ", elem, cellp
                  lgConv=.false.
               end if
            end do

             ! this was added to help MPI implementation
             NeTemp(cellP) = NeUsed

           end subroutine ionBalance2

        recursive subroutine ionBalance(lgConv)
            implicit none

            real(kind=8)                   :: ionRatio(2,2), ionProd(2,2), denominator(2)
            real(kind=8)                   :: mat(nElements)
            real(kind=8)                   :: collIonH, expFact
            real(kind=8)                   :: sumMat
            real(kind=8)                   :: ratio
            real(kind=8)                   :: revRate       ! reverse charge exchange rate
            real(kind=8)                   :: in(nElements) ! pop creation
            real(kind=8)                   :: out(nElements)! pop destruction

            real, dimension(nElements, nstages) :: &
                 & deltaE_k      ! deltaE/k [K]

            real, dimension(nElements, nstages,4) :: &
                 & chex         ! ch exchange coeff in cm^3/s

            real                   :: deltaHI       ! delta(X(H0))
            real                   :: deltaHeI      ! delta(X(He0))
            real                   :: deltaHeII     ! delta(X(HeII))
            real, save             :: HIOld         ! X(H0) from last iteration
            real, save             :: HeIOld        ! X(He0) from last iteration
            real, save             :: HeIIOld       ! X(He+) from last iteration
            real, parameter        :: limit = 0.01  ! convergence limit

            real                   :: t4            ! TeUsed/10000.

            real                   :: correction    ! used in lgNeInput = .t.

            integer                :: nAugerelec
            integer                :: elem, ion     ! element and ion counters
            integer                :: i             ! counter
            integer                :: maxim
            integer                :: nElec         ! of  e's in the ion
            integer                :: outShell      ! byproduct of proc to get stat weights
            integer                :: g0,g1         ! stat weights

            logical, intent(out)   :: lgConv        ! did ion bal converge?
            logical                :: lgNegative


            ! initialize variables
            correction  = 0.
            ionRatio    = 0.
            ionProd     = 1.
            denominator = 1.
            lgConv      = .true.

            ! take into account collisional ionization of H
            ! Drake & Ulrich, ApJS42(1980)351
            expFact = 157893.94/TeUsed

            if (expFact > 75.) then
                ! prevents underflow of exponential factor in collIon
                expFact = 75.
            end if

!print*, '1a ',   expFact
            collIonH = 2.75E-16*TeUsed*sqrt(TeUsed)*&
                 & (157893.94/TeUsed+2.)*exp(-expFact)
!print*, '1a '

            ! get HIOld, HeIOld, HeIIOld from last iteration
            HIOld   = ionDenUsed(elementXref(1),1)
            HeIOld  = ionDenUSed(elementXref(2),1)
            HeIIOld = ionDenUsed(elementXref(2),2)

            ! the set of charge exchange coeffs is not complete; the following might need
            ! to be changed when a more complete set is available



            chex          = 0.
            deltaE_k      = 0.

            chex(2,1,:)  = (/7.47e-6, 2.06, 9.93,-3.89/)! He0
            chex(2,2,:)  = (/1.e-5  , 0.  , 0.  , 0./)  ! He+
            chex(3,2,:)  = (/1.26   , 0.96,3.02 ,-0.65/)! Li+
            chex(3,3,:)  = (/1.e-5  , 0.  , 0.  , 0. /) ! Li+2
            chex(4,2,:)  = (/1.e-5  , 0.  , 0.  , 0. /) ! Be+
            chex(4,3,:)  = (/1.e-5  , 0. , 0.  , 0. /) ! Be+2
            chex(4,4,:)  = (/5.17   , 0.82, -.69, -1.12 /)! Be+3
            chex(5,2,:)  = (/2.e-2  , 0.  , 0.  , 0. /) ! B+
            chex(5,3,:)  = (/1.e-5  , 0.  , 0.  , 0. /) ! B+2
            chex(5,4,:)  = (/5.27e-1, 0.76,-0.63,-1.17/)! B+3
            chex(6,1,:)  = (/1.76e-9, 8.33, 4278.78, -6.41/)! C0
            chex(6,2,:)  = (/1.67e-4, 2.79, 304.72, -4.07/)! C+
            chex(6,3,:)  = (/3.25   , 0.21, 0.19, -3.29/)! C+2
            chex(6,4,:)  = (/332.46 ,-0.11,-0.995,-1.58e-3/)! C+3
            chex(7,1,:)  = (/1.01e-3,-0.29,-0.92, -8.38/)! N0
            chex(7,2,:)  = (/3.05e-1, 0.60, 2.65, -0.93/)! N+
            chex(7,3,:)  = (/4.54   , 0.57,-0.65, -0.89/)! N2+
            chex(7,4,:)  = (/3.28   , 0.52,-0.52, -0.19/)! N3+
            chex(8,1,:)  = (/1.04   , 3.15e-2, -0.61, -9.73/)! O0
            chex(8,2,:)  = (/1.04   , 0.27, 2.02, -5.92/)! O+
            chex(8,3,:)  = (/3.98   , 0.26, 0.56, -2.62/)! O2+
            chex(8,4,:)  = (/2.52e-1, 0.63, 2.08, -4.16/)! O3+
            chex(9,2,:)  = (/1.e-5  , 0. , 0.  , 0./) ! F+
            chex(9,3,:)  = (/9.86   , 0.29,-0.21,-1.15/) ! F+2
            chex(9,4,:)  = (/7.15e-1, 1.21,-0.70,-0.85/) ! F3+
            chex(10,2,:) = (/1.e-5  , 0.  , 0.  , 0.  /) ! Ne+
            chex(10,3,:) = (/14.73  , 4.52e-2, -0.84, -0.31 /) ! Ne+2
            chex(10,4,:) = (/6.47   , 0.54 , 3.59 , -5.22 /) ! Ne+3
            chex(11,2,:) = (/1.e-5  , 0.   , 0.   , 0. /) ! Na+
            chex(11,3,:) = (/1.33   , 1.15 , 1.20 , -0.32 /)! Na+2
            chex(11,4,:) = (/1.01e-1, 1.34 , 10.05, -6.41 /)! Na+3
            chex(12,2,:) = (/8.58e-5, 2.49e-3, 2.93e-2, -4.33 /)! Mg+
            chex(12,3,:) = (/6.49   , 0.53 , 2.82, -7.63 /) ! Mg+2
            chex(12,4,:) = (/6.36   , 0.55 , 3.86, -5.19 /) ! Mg+3
            chex(13,2,:) = (/1.e-5  , 0.   , 0.  , 0./) ! Al+
            chex(13,3,:) = (/7.11e-5, 4.12 , 1.72e4, -22.24/)! Al+2
            chex(13,4,:) = (/7.52e-1, 0.77 , 6.24, -5.67/) ! Al+3
            chex(14,2,:) = (/1.23   , 0.24 , 3.17, 4.18e-3/) ! Si+
            chex(14,3,:) = (/4.900e-1, -8.74e-2, -0.36, -0.79/)! Si+2
            chex(14,4,:) = (/7.58   , 0.37 , 1.06, -4.09/)! Si+3
            chex(16,1,:) = (/3.82e-7, 11.10, 2.57e4, -8.22/)! S0
            chex(16,2,:) = (/1.e-5  , 0.   , 0.   ,0. /)! S+
            chex(16,3,:) = (/2.29   , 4.02e-2, 1.59, -6.06/)! S+2
            chex(16,4,:) = (/6.44   , 0.13 , 2.69 , -5.69/)! S+3
            chex(18,2,:) = (/1.e-5  , 0.   , 0.    , 0./) ! Ar+
            chex(18,3,:) = (/4.57   , 0.27 , -0.18 , -1.57/)! Ar+2
            chex(18,4,:) = (/6.37   , 2.12 , 10.21 , -6.22/)! Ar+3
            chex(18,3,:) = (/3.17e-2, 2.12 , 12.06 , -0.40/)! Ca+2
            chex(18,4,:) = (/2.68   , 0.69 , -0.68 , -4.47/)! Ca+3
            chex(26,2,:) = (/1.26   , 7.72e-2, -0.41, -7.31/)! Fe+
            chex(26,3,:) = (/3.42   , 0.51 , -2.06 , -8.99/)! Fe+2.


            deltaE_k(7,1) = 10863.
            deltaE_k(8,1) = 2205.

            chex(:,:,1) = chex(:,:,1)*1.e-9


            t4 = TeUsed/10000.


            ! calculate the X(i+1)/X(i) ratio
            do elem = 1, 2
                do ion = 1, min(elem, nstages-1)
                    if (.not.lgElementOn(elem)) exit
!print*, 'aa ', chex(elem,ion,4), t4, chex(elem,ion,4)*t4

                    if (TeUsed < 6000. .or. TeUsed>5.e4) then
                       chex(elem,ion,1) = 0.
                    else
                       chex(elem,ion,1) = chex(elem,ion,1)*(t4**chex(elem,ion,2))*&
                            & (1.+chex(elem,ion,3)*exp(chex(elem,ion,4)*t4))
                    end if

                    if (chex(elem,ion,1) < 0. ) chex(elem,ion,1) = 0.


                    ! find the number of electron in this ion
                    nElec = elem - ion +1

                    ! find the stat weights
                    call getOuterShell(elem, nElec, outShell, g0, g1)

                    ! calculate the reverse charge exchange rate (only if deltaE_k > 1.)
                    if ( deltaE_k(elem,ion) > 1.) then

!print*, 'bb ', -deltaE_k(elem,ion)/TeUsed
                       revRate = chex(elem,ion,1) * &
                            & (1.-grid%ionDen(cellP, &
                            & elementXref(1),1))*grid%Hden(cellP)*&
                            & 2. * exp(-deltaE_k(elem,ion)/TeUsed)/(real(g0)/real(g1))
                    else
                       revRate = 0.
                    end if

                    if ( (elem>1) .or. (ion>1) ) collIonH = 0.

                    ! calculate the X(i+1)/X(i) ratio

                       ionRatio(elem,ion) = (photoIon(1,outshell,elem,ion)+&
                            ! ADDED *NeUsed
                            & collIonH*NeUsed+revRate)/&
                            & (NeUsed*alphaTot(elem,ion)&
                            & +chex(elem,ion,1)*grid%Hden(grid%active(xP,yp,zP)) * &
                            & grid%ionDen(cellP,elementXref(1),1))

                 end do
              end do

              ! calculate the products of ionRatio
              do elem = 1, 2
                 do ion = 1, min(elem, nstages-1)
                    if (.not.lgElementOn(elem)) exit

                    ! generate the product
                    do i = 1, ion
                       ionProd(elem, ion) = ionProd(elem,ion)*ionRatio(elem, i)

                    end do

                 end do
              end do

              ! calculate denominators for final ion abundances
              do elem = 1, 2
                do ion = 1, min(elem, nstages-1)
                    if (.not.lgElementOn(elem)) exit

                    denominator(elem) = denominator(elem) + ionProd(elem, ion)

                 end do
              end do


              ! calculate the new abundances for all ions
              do elem = 1, 2
                 do ion = 1, min(elem+1, nstages)

                    if (.not.lgElementOn(elem)) exit


                    if (ion == 1) then
                       grid%ionDen(cellP,elementXref(elem),ion) = 1./denominator(elem)
                    else

                       grid%ionDen(cellP,elementXref(elem),ion) = &
                            & ionProd(elem, ion-1)/denominator(elem)

                    end if

                    ! take back within the limit
                    if (grid%ionDen(cellP,elementXref(elem),ion) > xMax) &
                         & grid%ionDen(cellP,elementXref(elem),ion) = xMax
                    if (grid%ionDen(cellP,elementXref(elem),ion) < 1.e-20) &
                         & grid%ionDen(cellP,elementXref(elem),ion) = 1.e-20


                    ionDenUsed(elementXref(elem), ion) = &
                         & grid%ionDen(cellP,elementXref(elem),ion)
                    ! this was added to help MPI communication
                    ionDenTemp(cellP,elementXref(elem),ion) = &
                         & grid%ionDen(cellP,elementXref(elem),ion)


                 end do
              end do

              ! Now do the Heavy elements

              ! calculate the X(i+1)/X(i) ratio
              do elem = 3, nElements
                 if (lgElementOn(elem)) then

                    out = 0.
                    in = 0.

                    do ion = 1, min(elem, nstages-1)
                       out(ion) = out(ion)+collIon(1,elem,ion)*NeUsed

                       do nshell = 1, nshells(elem,ion)
!print*, 'out', ion, nshell,    out(ion),    photoIon(1,nshell,elem,ion)
                          out(ion) = out(ion)+photoIon(1,nshell,elem,ion)
                          nAugerelec = nauger(elem,ion,nshell)

                          ! loop over electrons that come out of shell
                          ! with multiple electron ejection
                          do nelec=2,nAugerelec

                             ! this is highest possible stage of ionization -
                             ! do not want to ignore ionization that go beyond this
                             maxim = min( ion+nelec-1, min(nstages-1,elem))

                             if(grid%ionDen(cellP,elementXref(elem),maxim) > 1e-30 ) then
                                ratio = dble(grid%ionDen(cellP,elementXref(elem),ion)) / &
                                     & dble(grid%ionDen(cellP,elementXref(elem),maxim))
                             else
                                ratio = 1.
                             endif
                             ! yield here is fraction removing ion electrons
                             out(maxim) = out(maxim) +&
                                  & photoIon(1,nshell,elem,ion) * auger(elem,ion,nshell,nelec) * ratio

                          end do
                       end do

                    end do


                    do ion = 1, min(elem, nstages-1)
                       in(ion) = dble(NeUsed*alphaTot(elem,ion)&
                            & +chex(elem,ion,1)*grid%Hden(grid%active(xP,yp,zP)) * &
                            & grid%ionDen(cellP,elementXref(1),1))

                    end do

                    ! invert and solve bidiagonal matrix
                    mat(1) = 1.d0
                    sumMat = 1.d0

                    do ion = 2, min(nstages,elem+1)
!print*, ion,  mat(ion-1), out(ion-1), in(ion-1)
                       mat(ion) = mat(ion-1) * out(ion-1)/in(ion-1)
                       sumMat = sumMat+mat(ion)
                    end do

                    sumMat = 1.d0/sumMat

                    lgNegative = .false.

                    grid%ionDen(cellP,elementXref(elem),min(nstages,elem+1)) = &
                         &mat(min(nstages,elem+1))*sumMat
                    if (grid%ionDen(cellP,elementXref(elem),min(nstages,elem+1))<0.) lgNegative=.true.

                    do ion=1, min(elem,nstages-1)
                       if (nIterateX < 2 .and. mat(ion)*sumMat>0. ) then
                          if (grid%ionDen(cellP,elementXref(elem),ion)<=0.) &
                               &grid%ionDen(cellP,elementXref(elem),ion) = 1.e-5
!                          ionDenOld = grid%ionDen(cellP,elementXref(elem),ion)
!print*, 'd', grid%ionDen(cellP,elementXref(elem),ion), mat(ion)*sumMat,sumMat, mat(ion), elem, ion
                          grid%ionDen(cellP,elementXref(elem),ion) = 10.**( 0.5*( &
                               & log10(grid%ionDen(cellP,elementXref(elem),ion))+ &
                               & log10(mat(ion)*sumMat)))

!                         if (grid%ionDen(cellP,elementXref(elem),i)>1.e-15) then
!                             if (abs(ionDenOld/grid%ionDen(cellP,elementXref(elem),i)-1) > 0.2 ) lgConv = .false.
!                          end if
                       else
!                           ionDenOld = grid%ionDen(cellP,elementXref(elem),ion)
                           grid%ionDen(cellP,elementXref(elem),ion) = mat(ion)*sumMat
!                           if (grid%ionDen(cellP,elementXref(elem),i)>1.e-10) then
!                              if (abs( (ionDenOld/grid%ionDen(cellP,&
!                                   & elementXref(elem),ion))-1.) > 0.2) lgConvEl(elem) = .false.
!                           end if
                        end if
                        if (grid%ionDen(cellP,elementXref(elem),ion) < 0.) lgNegative = .true.
                        if (lgNegative) then
                           print*, '! ionBalance: negative populations [elem,ion]', elem, ion, mat(ion), summat
                        end if

                        if (grid%ionDen(cellP,elementXref(elem),ion) > xMax) &
                             & grid%ionDen(cellP,elementXref(elem),ion) = xMax
                        if (grid%ionDen(cellP,elementXref(elem),ion) < 1.e-20) &
                             & grid%ionDen(cellP,elementXref(elem),ion) = 1.e-20


                     end do


                     ionDenUsed(elementXref(elem), :) = &
                          & grid%ionDen(cellP,elementXref(elem),:)
                     ! this was added to help MPI communication
                     ionDenTemp(cellP,elementXref(elem),:) = &
                          & grid%ionDen(cellP,elementXref(elem),:)

                  end if
               end do



               ! calculate new Ne
               NeUsed = 0.
               do elem = 1, nElements
                  do ion = 2, min(elem+1, nstages)
                     if (lgElementOn(elem)) then
                        if( ionDenUsed(elementXref(elem),ion) >= 1.e-10) &
                             & NeUsed = NeUsed + (ion-1)*&
                             &grid%elemAbun(grid%abFileIndex(xP,yP,zP),elem)*&
                             &ionDenUsed(elementXref(elem), ion)
                     end if

                  end do
               end do

               NeUsed = NeUsed * grid%Hden(cellP)

               if (NeUsed==0. .and. lgVerbose) print*, '! ionBalance [warning]: cell ', xP,yP,zP, &
                    &'; NeUsed = ',  NeUsed

               if (NeUsed == 0.) then
                  NeUsed = 1.
               end if

               if (LgNeInput) then
                  correction = NeUsed/grid%NeInput(cellP)
                  grid%Hden(cellP) = grid%Hden(cellP)/correction
               end if

               ! this was added to help MPI implementation
              NeTemp(cellP) = NeUsed

              grid%Ne(cellP) = NeUsed


              ! calculate the residuals
              deltaHI   = (ionDenUsed(elementXref(1),1) - HIOld)   / HIOld
              deltaHeI  = (ionDenUsed(elementXref(2),1) - HeIOld)  / HeIOld
              deltaHeII = (ionDenUsed(elementXref(2),2) - HeIIOld) / HeIIOld

              ! check for convergence
              if ( ( (abs(deltaHI)>limit) .or. (abs(deltaHeI)>limit) .or. &
                   &(abs(deltaHeII)>limit) ) .and. (nIterateX<maxIterateX) ) then

                 ! stepa up X iteration
                 nIterateX = nIterateX + 1
                 call ionBalance(lgConv)
                 return

              else if (( (abs(deltaHI)>limit) .or. (abs(deltaHeI)>limit) .or. &
                   &(abs(deltaHeII)>limit) ) .and. (nIterateX == maxIterateX) ) then

                 if (lgTalk)  print*, "! ionBalance: [warning] convergence not reached after ", &
                      &maxIterateX, " steps. Finishing up..."

                 lgConv = .false.
              end if

              ! this was added to help MPI implementation
              NeTemp(cellP) = NeUsed

            end subroutine ionBalance

        recursive subroutine ionBalance3(lgConv)
            implicit none

            logical, intent(out)   :: lgConv        ! did ion bal converge?

            ! local variables
            real                   :: collIon       ! collisional ionization of H
            real                   :: correction    ! used in lgNeInput = .t.
            real                   :: deltaHI       ! delta(X(H0))
            real                   :: deltaHeI      ! delta(X(He0))
            real                   :: deltaHeII     ! delta(X(HeII))
            real                   :: expFact       ! exponential factor
            real                   :: t4            ! TeUsed/10000.
            real, save             :: HIOld         ! X(H0) from last iteration
            real, save             :: HeIOld        ! X(He0) from last iteration
            real, save             :: HeIIOld       ! X(He+) from last iteration
            real, parameter        :: limit = 0.01  ! convergence limit

            integer                :: elem, ion     ! element and ion counters
            integer                :: g0,g1         ! stat weights
            integer                :: i             ! counter
            integer                :: nElec         ! of  e's in the ion
            integer                :: outShell      ! byproduct of proc to get stat weights


            real                   :: revRate       ! reverse charge exchange rate

            double precision, dimension(nELements) :: &
                 & denominator   ! denominator of final ion abundance


            real, dimension(nElements, nstages) :: &
                 & deltaE_k      ! deltaE/k [K]

            double precision, dimension(nElements, nstages) :: &
                 & ionRatio, &   ! X(i+1)/X(+i)
                 & ionProd       ! ionRatio products
            real, dimension(nElements, nstages,4) :: &
                 & chex         ! ch exchange coeff in cm^3/s


            ! initialize variables
            correction  = 0.
            ionRatio    = 0.
            ionProd     = 1.
            denominator = 1.
            lgConv      = .true.

            ! take into account collisional ionization of H
            ! Drake & Ulrich, ApJS42(1980)351
            expFact = 157893.94/TeUsed

            if (expFact > 75.) then
                ! prevents underflow of exponential factor in collIon
                expFact = 75.
            end if

!print*, '2a ', expFact

            collIon = 2.75E-16*TeUsed*sqrt(TeUsed)*&
                 & (157893.94/TeUsed+2.)*exp(-expFact)
!print*, '2a'

            ! get HIOld, HeIOld, HeIIOld from last iteration
            HIOld   = ionDenUsed(elementXref(1),1)
            HeIOld  = ionDenUSed(elementXref(2),1)
            HeIIOld = ionDenUsed(elementXref(2),2)

            ! the set of charge exchange coeffs is not complete; the following might need
            ! to be changed when a more complete set is available



            chex          = 0.
            deltaE_k      = 0.

            chex(2,1,:)  = (/7.47e-6, 2.06, 9.93,-3.89/)! He0
            chex(2,2,:)  = (/1.e-5  , 0.  , 0.  , 0./)  ! He+
            chex(3,2,:)  = (/1.26   , 0.96,3.02 ,-0.65/)! Li+
            chex(3,3,:)  = (/1.e-5  , 0.  , 0.  , 0. /) ! Li+2
            chex(4,2,:)  = (/1.e-5  , 0.  , 0.  , 0. /) ! Be+
            chex(4,3,:)  = (/1.e-5  , 0. , 0.  , 0. /) ! Be+2
            chex(4,4,:)  = (/5.17   , 0.82, -.69, -1.12 /)! Be+3
            chex(5,2,:)  = (/2.e-2  , 0.  , 0.  , 0. /) ! B+
            chex(5,3,:)  = (/1.e-5  , 0.  , 0.  , 0. /) ! B+2
            chex(5,4,:)  = (/5.27e-1, 0.76,-0.63,-1.17/)! B+3
            chex(6,1,:)  = (/1.76e-9, 8.33, 4278.78, -6.41/)! C0
            chex(6,2,:)  = (/1.67e-4, 2.79, 304.72, -4.07/)! C+
            chex(6,3,:)  = (/3.25   , 0.21, 0.19, -3.29/)! C+2
            chex(6,4,:)  = (/332.46 ,-0.11,-0.995,-1.58e-3/)! C+3
            chex(7,1,:)  = (/1.01e-3,-0.29,-0.92, -8.38/)! N0
            chex(7,2,:)  = (/3.05e-1, 0.60, 2.65, -0.93/)! N+
            chex(7,3,:)  = (/4.54   , 0.57,-0.65, -0.89/)! N2+
            chex(7,4,:)  = (/3.28   , 0.52,-0.52, -0.19/)! N3+
            chex(8,1,:)  = (/1.04   , 3.15e-2, -0.61, -9.73/)! O0
            chex(8,2,:)  = (/1.04   , 0.27, 2.02, -5.92/)! O+
            chex(8,3,:)  = (/3.98   , 0.26, 0.56, -2.62/)! O2+
            chex(8,4,:)  = (/2.52e-1, 0.63, 2.08, -4.16/)! O3+
            chex(9,2,:)  = (/1.e-5  , 0. , 0.  , 0./) ! F+
            chex(9,3,:)  = (/9.86   , 0.29,-0.21,-1.15/) ! F+2
            chex(9,4,:)  = (/7.15e-1, 1.21,-0.70,-0.85/) ! F3+
            chex(10,2,:) = (/1.e-5  , 0.  , 0.  , 0.  /) ! Ne+
            chex(10,3,:) = (/14.73  , 4.52e-2, -0.84, -0.31 /) ! Ne+2
            chex(10,4,:) = (/6.47   , 0.54 , 3.59 , -5.22 /) ! Ne+3
            chex(11,2,:) = (/1.e-5  , 0.   , 0.   , 0. /) ! Na+
            chex(11,3,:) = (/1.33   , 1.15 , 1.20 , -0.32 /)! Na+2
            chex(11,4,:) = (/1.01e-1, 1.34 , 10.05, -6.41 /)! Na+3
            chex(12,2,:) = (/8.58e-5, 2.49e-3, 2.93e-2, -4.33 /)! Mg+
            chex(12,3,:) = (/6.49   , 0.53 , 2.82, -7.63 /) ! Mg+2
            chex(12,4,:) = (/6.36   , 0.55 , 3.86, -5.19 /) ! Mg+3
            chex(13,2,:) = (/1.e-5  , 0.   , 0.  , 0./) ! Al+
            chex(13,3,:) = (/7.11e-5, 4.12 , 1.72e4, -22.24/)! Al+2
            chex(13,4,:) = (/7.52e-1, 0.77 , 6.24, -5.67/) ! Al+3
            chex(14,2,:) = (/1.23   , 0.24 , 3.17, 4.18e-3/) ! Si+
            chex(14,3,:) = (/4.900e-1, -8.74e-2, -0.36, -0.79/)! Si+2
            chex(14,4,:) = (/7.58   , 0.37 , 1.06, -4.09/)! Si+3
            chex(16,1,:) = (/3.82e-7, 11.10, 2.57e4, -8.22/)! S0
            chex(16,2,:) = (/1.e-5  , 0.   , 0.   ,0. /)! S+
            chex(16,3,:) = (/2.29   , 4.02e-2, 1.59, -6.06/)! S+2
            chex(16,4,:) = (/6.44   , 0.13 , 2.69 , -5.69/)! S+3
            chex(18,2,:) = (/1.e-5  , 0.   , 0.    , 0./) ! Ar+
            chex(18,3,:) = (/4.57   , 0.27 , -0.18 , -1.57/)! Ar+2
            chex(18,4,:) = (/6.37   , 2.12 , 10.21 , -6.22/)! Ar+3
            chex(18,3,:) = (/3.17e-2, 2.12 , 12.06 , -0.40/)! Ca+2
            chex(18,4,:) = (/2.68   , 0.69 , -0.68 , -4.47/)! Ca+3
            chex(26,2,:) = (/1.26   , 7.72e-2, -0.41, -7.31/)! Fe+
            chex(26,3,:) = (/3.42   , 0.51 , -2.06 , -8.99/)! Fe+2.


            deltaE_k(7,1) = 10863.
            deltaE_k(8,1) = 2205.

            chex(:,:,1) = chex(:,:,1)*1.e-9


            t4 = TeUsed/10000.


            ! calculate the X(i+1)/X(i) ratio
            do elem = 1, nElements
                do ion = 1, min(elem, nstages-1)
                    if (.not.lgElementOn(elem)) exit

!print*, 'cc ', chex(elem,ion,4)*t4

                    chex(elem,ion,1) = chex(elem,ion,1)*(t4**chex(elem,ion,2))*&
                         & (1.+chex(elem,ion,3)*exp(chex(elem,ion,4)*t4))

                    if (chex(elem,ion,1) < 0. ) chex(elem,ion,1) = 0.
                    if (TeUsed < 6000. .or. TeUsed>5.e4) chex(elem,ion,1) = 0.

                    ! find the number of electron in this ion
                    nElec = elem - ion +1

                    ! find the stat weights
                    call getOuterShell(elem, nElec, outShell, g0, g1)

                    ! calculate the reverse charge exchange rate (only if deltaE_k > 1.)
                    if ( deltaE_k(elem,ion) > 1.) then

!print*, 'dd ', -deltaE_k(elem,ion)/TeUsed

                       revRate = chex(elem,ion,1) * &
                            & (1.-grid%ionDen(cellP, &
                            & elementXref(1),1))*grid%Hden(cellP)*&
                            & 2. * exp(-deltaE_k(elem,ion)/TeUsed)/(real(g0)/real(g1))
                    else
                       revRate = 0.
                    end if

                    if ( (elem>1) .or. (ion>1) ) collIon = 0.

                    ! calculate the X(i+1)/X(i) ratio

                       ionRatio(elem,ion) = (photoIon(1,outshell,elem,ion)+&
                            & collIon+revRate)/&
                            & (NeUsed*alphaTot(elem,ion)&
                            & +chex(elem,ion,1)*grid%Hden(grid%active(xP,yp,zP)) * &
                            & grid%ionDen(cellP,elementXref(1),1))

                 end do
              end do

            ! calculate the products of ionRatio
            do elem = 1, nElements
                do ion = 1, min(elem, nstages-1)
                   if (.not.lgElementOn(elem)) exit

                   ! generate the product
                   do i = 1, ion
                       ionProd(elem, ion) = ionProd(elem,ion)*ionRatio(elem, i)

                   end do

                end do
            end do

            ! calculate denominators for final ion abundances
            do elem = 1, nElements
                do ion = 1, min(elem, nstages-1)
                    if (.not.lgElementOn(elem)) exit

                    denominator(elem) = denominator(elem) + ionProd(elem, ion)

               end do
            end do


            ! calculate the new abundances for all ions
            do elem = 1, nElements
                do ion = 1, min(elem+1, nstages)

                    if (.not.lgElementOn(elem)) exit


                    if (ion == 1) then
                        grid%ionDen(cellP,elementXref(elem),ion) = 1./denominator(elem)
                    else

                        grid%ionDen(cellP,elementXref(elem),ion) = &
                             & ionProd(elem, ion-1)/denominator(elem)

                   end if

                   ! take back within the limit
                   if (grid%ionDen(cellP,elementXref(elem),ion) > xMax) &
                        & grid%ionDen(cellP,elementXref(elem),ion) = xMax
                   if (grid%ionDen(cellP,elementXref(elem),ion) < 1.e-20) &
                        & grid%ionDen(cellP,elementXref(elem),ion) = 1.e-20


                   ionDenUsed(elementXref(elem), ion) = &
                        & grid%ionDen(cellP,elementXref(elem),ion)
                   ! this was added to help MPI communication
                   ionDenTemp(cellP,elementXref(elem),ion) = &
                        & grid%ionDen(cellP,elementXref(elem),ion)


                end do
            end do

            ! calculate new Ne
            NeUsed = 0.
            do elem = 1, nElements
                do ion = 2, min(elem+1, nstages)
                    if (lgElementOn(elem)) then
                      if( ionDenUsed(elementXref(elem),ion) >= 1.e-10) &
                           & NeUsed = NeUsed + (ion-1)*&
                           &grid%elemAbun(grid%abFileIndex(xP,yP,zP),elem)*&
                           &ionDenUsed(elementXref(elem), ion)
                    end if

                end do
            end do

            NeUsed = NeUsed * grid%Hden(cellP)

            if (NeUsed==0.) print*, '! ionBalance [warning]: cell ', xP,yP,zP, &
                 &'; NeUsed = ',  NeUsed

            if (NeUsed == 0.) then
               NeUsed = 1.
            end if

            if (LgNeInput) then
               correction = NeUsed/grid%NeInput(cellP)
               grid%Hden(cellP) = grid%Hden(cellP)/correction
            end if

            ! this was added to help MPI implementation
            NeTemp(cellP) = NeUsed

            grid%Ne(cellP) = NeUsed


            ! calculate the residuals
            deltaHI   = (ionDenUsed(elementXref(1),1) - HIOld)   / HIOld
            deltaHeI  = (ionDenUsed(elementXref(2),1) - HeIOld)  / HeIOld
            deltaHeII = (ionDenUsed(elementXref(2),2) - HeIIOld) / HeIIOld

            ! check for convergence
            if ( ( (abs(deltaHI)>limit) .or. (abs(deltaHeI)>limit) .or. &
                 &(abs(deltaHeII)>limit) ) .and. (nIterateX<maxIterateX) ) then

                ! stepa up X iteration
                nIterateX = nIterateX + 1
                call ionBalance(lgConv)
                return

            else if (( (abs(deltaHI)>limit) .or. (abs(deltaHeI)>limit) .or. &
                 &(abs(deltaHeII)>limit) ) .and. (nIterateX == maxIterateX) ) then

               if (lgTalk)  print*, "! ionBalance: [warning] convergence not reached after ", &
                    &maxIterateX, " steps. Finishing up..."

                lgConv = .false.
            end if

            ! this was added to help MPI implementation
            NeTemp(cellP) = NeUsed

        end subroutine ionBalance3

        subroutine calcAlpha()
            implicit none

            real, dimension(nElements, nstages) &
                 &:: diRec  ! total recombination coeffs

            integer :: nelectrons

            ! zero out alphaTot and contributors
            alphaTot = 0.
            diRec    = 0.

            call dielectronic(diRec)

            ! calculate radiative recombination part first

            if (lgBadnell) then

              ! use Badnell's data where available (up to Na-like)
               do elem = 1, nElements
                  if (lgElementOn(elem)) then
                     do ion = 1, min(nstages-1, elem)
                        nelectrons = elem-ion

                        if (nelectrons <= 11) then
                           alphaTot(elem, ion) = RRbadnell(elem, nelectrons)
                           alphaTot(elem, ion) = alphaTot(elem, ion)+DRbadnell(elem, nelectrons)
                        else
                           if ( (elem == 14   .and. ion==1) .or. &
                                & (elem == 16 .and. ion==3) .or. &
                                & (elem == 18 .and. ion==5) .or. &
                                & (elem == 20 .and. ion==7) .or. &
                                & (elem == 26 .and. ion==12) .or.&
                                & (elem == 6  .and. ion>=1 .and. ion<=6) .or.&
                                & (elem == 7  .and. ion>=1 .and. ion<=7) .or.&
                                & (elem == 8  .and. ion>=1 .and. ion<=8) ) then

                              if (lgNahar) then
                                 alphaTot(elem, ion) = nahar(elem,ion,TeUsed)
                              else
                                 alphaTot(elem, ion) = radRecFit(elem, elem-ion+1)
                              end if
                           else
                              alphaTot(elem, ion) = radRecFit(elem, elem-ion+1)
                           end if

                           if ( .not.lgNahar .or. ( .not.((elem == 6 .and. ion>=1 .and. ion<=6) .or.&
                                & (elem == 7 .and. ion>=1 .and. ion<=7) .or.&
                                & (elem == 8 .and. ion>=1 .and. ion<=8).or. &
                                & (elem == 14 .and. ion==1) .or. &
                                & (elem == 16 .and. ion==3) .or. &
                                & (elem == 18 .and. ion==5) .or. &
                                & (elem == 20 .and. ion==7) .or. &
                                & (elem == 26 .and. ion==12))) ) then
                              if (diRec(elem,ion) == 0..and. ion < 7 .and. lgElementOn(8)) &
                                   & diRec(elem,ion) = diRec(8,ion)
                              alphaTot(elem, ion) = alphaTot(elem, ion) + &
                                   & max(0.,diRec(elem, ion))

                           end if
                        end if
                     end do
                  end if
               end do

            else

               do elem = 1, nElements
                  if (lgElementOn(elem)) then
                     do ion = 1, min(nstages-1, elem)

                        if ( (elem == 14 .and. ion==1) .or. &
                             & (elem == 16 .and. ion==3) .or. &
                             & (elem == 18 .and. ion==5) .or. &
                             & (elem == 20 .and. ion==7) .or. &
                             & (elem == 26 .and. ion==12) .or.&
                             & (elem == 6 .and. ion>=1 .and. ion<=6) .or.&
                             & (elem == 7 .and. ion>=1 .and. ion<=7) .or.&
                             & (elem == 8 .and. ion>=1 .and. ion<=8) ) then

                           if (lgNahar) then
                              alphaTot(elem, ion) = nahar(elem,ion,TeUsed)
                           else
                              alphaTot(elem, ion) = radRecFit(elem, elem-ion+1)
                           end if
                        else
                           alphaTot(elem, ion) = radRecFit(elem, elem-ion+1)
                        end if
                     end do
                  end if
               end do


               ! calculate dielectronic recombination part
               do elem = 3, nElements
                  if (lgElementOn(elem)) then
                     do ion = 1, min(nstages-1, elem)
                        if ( .not.lgNahar .or. ( .not.((elem == 6 .and. ion>=1 .and. ion<=6) .or.&
                             & (elem == 7 .and. ion>=1 .and. ion<=7) .or.&
                             & (elem == 8 .and. ion>=1 .and. ion<=8).or. &
                             & (elem == 14 .and. ion==1) .or. &
                             & (elem == 16 .and. ion==3) .or. &
                             & (elem == 18 .and. ion==5) .or. &
                             & (elem == 20 .and. ion==7) .or. &
                             & (elem == 26 .and. ion==12))) ) then
                           if (diRec(elem,ion) == 0..and. ion < 7) &
                                & diRec(elem,ion) = diRec(8,ion)
                           alphaTot(elem, ion) = alphaTot(elem, ion) + &
                                & max(0.,diRec(elem, ion))
                        end if
                     end do
                  end if
               end do
            end if

       end subroutine calcAlpha

       function RRbadnell(elem, nel)
         implicit none

         real :: RRbadnell
         real :: bbad

         integer, intent(in) :: elem, nel


         bbad = RRbB(elem, nel)+RRbC(elem, nel)*exp(-RRbT2(elem, nel)/TeUsed)


         RRbadnell = RRbA(elem, nel) / ( &
              & ((TeUsed/RRbT0(elem, nel))**(0.5)) * &
              & ((1.+((TeUsed/RRbT0(elem, nel))**(0.5)))**(1.-bbad)) * &
              & ((1.+((TeUsed/RRbT1(elem, nel))**(0.5)))**(1.+bbad)) )


       end function RRbadnell

       function DRbadnell(elem, nel)
         implicit none

         real :: DRbadnell

         integer, intent(in) :: elem, nel
         integer :: idr

         DRbadnell = 0.
         do idr = 1, 9
            if ((DRbE(idr,elem,nel)/TeUsed) < 87.) &
                 & DRbadnell = DRbadnell+DRbC(idr,elem,nel)*exp(-DRbE(idr,elem,nel)/TeUsed)
         end do
         DRbadnell = DRbadnell*TeUsed**(-3./2.)

       end function DRbadnell

       ! total recombination coefficients from
       ! C and N Nahar 1997 ApJS 11, 339
       ! O Nahar 1999 ApJS 120, 131
       ! others: Nahar 2000, ApJS 126:537
       function nahar(iel,ist,tein)
         implicit none

         real, intent(in)    :: tein     ! e- temperature
         real                :: nahar    ! total recombination coeff [cm^3 s^-1]

         real :: lgTein      ! input log10Te
         real :: lgTeC(81)    ! array of log10(Te) for C ions
         real :: lgTeN(81)    ! array of log10(Te) for N ions
         real :: lgTeO(81)    ! array of log10(Te) for O ions
         real :: lgTe(71)    ! array of log10(Te)
         real :: c(81,6)     ! array of alpha total for O ions
         real :: n(81,7)     ! array of alpha total for O ions
         real :: o(81,8)     ! array of alpha total for O ions
         real :: si1(71)     ! array of alpha total for SiI
         real :: s3(71)      ! array of alpha total for SIII
         real :: ar5(71)     ! array of alpha total for ArV
         real :: ca7(71)     ! array of alpha total for CaVII
         real :: fe12(71)    ! array of alpha total for FeXII

         integer, intent(in) :: iel, ist ! atom number and ionic stage
         integer             :: ina

         lgTeC = (/1.0,1.1,1.2,1.3,1.4,1.5,&
              &1.6,1.7,1.8,1.9,2.0,2.1,&
              &2.2,2.3,2.4,2.5,2.6,2.7,&
              &2.8,2.9,3.0,3.1,3.2,3.3,&
              &3.4,3.5,3.6,3.7,3.8,3.9,&
              &4.0,4.1,4.2,4.3,4.4,4.5,&
              &4.6,4.7,4.8,4.9,5.0,5.1,&
              &5.2,5.3,5.4,5.5,5.6,5.7,&
              &5.8,5.9,6.0,6.1,6.2,6.3,&
              &6.4,6.5,6.6,6.7,6.8,6.9,&
              &7.0,7.1,7.2,7.3,7.4,7.5,&
              &7.6,7.7,7.8,7.9,8.0,8.1,&
              &8.2,8.3,8.4,8.5,8.6,8.7,&
              &8.8,8.9,9.0/)

         lgTeN = (/1.0,1.1,1.2,1.3,1.4,1.5,&
              &1.6,1.7,1.8,1.9,2.0,2.1,&
              &2.2,2.3,2.4,2.5,2.6,2.7,&
              &2.8,2.9,3.0,3.1,3.2,3.3,&
              &3.4,3.5,3.6,3.7,3.8,3.9,&
              &4.0,4.1,4.2,4.3,4.4,4.5,&
              &4.6,4.7,4.8,4.9,5.0,5.1,&
              &5.2,5.3,5.4,5.5,5.6,5.7,&
              &5.8,5.9,6.0,6.1,6.2,6.3,&
              &6.4,6.5,6.6,6.7,6.8,6.9,&
              &7.0,7.1,7.2,7.3,7.4,7.5,&
              &7.6,7.7,7.8,7.9,8.0,8.1,&
              &8.2,8.3,8.4,8.5,8.6,8.7,&
              &8.8,8.9,9.0/)

         lgTeO = (/1.0,1.1,1.2,1.3,1.4,&
              &1.5,1.6,1.7,1.8,1.9,2.0,&
              &2.1,2.2,2.3,2.4,2.5,2.6,&
              &2.7,2.8,2.9,3.0,3.1,3.2,&
              &3.3,3.4,3.5,3.6,3.7,3.8,&
              &3.9,4.0,4.1,4.2,4.3,4.4,&
              &4.5,4.6,4.7,4.8,4.9,5.0,&
              &5.1,5.2,5.3,5.4,5.5,5.6,&
              &5.7,5.8,5.9,6.0,6.1,6.2,&
              &6.3,6.4,6.5,6.6,6.7,6.8,&
              &6.9,7.0,7.1,7.2,7.3,7.4,&
              &7.5,7.6,7.7,7.8,7.9,8.0,&
              &8.1,8.2,8.3,8.4,8.5,8.6,&
              &8.7,8.8,8.9,9.0/)

         lgTe = (/1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,&
              & 2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,&
              & 3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.1,5.2,5.3,&
              & 5.4,5.5,5.6,5.7,5.8,5.9,6.0,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,&
              & 6.9,7.0,7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9,8.0/)

         c(:,1) = (/3.61E-11,3.15E-11,2.75E-11,2.40E-11,2.10E-11,1.83E-11,&
              &1.59E-11,1.39E-11,1.21E-11,1.05E-11,9.13E-12,7.94E-12,&
              &6.89E-12,5.99E-12,5.20E-12,4.51E-12,3.91E-12,3.39E-12,&
              &2.94E-12,2.55E-12,2.21E-12,1.91E-12,1.66E-12,1.44E-12,&
              &1.26E-12,1.10E-12,9.72E-13,8.67E-13,7.85E-13,7.25E-13,&
              &6.90E-13,7.05E-13,8.50E-13,1.30E-12,2.26E-12,3.84E-12,&
              &5.90E-12,8.02E-12,9.72E-12,1.07E-11,1.08E-11,1.01E-11,&
              &9.03E-12,7.69E-12,6.32E-12,5.05E-12,3.95E-12,3.04E-12,&
              &2.31E-12,1.75E-12,1.32E-12,9.93E-13,7.53E-13,5.74E-13,&
              &4.36E-13,3.35E-13,2.59E-13,2.00E-13,1.55E-13,1.20E-13,&
              &9.28E-14,7.15E-14,5.49E-14,4.20E-14,3.21E-14,2.44E-14,&
              &1.85E-14,1.39E-14,1.05E-14,7.89E-15,5.91E-15,4.42E-15,&
              &3.30E-15,2.46E-15,1.83E-15,1.36E-15,1.01E-15,7.53E-16,&
              &5.60E-16,4.17E-16,3.10E-16/)

         c(:,2) = (/1.63E-10,1.42E-10,1.25E-10,1.09E-10,9.56E-11,8.35E-11,&
              &7.30E-11,6.38E-11,5.56E-11,4.86E-11,4.24E-11,3.70E-11,&
              &3.22E-11,2.81E-11,2.45E-11,2.13E-11,1.85E-11,1.61E-11,&
              &1.41E-11,1.26E-11,1.16E-11,1.12E-11,1.12E-11,1.13E-11,&
              &1.12E-11,1.08E-11,1.00E-11,9.03E-12,7.91E-12,6.86E-12,&
              &6.02E-12,5.49E-12,5.40E-12,6.02E-12,7.86E-12,1.14E-11,&
              &1.63E-11,2.17E-11,2.62E-11,2.88E-11,2.91E-11,2.74E-11,&
              &2.45E-11,2.08E-11,1.71E-11,1.37E-11,1.06E-11,8.14E-12,&
              &6.13E-12,4.56E-12,3.36E-12,2.46E-12,1.79E-12,1.30E-12,&
              &9.36E-13,6.73E-13,4.79E-13,3.43E-13,2.45E-13,1.75E-13,&
              &1.25E-13,8.88E-14,6.33E-14,4.51E-14,3.21E-14,2.29E-14,&
              &1.63E-14,1.16E-14,8.28E-15,5.91E-15,4.22E-15,3.02E-15,&
              &2.16E-15,1.55E-15,1.12E-15,8.05E-16,5.83E-16,4.24E-16,&
              &3.10E-16,2.28E-16,1.69E-16/)

         c(:,3) = (/9.01E-10,8.12E-10,7.35E-10,6.68E-10,6.09E-10,5.57E-10,&
              &5.09E-10,4.62E-10,4.16E-10,3.70E-10,3.23E-10,2.78E-10,&
              &2.35E-10,1.96E-10,1.61E-10,1.32E-10,1.07E-10,8.69E-11,&
              &7.09E-11,5.86E-11,4.97E-11,4.37E-11,3.97E-11,3.69E-11,&
              &3.43E-11,3.15E-11,2.83E-11,2.48E-11,2.13E-11,1.80E-11,&
              &1.53E-11,1.38E-11,1.42E-11,1.73E-11,2.28E-11,2.91E-11,&
              &3.46E-11,3.77E-11,3.79E-11,3.58E-11,3.20E-11,2.74E-11,&
              &2.27E-11,1.83E-11,1.45E-11,1.13E-11,8.65E-12,6.56E-12,&
              &4.92E-12,3.66E-12,2.71E-12,1.99E-12,1.46E-12,1.06E-12,&
              &7.72E-13,5.60E-13,4.07E-13,2.96E-13,2.10E-13,1.52E-13,&
              &1.10E-13,7.97E-14,5.78E-14,4.19E-14,3.04E-14,2.21E-14,&
              &1.61E-14,1.17E-14,8.48E-15,6.17E-15,4.49E-15,3.27E-15,&
              &2.39E-15,1.74E-15,1.28E-15,9.36E-16,6.89E-16,5.09E-16,&
              &3.77E-16,2.81E-16,2.11E-16/)

         c(:,4) = (/6.03E-10,5.29E-10,4.65E-10,4.07E-10,3.56E-10,3.12E-10,&
              &2.72E-10,2.38E-10,2.07E-10,1.81E-10,1.58E-10,1.37E-10,&
              &1.20E-10,1.04E-10,9.03E-11,7.86E-11,6.82E-11,5.91E-11,&
              &5.13E-11,4.45E-11,3.85E-11,3.33E-11,2.89E-11,2.50E-11,&
              &2.16E-11,1.87E-11,1.61E-11,1.39E-11,1.20E-11,1.03E-11,&
              &8.83E-12,7.55E-12,6.44E-12,5.48E-12,4.65E-12,3.94E-12,&
              &3.33E-12,2.80E-12,2.35E-12,1.97E-12,1.65E-12,1.38E-12,&
              &1.14E-12,9.50E-13,7.86E-13,6.52E-13,5.59E-13,5.41E-13,&
              &6.49E-13,9.10E-13,1.29E-12,1.68E-12,1.98E-12,2.13E-12,&
              &2.12E-12,1.96E-12,1.73E-12,1.45E-12,1.18E-12,9.34E-13,&
              &7.22E-13,5.49E-13,4.12E-13,3.05E-13,2.24E-13,1.64E-13,&
              &1.18E-13,8.51E-14,6.13E-14,4.40E-14,3.14E-14,2.25E-14,&
              &1.61E-14,1.15E-14,8.17E-15,5.82E-15,4.15E-15,2.96E-15,&
              &2.11E-15,1.50E-15,1.07E-15/)

         c(:,5) = (/1.05E-09,9.25E-10,8.13E-10,7.14E-10,6.28E-10,5.51E-10,&
              &4.83E-10,4.23E-10,3.70E-10,3.24E-10,2.83E-10,2.47E-10,&
              &2.16E-10,1.88E-10,1.64E-10,1.43E-10,1.25E-10,1.09E-10,&
              &9.46E-11,8.22E-11,7.15E-11,6.21E-11,5.39E-11,4.67E-11,&
              &4.05E-11,3.51E-11,3.03E-11,2.63E-11,2.27E-11,1.96E-11,&
              &1.69E-11,1.45E-11,1.25E-11,1.07E-11,9.23E-12,7.92E-12,&
              &6.78E-12,5.80E-12,4.95E-12,4.22E-12,3.60E-12,3.06E-12,&
              &2.60E-12,2.20E-12,1.86E-12,1.57E-12,1.34E-12,1.19E-12,&
              &1.19E-12,1.42E-12,1.90E-12,2.53E-12,3.12E-12,3.53E-12,&
              &3.67E-12,3.54E-12,3.22E-12,2.79E-12,2.32E-12,1.87E-12,&
              &1.47E-12,1.13E-12,8.59E-13,6.43E-13,4.77E-13,3.51E-13,&
              &2.57E-13,1.87E-13,1.33E-13,9.64E-14,6.95E-14,4.98E-14,&
              &3.57E-14,2.56E-14,1.84E-14,1.32E-14,9.43E-15,6.76E-15,&
              &4.84E-15,3.48E-15,2.50E-15/)

         c(:,6) = (/1.60E-09,1.41E-09,1.24E-09,1.09E-09,9.60E-10,8.44E-10,&
              &7.41E-10,6.50E-10,5.71E-10,5.00E-10,4.38E-10,3.83E-10,&
              &3.35E-10,2.93E-10,2.56E-10,2.24E-10,1.95E-10,1.70E-10,&
              &1.48E-10,1.29E-10,1.13E-10,9.80E-11,8.53E-11,7.42E-11,&
              &6.44E-11,5.60E-11,4.86E-11,4.21E-11,3.65E-11,3.16E-11,&
              &2.73E-11,2.36E-11,2.04E-11,1.76E-11,1.52E-11,1.31E-11,&
              &1.13E-11,9.67E-12,8.30E-12,7.12E-12,6.09E-12,5.21E-12,&
              &4.45E-12,3.79E-12,3.23E-12,2.74E-12,2.32E-12,1.96E-12,&
              &1.66E-12,1.40E-12,1.17E-12,9.81E-13,8.19E-13,6.81E-13,&
              &5.64E-13,4.65E-13,3.82E-13,3.12E-13,2.54E-13,2.05E-13,&
              &1.65E-13,1.32E-13,1.05E-13,8.25E-14,6.47E-14,5.05E-14,&
              &3.92E-14,3.01E-14,2.31E-14,1.76E-14,1.33E-14,1.01E-14,&
              &7.56E-15,5.65E-15,4.21E-15,3.12E-15,2.30E-15,1.69E-15,&
              &1.24E-15,9.11E-16,6.63E-16/)

         n(:,1) = (/3.06E-11,2.66E-11,2.31E-11,2.01E-11,1.75E-11,1.52E-11,&
              &1.32E-11,1.14E-11,9.87E-12,8.54E-12,7.38E-12,6.38E-12,&
              &5.50E-12,4.75E-12,4.09E-12,3.54E-12,3.07E-12,2.69E-12,&
              &2.37E-12,2.11E-12,1.89E-12,1.69E-12,1.50E-12,1.34E-12,&
              &1.19E-12,1.05E-12,9.26E-13,8.10E-13,7.05E-13,6.11E-13,&
              &5.28E-13,4.62E-13,4.32E-13,4.93E-13,7.40E-13,1.27E-12,&
              &2.11E-12,3.16E-12,4.22E-12,5.07E-12,5.55E-12,5.60E-12,&
              &5.29E-12,4.72E-12,4.04E-12,3.33E-12,2.66E-12,2.08E-12,&
              &1.60E-12,1.21E-12,9.01E-13,6.68E-13,4.91E-13,3.60E-13,&
              &2.57E-13,1.86E-13,1.34E-13,9.64E-14,6.93E-14,4.98E-14,&
              &3.57E-14,2.56E-14,1.83E-14,1.32E-14,9.43E-15,6.77E-15,&
              &4.86E-15,3.50E-15,2.52E-15,1.82E-15,1.31E-15,9.54E-16,&
              &6.95E-16,5.09E-16,3.75E-16,2.78E-16,2.07E-16,1.56E-16,&
              &1.19E-16,9.14E-17,7.11E-17/)

         n(:,2) = (/1.61E-10,1.41E-10,1.23E-10,1.08E-10,9.43E-11,8.24E-11,&
              &7.20E-11,6.29E-11,5.49E-11,4.80E-11,4.18E-11,3.64E-11,&
              &3.17E-11,2.77E-11,2.41E-11,2.09E-11,1.82E-11,1.58E-11,&
              &1.38E-11,1.20E-11,1.04E-11,9.02E-12,7.84E-12,6.83E-12,&
              &5.98E-12,5.27E-12,4.69E-12,4.21E-12,3.79E-12,3.42E-12,&
              &3.08E-12,2.77E-12,2.54E-12,2.49E-12,2.90E-12,4.20E-12,&
              &6.73E-12,1.04E-11,1.45E-11,1.81E-11,2.05E-11,2.11E-11,&
              &2.03E-11,1.83E-11,1.57E-11,1.30E-11,1.04E-11,8.17E-12,&
              &6.26E-12,4.72E-12,3.52E-12,2.59E-12,1.90E-12,1.38E-12,&
              &1.00E-12,7.16E-13,5.14E-13,3.69E-13,2.64E-13,1.88E-13,&
              &1.34E-13,9.59E-14,6.84E-14,4.87E-14,3.47E-14,2.47E-14,&
              &1.76E-14,1.26E-14,8.97E-15,6.40E-15,4.57E-15,3.27E-15,&
              &2.34E-15,1.68E-15,1.21E-15,8.73E-16,6.33E-16,4.60E-16,&
              &3.36E-16,2.48E-16,1.83E-16/)

         n(:,3) = (/3.65E-10,3.20E-10,2.81E-10,2.46E-10,2.15E-10,1.88E-10,&
              &1.65E-10,1.44E-10,1.26E-10,1.10E-10,9.59E-11,8.36E-11,&
              &7.29E-11,6.36E-11,5.54E-11,4.85E-11,4.28E-11,3.84E-11,&
              &3.53E-11,3.32E-11,3.16E-11,3.00E-11,2.83E-11,2.65E-11,&
              &2.47E-11,2.29E-11,2.10E-11,1.91E-11,1.73E-11,1.58E-11,&
              &1.46E-11,1.36E-11,1.27E-11,1.19E-11,1.18E-11,1.32E-11,&
              &1.69E-11,2.28E-11,2.96E-11,3.55E-11,3.89E-11,3.94E-11,&
              &3.73E-11,3.34E-11,2.85E-11,2.35E-11,1.88E-11,1.47E-11,&
              &1.13E-11,8.53E-12,6.36E-12,4.70E-12,3.45E-12,2.52E-12,&
              &1.83E-12,1.32E-12,9.55E-13,6.88E-13,4.88E-13,3.50E-13,&
              &2.51E-13,1.79E-13,1.28E-13,9.13E-14,6.53E-14,4.66E-14,&
              &3.33E-14,2.38E-14,1.70E-14,1.22E-14,8.71E-15,6.24E-15,&
              &4.48E-15,3.22E-15,2.32E-15,1.67E-15,1.21E-15,8.82E-16,&
              &6.44E-16,4.73E-16,3.49E-16/)

         n(:,4) = (/6.50E-10,5.71E-10,5.02E-10,4.40E-10,3.86E-10,3.38E-10,&
              &2.96E-10,2.59E-10,2.26E-10,1.98E-10,1.73E-10,1.51E-10,&
              &1.31E-10,1.15E-10,9.98E-11,8.72E-11,7.63E-11,6.72E-11,&
              &5.99E-11,5.39E-11,4.90E-11,4.45E-11,4.05E-11,3.70E-11,&
              &3.40E-11,3.18E-11,3.00E-11,2.86E-11,2.70E-11,2.51E-11,&
              &2.29E-11,2.08E-11,1.96E-11,2.05E-11,2.45E-11,3.12E-11,&
              &3.88E-11,4.52E-11,4.86E-11,4.84E-11,4.54E-11,4.03E-11,&
              &3.44E-11,2.85E-11,2.30E-11,1.82E-11,1.42E-11,1.10E-11,&
              &8.36E-12,6.30E-12,4.72E-12,3.51E-12,2.59E-12,1.91E-12,&
              &1.40E-12,1.02E-12,7.47E-13,5.46E-13,4.00E-13,2.93E-13,&
              &2.09E-13,1.53E-13,1.12E-13,8.19E-14,6.00E-14,4.40E-14,&
              &3.23E-14,2.37E-14,1.74E-14,1.28E-14,9.40E-15,6.91E-15,&
              &5.08E-15,3.73E-15,2.75E-15,2.03E-15,1.50E-15,1.11E-15,&
              &8.21E-16,6.11E-16,4.56E-16/)

         n(:,5) = (/9.71E-10,8.54E-10,7.49E-10,6.58E-10,5.77E-10,5.06E-10,&
              &4.43E-10,3.87E-10,3.38E-10,2.96E-10,2.58E-10,2.25E-10,&
              &1.96E-10,1.71E-10,1.49E-10,1.29E-10,1.12E-10,9.76E-11,&
              &8.47E-11,7.34E-11,6.37E-11,5.52E-11,4.78E-11,4.13E-11,&
              &3.58E-11,3.09E-11,2.66E-11,2.30E-11,1.98E-11,1.70E-11,&
              &1.46E-11,1.25E-11,1.07E-11,9.14E-12,7.79E-12,6.62E-12,&
              &5.61E-12,4.76E-12,4.02E-12,3.38E-12,2.85E-12,2.39E-12,&
              &2.00E-12,1.67E-12,1.39E-12,1.15E-12,9.54E-13,7.97E-13,&
              &7.01E-13,7.11E-13,8.68E-13,1.17E-12,1.54E-12,1.87E-12,&
              &2.08E-12,2.13E-12,2.04E-12,1.83E-12,1.57E-12,1.30E-12,&
              &1.04E-12,8.14E-13,6.24E-13,4.72E-13,3.52E-13,2.61E-13,&
              &1.91E-13,1.40E-13,1.00E-13,7.27E-14,5.25E-14,3.77E-14,&
              &2.71E-14,1.95E-14,1.40E-14,1.00E-14,7.16E-15,5.13E-15,&
              &3.68E-15,2.64E-15,1.89E-15/)

         n(:,6) = (/1.53E-09,1.35E-09,1.19E-09,1.05E-09,9.21E-10,8.09E-10,&
              &7.10E-10,6.23E-10,5.46E-10,4.78E-10,4.18E-10,3.66E-10,&
              &3.20E-10,2.79E-10,2.44E-10,2.13E-10,1.85E-10,1.62E-10,&
              &1.41E-10,1.22E-10,1.06E-10,9.26E-11,8.04E-11,6.98E-11,&
              &6.06E-11,5.25E-11,4.55E-11,3.94E-11,3.40E-11,2.94E-11,&
              &2.54E-11,2.19E-11,1.89E-11,1.63E-11,1.40E-11,1.20E-11,&
              &1.03E-11,8.82E-12,7.55E-12,6.45E-12,5.50E-12,4.69E-12,&
              &3.99E-12,3.39E-12,2.87E-12,2.43E-12,2.05E-12,1.73E-12,&
              &1.48E-12,1.31E-12,1.29E-12,1.45E-12,1.79E-12,2.19E-12,&
              &2.53E-12,2.72E-12,2.72E-12,2.56E-12,2.27E-12,1.94E-12,&
              &1.60E-12,1.27E-12,9.95E-13,7.64E-13,5.78E-13,4.33E-13,&
              &3.21E-13,2.37E-13,1.74E-13,1.24E-13,9.00E-14,6.52E-14,&
              &4.69E-14,3.38E-14,2.44E-14,1.76E-14,1.27E-14,9.11E-15,&
              &6.56E-15,4.73E-15,3.41E-15/)

         n(:,7) = (/2.20E-09,1.94E-09,1.71E-09,1.51E-09,1.33E-09,1.17E-09,&
              &1.03E-09,9.02E-10,7.94E-10,6.95E-10,6.10E-10,5.34E-10,&
              &4.68E-10,4.09E-10,3.58E-10,3.13E-10,2.73E-10,2.38E-10,&
              &2.08E-10,1.81E-10,1.58E-10,1.38E-10,1.20E-10,1.04E-10,&
              &9.07E-11,7.89E-11,6.85E-11,5.95E-11,5.16E-11,4.46E-11,&
              &3.88E-11,3.35E-11,2.90E-11,2.50E-11,2.16E-11,1.86E-11,&
              &1.61E-11,1.38E-11,1.19E-11,1.02E-11,8.76E-12,7.49E-12,&
              &6.41E-12,5.48E-12,4.67E-12,3.98E-12,3.38E-12,2.87E-12,&
              &2.43E-12,2.05E-12,1.73E-12,1.45E-12,1.22E-12,1.02E-12,&
              &8.46E-13,7.02E-13,5.80E-13,4.77E-13,3.91E-13,3.17E-13,&
              &2.58E-13,2.07E-13,1.66E-13,1.32E-13,1.04E-13,8.21E-14,&
              &6.41E-14,4.99E-14,3.85E-14,2.95E-14,2.26E-14,1.71E-14,&
              &1.29E-14,9.73E-15,7.28E-15,5.43E-15,4.03E-15,2.98E-15,&
              &2.20E-15,1.61E-15,1.18E-15/)

         o(:,1) = (/2.91E-11,2.53E-11,2.20E-11,1.91E-11,1.65E-11,1.43E-11,&
              &1.24E-11,1.07E-11,9.28E-12,8.01E-12,6.91E-12,5.96E-12,&
              &5.12E-12,4.41E-12,3.79E-12,3.25E-12,2.78E-12,2.38E-12,&
              &2.03E-12,1.73E-12,1.48E-12,1.26E-12,1.07E-12,9.06E-13,&
              &7.68E-13,6.54E-13,5.58E-13,4.79E-13,4.14E-13,3.60E-13,&
              &3.14E-13,2.75E-13,2.46E-13,2.49E-13,3.45E-13,6.26E-13,&
              &1.16E-12,1.90E-12,2.72E-12,3.41E-12,3.84E-12,3.96E-12,&
              &3.80E-12,3.44E-12,2.96E-12,2.46E-12,1.98E-12,1.55E-12,&
              &1.20E-12,9.06E-13,6.79E-13,5.04E-13,3.71E-13,2.72E-13,&
              &1.98E-13,1.41E-13,1.02E-13,7.35E-14,5.27E-14,3.78E-14,&
              &2.72E-14,1.95E-14,1.40E-14,1.00E-14,7.19E-15,5.16E-15,&
              &3.71E-15,2.67E-15,1.92E-15,1.39E-15,1.00E-15,7.27E-16,&
              &5.29E-16,3.87E-16,2.85E-16,2.11E-16,1.57E-16,1.18E-16,&
              &8.97E-17,6.88E-17,5.33E-17/)

         o(:,2) = (/1.45E-10,1.27E-10,1.11E-10,9.71E-11,8.48E-11,7.39E-11,&
              &6.45E-11,5.62E-11,4.88E-11,4.26E-11,3.70E-11,3.21E-11,&
              &2.79E-11,2.42E-11,2.10E-11,1.82E-11,1.57E-11,1.36E-11,&
              &1.18E-11,1.03E-11,9.06E-12,8.00E-12,7.09E-12,6.28E-12,&
              &5.55E-12,4.89E-12,4.29E-12,3.75E-12,3.28E-12,2.87E-12,&
              &2.52E-12,2.21E-12,1.97E-12,1.83E-12,1.90E-12,2.42E-12,&
              &3.63E-12,5.62E-12,8.17E-12,1.08E-11,1.30E-11,1.43E-11,&
              &1.45E-11,1.38E-11,1.24E-11,1.07E-11,8.82E-12,7.10E-12,&
              &5.57E-12,4.29E-12,3.25E-12,2.43E-12,1.80E-12,1.33E-12,&
              &9.70E-13,7.07E-13,5.05E-13,3.64E-13,2.63E-13,1.89E-13,&
              &1.35E-13,9.69E-14,6.95E-14,4.98E-14,3.57E-14,2.55E-14,&
              &1.83E-14,1.31E-14,9.41E-15,6.76E-15,4.86E-15,3.50E-15,&
              &2.53E-15,1.83E-15,1.33E-15,9.71E-16,7.12E-16,5.25E-16,&
              &3.89E-16,2.91E-16,2.20E-16/)

         o(:,3) = (/4.49E-10,3.96E-10,3.48E-10,3.06E-10,2.69E-10,2.37E-10,&
              &2.08E-10,1.83E-10,1.60E-10,1.41E-10,1.23E-10,1.08E-10,&
              &9.49E-11,8.33E-11,7.30E-11,6.40E-11,5.60E-11,4.89E-11,&
              &4.24E-11,3.67E-11,3.15E-11,2.70E-11,2.31E-11,1.99E-11,&
              &1.74E-11,1.56E-11,1.46E-11,1.41E-11,1.39E-11,1.37E-11,&
              &1.33E-11,1.26E-11,1.17E-11,1.07E-11,9.98E-12,1.01E-11,&
              &1.18E-11,1.56E-11,2.13E-11,2.76E-11,3.32E-11,3.69E-11,&
              &3.83E-11,3.74E-11,3.47E-11,3.08E-11,2.63E-11,2.18E-11,&
              &1.76E-11,1.38E-11,1.07E-11,8.11E-12,6.08E-12,4.52E-12,&
              &3.32E-12,2.43E-12,1.77E-12,1.28E-12,9.24E-13,6.59E-13,&
              &4.73E-13,3.39E-13,2.42E-13,1.73E-13,1.23E-13,8.81E-14,&
              &6.29E-14,4.48E-14,3.20E-14,2.28E-14,1.63E-14,1.16E-14,&
              &8.30E-15,5.93E-15,4.25E-15,3.04E-15,2.18E-15,1.57E-15,&
              &1.13E-15,8.20E-16,5.96E-16/)

         o(:,4) = (/6.91E-10,6.08E-10,5.35E-10,4.69E-10,4.12E-10,3.61E-10,&
              &3.16E-10,2.77E-10,2.43E-10,2.13E-10,1.86E-10,1.63E-10,&
              &1.42E-10,1.25E-10,1.09E-10,9.66E-11,8.62E-11,7.79E-11,&
              &7.12E-11,6.55E-11,6.05E-11,5.64E-11,5.32E-11,5.06E-11,&
              &4.84E-11,4.62E-11,4.38E-11,4.11E-11,3.82E-11,3.52E-11,&
              &3.24E-11,2.98E-11,2.76E-11,2.58E-11,2.47E-11,2.59E-11,&
              &3.20E-11,4.50E-11,6.39E-11,8.43E-11,1.01E-10,1.10E-10,&
              &1.10E-10,1.03E-10,9.16E-11,7.77E-11,6.36E-11,5.06E-11,&
              &3.94E-11,3.01E-11,2.26E-11,1.68E-11,1.24E-11,9.08E-12,&
              &6.61E-12,4.79E-12,3.46E-12,2.50E-12,1.80E-12,1.30E-12,&
              &9.00E-13,6.43E-13,4.60E-13,3.26E-13,2.33E-13,1.66E-13,&
              &1.19E-13,8.48E-14,6.06E-14,4.34E-14,3.11E-14,2.23E-14,&
              &1.61E-14,1.16E-14,8.40E-15,6.11E-15,4.47E-15,3.29E-15,&
              &2.43E-15,1.82E-15,1.37E-15/)

         o(:,5) = (/1.02E-09,9.01E-10,7.91E-10,6.95E-10,6.11E-10,5.36E-10,&
              &4.69E-10,4.11E-10,3.59E-10,3.14E-10,2.74E-10,2.40E-10,&
              &2.09E-10,1.82E-10,1.59E-10,1.38E-10,1.20E-10,1.05E-10,&
              &9.11E-11,7.91E-11,6.87E-11,5.96E-11,5.16E-11,4.47E-11,&
              &3.87E-11,3.35E-11,2.91E-11,2.57E-11,2.30E-11,2.11E-11,&
              &1.96E-11,1.85E-11,1.77E-11,1.81E-11,2.06E-11,2.60E-11,&
              &3.36E-11,4.15E-11,4.74E-11,5.00E-11,4.91E-11,4.54E-11,&
              &3.99E-11,3.39E-11,2.79E-11,2.25E-11,1.80E-11,1.42E-11,&
              &1.10E-11,8.53E-12,6.55E-12,4.97E-12,3.76E-12,2.81E-12,&
              &2.09E-12,1.56E-12,1.15E-12,8.58E-13,6.33E-13,4.72E-13,&
              &3.53E-13,2.36E-13,1.75E-13,1.30E-13,9.50E-14,7.10E-14,&
              &5.31E-14,3.99E-14,3.00E-14,2.26E-14,1.71E-14,1.29E-14,&
              &9.76E-15,7.40E-15,5.62E-15,4.28E-15,3.27E-15,2.51E-15,&
              &1.93E-15,1.50E-15,1.17E-15/)

         o(:,6) = (/1.41E-09,1.25E-09,1.10E-09,9.63E-10,8.45E-10,7.42E-10,&
              &6.50E-10,5.69E-10,4.99E-10,4.36E-10,3.81E-10,3.32E-10,&
              &2.90E-10,2.53E-10,2.20E-10,1.92E-10,1.67E-10,1.45E-10,&
              &1.26E-10,1.09E-10,9.47E-11,8.21E-11,7.12E-11,6.16E-11,&
              &5.33E-11,4.61E-11,3.98E-11,3.44E-11,2.96E-11,2.55E-11,&
              &2.20E-11,1.89E-11,1.62E-11,1.38E-11,1.18E-11,1.01E-11,&
              &8.60E-12,7.30E-12,6.19E-12,5.24E-12,4.43E-12,3.73E-12,&
              &3.13E-12,2.63E-12,2.20E-12,1.83E-12,1.52E-12,1.26E-12,&
              &1.05E-12,9.07E-13,8.77E-13,1.02E-12,1.35E-12,1.79E-12,&
              &2.22E-12,2.52E-12,2.63E-12,2.55E-12,2.32E-12,2.01E-12,&
              &1.68E-12,1.36E-12,1.07E-12,8.22E-13,6.24E-13,4.68E-13,&
              &3.47E-13,2.56E-13,1.88E-13,1.32E-13,9.59E-14,6.93E-14,&
              &4.96E-14,3.57E-14,2.57E-14,1.85E-14,1.33E-14,9.55E-15,&
              &6.87E-15,4.95E-15,3.58E-15/)

         o(:,7) = (/2.11E-09,1.86E-09,1.64E-09,1.44E-09,1.27E-09,1.12E-09,&
              &9.80E-10,8.60E-10,7.55E-10,6.61E-10,5.79E-10,5.07E-10,&
              &4.44E-10,3.87E-10,3.39E-10,2.96E-10,2.58E-10,2.25E-10,&
              &1.96E-10,1.70E-10,1.49E-10,1.29E-10,1.12E-10,9.75E-11,&
              &8.47E-11,7.35E-11,6.37E-11,5.52E-11,4.78E-11,4.13E-11,&
              &3.57E-11,3.08E-11,2.66E-11,2.29E-11,1.97E-11,1.70E-11,&
              &1.46E-11,1.25E-11,1.07E-11,9.18E-12,7.85E-12,6.70E-12,&
              &5.71E-12,4.86E-12,4.13E-12,3.50E-12,2.96E-12,2.50E-12,&
              &2.11E-12,1.80E-12,1.57E-12,1.47E-12,1.52E-12,1.69E-12,&
              &1.92E-12,2.10E-12,2.18E-12,2.12E-12,1.96E-12,1.72E-12,&
              &1.46E-12,1.20E-12,9.52E-13,7.43E-13,5.72E-13,4.33E-13,&
              &3.26E-13,2.43E-13,1.80E-13,1.33E-13,9.42E-14,6.89E-14,&
              &5.03E-14,3.63E-14,2.64E-14,1.91E-14,1.39E-14,1.01E-14,&
              &7.31E-15,5.30E-15,3.85E-15/)

         o(:,8) = (/2.90E-09,2.56E-09,2.26E-09,2.00E-09,1.76E-09,1.55E-09,&
              &1.36E-09,1.20E-09,1.05E-09,9.26E-10,8.12E-10,7.12E-10,&
              &6.23E-10,5.46E-10,4.78E-10,4.18E-10,3.65E-10,3.19E-10,&
              &2.78E-10,2.43E-10,2.12E-10,1.85E-10,1.61E-10,1.40E-10,&
              &1.22E-10,1.06E-10,9.22E-11,8.01E-11,6.95E-11,6.03E-11,&
              &5.22E-11,4.53E-11,3.92E-11,3.39E-11,2.93E-11,2.53E-11,&
              &2.18E-11,1.88E-11,1.62E-11,1.39E-11,1.19E-11,1.03E-11,&
              &8.78E-12,7.51E-12,6.42E-12,5.48E-12,4.67E-12,3.97E-12,&
              &3.37E-12,2.85E-12,2.40E-12,2.03E-12,1.71E-12,1.43E-12,&
              &1.20E-12,9.96E-13,8.27E-13,6.84E-13,5.63E-13,4.61E-13,&
              &3.75E-13,3.05E-13,2.45E-13,1.96E-13,1.57E-13,1.24E-13,&
              &9.76E-14,7.63E-14,5.93E-14,4.59E-14,3.51E-14,2.70E-14,&
              &2.05E-14,1.55E-14,1.16E-14,8.72E-15,6.51E-15,4.83E-15,&
              &3.58E-15,2.64E-15,1.93E-15/)

         si1 = (/5.02E-11, 4.41E-11, 3.87E-11, 3.40E-11, 2.99E-11, 2.62E-11, &
              & 2.30E-11, 2.02E-11, 1.77E-11, 1.55E-11, 1.36E-11, 1.19E-11, &
              & 1.05E-11, 9.19E-12, 8.06E-12, 7.07E-12, 6.20E-12, 5.45E-12, &
              & 4.79E-12, 4.21E-12, 3.72E-12, 3.30E-12, 2.94E-12, 2.65E-12, &
              & 2.40E-12, 2.19E-12, 2.01E-12, 1.87E-12, 1.76E-12, 1.67E-12, &
              & 1.60E-12, 1.62E-12, 1.83E-12, 2.44E-12, 3.54E-12, 5.06E-12, &
              & 6.67E-12, 7.98E-12, 8.66E-12, 8.72E-12, 8.24E-12, 7.33E-12, &
              & 6.28E-12, 5.18E-12, 4.16E-12, 3.29E-12, 2.54E-12, 1.95E-12, &
              & 1.48E-12, 1.12E-12, 8.44E-13, 5.73E-13, 4.24E-13, 3.13E-13, &
              & 2.26E-13, 1.66E-13, 1.22E-13, 8.95E-14, 6.57E-14, 4.83E-14, &
              & 3.55E-14, 2.62E-14, 1.94E-14, 1.44E-14, 1.07E-14, 8.01E-15, &
              & 6.03E-15, 4.57E-15, 3.49E-15, 2.69E-15, 2.09E-15/)

         s3 = (/8.37E-10, 7.40E-10, 6.54E-10, 5.77E-10, 5.09E-10,&
              & 4.49E-10, 3.95E-10, 3.48E-10, 3.05E-10, 2.68E-10,&
              & 2.34E-10, 2.04E-10, 1.77E-10, 1.54E-10, 1.33E-10,&
              & 1.14E-10, 9.79E-11, 8.37E-11, 7.15E-11, 6.10E-11,&
              & 5.23E-11, 4.51E-11, 3.92E-11, 3.43E-11, 3.04E-11,&
              & 2.73E-11, 2.47E-11, 2.26E-11, 2.09E-11, 1.94E-11,&
              & 1.81E-11, 1.69E-11, 1.59E-11, 1.54E-11, 1.68E-11,&
              & 2.23E-11, 3.42E-11, 5.28E-11, 7.48E-11, 9.51E-11,&
              & 1.09E-10, 1.14E-10, 1.10E-10, 1.00E-10, 8.66E-11,&
              & 7.20E-11, 5.79E-11, 4.54E-11, 3.48E-11, 2.63E-11,&
              & 1.96E-11, 1.45E-11, 1.06E-11, 7.71E-12, 5.57E-12,&
              & 4.02E-12, 2.88E-12, 2.06E-12, 1.47E-12, 1.05E-12,&
              & 7.50E-13, 5.34E-13, 3.80E-13, 2.70E-13, 1.92E-13,&
              & 1.36E-13, 9.69E-14, 6.89E-14, 4.89E-14, 3.48E-14,&
              & 2.47E-14/)

         ar5 = (/1.37E-09, 1.21E-09, 1.07E-09, 9.40E-10, 8.29E-10,&
              & 7.31E-10, 6.43E-10, 5.66E-10, 4.97E-10, 4.38E-10,&
              & 3.85E-10, 3.38E-10, 2.97E-10, 2.62E-10, 2.31E-10,&
              & 2.05E-10, 1.83E-10, 1.65E-10, 1.50E-10, 1.39E-10,&
              & 1.32E-10, 1.28E-10, 1.28E-10, 1.28E-10, 1.27E-10,&
              & 1.25E-10, 1.20E-10, 1.14E-10, 1.07E-10, 9.94E-11,&
              & 9.18E-11, 8.42E-11, 7.68E-11, 6.97E-11, 6.37E-11,&
              & 5.99E-11, 6.13E-11, 7.11E-11, 9.08E-11, 1.18E-10,&
              & 1.46E-10, 1.66E-10, 1.76E-10, 1.73E-10, 1.60E-10,&
              & 1.40E-10, 1.18E-10, 9.58E-11, 7.58E-11, 5.87E-11,&
              & 4.46E-11, 3.34E-11, 2.48E-11, 1.82E-11, 1.33E-11,&
              & 9.64E-12, 6.96E-12, 5.02E-12, 3.59E-12, 2.57E-12,&
              & 1.84E-12, 1.31E-12, 9.35E-13, 6.67E-13, 4.75E-13,&
              & 3.38E-13, 2.41E-13, 1.72E-13, 1.22E-13, 8.70E-14,&
              & 6.19E-14/)

         ca7 = (/3.39E-09, 3.00E-09, 2.65E-09, 2.34E-09, 2.06E-09,&
              & 1.81E-09, 1.60E-09, 1.40E-09, 1.23E-09, 1.08E-09,&
              & 9.60E-10, 8.66E-10, 8.02E-10, 7.62E-10, 7.35E-10,&
              & 7.11E-10, 6.88E-10, 6.68E-10, 6.47E-10, 6.23E-10,&
              & 5.90E-10, 5.45E-10, 4.91E-10, 4.34E-10, 3.79E-10,&
              & 3.29E-10, 2.87E-10, 2.53E-10, 2.27E-10, 2.06E-10,&
              & 1.90E-10, 1.76E-10, 1.63E-10, 1.51E-10, 1.40E-10,&
              & 1.30E-10, 1.23E-10, 1.22E-10, 1.31E-10, 1.53E-10,&
              & 1.85E-10, 2.17E-10, 2.41E-10, 2.50E-10, 2.43E-10,&
              & 2.23E-10, 1.95E-10, 1.64E-10, 1.33E-10, 1.05E-10,&
              & 8.16E-11, 6.21E-11, 4.66E-11, 3.46E-11, 2.55E-11,&
              & 1.86E-11, 1.35E-11, 9.80E-12, 7.08E-12, 5.04E-12,&
              & 3.62E-12, 2.59E-12, 1.85E-12, 1.32E-12, 9.46E-13,&
              & 6.76E-13, 4.82E-13, 3.44E-13, 2.46E-13, 1.75E-13,&
              & 1.25E-13/)

         fe12 = (/1.84E-08, 1.63E-08, 1.44E-08, 1.27E-08, 1.12E-08,&
              & 9.81E-09, 8.59E-09, 7.50E-09, 6.53E-09, 5.67E-09,&
              & 4.89E-09, 4.21E-09, 3.60E-09, 3.07E-09, 2.63E-09,&
              & 2.26E-09, 1.97E-09, 1.77E-09, 1.61E-09, 1.49E-09,&
              & 1.38E-09, 1.28E-09, 1.19E-09, 1.12E-09, 1.08E-09,&
              & 1.06E-09, 1.06E-09, 1.06E-09, 1.05E-09, 1.02E-09,&
              & 9.68E-10, 8.92E-10, 8.00E-10, 7.02E-10, 6.07E-10,&
              & 5.22E-10, 4.49E-10, 3.87E-10, 3.37E-10, 2.94E-10,&
              & 2.59E-10, 2.30E-10, 2.06E-10, 1.84E-10, 1.64E-10,&
              & 1.43E-10, 1.23E-10, 1.03E-10, 8.39E-11, 6.72E-11,&
              & 5.28E-11, 4.09E-11, 3.12E-11, 2.36E-11, 1.77E-11,&
              & 1.32E-11, 9.80E-12, 7.23E-12, 5.34E-12, 3.94E-12,&
              & 2.90E-12, 2.04E-12, 1.50E-12, 1.09E-12, 7.93E-13,&
              & 5.79E-13, 4.23E-13, 3.08E-13, 2.25E-13, 1.65E-13,&
              & 1.20E-13/)

         if (Tein>0.) then
            lgTein = log10(tein)
         else
            print*, '! nahar: Tein <= 0.', tein
            stop
         end if

         ! C ions
         if (iel == 6 .and. ist>=1 .and. ist<=6) then

            if (lgTein <= lgTeC(1)) then
               nahar = c(1,ist)
            elseif (lgTein >= lgTeC(81)) then
               nahar = c(81,ist)
            else
               do ina = 2, 81
                  if (lgTein<lgTeC(ina)) then

                     nahar = c(ina-1,ist) + (c(ina,ist)-c(ina-1,ist))*&
                          & (lgTein-lgTeC(ina-1))/(lgTeC(ina)-lgTeC(ina-1))
                     exit

                  end if
               end do
            end if


         ! N ions
         elseif (iel == 7 .and. ist>=1 .and. ist<=7) then

            if (lgTein <= lgTeN(1)) then
               nahar = n(1,ist)
            elseif (lgTein >= lgTeN(81)) then
               nahar = n(81,ist)
            else
               do ina = 2, 81
                  if (lgTein<lgTeN(ina)) then
                     nahar = n(ina-1,ist) + (n(ina,ist)-n(ina-1,ist))*&
                          & (lgTein-lgTeN(ina-1))/(lgTeN(ina)-lgTeN(ina-1))
                     exit
                  end if
               end do
            end if

         ! O ions
         elseif (iel == 8 .and. ist>=1 .and. ist<=8) then

            if (lgTein <= lgTeO(1)) then
               nahar = o(1,ist)
            elseif (lgTein >= lgTeO(81)) then
               nahar = o(81,ist)
            else
               do ina = 2, 81
                  if (lgTein<lgTeO(ina)) then
                     nahar = o(ina-1,ist) + (o(ina,ist)-o(ina-1,ist))*&
                          & (lgTein-lgTeO(ina-1))/(lgTeO(ina)-lgTeO(ina-1))
                     exit
                  end if
               end do
            end if

         !Si I
         elseif (iel == 14 .and. ist ==1) then

            if (lgTein <= lgTe(1)) then
               nahar = si1(1)
            elseif (lgTein >= lgTe(71)) then
               nahar = si1(71)
            else
               do ina = 2, 71
                  if (lgTein<lgTe(ina)) then
                     nahar = si1(ina-1) + (si1(ina)-si1(ina-1))*&
                          & (lgTein-lgTe(ina-1))/(lgTe(ina)-lgTe(ina-1))
                     exit
                  end if
               end do
            end if

         ! SIII
         else if (iel ==16 .and. ist==3) then

            if (lgTein <= lgTe(1)) then
               nahar = s3(1)
            elseif (lgTein >= lgTe(71)) then
               nahar = s3(71)
            else
               do ina = 1, 71
                  if (lgTein<lgTe(ina)) then
                     nahar = s3(ina-1) + (s3(ina)-s3(ina-1))*&
                          & (lgTein-lgTe(ina-1))/(lgTe(ina)-lgTe(ina-1))
                     exit
                  end if
               end do
            end if

         ! Ar V
         else if (iel ==18 .and. ist==5) then

            if (lgTein <= lgTe(1)) then
               nahar = ar5(1)
            elseif (lgTein >= lgTe(71)) then
               nahar = ar5(71)
            else
               do ina = 1, 71
                  if (lgTein<lgTe(ina)) then
                     nahar = ar5(ina-1) + (ar5(ina)-ar5(ina-1))*&
                          & (lgTein-lgTe(ina-1))/(lgTe(ina)-lgTe(ina-1))
                     exit
                  end if
               end do
            end if

         ! Ca VII
         else if (iel ==20 .and. ist==7) then

            if (lgTein <= lgTe(1)) then
               nahar = ca7(1)
            elseif (lgTein >= lgTe(71)) then
               nahar = ca7(71)
            else
               do ina = 1, 71
                  if (lgTein<lgTe(ina)) then
                     nahar = ca7(ina-1) + (ca7(ina)-ca7(ina-1))*&
                          & (lgTein-lgTe(ina-1))/(lgTe(ina)-lgTe(ina-1))
                     exit
                  end if
               end do
            end if

         ! Fe XII
         else if (iel ==26 .and. ist==12) then

            if (lgTein <= lgTe(1)) then
               nahar = fe12(1)
            elseif (lgTein >= lgTe(71)) then
               nahar = fe12(71)
            else
               do ina = 1, 71
                  if (lgTein<lgTe(ina)) then
                     nahar = fe12(ina-1) + (fe12(ina)-fe12(ina-1))*&
                          & (lgTein-lgTe(ina-1))/(lgTe(ina)-lgTe(ina-1))
                     exit
                  end if
               end do
            end if

         else
            print*, '! nahar: unsopported ion', iel, ist
            stop
         end if


       end function nahar


       ! This subroutine calculates rates of radiative recombination for all ions
       ! of all elements from H through Zn by use of the following fits:
       ! H-like, He-like, Li-like, Na-like - Verner & Ferland, 1996, ApJS, 103, 467
       ! Other ions of C, N, O, Ne - Pequignot et al. 1991, A&A, 251, 680,
       ! refitted by Verner & Ferland formula to ensure correct asymptotes
       ! Fe XVII-XXIII - Arnaud & Raymond, 1992, ApJ, 398, 394
       ! Fe I-XV - refitted by Verner & Ferland formula to ensure correct asymptotes
       ! Other ions of Mg, Si, S, Ar, Ca, Fe, Ni -
       !                     - Shull & Van Steenberg, 1982, ApJS, 48, 95
       ! Other ions of Na, Al - Landini & Monsignori Fossi, 1990, A&AS, 82, 229
       ! Other ions of F, P, Cl, K, Ti, Cr, Mn, Co (excluding Ti I-II, Cr I-IV,
       ! Mn I-V, Co I)        - Landini & Monsignori Fossi, 1991, A&AS, 91, 183
       ! All other species    - interpolations of the power-law fits
       !  Input parameters: z - atomic number
       !                    n - number of electrons from 1 to z
       ! Output parameter:  radRecFit  - rate coefficient, cm^3 s^(-1)

       function radRecFit(z, n)
            implicit none

            integer, intent(in)             :: z     ! atomic weight of the element
            integer, intent(in)             :: n     ! number of electrons (from 1 to z)


            real                            :: radRecFit    ! rad rate coeff [cm^3/s]

            ! local variables
            integer                              :: elem    ! element counter
            integer                              :: i       ! counter
            integer                              :: ion     ! ion stage counter
            integer                              :: ios     ! I/O error status


            logical, save                        :: lgFirst = .true.
                                                            ! first time this is evaluated?

            real                                 :: tt      ! temp dep fact in interpolation

            real, dimension(2, nElements, nElements),&      ! coefficients for the
                 & save  :: rrec    ! calculation of the
            real, dimension(4, nElements, nElements),&      ! radiative rates
                 & save  :: rnew    !
            real, dimension(3, 4:13), save         :: fe    !

            ! if this is the first time this procedure is called
            ! read in radiative recombination coefficient file
            ! and the dielectronic recombination data
            if (lgFirst) then
                close(17)
                open (unit=17, file=PREFIX//'/share/mocassin/data/radrec.dat', status='old',position='rewind', &
                     & iostat = ios, action="read")

                do ion = 4, 30
                    if (ion /= 11) then
                        do elem = ion, 30
                            if ( (elem /= 26) .or. (ion >= 14) ) then
                                read(unit=17, fmt=*, iostat=ios) (rrec(i,elem,ion), i=1,2)
                            end if
                        end do
                    end if
                end do

                do ion = 1, 3
                    do elem = ion, 30
                        read(unit=17, fmt=*, iostat=ios) (rnew(i,elem,ion), i=1,4)
                    end do
                end do

                do elem = 11, 30
                   read(unit=17, fmt=*, iostat=ios) (rnew(i,elem,11), i=1,4)
                end do

                do ion = 4, 10
                    do elem = max(6,ion), 8
                        read(unit=17, fmt=*, iostat=ios) (rnew(i,elem,ion), i=1,4)
                    end do
                    read(unit=17, fmt=*, iostat=ios) (rnew(i,10,ion), i=1,4)
                end do

                do ion = 12, 26
                    read(unit=17, fmt=*, iostat=ios) (rnew(i,26,ion), i=1,4)
                end do

                ! close the files
                close(17)

                ! set data for Fe
                fe(1,:) = (/4.33e-10,3.91e-10,3.49e-10,3.16e-10,&
                     &2.96e-10,2.59e-10,2.24e-10,1.91e-10,1.68e-10,1.46e-10/)
                fe(2,:) = (/0.531,0.523,0.521,0.534, &
                     &0.557,0.567,0.579,0.601,0.602,0.597/)
                fe(3,:) = (/5.77e-02,6.15e-02,6.22e-02,6.02e-02,&
                     &5.79e-02,5.65e-02,5.49e-02,5.10e-02,5.07e-02,5.22e-02/)

                ! set lgFirst to .false.
                lgFirst = .false.
            end if

            ! check right element and number of electron reference
            if ( (z<1) .or. (z>30) ) then
                print*, "! radRecFit: insane atomic number", z
                stop
            end if
            if ( (n<1) .or. (n>z) ) then
                print*, "! radRecFit: insane number of electrons", n
                stop
            end if

            ! calculate the rates
            if ( (n<=3) .or. (n==11) .or. ((z>5) .and. (z<9)) .or. (z==10) .or.&
                 &((z==26) .and. (n>11)) ) then

                tt = sqrt(TeUsed/rnew(3,z,n))

                radRecFit = rnew(1,z,n) / (tt*(tt+1.)**(1.-rnew(2,z,n))*&
                     &(1.+sqrt(TeUsed/rnew(4,z,n)))**(1.+rnew(2,z,n)))
            else
               tt = TeUsed*1.e-4
               if ( (z==26) .and. (n<=13) ) then
                   radRecFit = fe(1,n)/tt**(fe(2,n)+fe(3,n)*log10(tt))
               else
                   radRecFit = rrec(1,z,n)/tt**rrec(2,z,n)
               end if
            end if

        end function radRecFit

        subroutine dielectronic(diRec)
            implicit none

            real, intent(out), dimension(nElements, nstages) &
                 &:: diRec   ! total recombination coeffs

            ! local variables

            integer              :: elem      ! atomic number
            integer              :: ion       ! ionization stage
            integer              :: n         ! number of electrons
            integer              :: ios       ! I/O error status
            integer              :: g         ! temperature flag

            real, dimension(nElements, nstages) :: aldroPequi! high T dielec rec coeff by A&P73

            real                 :: a,b,c,d,f ! fitting coefficients
            real                 :: t,t0,t1   ! t = TeUsed/10000., t0,t1 are fitting par
            real                 :: alpha

! Change this! CHECK WHAT THE REAL RANGE IS!!!!
            if (TeUsed>1000.) then
               t = TeUsed/10000.
            else
               t = 0.1
            end if

            alpha = 0.
            aldroPequi=0.

            close(18)
            open (unit=18, file=PREFIX//'/share/mocassin/data/dielectronic.dat', status='old',position='rewind', &
                 &iostat = ios, action="read")
            do i = 1, 10000
               read(unit=18, fmt=*, iostat=ios) elem, n, a, b, c, d, f, g
               if (ios < 0) exit ! end of file reached
               ion = elem + 1 - n

               if (ion <= nstages) then

                  if (g == 0) then
!print*, 'ee ', f/t, f, t
                     diRec(elem, ion) = (10.**(-12))*(a/t+b+c*t+d*t**2.)*t**(-3./2.)*exp(-f/t)
                  else if (g == 1) then
                     if (TeUsed < 20000.) then
!print*, 'ff ', f/t, f, t
                        diRec(elem, ion) = (10.**(-12))*(a/t+b+c*t+d*t**2.)*t**(-3./2.)*exp(-f/t)
                     end if
                  else if (g == 2) then
                     if (TeUsed >= 20000.) then
!print*, 'gg ', f/t, f, t
                        diRec(elem, ion) = (10.**(-12))*(a/t+b+c*t+d*t**2.)*t**(-3./2.)*exp(-f/t)
                     end if
                  end if
               end if

            end do
            close(18)

            ! calculate the high temperatures dielectronic recombination coeficients of
            ! Aldrovandi and Pequignot 1973
            t = TeUsed

            close(17)
            open (unit=17, file=PREFIX//'/share/mocassin/data/aldrovandi.dat', status='old',&
                 &position='rewind', iostat = ios, action="read")
            if (ios /= 0) then
               print*, "! dielectronic: can't open file ",PREFIX,"/share/mocassin/data/alrovandi.dat"
               stop
            end if

            do i = 1, 100000
               read(unit=17, fmt=*, iostat=ios) elem, n, a, b, t0, t1
               if (ios < 0) exit ! end of file reached
               ion = elem + 1 - n

!print*, 'gg ', t0/t, t1/t, t1, t0,t
               alpha = a*t**(-3./2.)*exp(-t0/t)*(1.+b*exp(-t1/t))

               ion = elem-n+1

               if (ion <= nstages) aldroPequi(elem, ion) = alpha

            end do

            close(17)


            if (TeUsed>60000.) then

               diRec = aldroPequi

!            else
!               do elem = 1, nElements
!                  do ion = 1, nstages
!                     if (diRec(elem,ion) == 0.) diRec(elem,ion) = aldroPequi(elem,ion)
!                  end do
!               end do
            end if


          end subroutine dielectronic


          subroutine getRadAcc()
            implicit none

!            real :: radius2

            integer :: i, ns, ai

!            radius2 = 1.e20*((grid%xAxis(xP)/1.e10)*(grid%xAxis(xP)/1.e10)+&
!                 & (grid%yAxis(yP)/1.e10)*(grid%yAxis(yP)/1.e10)+&
!                 & (grid%zAxis(zP)/1.e10)*(grid%zAxis(zP)/1.e10))

            grid%arad(cellP) = 0.

            do ai = 1, nSizes
               do nS = 1, nSpecies
                  if (grid%Tdust(nS,ai,cellP)<TdustSublime(nS)) then
                     do i = 1, nbins

                        grid%arad(cellP) = grid%arad(cellP)+Cpr(nS,ai,i)*grid%H(cellP,i)*&
                             & grainAbun(nS)*grainWeight(ai)

                     end do
                  end if
               end do
            end do


            ! arad(r) = Integral[kappa_pr(nu)*Lnu(r)dnu, {nu, 1, nbin}]/(4 pi r^2 c)
            ! the extra 4Pi is because Jste = 4Pi Jnu



            grid%arad(cellP) = grid%arad(cellP)/(4.*Pi*c)

          end subroutine getRadAcc



         subroutine getDustT()
            implicit none

            real                   :: dustAbsIntegral   ! dust absorption integral
            real                   :: dabs
            real                   :: resLineHeat       ! resonance line heating
            real, dimension(nbins) :: radField          ! radiation field

            integer :: nS, i, ai ! counters
            integer :: iT        ! pointer to dust temp in dust temp array

            ! radiation field at this location
            if (lgDebug) then
               radField = ((grid%Jste(cellP,:) + grid%Jdif(cellP,:)))/Pi
            else
               radField = grid%Jste(cellP,:)/Pi
            end if

            ! zero out dust temperature arrays
            grid%Tdust(:,:,cellP) = 0.


            ! calculate absorption integrals for each species
            do nS = 1, nSpecies

               do ai = 1, nSizes

                  dustAbsIntegral=0.
                  if (lgTraceHeating.and.taskid==0) dabs=0.

                  do i = 1, nbins
                     dustAbsIntegral = dustAbsIntegral+xSecArray(dustAbsXsecP(nS,ai)+i-1)*radField(i)
                     if (lgTraceHeating.and.taskid==0) then
                        dabs = dabs+xSecArray(dustAbsXsecP(nS,ai)+i-1)*radField(i)
                     end if
                  end do

                  if (lgGas .and. convPercent>=resLinesTransfer .and. (.not.lgResLinesFirst) .and. &
                       & (.not.nIterateMC==1) ) then
                     dustHeatingBudget(grid%abFileIndex(xp,yp,zp),0) = &
                          &dustHeatingBudget(grid%abFileIndex(xp,yp,zp),0)+&
                          & dustAbsIntegral*grainWeight(ai)*grainAbun(nS)*grid%Ndust(cellP)
                     dustHeatingBudget(0,0) = dustHeatingBudget(0,0)+&
                          & dustAbsIntegral*grainWeight(ai)*grainAbun(nS)*grid%Ndust(cellP)
                     resLineHeat = resLineHeating(ai,ns)
                     dustAbsIntegral = dustAbsIntegral+resLineHeat
                  end if

                  if (lgGas .and. lgPhotoelectric) then
                     dustAbsIntegral = dustAbsIntegral+gasDustColl_d(nS,ai)-&
                          & photoelHeat_d(nS,ai)
                  end if

                  if (lgTraceHeating.and.taskid==0) then
                     write(57,*) 'Dust Species: ', ns, ai
                     write(57,*) 'Abs of cont rad field  ', dabs
                     write(57,*) 'Res Lines Heating: ', reslineheat
                     if (lgGas .and. lgPhotoelectric) then
                        write(57,*) 'Grain potential ', grainPot(ns,ai)
                        write(57,*) 'Gas-dust collision heat', gasDustColl_d(nS,ai)
                        write(57,*) 'Photoelctric cooling: ', photoelHeat_d(nS,ai)
                     end if
                  end if

                  call locate(dustEmIntegral(nS,ai,:), dustAbsIntegral, iT)

                  if (iT<=0) then
                     print*, "getDustT: [warning] temperature of grain = 0. K!!!!"
                     print*, cellP
                     print*, nS, dustAbsIntegral
                     print*, dustEmIntegral(nS,ai,1)
                     !stop
                     iT=1
                     grid%Tdust(nS,ai,cellP) = 1.
                  else if (iT>=nTemps) then
                     grid%Tdust(nS,ai,cellP) = real(nTemps)
                  else
                     grid%Tdust(nS,ai,cellP) = real(iT) + &
                          & (dustAbsIntegral-dustEmIntegral(nS,ai,iT))*&
                          & (real(iT+1)-real(iT))/(dustEmIntegral(nS,ai,iT+1)-&
                          & dustEmIntegral(nS,ai,iT))
                  end if

                  if (lgTraceHeating.and.taskid==0) then
                     write(57,*) 'Radiative cooling: ', dustEmIntegral(ns,ai,iT)
                     write(57,*) 'Grain temperature: ', grid%Tdust(nS,ai,cellP), &
                          & " species ", grainLabel(nS), " size:", grainRadius(ai)

                  end if

                  if (lgTalk) &
                       & print*, "! getDustT: [talk] cell ", xP,yP,zP,"; Grain temperature: "&
                       &, grid%Tdust(nS,ai,cellP), " species ", grainLabel(nS), " size:", grainRadius(ai)

                  ! find weighted mean
                  grid%Tdust(nS,0,cellP) = grid%Tdust(nS,0,cellP)+&
                       & grid%Tdust(nS,ai,cellP)*grainWeight(ai)

               end do

               grid%Tdust(0,0,cellP) = grid%Tdust(0,0,cellP)+&
                    & grid%Tdust(nS,0,cellP)*grainAbun(nS)

            end do


          end subroutine getDustT

          function resLineHeating(sizeP,speciesP)
            implicit none

            real                :: resLineHeating ! dust heating due to res lines

            real                :: Gline   ! energy is the line [erg sec^-1]
            real                :: heat    ! sub calculation var

            integer, intent(in) :: sizeP   ! size pointer
            integer, intent(in) :: speciesP! species pointer

            integer             :: iL      ! line counter
            integer             :: imul    ! multiplet counter

            resLineHeating = 0.
            do iL = 1, nResLines

               Gline = 0.
               do imul = 1, resLine(iL)%nmul
                  if (resLine(iL)%elem==1) then

                     if ( resLine(iL)%ion == 1 .and. resLine(iL)%moclow(imul)==1 &
                          &.and. resLine(iL)%mochigh(imul)==2 ) then

                        ! fits to Storey and Hummer MNRAS 272(1995)41
!print*, 'f', TeUsed

                        Gline = Gline + 10**(-0.897*log10(TeUsed) + 5.05)*&
                             & grid%elemAbun(grid%abFileIndex(xP,yP,zP),1)&
                             & *1.e-25*ionDenUsed(elementXref(1),2)*&
                             & NeUsed

                     else

                        print*, "! resLineHeating: [warning] only dust heating &
                             &from H Lyman alpha and resonance lines from heavy"
                        print*, "elements is implemented in this version. &
                             &please contact author B. Ercolano -1-", &
                             &iL, resLine(iL)%elem, &
                             & resLine(iL)%ion

                     end if

                  else if (resLine(iL)%elem>1 .and. resLine(iL)%elem &
                       &== resLine(iL)%ion) then

                     Gline = Gline+hydrolines(resLine(iL)%elem,&
                          &resLine(iL)%mochigh(imul),resLine(iL)%moclow(imul))


                  else if (resLine(iL)%elem>2 .and. resLine(iL)%elem &
                       &/= resLine(iL)%ion)  then

                     if (forbiddenLines(resLine(iL)%elem,resLine(iL)%ion, &
                          &resline(iL)%moclow(imul), resline(iL)%mochigh(imul)) > 1.e-35) &
                          & Gline = Gline+forbiddenLines(resLine(iL)%elem,resLine(iL)%ion, &
                          &resline(iL)%moclow(imul), resline(iL)%mochigh(imul))

                  else

                     print*, "! resLineHeating: [warning] only dust heating from H ",&
                          &"Lyman alpha and resonance lines from heavy"
                     print*, "elements is implemented in this version. please contact &
                          &author B. Ercolano -2-", &
                          &iL, resLine(iL)%elem, &
                             & resLine(iL)%ion

                  end if


               end do

               Gline=Gline/Pi

               heat = Gline*grid%Hden(cellP)* &
                    & (1.-grid%fEscapeResPhotons(cellP, iL))*&
                    & xSecArray(dustAbsXSecP(speciesP,sizeP)+resLine(iL)%nuP-1)/&
                    & (grid%Ndust(cellP)*&
                    & absOpacSpecies(speciesP,resLine(iL)%nuP))

               ! Harrington Monk and Clegg 1988 (section 3.2)
               resLineHeating = resLineHeating + heat

               dustHeatingBudget(grid%abFileIndex(xP,yP,zP),iL) = &
                    & dustHeatingBudget(grid%abFileIndex(xP,yP,zP),iL) + &
                    & heat*grainWeight(sizeP)*grainAbun(speciesP)*grid%Ndust(cellP)
               dustHeatingBudget(0,iL) = dustHeatingBudget(0,iL) + &
                    & heat*grainWeight(sizeP)*grainAbun(speciesP)*grid%Ndust(cellP)


            end do

          end function resLineHeating

          ! calculate collisional ionisation and energy rates
          ! Voronov 1997
          ! adapted from CLOUDY
          subroutine collisionalIonisation
            implicit none

            real    :: te, u
            integer :: elem,ion, elec

            collIon = 0.

            te=TeUsed*8.617385e-05
            do elem = 1, nElements
               if (lgElementOn(elem)) then
                  do ion = 1, min(nstages-1, elem -1)
                     elec=elem-ion+1
                     if(elem < 1.or. elem>30)return
                     if(elec<1   .or. elec>elem)return
                     u=CF(1,elem,elec)/te
                     if(u>80.0) return

!print*, 'hh ', u
                     collIon(1,elem,ion)=CF(3,elem,elec)*(1.0+CF(2,elem,elec)*&
                          & sqrt(u))/(CF(4,elem,elec)+u)*&
                          & u**CF(5,elem,elec)*exp(-u)

                     collIon(2,elem,ion)=cionTherm(elem,ion)*collIon(1,elem,ion)*kBoltzmann

                  end do
               end if
            end do
!            collIon = collIon*NeUsed

          end subroutine collisionalIonisation

          subroutine photoionisation()
            implicit none

            real, dimension(nElements, nstages) :: &
                 & deltaE_k      ! deltaE/k [K]

            real :: radField(nbins)
            real, dimension(10,nElements, nStages) :: &
                 & nPhotoSte , &    ! # of stellar photoionizations
                 & nPhotoDif , &    ! # of diffuse photoionizations
                 & heatIonSte,&   ! # of stellar photoionizations
                 & heatIonDif      ! # of diffuse photoionizations

!           real, dimension(nelements, nstages) :: heleion

            real :: revRate, phXSec, heatef, xe

            real :: t4

            integer :: elem,ion,nshell,outshell,g0,g1,nElec, ipNuP, highNuP

            nPhotoSte =1.e-20
            nPhotoDif =1.e-20
            heatIonSte=0.
            heatIonDif=0.

            xe = 0.
            do elem = 1, nelements
               if (lgElementOn(elem)) then
                  xe = xe + ionDenUsed(elementXref(elem),1)*elemAbunUsed(elem)*HdenUsed
               end if
            end do
            xe = NeUsed / xe

            if (xe > 0.95) then
               heatef = 1.
            else
               ! Xu and McCray 1991, Ap.J. 375, 190.
               !  everything goes to asymptote that is not present in Shull and
               !  Van Steenberg
               xe = max(xe, 1.e-4)
               heatef = 0.9971 * (1. - (1.-xe**0.2663)**1.3163)
            end if

            ! calculate the number of stellar and diffuse photoionizations for each species
            ! NOTE: no need to time by the frequency bin width (widflx(j)) because JSte and JDif
            !       were calculated for each individual bin

            do elem = 1, nElements  ! begin element loop
               do ion = 1, min(elem, nStages-1) ! begin ion loop
                  if(.not.lgElementOn(elem)) exit

                  do nShell = 1, nShells(elem, ion)
                     if (elem > 2) then

                        IPnuP = elementP(elem, ion, nShell, 1)
                        highNuP = elementP(elem, ion, nShell, 2)

                     else if (elem == 1) then ! HI

                        IPNuP   = HlevNuP(1)
                        highNuP = nbins

                     else if ( (elem == 2) .and. (ion == 1) ) then ! HeI

                        IPNuP   = HeIlevNuP(1)
                        highNuP = nbins

                     else if ( (elem == 2) .and. (ion == 2) ) then ! HeII

                        IPNuP   = HeIIlevNuP(1)
                        highNuP = nbins

                     end if

                     if (elem < 3) then


                        do j = IPnuP, highNuP
                           if ( elem == 1 ) then ! HI
                              phXSec  = xSecArray(j-IPnuP+1+HlevXSecP(1)-1)
                           else if ( (elem == 2) .and. (ion == 1) ) then ! HeI
                              phXSec  = xSecArray(j-IPnuP+1+HeISingXSecP(1)-1)
                           else if ( (elem == 2) .and. (ion == 2) ) then ! HeII
                              phXSec  = xSecArray(j-IPnuP+1+HeIIXSecP(1)-1)
                           end if

                           if ((phXSec < 1.e-35) ) then
                              phXSec = 0.
                           end if

                           if ( grid%JSte(cellP,j) >0.) then
                              nPhotoSte(nshell,elem,ion) = nPhotoSte(nshell,elem,ion) + &
                                   & grid%JSte(cellP,j)*phXSec/(hcRyd*nuArray(j))

                              heatIonSte(nshell,elem,ion) = heatIonSte(nshell,elem,ion) + &
                                   & phXSec*grid%Jste(cellP,j)*&
                                   & (nuArray(j)-nuArray(IPNuP)) / (nuArray(j))
!print*, elem, ion, nshell, nPhotoSte(nshell,elem,ion) , phXSec, grid%Jste(cellP,j), (hcRyd*nuArray(j))
                           end if
                           if ( lgDebug) then
                              if (grid%JDif(cellP,j) >0.) then
                                 nPhotoDif(nshell,elem,ion) = nPhotoDif(nshell,elem,ion) + &
                                      & grid%JDif(cellP,j)*phXSec/(hcRyd*nuArray(j))
                                 heatIonDif(nshell,elem,ion) = heatIonDif(nshell,elem,ion) + &
                                      & phXSec*grid%JDif(cellP,j)*&
                                      & (nuArray(j)-nuArray(IPNuP)) / (nuArray(j))
                              end if
                           end if

                        end do

                     else

                        radField = grid%JSte(cellP,:)
                        call getPhotoRates(nPhotoSte(nshell,elem,ion), heatIonSte(nshell,elem,ion),&
                             & elem, ion, nshell, IPnuP, highNuP, auger(elem,ion,nshell,1),radField, heatef)

                        if (lgDebug) then
                           radField = grid%JDif(cellP,:)
                           call getPhotoRates(nPhotoDif(nshell,elem,ion), heatIonDif(nshell,elem,ion),&
                                & elem, ion, nshell, IPnuP, highNuP, auger(elem,ion,nshell,1),radField, heatef)
                        end if

                     end if


                  end do ! end shell loop


               end do ! end ion loop
            end do ! end element loop

            photoIon(1, 1:10,1:30,1:nstages) = nPhotoSte

            if (lgDebug) photoIon(1, 1:10,1:30,1:nstages) = nPhotoSte+nPhotoDif
            photoIon(2, 1:10,1:30,1:nstages) = heatIonSte
            if (lgDebug) photoIon(2, 1:10,1:30,1:nstages) = heatIonSte+heatIonDif

            ! the set of charge exchange coeffs is not complete; the following might need
            ! to be changed when a more complete set is available

            chex          = 0.
            deltaE_k      = 0.

            chex(2,1,:)  = (/7.47e-6, 2.06, 9.93,-3.89/)! He0
            chex(2,2,:)  = (/1.e-5  , 0.  , 0.  , 0./)  ! He+
            chex(3,2,:)  = (/1.26   , 0.96,3.02 ,-0.65/)! Li+
            chex(3,3,:)  = (/1.e-5  , 0.  , 0.  , 0. /) ! Li+2
            chex(4,2,:)  = (/1.e-5  , 0.  , 0.  , 0. /) ! Be+
            chex(4,3,:)  = (/1.e-5  , 0. , 0.  , 0. /) ! Be+2
            chex(4,4,:)  = (/5.17   , 0.82, -.69, -1.12 /)! Be+3
            chex(5,2,:)  = (/2.e-2  , 0.  , 0.  , 0. /) ! B+
            chex(5,3,:)  = (/1.e-5  , 0.  , 0.  , 0. /) ! B+2
            chex(5,4,:)  = (/5.27e-1, 0.76,-0.63,-1.17/)! B+3
            chex(6,1,:)  = (/1.76e-9, 8.33, 4278.78, -6.41/)! C0
            chex(6,2,:)  = (/1.67e-4, 2.79, 304.72, -4.07/)! C+
            chex(6,3,:)  = (/3.25   , 0.21, 0.19, -3.29/)! C+2
            chex(6,4,:)  = (/332.46 ,-0.11,-0.995,-1.58e-3/)! C+3
            chex(7,1,:)  = (/1.01e-3,-0.29,-0.92, -8.38/)! N0
            chex(7,2,:)  = (/3.05e-1, 0.60, 2.65, -0.93/)! N+
            chex(7,3,:)  = (/4.54   , 0.57,-0.65, -0.89/)! N2+
            chex(7,4,:)  = (/3.28   , 0.52,-0.52, -0.19/)! N3+
            chex(8,1,:)  = (/1.04   , 3.15e-2, -0.61, -9.73/)! O0
            chex(8,2,:)  = (/1.04   , 0.27, 2.02, -5.92/)! O+
            chex(8,3,:)  = (/3.98   , 0.26, 0.56, -2.62/)! O2+
            chex(8,4,:)  = (/2.52e-1, 0.63, 2.08, -4.16/)! O3+
            chex(9,2,:)  = (/1.e-5  , 0. , 0.  , 0./) ! F+
            chex(9,3,:)  = (/9.86   , 0.29,-0.21,-1.15/) ! F+2
            chex(9,4,:)  = (/7.15e-1, 1.21,-0.70,-0.85/) ! F3+
            chex(10,2,:) = (/1.e-5  , 0.  , 0.  , 0.  /) ! Ne+
            chex(10,3,:) = (/14.73  , 4.52e-2, -0.84, -0.31 /) ! Ne+2
            chex(10,4,:) = (/6.47   , 0.54 , 3.59 , -5.22 /) ! Ne+3
            chex(11,2,:) = (/1.e-5  , 0.   , 0.   , 0. /) ! Na+
            chex(11,3,:) = (/1.33   , 1.15 , 1.20 , -0.32 /)! Na+2
            chex(11,4,:) = (/1.01e-1, 1.34 , 10.05, -6.41 /)! Na+3
            chex(12,2,:) = (/8.58e-5, 2.49e-3, 2.93e-2, -4.33 /)! Mg+
            chex(12,3,:) = (/6.49   , 0.53 , 2.82, -7.63 /) ! Mg+2
            chex(12,4,:) = (/6.36   , 0.55 , 3.86, -5.19 /) ! Mg+3
            chex(13,2,:) = (/1.e-5  , 0.   , 0.  , 0./) ! Al+
            chex(13,3,:) = (/7.11e-5, 4.12 , 1.72e4, -22.24/)! Al+2
            chex(13,4,:) = (/7.52e-1, 0.77 , 6.24, -5.67/) ! Al+3
            chex(14,2,:) = (/1.23   , 0.24 , 3.17, 4.18e-3/) ! Si+
            chex(14,3,:) = (/4.900e-1, -8.74e-2, -0.36, -0.79/)! Si+2
            chex(14,4,:) = (/7.58   , 0.37 , 1.06, -4.09/)! Si+3
            chex(16,1,:) = (/3.82e-7, 11.10, 2.57e4, -8.22/)! S0
            chex(16,2,:) = (/1.e-5  , 0.   , 0.   ,0. /)! S+
            chex(16,3,:) = (/2.29   , 4.02e-2, 1.59, -6.06/)! S+2
            chex(16,4,:) = (/6.44   , 0.13 , 2.69 , -5.69/)! S+3
            chex(18,2,:) = (/1.e-5  , 0.   , 0.    , 0./) ! Ar+
            chex(18,3,:) = (/4.57   , 0.27 , -0.18 , -1.57/)! Ar+2
            chex(18,4,:) = (/6.37   , 2.12 , 10.21 , -6.22/)! Ar+3
            chex(18,3,:) = (/3.17e-2, 2.12 , 12.06 , -0.40/)! Ca+2
            chex(18,4,:) = (/2.68   , 0.69 , -0.68 , -4.47/)! Ca+3
            chex(26,2,:) = (/1.26   , 7.72e-2, -0.41, -7.31/)! Fe+
            chex(26,3,:) = (/3.42   , 0.51 , -2.06 , -8.99/)! Fe+2.


            deltaE_k(7,1) = 10863.
            deltaE_k(8,1) = 2205.

            chex(:,:,1) = chex(:,:,1)*1.e-9

            t4 = TeUsed/10000.

            do elem = 1, nElements
               do ion = 1, min(elem, nstages-1)
                  if (.not.lgElementOn(elem)) exit

!print*, 'iil ', chex(elem,ion,4)*t4, t4
                  if (TeUsed < 6000. .or. TeUsed>5.e4) then
                     chex(elem,ion,1) = 0.
                  else
                     chex(elem,ion,1) = chex(elem,ion,1)*(t4**chex(elem,ion,2))*&
                          & (1.+chex(elem,ion,3)*exp(chex(elem,ion,4)*t4))
                  end if

                  if (chex(elem,ion,1) < 0. ) chex(elem,ion,1) = 0.

                  ! find the number of electron in this ion
                  nElec = elem - ion +1

                  ! find the stat weights
                  call getOuterShell(elem, nElec, outShell, g0, g1)

                  ! calculate the reverse charge exchange rate (only if deltaE_k > 1.)
                  if ( deltaE_k(elem,ion) > 1.) then

!print*, 'lli ', deltaE_k(elem,ion)/TeUsed, TeUsed

                     revRate = chex(elem,ion,1) * &
                          & (1.-grid%ionDen(cellP, &
                          & elementXref(1),1))*grid%Hden(cellP)*&
                          & 2. * exp(-deltaE_k(elem,ion)/TeUsed)/(real(g0)/real(g1))
                  else
                     revRate = 0.
                  end if
!print*, elem, ion, outshell, photoIon(1, outShell,elem, ion), revRate

                  photoIon(1, outShell,elem, ion) = photoIon(1, outShell,elem, ion)+revRate


               end do
            end do



          end subroutine photoionisation

          subroutine compton
            implicit none

            real ::  s1, s2, radfield, comptonRecoilHeatH, comptonRecoilHeatHe, energy

            integer :: i


            comptonCool = 0.
            comptonHeat = 0.
            comptonRecoilHeat = 0.
            comptonRecoilHeatH = 0.
            comptonRecoilHeatHe = 0.
            comptonRecoilIonH = 0.

            s1=0.
            s2=0.
            do i=1,nbins

               radfield = grid%Jste(cellP,i)
               if (lgDebug) radfield = radfield+grid%Jdif(cellP,i)

               s1 = s1+radField
               s2 = s2+radField/nuArray(i)


               if (i > cRecoilHP) then
                  ! bound electron scattering of >2.3 kev photons if neutral
                  ! recoil starts at 194 Ryd = 2.6
                  ! heating modified for suprathermal secondaries below; nuArray(i)^2
                  energy = 2.66e-5 * nuArray(i)*nuArray(i) - 1.
                  comptonRecoilHeatH = comptonRecoilHeatH+xSecRecoil(i)*radfield*energy

                  ! direct H ionisation
                  comptonRecoilIonH = comptonRecoilIonH+xSecRecoil(i)*radfield

               end if

               if (i > cRecoilHeP) then
                  ! bound electron scattering of >2.3 kev photons if neutral
                  ! recoil starts at 194 Ryd = 2.6
                  ! heating modified for suprathermal secondaries below; nuArray(i)^2
                  energy = 2.66e-5 * nuArray(i)*nuArray(i) - 1.8
                  comptonRecoilHeatHe = comptonRecoilHeatHe+xSecRecoil(i)*radfield*energy

               end if

            end do


            s1 = s1*Ryd2erg
!            ebar = s1/s2


            ! Compton cooling and heating
            comptonCool = (6.65e-25*6.75e-10*TeUsed*s2)/(grid%Hden(cellP))
            comptonHeat = (6.65e-25*s1/8.184e-7)/(grid%Hden(cellP))

            comptonRecoilHeat=(comptonRecoilHeatH*grid%ionDen(cellP,elementXref(1),1)*&
                        & grid%elemAbun(grid%abFileIndex(xP,yP,zP),1)+comptonRecoilHeatHe* &
                        & (grid%ionDen(cellP,elementXref(2),1)*&
                        & grid%elemAbun(grid%abFileIndex(xP,yP,zP),2))*2.)*Ryd2erg

          end subroutine compton

          subroutine compton2
            implicit none

            real ::  heatin,hin, radfield, comptonRecoilHeatH, comptonRecoilHeatHe, energy

            integer :: i


            comptonCool = 0.
            comptonHeat = 0.
            comptonRecoilHeat = 0.
            comptonRecoilHeatH = 0.
            comptonRecoilHeatHe = 0.
            comptonRecoilIonH = 0.

            heatin = 0.

            do i=1,nbins

               radfield = grid%Jste(cellP,i)
               if (lgDebug) radfield = radfield+grid%Jdif(cellP,i)

               ! Compton cooling
               ! comXSecC isTarter expression times nuArray(i)*3.858E-25
               ! 6.338E-6 is k in inf mass Rydbergs, still needs factor of TE
                comptonCool = comptonCool + radfield * comXSecC(i) *&
                    & (4.d0*6.338e-6*1e-15)

               ! Compton heating
               ! comXSecC is  Tarter expression times nuArray(I)**2 * 3.858E-25
               ! CMHEAT is just spontaneous, HEATIN is just induced
                comptonHeat= comptonHeat+radfield * comXSecH(i)*1.e-15

               ! induced Compton heating
               hin = comXSecH(i) *(radfield/(hcRyd*nuArray(i)))*1.e-15

               comptonHeat= comptonHeat+hin

               if (i > cRecoilHP) then
                  ! bound electron scattering of >2.3 kev photons if neutral
                  ! recoil starts at 194 Ryd = 2.6
                  ! heating modified for suprathermal secondaries below; nuArray(i)^2
                  energy = 2.66e-5 * nuArray(i)*nuArray(i) - 1.
                  comptonRecoilHeatH = comptonRecoilHeatH+xSecRecoil(i)*radfield*energy

                  ! direct H ionisation
                  comptonRecoilIonH = comptonRecoilIonH+xSecRecoil(i)*radfield

               end if

               if (i > cRecoilHeP) then
                  ! bound electron scattering of >2.3 kev photons if neutral
                  ! recoil starts at 194 Ryd = 2.6
                  ! heating modified for suprathermal secondaries below; nuArray(i)^2
                  energy = 2.66e-5 * nuArray(i)*nuArray(i) - 1.8
                  comptonRecoilHeatHe = comptonRecoilHeatHe+xSecRecoil(i)*radfield*energy

               end if



            end do
            comptonRecoilHeat=(comptonRecoilHeatH*grid%ionDen(cellP,elementXref(1),1)*&
                        & grid%elemAbun(grid%abFileIndex(xP,yP,zP),1)+comptonRecoilHeatHe* &
                        & (grid%ionDen(cellP,elementXref(2),1)*&
                        & grid%elemAbun(grid%abFileIndex(xP,yP,zP),2))*2.)*Ryd2erg

            comptonCool = comptonCool*TeUsed


          end subroutine compton2

          subroutine getPhotoRates(gamma,heat,elem,ion,nshell,nu1P,nu2P,&
               & yield1,radField, heatef)
            implicit none

            real, intent(in)    :: radField(nbins), yield1, heatef
            real, intent(out)   :: gamma, heat

            real                :: gammaHigh, heatHigh, augerE

            integer, intent(in) :: elem, ion, nu1P, nShell, nu2P

            integer             :: iup, ilow, IPnuP, highNuP, xSecP

            heat  = 0.
            gamma = 0.

            iup = min(nu2P, nbins)

            augerE = hcRyd*nuArray(nu1P)*yield1

            if(.not.lgElementOn(elem)) return

            IPnuP = nu1P
            highNuP = nu2P
            xSecP = elementP(elem, ion, nShell, 3)
!print*, elem,ion,nShell , nuArray(nu1p), nuArray(nu2p)
            ! low energies - no secondary ionisation
            do i = nu1P+1, min(iup, secIonP-1)
               if (radField(i) > 0.) then

                  phXsec = xSecArray(i+xSecP-IPnuP)
                  if (phXsec < 1.e-35) then
                     phXsec = 0.
                  end if
                  gamma = gamma + radField(i)*phXsec/(hcRyd*nuArray(i))
                  heat = heat + radField(i)*phXsec
!print*, nuArray(i), xSecArray(i+xSecP-IPnuP), gamma, heat
               end if

            end do

            ! heating is now corrected for the work function
            heat = heat -gamma*augerE
!print*, heat, augerE

            ! threshold
            phXsec = xSecArray(xSecP)
            gamma = gamma + radField(IPnuP)*phXsec/(hcRyd*nuArray(IPnuP))
            heat = heat + radField(IPnuP)*phXsec*(hcRyd*nuArray(IPnuP)-augerE)/(hcRyd*nuArray(IPnuP))

!print*, gamma, heat, radField(IPnuP)*phXsec*(hcRyd*nuArray(IPnuP)-augerE)/(hcRyd*nuArray(IPnuP))
            ! higher energy end - secondary ionisation may occur -
            heatHigh  = 0.
            gammaHigh = 0.
!print*, 'secion', xsecArray(secionP)
            ilow = max(IPnuP+1, secIonP)
            do i = ilow, iup

               phXsec = xSecArray(i+xSecP-IPnuP)
               gammaHigh = gammaHigh + radField(IPnuP)*phXsec/(hcRyd*nuArray(i))
               HeatHigh = heatHigh +  radField(IPnuP)*phXsec
!print*, nuArray(i), phXsec, gammaHigh, HeatHigh
            end do

            gamma = gamma+gammaHigh
            heatHigh = heatHigh - gammaHigh*augerE
!print*, gamma, gammaHigh*augerE, heatEf
            heat = heat + heatHigh*heatEf
!print*, heat
          end subroutine getPhotoRates


          subroutine RecLinesEmission()
            implicit none

            ! local variables
            integer                    :: itemp, iden, idenp, izp
            integer                    :: ios         ! I/O error status
            integer                    :: i, denint   ! counters
            integer                    :: ilow,&      ! pointer to lower level
                 &iup                                 ! pointer to upper level
            integer                    :: elUp
            integer                    :: ix,iy,iz
            real                       :: Afit,Bfit,zFit ! fit coeffs
            real                       :: Hbeta(30)    ! Hbeta emission
            real                       :: T4,T4z,Nez  ! TeUsed/10000., scaled, Ne scaled
            real                       :: log10NeZ    ! log10 NeUsed scaled
            real                       :: log10TeZ    ! log10(6^2*Te/Z^2)
            real                       :: x1, x2
            real                       :: dens(14)    !
            real                       :: hydrolinesloc(1:30,1:13,2:15,1:8)
            real                       :: x,y1,y2

            T4 = TeUsed / 10000.
            ix = xP
            iy = yP
            iz = zP

            ! find the nearest temp bin
            itemp = 1
            do i = 1, 12
               if (T4>hydroLinesTemps(i)) itemp = i
            end do
            if (itemp<12) then
               if (T4 > ( hydroLinesTemps(itemp+1)+hydroLinesTemps(itemp))/2.) &
                    & itemp = itemp+1
            end if

            ! read in HI recombination lines [e-25 ergs*cm^3/s] normalised to Hbeta
            ! (subset from Storey and Hummer MNRAS 272(1995)41)
            ! for Z>8 apply hydrogenic T-scaling
            if (T4 > hydroLinesTemps(12)) then
               elUp = 2
            else
               elUp = 8
            end if

            do izp = 1, elUp
               if (lgElementOn(izp) .and. nstages > izp) then

                  close(94)
                  open(unit = 94,  action="read", file = PREFIX//"/share/mocassin/"//hydroLinesFile(izp,itemp), &
                       status = "old", position = "rewind", iostat=ios)
                  if (ios /= 0) then
                     print*, "! RecLinesEmission: can't open file: ",PREFIX,"/share/mocassin/", hydroLinesFile(izp,itemp)
                     stop
                  end if
                  dens = 0.
                  hydrolinesloc = 0.
                  do iden = 1, 100
                     read(unit=94, fmt=*, iostat=ios) dens(iden)
                     if (ios<0) exit

                     do iup = 15, 2, -1
                        read(94, fmt=*) (hydrolinesloc(izp, iden,iup, ilow), ilow = 1, min(8, iup-1))
                     end do
                  end do
                  close(94)

                  ! look at density
                  idenp = 1
                  do iden = 1, 13
                     if (NeUsed>dens(iden) .and. dens(iden)> 0.) idenp = iden
                  end do

                  if ((dens(13)>0. .and. idenp<13) .or. (dens(13)==0. &
                       & .and. idenp<9)) then
                     ! interpolate
                     do iup = 2, 15
                        do ilow = 1, min(8,iup-1)
                           hydroLines(izp,iup,ilow) = &
                                &(hydrolinesloc(izp, idenp,iup, ilow)+&
                                & (hydrolinesloc(izp, idenp+1,iup, ilow)-&
                                & hydrolinesloc(izp, idenp,iup, ilow))*&
                                & (NeUsed-dens(idenp))/(dens(idenp+1)-dens(idenp)))
                        end do
                     end do

                  else
                     hydroLines(izp,:,:) = hydrolinesloc(izp, idenp,:,:)
                  end if

                  ! calculate Hbeta
                  ! fits to Storey and Hummer MNRAS 272(1995)41!
                  zfit = (log10Ne-HbACoeff(izp,2))/HbACoeff(izp,3)
                  Afit = HbACoeff(izp,1)*exp(-zfit*zfit/2.)+HbACoeff(izp,4)
                  zfit = (log10Ne-HbBCoeff(izp,2))/HbBCoeff(izp,3)
                  Bfit = HbBCoeff(izp,1)*exp(-zfit*zfit/2.)+HbBCoeff(izp,4)

                  if ((izp == 1) .or. &
                       &(izp >1 .and. T4<=hydrolinesTemps(12))) then
                     Hbeta(izp) = 10.**(Afit+Bfit*log10Te)
                  else
                     if (izp > 2) then
                        print*, "recLinesEmission [emission] : izp insanity"
                        stop
                     end if

                     if (idenp<13) then
                        ! interpolate in density
                        y1 = rbEdge(izp,3,idenp)+(rbEdge(izp,3,idenp)-rbEdge(izp,3,idenp+1))*&
                             & (NeUsed-dens(idenp))/(dens(idenp+1)-dens(idenp))
                        y2 = rbEdge(izp,3,9+idenp)+(rbEdge(izp,3,9+idenp)-rbEdge(izp,3,9+idenp+1))*&
                             & (NeUsed-dens(idenp))/(dens(idenp+1)-dens(idenp))
                     else
                        y1 = rbEdge(izp,3,idenp)
                        y2 = rbEdge(izp,3,13+idenp)
                     end if
                     x1 =  rbEdge(izp,1,1)
                     x2 =  rbEdge(izp,1,14)

                     x  = log10(T4*1.e4)
                     Hbeta(izp) = y2+(y1-y2)*(x-x2)/(x1-x2)
                     Hbeta(izp) = 10.**Hbeta(izp)
                  end if

                  Hbeta(izp) = Hbeta(izp)*NeUsed*&
                       &ionDenUsed(elementXref(izp),izp+1)*&
                       &grid%elemAbun(grid%abFileIndex(ix,iy,iz),izp)
                  hydroLines(izp,:,:) = hydroLines(izp,:,:)*Hbeta(izp)
               end if


            end do

            do izp = elUp+1, nElements

               if (lgElementOn(izp) .and. nstages > izp) then

                  ! scale T4 to HeII (Z=2)
                  T4z = T4*4./real(izp*izp)
                  log10TeZ = log10(TeUsed*4./real(izp*izp))

                  Nez = NeUsed*4./real(izp*izp)
                  log10Nez = log10(NeUsed*4./real(izp*izp))

                  ! find the nearest temp bin
                  itemp = 1
                  do i = 1, 12
                     if (T4z>hydroLinesTemps(i)) itemp = i
                  end do
                  if (itemp<12) then
                     if (T4z > ( hydroLinesTemps(itemp+1)+&
                          & hydroLinesTemps(itemp))/2.) &
                          & itemp = itemp+1
                  end if

                  close(94)
                  ! this is the HeII case A
                  open(unit = 94,  action="read", file = PREFIX//"/share/mocassin/"//hydroLinesFile(9,itemp), &
                       status = "old", position = "rewind", iostat=ios)
                  if (ios /= 0) then
                     print*, "! RecLinesEmission: can't open file: ",PREFIX,"/share/mocassin/", hydroLinesFile(9,itemp)
                     stop
                  end if
                  dens = 0.
                  hydrolinesloc = 0.
                  do iden = 1, 100
                     read(unit=94, fmt=*, iostat=ios) dens(iden)
                     if (ios<0) exit

                     do iup = 15, 2, -1
                        read(94, fmt=*) (hydrolinesloc(izp, iden,iup, ilow), &
                             &ilow = 1, min(8, iup-1))
                     end do
                  end do
                  close(94)

                  ! look at density
                  idenp = 1
                  do iden = 1, 13
                     if (Nez>dens(iden) .and. dens(iden)> 0.) idenp = iden
                  end do

                  if ((dens(13)>0. .and. idenp<13) .or. (dens(13)==0. &
                       &.and. idenp<9)) then
                     ! interpolate
                     do iup = 2, 15
                        do ilow = 1, min(8,iup-1)
                           hydroLines(izp,iup,ilow) = (hydrolinesloc(izp, &
                                &idenp,iup, ilow)+&
                                & (hydrolinesloc(izp, idenp+1,iup, ilow)-&
                                & hydrolinesloc(izp, idenp,iup, ilow))*&
                                & (Nez-dens(idenp))/(dens(idenp+1)-dens(idenp)))
                        end do
                     end do

                  else
                     hydroLines(izp,:,:) = hydrolinesloc(izp, idenp,:,:)
                  end if

                  ! calculate Hbeta
                  ! fits to Storey and Hummer MNRAS 272(1995)41!
                  zfit = (log10Nez-HbACoeff(9,2))/HbACoeff(9,3)
                  Afit = HbACoeff(9,1)*exp(-zfit*zfit/2.)+HbACoeff(9,4)
                  zfit = (log10Nez-HbBCoeff(9,2))/HbBCoeff(9,3)
                  Bfit = HbBCoeff(9,1)*exp(-zfit*zfit/2.)+HbBCoeff(9,4)

                  if (T4Z<=hydrolinesTemps(12)) then
                     Hbeta(izp) = (real(izp**3)/8.)*10.**(Afit+Bfit*log10TeZ)
                  else
                     if (idenp<9) then
                        ! interpolate in density
                        y1 = r2aEdge(3,idenp)+(r2aEdge(3,idenp)-r2aEdge(3,idenp+1))*&
                             & (NeUsed-dens(idenp))/(dens(idenp+1)-dens(idenp))
                        y2 = r2aEdge(3,9+idenp)+(r2aEdge(3,9+idenp)-r2aEdge(3,9+idenp+1))*&
                             & (NeUsed-dens(idenp))/(dens(idenp+1)-dens(idenp))
                     else
                        y1 = r2aEdge(3,idenp)
                        y2 = r2aEdge(3,9+idenp)
                     end if
                     x1 =  r2aEdge(1,1)
                     x2 =  r2aEdge(1,10)

                     x  = log10(T4Z*1.e4)
                     Hbeta(izp) = y2+(y1-y2)*(x-x2)/(x1-x2)
                     Hbeta(izp) = (real(izp**3)/8.)*(10.**Hbeta(izp))
                  end if


                  Hbeta(izp) = (real(izp**3)/8.)*10.**(Afit+Bfit*log10TeZ)

                  Hbeta(izp) = Hbeta(izp)*NeUsed*ionDenUsed(elementXref(izp),izp+1)*&
                       &grid%elemAbun(grid%abFileIndex(ix,iy,iz),izp)

                  hydroLines(izp,:,:) = hydroLines(izp,:,:)*Hbeta(izp)

               end if

            end do

            ! add contribution of Lyman alpha
            ! fits to Storey and Hummer MNRAS 272(1995)41
!            Lalpha = 10**(-0.897*log10Te + 5.05)
            !print*, Lalpha, hydroLines(1,2,1)
!            hydroLines(1,15, 8) = hydroLines(1,15, 8) + &
!                 & grid%elemAbun(grid%abFileIndex(ix,iy,iz),1)*&
!                 & ionDenUsed(elementXref(1),2)*&
!                 & NeUsed*Lalpha

            ! now do HeI

            ! atomic data limits
            if (T4 < 0.5) T4 = 0.5
            if (T4 < 2.0) T4 = 2.0

            if (grid%Ne(grid%active(ix,iy,iz)) <= 100.) then
               denint=0
            elseif (grid%Ne(grid%active(ix,iy,iz)) > 100. .and. &
                 &grid%Ne(grid%active(ix,iy,iz)) <= 1.e4) then
               denint=1
            elseif (grid%Ne(grid%active(ix,iy,iz)) > 1.e4 .and. &
                 &grid%Ne(grid%active(ix,iy,iz)) <= 1.e6) then
               denint=2
            elseif (grid%Ne(grid%active(ix,iy,iz)) > 1.e6) then
               denint=3
            end if

            ! data from Benjamin, Skillman and Smits ApJ514(1999)307 [e-25 ergs*cm^3/s]
            if (denint>0.and.denint<3) then
               do i = 1, 34
                  x1=HeIrecLineCoeff(i,denint,1)*&
                       &(T4**(HeIrecLineCoeff(i,denint,2)))*&
                       &exp(HeIrecLineCoeff(i,denint,3)/T4)
                  x2=HeIrecLineCoeff(i,denint+1,1)*&
                       &(T4**(HeIrecLineCoeff(i,denint+1,2)))*&
                       &exp(HeIrecLineCoeff(i,denint+1,3)/T4)

                  HeIRecLines(i) = x1+((x2-x1)*(NeUsed-100.**denint)/&
                       &(100.**(denint+1)-100.**(denint)))

               end do
            elseif(denint==0) then
               do i = 1, 34
                  HeIRecLines(i) = HeIrecLineCoeff(i,1,1)*&
                       & (T4**(HeIrecLineCoeff(i,1,2)))*&
                       &exp(HeIrecLineCoeff(i,1,3)/T4)
               end do
            elseif(denint==3) then
               do i = 1, 34
                  HeIRecLines(i) = HeIrecLineCoeff(i,3,1)*&
                       &(T4**(HeIrecLineCoeff(i,3,2)))*&
                       &exp(HeIrecLineCoeff(i,3,3)/T4)
               end do
            end if
            HeIRecLines=HeIRecLines*NeUsed*&
                 &grid%elemAbun(grid%abFileIndex(ix,iy,iz),2)*&
                 &ionDenUsed(elementXref(2),2)


          end subroutine RecLinesEmission

        end subroutine updateCell

end module update_mod
