! Copyright (C) 2008 Barbara Ercolano
!
! Version 3.04
module photon_mod
    use common_mod
    use constants_mod
    use continuum_mod
    use grid_mod
    use interpolation_mod
    use pathIntegration_mod
    use vector_mod

    ! common variables

    real :: Qphot = 0.
    real :: QfluoBins(14) = 0.
    real :: QFluo(14) = 0.
    real, dimension(14) :: fluorescenceTot


    type(vector), parameter :: origin=vector(0.,0.,0.)  ! origin of the cartesian grid axes

    integer     , parameter :: safeLim = 10000          ! safety limit for the loops
    integer                        :: totalEscaped

    contains

    subroutine energyPacketDriver(iStar, n, grid, gpLoc, cellLoc)
        implicit none

        real, pointer                  :: JnuStart(:,:,:) ! impinging flux (plane parallel only)

        type(vector)                   :: posDiff     ! initial position vector for diff ext
        type(vector)                   :: posVector   ! initial position vector for dust emi
        type(vector)                   :: positionIn  ! position of packet absorption

        integer, intent(in)            :: n           ! number of energy packets
        integer, intent(in)            :: iStar       ! central star index

        integer                        :: gPIn        !
        integer                        :: igp         ! 1= mother 2 =sub
        integer                        :: i,j,k       ! counters
        integer                        :: iCell       ! cell counter
        integer                        :: igrid,ix,iy,iz ! location indeces
        integer                        :: ierr        ! allocation error status
        integer                        :: iPhot       ! counter
        integer                        :: seedSize    ! pseudo random number generator seed
        integer, dimension(2)          :: inX,inY,inZ ! initial position indeces
        integer, pointer               :: seed(:)     ! seed array
        integer                        :: msec        ! millisecs of the sec
        integer                        :: dt(8)       ! date and time values
        integer                        :: trapped
        integer                        :: reRun


        integer, intent(inout) &
             & :: gpLoc                               ! local grid (only used for extra diffuse sources)
        integer, intent(inout) &
             & :: cellLoc(3)                          ! local cell (only used for extra diffuse sources)

        type(grid_type), dimension(:), intent(inout) :: grid        ! the grid(s)

        character(len=7)               :: chTypeD     ! character type for driver
        character(len=7)               :: chTypeIn    ! character type




        if (iStar == 0) then
!           deltaE(0) = grid(gpLoc)%LdiffuseLoc(grid(gpLoc)%active(cellLoc(1),cellLoc(2),cellLoc(3)))/NphotonsDiffuseLoc
           deltaEUsed = grid(gpLoc)%LdiffuseLoc(grid(gpLoc)%active(cellLoc(1),cellLoc(2),cellLoc(3)))/NphotonsDiffuseLoc

        else if (iStar >= 1) then

           deltaEUsed = Lstar(istar)/nPhotons(istar)

        end if

        if (deltaEUsed< 0.) then
           print*, "! energyPacketDriver: deltaEUsed < 0! - 1"
           stop
        end if

        call date_and_time(values=dt)
        msec=dt(8)

        call random_seed(seedSize)

        allocate(seed(1:seedSize), stat= ierr)
        if (ierr /= 0) then
            print*, "energyPacketDriver: can't allocate array memory: seed"
            stop
        end if

        seed = 0

        call random_seed(get = seed)

        seed = seed + msec + taskid

        call random_seed(put = seed)

        if (associated(seed)) deallocate(seed)

        Qphot = 0.
        QfluoBins = 0.
        QFluo = 0.

        fluorescenceTot = 0.


        if (taskid ==0 .and. lgPlaneIonization .and. nIterateMC==1) then
           allocate(JnuStart(grid(1)%nx, grid(1)%nz, nbins))
           JnuStart = 0.
        end if

        trapped = 0

        do iPhot = 1, n

           if (iStar>=1) then

              chTypeD = "stellar"

              if (starIndeces(iStar,4)==1) then
                 igp = 1
              else if (starIndeces(iStar,4)>1) then
                 igp = 2
              else
                 print*, "! energyPacketDriver: insane grid&
                      & pointer, starIndeces(iStar,4)", starIndeces(iStar,4)
                 stop
              end if

              inX=-1
              inY=-1
              inZ=-1

              inX(igp) = starIndeces(iStar,1)
              inY(igp) = starIndeces(iStar,2)
              inZ(igp) = starIndeces(iStar,3)


              chTypeIn = chTypeD
              positionIn = starPosition(iStar)
              gPIn = starIndeces(iStar,4)
              reRun = 0
              do i = 1, recursionLimit
                 call energyPacketRun(chTypeIn, positionIn, inX, inY, &
                      &inZ, gPIn, reRun)
                 if (rerun == 0) exit
!                 elseif (reRun == 1) then
!                    call energyPacketRun(chTypeIn, positionIn, inX, inY, &
!                      &inZ, gPIn, reRun)
!                 else
!                    print*, "! energyPacketDriver: insane reRun value ", reRun
!                    stop
!                 end if
              end do
              if (i>=recursionLimit) trapped = trapped+1

           else if (iStar==0) then

              chTypeD = "diffExt"
              inX=-1
              inY=-1
              inZ=-1

              if (gpLoc==1) then
                 igp = 1
              else if (gpLoc>1) then
                 igp = 2
              else
                 print*, "! energyPacketDriver: insane grid pointer"
                 stop
              end if

              inX(igp) = cellLoc(1)
              inY(igp) = cellLoc(2)
              inZ(igp) = cellLoc(3)

              posDiff%x = grid(gpLoc)%xAxis(cellLoc(1))
              posDiff%y = grid(gpLoc)%yAxis(cellLoc(2))
              posDiff%z = grid(gpLoc)%zAxis(cellLoc(2))

!              countRecursive = 0


              chTypeIn = chTypeD
              positionIn = posDiff
              gPIn = gPLoc
              reRun = 0
              do i = 1, recursionLimit
                 call energyPacketRun(chTypeIn, positionIn, inX, inY, &
                      &inZ, gPIn, reRun)
                 if (rerun == 0) exit
!                 elseif (reRun == 1) then
!                    call energyPacketRun(chTypeIn, positionIn, inX, inY, &
!                      &inZ, gPIn, reRun)
!                 else
!                    print*, "! energyPacketDriver: insane reRun value ", reRun
!                    stop
!                 end if
              end do
              if (i>=recursionLimit) trapped = trapped+1


!              call energyPacketRun(chType=chTypeD, position=posDiff, xp=inX, &
!                   & yp=inY, zp=inZ, gp=gpLoc)

           else

              print*, '! energyPacketDriver: insanity in iStar value'
              stop

           end if


        end do



        if (iStar>0) then
           print*, 'Star: ', iStar
           print*, 'Qphot = ', Qphot
           if (lgFluorescence) then
              print*, 'Label, Threshold [KeV], Qfluo , QfluoBins [erg/sec], nuFluo'

              print*, fluoLabArray(1), fluoThreshArray(1)*13.6/1.e3, Qfluo(1), &
                   &  QfluoBins(1)/(widflx(FeKaColdP)*ryd2erg), nuArray(FeKaColdP)
              print*, fluoLabArray(2), fluoThreshArray(2)*13.6/1.e3, Qfluo(2), &
                   & QfluoBins(2)/(widflx(FeL1P)*ryd2erg), nuArray(FeL1P)
              print*, fluoLabArray(3), fluoThreshArray(3)*13.6/1.e3, Qfluo(3), &
                   & QfluoBins(3)/(widflx(FeL2P)*ryd2erg), nuArray(FeL2P)
              print*, fluoLabArray(4), fluoThreshArray(4)*13.6/1.e3, Qfluo(4), &
                   & QfluoBins(4)/(widflx(CKaP)*ryd2erg), nuArray(CKaP)
              print*, fluoLabArray(5), fluoThreshArray(5)*13.6/1.e3, Qfluo(5), &
                   & QfluoBins(5)/(widflx(NKaP)*ryd2erg), nuArray(NKaP)
              print*, fluoLabArray(6), fluoThreshArray(6)*13.6/1.e3, Qfluo(6), &
                   & QfluoBins(6)/(widflx(OKaP)*ryd2erg), nuArray(OKaP)
              print*, fluoLabArray(7), fluoThreshArray(7)*13.6/1.e3, Qfluo(7), &
                   & QfluoBins(7)/(widflx(NeKaP)*ryd2erg), nuArray(NeKaP)
              print*, fluoLabArray(8), fluoThreshArray(8)*13.6/1.e3, Qfluo(8), &
                   & QfluoBins(8)/(widflx(MgKaP)*ryd2erg), nuArray(MgKaP)
              print*, fluoLabArray(9), fluoThreshArray(9)*13.6/1.e3, Qfluo(9), &
                   & QfluoBins(9)/(widflx(AlKaP)*ryd2erg), nuArray(AlKaP)
              print*, fluoLabArray(10),fluoThreshArray(10)*13.6/1.e3, Qfluo(10),&
                   & QfluoBins(10)/(widflx(SiKaP)*ryd2erg), nuArray(SiKaP)
              print*, fluoLabArray(11),fluoThreshArray(11)*13.6/1.e3, Qfluo(11), &
                   &QfluoBins(11)/(widflx(SKaP)*ryd2erg), nuArray(SKaP)
              print*, fluoLabArray(12),fluoThreshArray(12)*13.6/1.e3, Qfluo(12), &
                   &QfluoBins(12)/(widflx(ArKaP)*ryd2erg), nuArray(ArKaP)
              print*, fluoLabArray(13),fluoThreshArray(13)*13.6/1.e3, Qfluo(13), &
                   &QfluoBins(13)/(widflx(CaKaP)*ryd2erg), nuArray(CaKaP)
              print*, fluoLabArray(14),fluoThreshArray(14)*13.6/1.e3, Qfluo(14), &
                   &QfluoBins(14)/(widflx(TiKaP)*ryd2erg), nuArray(TiKaP)

           endif
        end if

        if (lgFluorescence) then
           print*, 'Total number of fluorescence photons EMITTED -before transfer-: '
           do i = 1, nfluo
              print*, fluorescenceLabel(i), fluorescenceTot(i)
           end do
        end if


        if (lgGas .and. convPercent>=resLinesTransfer .and.&
             & .not.lgResLinesFirst&
             & .and. (.not.nIterateMC==1)) then

           print*, "! energyPacketDriver: starting resonance line packets transfer"


           iCell = 0
           do igrid = 1, nGrids

              if (igrid==1) then
                 igp = 1
              else if (igrid>1) then
                 igp = 2
              else
                 print*, "! energyPacketDriver: insane grid pointer"
                 stop
              end if


              do ix = 1, grid(igrid)%nx
                 do iy = 1, grid(igrid)%ny
                    do iz = 1, grid(igrid)%nz
                       iCell = iCell+1

                       if (mod(iCell-(taskid+1),numtasks)==0) then
                          if (grid(igrid)%active(ix,iy,iz)>0) then

                             do iPhot = 1, grid(igrid)%resLinePackets(grid(igrid)%active(ix,iy,iz))

                                chTypeD = "diffuse"
                                posVector%x = grid(igrid)%xAxis(ix)
                                posVector%y = grid(igrid)%yAxis(iy)
                                posVector%z = grid(igrid)%zAxis(iz)

                                inX=-1
                                inY=-1
                                inZ=-1
                                inX(igp)=ix
                                inY(igp)=iy
                                inZ(igp)=iz

                                if (igrid>1) then
                                   ! check location on mother grid
                                   call locate(grid(grid(igrid)%motherP)%xAxis, &
                                        & posVector%x,inX(1))
                                   if (posVector%x > (grid(grid(igrid)%motherP)%xAxis(inX(1))+&
                                        & grid(grid(igrid)%motherP)%xAxis(inX(1)+1) )/2.) &
                                        & inX(1) = inX(1)+1
                                   call locate(grid(grid(igrid)%motherP)%yAxis, &
                                        & posVector%y,inY(1))
                                   if (posVector%y > (grid(grid(igrid)%motherP)%yAxis(inY(1))+&
                                        & grid(grid(igrid)%motherP)%yAxis(inY(1)+1) )/2.) &
                                        & inY(1) = inY(1)+1
                                   call locate(grid(grid(igrid)%motherP)%zAxis, &
                                        & posVector%z,inZ(1))
                                   if (posVector%z > (grid(grid(igrid)%motherP)%zAxis(inZ(1))+&
                                        & grid(grid(igrid)%motherP)%zAxis(inZ(1)+1) )/2.) &
                                        & inZ(1) = inZ(1)+1
                                end if



                                chTypeIn = chTypeD
                                positionIn = posVector
                                gPIn = igrid
                                reRun = 0

                                do i = 1, recursionLimit
                                   call energyPacketRun(chTypeIn, positionIn, inX, inY, &
                                        &inZ, gPIn, reRun)
                                   if (rerun == 0) exit
!                                   elseif (reRun == 1) then
!                                      call energyPacketRun(chTypeIn, positionIn, inX, inY, &
!                                           &inZ, gPIn, reRun)
!                                   else
!                                      print*, "! energyPacketDriver: insane reRun value ", reRun
!                                      stop
!                                   end if
                                end do
                                if (i>=recursionLimit) trapped = trapped+1

!                                call energyPacketRun(chTypeD,posVector,inX,inY,inZ,igrid)

                             end do
                          end if
                       end if

                    end do
                 end do
              end do

           end do

           print*, "! energyPacketDriver: ending resonance line packets transfer"

        end if


        if (iStar>0) print*, 'Qphot = ', Qphot, taskid
        print*, 'Packets trapped = ', trapped, taskid


        ! evaluate Jste and Jdif
        ! NOTE : deltaE is in units of [E36 erg/s] however we also need a factor of
        ! 1.e-45 from the calculations of the volume of the cell hence these
        ! two factors cancel each other out giving units of [E-9erg/s] so we need to
        ! multiply by 1.E-9
        ! NOTE : Jste and Jdif calculated this way are in units of
        ! [erg sec^-1 cm^-2] -> no Hz^-1 as they are summed over separate bins (see
        ! Lucy A&A (1999)

        if(iStar>0.) then
           print*, 'Lstar', Lstar(iStar)
        end if

        if (taskid==0 .and. lgPlaneIonization .and. niterateMC ==1) then

           if (taskid == 0 ) then
              close(29)
              open(unit=29, status='unknown', position='rewind', file='output/radField0.out',  action="write")
              do i = 1, grid(1)%nx
                 do j = 1, grid(1)%nz
                    do k = 1, nbins
                       write(29,*) nuArray(k), JnuStart(i,j,k)
                    end do
                 end do
              end do
              close(29)
              deallocate(Jnustart)
           end if

        end if

        contains

          subroutine energyPacketRun(chType, position, xP, yP, zP, gP, rR)
            implicit none

            type(vector),intent(inout) :: position         ! the position of the photon

            integer, dimension(2), intent(inout)    :: xP, yP, &
                 & zP                                            ! cartesian axes indeces
                                                                 ! 1= mother; 2=sub


            integer, intent(inout) :: rR               ! rerun?
            integer, intent(inout) :: gP               ! grid index
            integer                          :: igpr             ! grid pointer 1= mother 2=sub
            integer                          :: difSourceL(3)    ! cell indeces

            integer                          :: idirP, idirT     ! direction cosines
            integer                          :: nP

            type(photon_packet)              :: enPacket         ! the energu packet

            character(len=7), intent(inout)     :: chType           ! stellar or diffuse?

!            countRecursive = countRecursive+1
!            if (countRecursive > recursionLimit) then
!               trapped = trapped+1
!               return
!            end if

            rR = 0

            if (gP==1) then
               igpr = 1
            else if (gP>1) then
               igpr = 2
            else
               print*,  "! energyPacketRun: insane grid index"
               stop
            end if

            ! create a new photon packet
            select case (chType)

            ! if the energy packet is stellar
            case ("stellar")
               ! check for errors in the sources position
               if( position /= starPosition(iStar) ) then
                  print*, "! energyPacketRun: stellar energy packet must&
                       & start at the stellar position"
                  stop
               end if

               ! create the packet
               enPacket = newPhotonPacket(chType,position, xp, yp, zp, gp, noCellLoc,nP)

               if (nP == 1) then

                   ! the packet escapes without further interaction

                   idirT = int(acos(enPacket%direction%z)/dTheta)+1
                   if (idirT>totangleBinsTheta) then
                      idirT=totangleBinsTheta
                   end if
                   if (idirT<1 .or. idirT>totAngleBinsTheta) then
                      print*, '! energyPacketRun: error in theta direction cosine assignment',&
                           &  idirT, enPacket, dTheta, totAngleBinsTheta
                      stop
                   end if

                   if (enPacket%direction%x<1.e-35) then
                      idirP = 0
                   else
                      idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
                   end if
                   if (idirP<0) idirP=totAngleBinsPhi+idirP
                   idirP=idirP+1
                   if (idirP>totangleBinsPhi) then
                      idirP=totangleBinsPhi
                   end if

                   if (idirP<1 .or. idirP>totAngleBinsPhi) then
                      print*, '! energyPacketRun: error in Phi direction cosine assignment',&
                           &  idirP, enPacket, dPhi, totAngleBinsPhi
                      stop
                   end if


                   if (nAngleBins>0) then
                      if ( (totangleBinsPhi>1 .and. viewPointPtheta(idirT) == viewPointPphi(idirP)) .or. &
                           & (totangleBinsPhi==1 .and. viewPointPtheta(idirT)>0) ) then
                         grid(1)%escapedPackets(1, enPacket%nuP,&
                              & viewPointPtheta(idirT)) = &
                              &grid(1)%escapedPackets(1, &
                              & enPacket%nuP,viewPointPtheta(idirT)) +  deltaEUsed
                         grid(1)%escapedPackets(1, enPacket%nuP,0)=&
                              &grid(1)%escapedPackets(1, &
                              & enPacket%nuP,0) +  deltaEUsed
                         if (lgSeparateSED) then
                            grid(1)%escapedPacketsComponents(1, enPacket%nuP,&
                                 & viewPointPtheta(idirT),0) = &
                                 &grid(1)%escapedPacketsComponents(1, &
                                 & enPacket%nuP,viewPointPtheta(idirT),0) +  deltaEUsed
                            grid(1)%escapedPacketsComponents(1, enPacket%nuP,0,0)=&
                                 &grid(1)%escapedPacketsComponents(1, &
                                 & enPacket%nuP,0,0) +  deltaEUsed
                         end if
                      else
                         grid(1)%escapedPackets(1,enPacket%nuP,0) = &
                              & grid(1)%escapedPackets(1, &
                              & enPacket%nuP,0) +  deltaEUsed
                         if (lgSeparateSED) then
                            grid(1)%escapedPacketsComponents&
                                 &(1,enPacket%nuP,0,enPacket%SEDtype) = &
                                 & grid(1)%escapedPacketsComponents&
                                 &(1,  enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                         endif
                      end if
                   else
                      grid(1)%escapedPackets(1, &
                           enPacket%nuP,0) = &
                           & grid(1)%escapedPackets(1, &
                           & enPacket%nuP,0) +  deltaEUsed
                      if (lgSeparateSED) then
                         grid(1)%escapedPacketsComponents(1, &
                              enPacket%nuP,0,enPacket%SEDtype) = &
                              & grid(1)%escapedPacketsComponents(1, &
                              & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                      end if

                   end if


                   return
                end if


            ! if the photon is from an extra source of diffuse radiation
             case ("diffExt")

                ! check that the grid and cell have been specified
!                if (.not.(present(gp).and.present(xP).and.present(yP).and.present(zP))) then
!                    print*, "! energyPacketRun: gp and xp,yp and zp must be specified if iStar=0"
!                    stop
!                end if

                difSourceL(1) = xP(igpr)
                difSourceL(2) = yP(igpr)
                difSourceL(3) = zP(igpr)


                ! create the packet
!                if (.not.present(gP)) gP=1

                enPacket = newPhotonPacket(chType, position, xp, yp, zp, gP, difSourceL, nP)

            ! if the photon is diffuse

            case ("diffuse")

                ! check that the position has been specified
!                if (.not.present(position)) then
!                    print*, "! energyPacketRun: position of the new diffuse&
!                         & energy packet has not been specified"
!                    stop
!                end if

                ! check also that axes indeces have been carried through (save time)
!                if (.not.(present(xP).and.present(yP).and.present(zP))) then
!                    print*, "! energyPacketRun: cartesian axes indeces of the new diffuse &
!                         & energy packet has not been specified"
!                    stop
!                end if
                ! check also that grid index have been carried through
!                if (.not.(present(gP))) then
!                    print*, "! energyPacketRun: cartesian axes indeces of the new diffuse &
!                         & energy packet has not been specified"
!                    stop
!                end if

                ! create the packet

                enPacket = newPhotonPacket(chType=chType, position=position, xP=xP, yP=yP, &
                     & zP=zP, gP=gP, difSource=noCellLoc, nextPacket=nP)

            case ("dustEmi")

                 ! check that the position has been specified
!                if (.not.present(position)) then
!                    print*, "! energyPacketRun: position of the new dust emitted &
!                         & energy packet has not been specified"
!                    stop
!                end if

                ! check also that axes indeces have been carried through (save time)
!                if (.not.(present(xP).and.present(yP).and.present(zP))) then
!                    print*, "! energyPacketRun: cartesian axes indeces of the new dust emitted &
!                         &energy packet has not been specified"
!                    stop
!                end if
                ! check also that grud index have been carried through
!                if (.not.(present(gP))) then
!                    print*, "! energyPacketRun: cartesian axes indeces of the new dust emitted &
!                         &energy packet has not been specified"
!                    stop
!                end if

                ! crate the packet
                enPacket = newPhotonPacket(chType=chType, position=position, xP=xP, yP=yP, &
                     &zP=zP, gP=gP, difSource = noCellLoc,  nextPacket=nP)

            end select

            if (.not.lgDust .and. enPacket%nu < ionEdge(1) .and. .not.enPacket%lgLine) then

               ! the packet escapes without further interaction
               idirT = int(acos(enPacket%direction%z)/dTheta)+1
               if (idirT>totangleBinsTheta) then
                  idirT=totangleBinsTheta
               end if
               if (idirT<1 .or. idirT>totAngleBinsTheta) then
                  print*, '! energyPacketRun: error in theta direction cosine assignment',&
                       &  idirT, enPacket, dTheta, totAngleBinsTheta
                  stop
               end if


               if (enPacket%direction%x<1.e-35) then
                  idirP = 0
               else
                  idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
               end if
               if (idirP<0) idirP=totAngleBinsPhi+idirP
               idirP=idirP+1

               if (idirP>totangleBinsPhi) then
                  idirP=totangleBinsPhi
!                  print*, '! energyPacketRun: idir>totanglebins - error corrected', &
!                       & idir, totanglebins, enPacket%direction, dtheta
               end if

               if (idirP<1 .or. idirP>totAngleBinsPhi) then
                  print*, '! energyPacketRun: error in phi direction cosine assignment',&
                       &  idirP, enPacket, dPhi, totAngleBinsPhi
                  stop
               end if


               if (nAngleBins>0) then
                  if ( (totangleBinsPhi>1 .and. viewPointPtheta(idirT) == viewPointPphi(idirP)) .or. &
                       & (totangleBinsPhi==1 .and. viewPointPtheta(idirT)>0) ) then
                     grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,&
                          & viewPointPtheta(idirT)) = &
                          &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                          & enPacket%nuP,viewPointPtheta(idirT)) +  deltaEUsed
                     grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,0) = &
                          &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                          & enPacket%nuP,0) +  deltaEUsed

                     if (lgSeparateSED) then
                        grid(enPacket%origin(1))%escapedPacketsComponents&
                             & (enPacket%origin(2), enPacket%nuP,&
                             & viewPointPtheta(idirT),enPacket%SEDtype) = &
                             &grid(enPacket%origin(1))%escapedPacketsComponents&
                             &(enPacket%origin(2), &
                             & enPacket%nuP,viewPointPtheta(idirT),enPacket%SEDtype) +  deltaEUsed
                        grid(enPacket%origin(1))%escapedPacketsComponents&
                             &(enPacket%origin(2), enPacket%nuP,0,enPacket%SEDtype) = &
                             &grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                             & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                     endif
                  else
                     grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2),enPacket%nuP,0) = &
                          & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                          & enPacket%nuP,0) +  deltaEUsed
                     if (lgSeparateSED) then
                        grid(enPacket%origin(1))%escapedPacketsComponents&
                             &(enPacket%origin(2),enPacket%nuP,0,enPacket%SEDtype) = &
                             & grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                             & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                     endif

                  end if
               else
                  grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                       & enPacket%nuP,0) = &
                       & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                       & enPacket%nuP,0) +  deltaEUsed
                  if (lgSeparateSED) then
                     grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                          & enPacket%nuP,0,enPacket%SEDtype) = &
                          & grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                          & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                  endif

               end if

               return
            end if

            ! if the new packet is capable of re-ionizing or we have dust
            if (.not.enPacket%lgLine) then

                ! compute the next segment of trajectory
                call pathSegment(enPacket)

                return

            else ! if the packet is a line packet
                ! add to respective line packet bin

               if (lgDebug) &
                    & grid(gP)%linePackets(grid(gP)%active(enPacket%xP(igpr), &
                    & enPacket%yP(igpr), enPacket%zP(igpr)), enPacket%nuP) = &
                    & grid(gP)%linePackets(grid(gP)%active(enPacket%xP(igpr), &
                    & enPacket%yP(igpr), enPacket%zP(igpr)), enPacket%nuP) +deltaEUsed


            end if

        end subroutine energyPacketRun

        ! this function initializes a photon packet
        function initPhotonPacket(nuP,  position, lgLine, lgStellar, xP, yP, zP, gP)
            implicit none

            type(vector), intent(in) :: position          ! the position at which the photon
                                                          ! packet is created
            real                     :: random            ! random number

            integer, intent(in)      :: nuP               ! the frequency of the photon packet
            integer, intent(in),dimension(2) :: xP, yP, &
                 & zP                                     ! indeces of position on the x, y and z axes
            integer, intent(in)      :: gP                ! grid index
            integer                  :: igpi              ! grid pointer 1=mother, 2=sub
            integer                  :: irepeat           ! counter

            type(photon_packet)      :: initPhotonPacket  ! the photon packet

            logical, intent(in)      :: lgLine, lgStellar ! line, stellar packet?


            initPhotonPacket%position = position

            initPhotonPacket%iG  = gP

            if (gP==1) then
               igpi=1
            else if (gp>1) then
               igpi=2
            else
               print*, "! initPhotonPacket: insane gridp pointer"
               stop
            end if

            initPhotonPacket%nuP      = nuP

            initPhotonPacket%lgStellar = lgStellar

            ! check if photon packen is line or continuum photon
            if ( lgLine ) then
                ! line photon
                initPhotonPacket%nu       = 0.
                initPhotonPacket%lgLine   = .true.
            else
                ! continuum photon
                initPhotonPacket%nu       = nuArray(nuP)
                initPhotonPacket%lgLine   = .false.
            end if

            initPhotonPacket%xP  = xP
            initPhotonPacket%yP  = yP
            initPhotonPacket%zP  = zP

            ! cater for plane parallel ionization case
            if (initPhotonPacket%lgStellar .and. lgPlaneIonization) then

               ! get position

               ! x-direction
               call random_number(random)
               random = 1. - random
               initPhotonPacket%position%x = &
                    & grid(gP)%xAxis(1) + random*( &
                    & grid(gP)%xAxis(grid(gP)%nx)-grid(gP)%xAxis(1))
               if (initPhotonPacket%position%x<grid(gP)%xAxis(1) .or.&
                    &initPhotonPacket%position%x>grid(gP)%xAxis(grid(gP)%nx) ) then
                  print*, ' ! initPhotonPacket: planeParallelIonisation &
                       &photon born at insane x-location', random
                  stop
               end if

               call locate(grid(gP)%xAxis, initPhotonPacket%position%x, initPhotonPacket%xP(igpi))
               if (initPhotonPacket%xP(igpi) < grid(gP)%nx) then
                  if (initPhotonPacket%position%x >= (grid(gP)%xAxis(initPhotonPacket%xP(igpi))+&
                       & grid(gP)%xAxis(initPhotonPacket%xP(igpi)+1))/2.) &
                       & initPhotonPacket%xP(igpi) = initPhotonPacket%xP(igpi)+1
               end if

               ! y-direction
               initPhotonPacket%position%y = 0.
               initPhotonPacket%yP(igpi) = 1

               ! z-direction
               call random_number(random)
               random = 1. - random
               initPhotonPacket%position%z = &
                    & grid(gP)%zAxis(1) + random*( &
                    & grid(gP)%zAxis(grid(gP)%nz)-grid(gP)%zAxis(1))
               if (initPhotonPacket%position%z<grid(gP)%zAxis(1) .or.&
                    &initPhotonPacket%position%z>grid(gP)%zAxis(grid(gP)%nz) ) then
                  print*, ' ! initPhotonPacket: planeParallelIonisation &
                       &photon born at insane z-location', random
                  stop
               end if
               call locate(grid(gP)%zAxis, initPhotonPacket%position%z, initPhotonPacket%zP(igpi))
               if (initPhotonPacket%zP(igpi) < grid(gP)%nz) then
                  if (initPhotonPacket%position%z >= (grid(gP)%zAxis(initPhotonPacket%zP(igpi))+&
                       & grid(gP)%zAxis(initPhotonPacket%zP(igpi)+1))/2.) initPhotonPacket%zP(igpi) =&
                       & initPhotonPacket%zP(igpi)+1
               end if

               if (initPhotonPacket%xP(igpi)<1) initPhotonPacket%xP(igpi)=1
               if (initPhotonPacket%zP(igpi)<1) initPhotonPacket%zP(igpi)=1

               ! direction is parallel to y-axis direction
               initPhotonPacket%direction%x = 0.
               initPhotonPacket%direction%y = 1.
               initPhotonPacket%direction%z = 0.

               if (initPhotonPacket%xP(igpi) >  grid(gP)%xAxis(grid(gP)%nx) .or. &
                    & initPhotonPacket%zP(igpi) >  grid(gP)%zAxis(grid(gP)%nz)) then
                  print*, "! initPhotonPacket: insanity in planeIonisation init"
                  print*, igpi, initPhotonPacket%xP(igpi),  grid(gP)%xAxis(grid(gP)%nx), &
                       & initPhotonPacket%zP(igpi), grid(gP)%zAxis(grid(gP)%nz),  random, &
                       &initPhotonPacket%position%z

                  stop
               end if

                planeIonDistribution(initPhotonPacket%xP(igpi),initPhotonPacket%zP(igpi)) = &
                     & planeIonDistribution(initPhotonPacket%xP(igpi),initPhotonPacket%zP(igpi)) + 1

             else

                do irepeat = 1, 1000000
                   ! get a random direction
                   initPhotonPacket%direction = randomUnitVector()
                   if (initPhotonPacket%direction%x/=0. .and. &
                        & initPhotonPacket%direction%y/=0. .and. &
                        & initPhotonPacket%direction%z/=0.) exit
                end do
            end if

            if ((lgSymmetricXYZ) .and. initPhotonPacket%lgStellar .and. .not.lgMultistars) then
                if (initPhotonPacket%direction%x<0.) &
                     & initPhotonPacket%direction%x = -initPhotonPacket%direction%x
                if (initPhotonPacket%direction%y<0.) &
                     & initPhotonPacket%direction%y = -initPhotonPacket%direction%y
                if (initPhotonPacket%direction%z<0.) &
                     & initPhotonPacket%direction%z = -initPhotonPacket%direction%z
            end if

            initPhotonPacket%origin(1) = gP
            if (initPhotonPacket%xP(igpi)>0 .and. initPhotonPacket%yP(igpi) > 0. .and. &
                 & initPhotonPacket%zP(igpi) > 0.) then
               initPhotonPacket%origin(2) = grid(gP)%active(initPhotonPacket%xP(igpi),&
                    & initPhotonPacket%yP(igpi), initPhotonPacket%zP(igpi))
            else
               initPhotonPacket%origin(2) = 0
            end if


            if (lgStellar) initPhotonPacket%SEDtype = 0

        end function initPhotonPacket


        ! this subroutine determines the frequency of a newly created photon packet
        ! according to the given probability density
        subroutine getNu(probDen, nuP)

            real, dimension(:), intent(in) :: probDen    ! probability density function

            integer, intent(out)           :: nuP         ! frequency index of the new

            ! local variables
            real                           :: random     ! random number

            ! get a random number
            call random_number(random)

            random = 1.-random

            ! see what frequency random corresponds to
            call locate(probDen, random, nuP)
             if (nuP <= 0) nuP = 1

             if (nuP<nbins) then
                nuP=nuP+1
             end if

        end subroutine getNu

        ! this subroutine determines the frequency of a newly created photon packet
        ! according to the given probability density
        ! does not use bisection to locate nu on array
        subroutine getNu2(probDen, nuP)

            real, dimension(:), intent(in) :: probDen    ! probability density function

            real                           :: random     ! random number

            integer, intent(out)           :: nuP        ! frequency index of the new

            integer                        :: isearch,i  !

            ! get a random number
            call random_number(random)

            do i = 1, 10000
               if (random==0 .or. random==1.) then
                  call random_number(random)
               else
                  exit
               end if
            end do
            if (i>=10000) then
               print*, '! getNu2: problem with random number generator', random, i
               stop
            end if

            ! see what frequency random corresponds to
            nuP=1
            do isearch = 1, nbins
               if (random>=probDen(isearch)) then
                  nuP=isearch
               else
                  exit
               end if
            end do

            if (nuP<nbins) then
               nuP=nuP+1
            end if

            if (nuP>nbins) then
               print*, 'random: ', random
               print*, 'probDen: ', probDen
            end if

          end subroutine getNu2


        ! this function creates a new photon packet
        function newPhotonPacket(chType, position, xP, yP, zP, gP, difSource,nextPacket)

            type(photon_packet)                :: newPhotonPacket! the photon packet to be created
            type(vector), intent(in)           :: position       ! the position of the photon packet
            type(vector)                       :: positionLoc    ! the position of the photon packet

            real                               :: random         ! random number
            real                    :: xloc, yloc, zloc, dSx, dSy, dSz, dS, dSbis, dSbisbis

            integer                            :: nuP            ! the frequency index of the photon packet
            integer, dimension(2)         :: orX,orY,orZ    ! dummy
            integer, dimension(2),intent(in) :: xP, yP, &
                 & zP                                            ! cartesian axes indeces
            integer, intent(in)      :: difSource(3)  ! grid and cell indeces
            integer                            :: ids            ! counter

            integer, intent(inout)   :: gP
            integer                            :: igpn           ! grid pointe 1=motehr, 2=sub
            integer, intent(out)               :: nextPacket

            character(len=7), intent(in)       :: chType         ! stellar or diffuse?

            nextPacket = 0

            if (gP==1) then
               igpn = 1
            else if (gp>1) then
               igpn = 2
            else
               print*,  "! newPhotonPacket: insane grid pointer"
               stop
            end if

            select case (chType)

            ! if the photon is stellar
            case ("stellar")

                ! check for errors in the sources position
               if( position /= starPosition(iStar) ) then
                  print*, "! newPhotonPacket: stellar photon packet must&
                       & start at the stellar position"
                  stop
               end if


               !gP = starIndeces(iStar,4)

               if (starIndeces(iStar,4) == 1) then
                  igpn = 1
               else if (starIndeces(iStar,4) > 1) then
                  igpn = 2
               else
                  print*,  "! newPhotonPacket: insane grid pointer -star position- "
                  stop
               end if

               ! determine the frequency of the newly created photon packet
               call getNu2(inSpectrumProbDen(iStar,1:nbins), nuP)

                if (nuP>nbins) then
                   print*, "! newPhotonPacket: insanity occured in stellar photon &
                        &nuP assignment (nuP,xP,yP,zP,activeP)", nuP, xP(igpn),yP(igpn),zP(igpn), &
                        & grid(starIndeces(iStar,4))%active(xP(igpn),yP(igpn),zP(igpn))
                   print*, "inSpectrumProbDen: ",iStar,inSpectrumProbDen(iStar,:), nuP
                   stop
                end if

                if (nuP < 1) then
                    print*, "! newPhotonPacket: insanity occured in stellar photon &
                         &nuP assignment"
                    stop
                 end if


                ! initialize the new photon packet
                if (starIndeces(iStar,1) > 0 .and. starIndeces(iStar,2) > 0 &
                     & .and. starIndeces(iStar,3) > 0 ) then
                   orX(igpn) = starIndeces(iStar,1)
                   orY(igpn) = starIndeces(iStar,2)
                   orZ(igpn) = starIndeces(iStar,3)

                   if (grid(starIndeces(iStar,4))%active(orX(igpn),&
                        &  orY(igpn), orZ(igpn)) < 0.) then
                      print*, "! newPhotonPacket: new packet cannot be&
                           & emitted from re-mapped cell -1-"
                      print*, "chType, nuP, starPosition(iStar), .false., &
                           & .true., orX,orY,orZ, gp"
                      print*, chType, nuP, starPosition(iStar), .false., .true., orX,orY,orZ, gp
                      stop
                   end if

                else
                   orX=-1
                   orY=-1
                   orZ=-1
                end if
!print*, orx, ory, orz
                newPhotonPacket = initPhotonPacket(nuP, starPosition(iStar), .false., &
                     &.true., orX,orY,orZ, starIndeces(iStar,4))

                newPhotonPacket%SEDtype = 0

                if (starIndeces(iStar,1) < 0 .or. starIndeces(iStar,2) < 0 &
                     & .or. starIndeces(iStar,3) < 0) then

                   ! use direction info and grid extents to find the nearest entry to the grid
                   xloc = starPosition(iStar)%x
                   yloc = starPosition(iStar)%y
                   zloc = starPosition(iStar)%z

                   dSx=0.
                   dSy=0.
                   dSz=0.
!print*, xloc, yloc, zloc

!newPhotonPacket%direction%x = 0.
!newPhotonPacket%direction%y = 0.
!newPhotonPacket%direction%z = -1.

!print*, newPhotonPacket%direction%x, newPhotonPacket%direction%y, newPhotonPacket%direction%z

                   if ( xloc > grid(1)%xAxis(grid(1)%nx) .or. &
                        & xloc < grid(1)%xAxis(1)) then

                      if (xloc < grid(1)%xAxis(1) .and. lgSymmetricXYZ) then
                         print*, ' ! newPhotonPacket: star cannot be placed &
                              &behind a mirror -x', xloc
                         stop
                      end if

                      if( newPhotonPacket%direction%x > 0.) then
                         if (xloc > grid(1)%xAxis(grid(1)%nx)) then
                            dSx = 0.
                         else if (xloc < grid(1)%xAxis(1)) then
                            dSx = ( xloc-grid(1)%xAxis(1))/&
                                 &newPhotonPacket%direction%x
                         end if
                      else if (newPhotonPacket%direction%x < 0.) then
                         if (xloc < grid(1)%xAxis(1)) then
                            dSx = 0.
                         elseif (xloc > grid(1)%xAxis(grid(1)%nx)) then
                            dSx = ( xloc-grid(1)%xAxis(grid(1)%nx))/&
                                 &newPhotonPacket%direction%x
                         end if
                      end if

                   endif


                   if ( yloc > grid(1)%yAxis(grid(1)%ny) .or. &
                        & yloc < grid(1)%yAxis(1)) then

                      if (yloc < grid(1)%yAxis(1) .and. lgSymmetricXYZ) then
                         print*, ' ! newPhotonPacket: star cannot be placed &
                              &behind a mirror -y', yloc
                         stop
                      end if

                      if( newPhotonPacket%direction%y > 0.) then
                         if (yloc > grid(1)%yAxis(grid(1)%ny)) then
                            dSy = 0.
                         else if (yloc < grid(1)%yAxis(1)) then
                            dSy = ( yloc-grid(1)%yAxis(1))/&
                                 &newPhotonPacket%direction%y
                         end if
                      else if (newPhotonPacket%direction%y < 0.) then
                         if (yloc < grid(1)%yAxis(1)) then
                            dSy = 0.
                         elseif (yloc > grid(1)%yAxis(grid(1)%ny)) then
                            dSy = ( yloc-grid(1)%yAxis(grid(1)%ny))/&
                                 &newPhotonPacket%direction%y
                         end if
                      end if

                   endif


                   if ( zloc > grid(1)%zAxis(grid(1)%nz) .or. &
                        & zloc < grid(1)%zAxis(1)) then

                      if (zloc < grid(1)%zAxis(1) .and. lgSymmetricXYZ) then
                         print*, ' ! newPhotonPacket: star cannot be placed &
                              &behind a mirror -z', zloc
                         stop
                      end if

                      if( newPhotonPacket%direction%z > 0.) then
                         if (zloc > grid(1)%zAxis(grid(1)%nz)) then
                            dSz = 0.
                         else if (zloc < grid(1)%zAxis(1)) then
                            dSz = ( zloc-grid(1)%zAxis(1))/&
                                 &newPhotonPacket%direction%z
                         end if
                      else if (newPhotonPacket%direction%z < 0.) then
                         if (zloc < grid(1)%zAxis(1)) then
                            dSz = 0.
                         elseif (zloc > grid(1)%zAxis(grid(1)%nz)) then
                            dSz = ( zloc-grid(1)%zAxis(grid(1)%nz))/&
                                 &newPhotonPacket%direction%z
                         end if
                      end if

                   endif



                   dSx=abs(dSx)
                   dSy=abs(dSy)
                   dSz=abs(dSz)
                   dS=0.
                   dSbis=0.
                   dSbisbis=0.

                   if (dSx /= 0.) then
                      dS = dSx
                   end if
                   if (dSy /= 0.) then
                      if (dS /= 0.) then
                         dSbis = max(dS,dSy)
                         dS = min(dS,dSy)
                      else
                         dS = dSy
                      endif
                   end if
                   if (dSz /= 0.) then
                      if (dS /= 0.) then
                         dSbisbis = dSbis
                         dSbis = max(dS,dSz)
                         dS = min(dS,dSz)
                      else
                         dS = dSz
                      endif
                   end if


                   if (dS == 0.) then ! this packet is not going to hit the grid
                      nextPacket = 1
                   else
                      do ids = 1, 3
                         if (ids == 2) then
                            if (dSbis /= 0.) then
                               dS = dSbis-dS
                            else
                               exit
                            end if
                         end if
                         if (ids == 3) then
                            if (dSbisbis /= 0.) then
                               dS = dSbisbis-dS
                            else
                               exit
                            end if
                         end if

                         if (ids ==1 ) then

                            if (dS == dSx .and. dSx /= 0.) then
                               if (xloc> grid(1)%xAxis(grid(1)%nx)) orX(igpn) = grid(1)%nx
                               if (xloc< grid(1)%xAxis(1)) orX(igpn) = 1

                               newPhotonPacket%position = newPhotonPacket%position + &
                                    &dS*newPhotonPacket%direction

                               if (newPhotonPacket%position%y <=  grid(1)%yAxis(grid(1)%ny)) then
                                  call locate( grid(1)%yAxis, newPhotonPacket%position%y,  orY(igpn))
                                  if (orY(igpn) < grid(1)%ny .and. igpn >0) then
                                     if (newPhotonPacket%position%y > (grid(1)%yAxis(orY(igpn)) + &
                                          & grid(1)%yAxis(orY(igpn)+1))/2.) then
                                        orY(igpn) = orY(igpn)+1
                                     end if
                                  end if
                               elseif (newPhotonPacket%position%y >  grid(1)%yAxis(grid(1)%ny).and. &
                                    & newPhotonPacket%position%y <= &
                                    &(3.*grid(1)%yAxis(grid(1)%ny)-grid(1)%yAxis(grid(1)%ny-1))/2.) then
                                  orY(igpn)  = grid(1)%ny
                               elseif (newPhotonPacket%position%y <  grid(1)%yAxis(1) .and. &
                                    & newPhotonPacket%position%y >= &
                                    &(3.*grid(1)%yAxis(1)-grid(1)%yAxis(2))/2.) then
                                  orY(igpn)  = 1
                               else
                                  orY(igpn) = -1
                               endif

                               if (newPhotonPacket%position%z <=  grid(1)%zAxis(grid(1)%nz)) then
                                  call locate( grid(1)%zAxis, newPhotonPacket%position%z,  orZ(igpn))
                                  if (orZ(igpn) < grid(1)%nz.and. igpn >0) then
                                     if (newPhotonPacket%position%z > (grid(1)%zAxis(orZ(igpn)) + &
                                          & grid(1)%zAxis(orZ(igpn)+1))/2.) then
                                        orZ(igpn) = orZ(igpn)+1
                                     end if
                                  end if
                               elseif (newPhotonPacket%position%z >  grid(1)%zAxis(grid(1)%nz).and. &
                                    & newPhotonPacket%position%z <= &
                                    &(3.*grid(1)%zAxis(grid(1)%nz)-grid(1)%zAxis(grid(1)%nz-1))/2.) then
                                  orZ(igpn)  = grid(1)%nz
                               elseif (newPhotonPacket%position%z <  grid(1)%zAxis(1) .and. &
                                    & newPhotonPacket%position%z >= &
                                    &(3.*grid(1)%zAxis(1)-grid(1)%zAxis(2))/2.) then
                                  orZ(igpn)  = 1
                               else
                                  orZ(igpn) = -1
                               endif

                            elseif  (dS == dSy .and. dSy /= 0.) then
                               if (yloc> grid(1)%yAxis(grid(1)%ny)) orY(igpn) = grid(1)%ny
                               if (yloc< grid(1)%yAxis(1)) orY(igpn) = 1
                               newPhotonPacket%position = newPhotonPacket%position + &
                                    &dS*newPhotonPacket%direction

                               if (newPhotonPacket%position%x <=  grid(1)%xAxis(grid(1)%nx)) then
                                  call locate( grid(1)%xAxis, newPhotonPacket%position%x,  orX(igpn))
                                  if (orX(igpn) < grid(1)%nx.and. igpn >0) then
                                     if (newPhotonPacket%position%x > (grid(1)%xAxis(orX(igpn)) + &
                                          & grid(1)%xAxis(orX(igpn)+1))/2.) then
                                        orX(igpn) = orX(igpn)+1
                                     end if
                                  end if
                               elseif (newPhotonPacket%position%x >  grid(1)%xAxis(grid(1)%nx).and. &
                                    & newPhotonPacket%position%x <= &
                                    &(3.*grid(1)%xAxis(grid(1)%nx)-grid(1)%xAxis(grid(1)%nx-1))/2.) then
                                  orX(igpn)  = grid(1)%nx
                               elseif (newPhotonPacket%position%x <  grid(1)%xAxis(1) .and. &
                                    & newPhotonPacket%position%x >= &
                                    &(3.*grid(1)%xAxis(1)-grid(1)%xAxis(2))/2.) then
                                  orX(igpn)  = 1
                               else
                                  orX(igpn) = -1
                               endif

                               if (newPhotonPacket%position%z <=  grid(1)%zAxis(grid(1)%nz)) then
                                  call locate( grid(1)%zAxis, newPhotonPacket%position%z,  orZ(igpn))
                                  if (orZ(igpn) < grid(1)%nz.and. igpn >0) then
                                     if (newPhotonPacket%position%z > (grid(1)%zAxis(orZ(igpn)) + &
                                          & grid(1)%zAxis(orZ(igpn)+1))/2.) then
                                        orZ(igpn) = orZ(igpn)+1
                                     end if
                                  end if
                               elseif (newPhotonPacket%position%z >  grid(1)%zAxis(grid(1)%nz).and. &
                                    & newPhotonPacket%position%z <= &
                                    &(3.*grid(1)%zAxis(grid(1)%nz)-grid(1)%zAxis(grid(1)%nz-1))/2.) then
                                  orZ(igpn)  = grid(1)%nz
                               elseif (newPhotonPacket%position%z <  grid(1)%zAxis(1) .and. &
                                    & newPhotonPacket%position%z >= &
                                    &(3.*grid(1)%zAxis(1)-grid(1)%zAxis(2))/2.) then
                                  orZ(igpn)  = 1
                               else
                                  orZ(igpn) = -1
                               endif

                            elseif  (dS == dSz .and. dSz /= 0.) then
                               if (zloc> grid(1)%zAxis(grid(1)%nz)) orZ(igpn) = grid(1)%nz
                               if (zloc< grid(1)%zAxis(1)) orZ(igpn) = 1
                               newPhotonPacket%position = newPhotonPacket%position + &
                                    &dS*newPhotonPacket%direction

                               if (newPhotonPacket%position%y <=  grid(1)%yAxis(grid(1)%ny)) then
                                  call locate( grid(1)%yAxis, newPhotonPacket%position%y,  orY(igpn))
                                  if (orY(igpn) < grid(1)%ny .and. igpn >0) then
                                     if (newPhotonPacket%position%y > (grid(1)%yAxis(orY(igpn)) + &
                                          & grid(1)%yAxis(orY(igpn)+1))/2.) then
                                        orY(igpn) = orY(igpn)+1
                                     end if
                                  end if

                               elseif (newPhotonPacket%position%y >  grid(1)%yAxis(grid(1)%ny).and. &
                                    & newPhotonPacket%position%y <= &
                                    &(3.*grid(1)%yAxis(grid(1)%ny)-grid(1)%yAxis(grid(1)%ny-1))/2.) then
                                  orY(igpn)  = grid(1)%ny

                               elseif (newPhotonPacket%position%y <  grid(1)%yAxis(1) .and. &
                                    & newPhotonPacket%position%y >= &
                                    &(3.*grid(1)%yAxis(1)-grid(1)%yAxis(2))/2.) then
                                  orY(igpn)  = 1


                               else
                                  orY(igpn) = -1
                               endif


                               if (newPhotonPacket%position%x <=  grid(1)%xAxis(grid(1)%nx)) then
                                  call locate( grid(1)%xAxis, newPhotonPacket%position%x,  orX(igpn))
                                  if (orX(igpn) < grid(1)%nx.and. igpn >0) then
                                     if (newPhotonPacket%position%x > (grid(1)%xAxis(orX(igpn)) + &
                                       & grid(1)%xAxis(orX(igpn)+1))/2.) then
                                        orX(igpn) = orX(igpn)+1
                                     end if
                                  end if
                               elseif (newPhotonPacket%position%x >  grid(1)%xAxis(grid(1)%nx).and. &
                                    & newPhotonPacket%position%x <= &
                                    &(3.*grid(1)%xAxis(grid(1)%nx)-grid(1)%xAxis(grid(1)%nx-1))/2.) then
                                  orX(igpn)  = grid(1)%nx
                               elseif (newPhotonPacket%position%x <  grid(1)%xAxis(1) .and. &
                                    & newPhotonPacket%position%x >= &
                                    &(3.*grid(1)%xAxis(1)-grid(1)%xAxis(2))/2.) then
                                  orX(igpn)  = 1
                               else
                                  orX(igpn) = -1
                               end if

                            else
                               print*, '! newPhotonPacket: insane dS', dS, dSx, dSy, dSz
                               stop
                            end if
                         else
                            newPhotonPacket%position = newPhotonPacket%position + &
                                 &dS*newPhotonPacket%direction

                            if (newPhotonPacket%position%y <=  grid(1)%yAxis(grid(1)%ny)) then
                               call locate( grid(1)%yAxis, newPhotonPacket%position%y,  orY(igpn))
                               if (orY(igpn) < grid(1)%ny .and. igpn >0) then
                                  if (newPhotonPacket%position%y > (grid(1)%yAxis(orY(igpn)) + &
                                       & grid(1)%yAxis(orY(igpn)+1))/2.) then
                                     orY(igpn) = orY(igpn)+1
                                  end if
                               end if
                            elseif (newPhotonPacket%position%y >  grid(1)%yAxis(grid(1)%ny).and. &
                                 & newPhotonPacket%position%y <= &
                                 &(3.*grid(1)%yAxis(grid(1)%ny)-grid(1)%yAxis(grid(1)%ny-1))/2.) then
                               orY(igpn)  = grid(1)%ny
                            elseif (newPhotonPacket%position%y <  grid(1)%yAxis(1) .and. &
                                 & newPhotonPacket%position%y >= &
                                 &(3.*grid(1)%yAxis(1)-grid(1)%yAxis(2))/2.) then
                               orY(igpn)  = 1
                            else
                               orY(igpn) = -1
                            endif

                            if (newPhotonPacket%position%x <=  grid(1)%xAxis(grid(1)%nx)) then
                               call locate( grid(1)%xAxis, newPhotonPacket%position%x,  orX(igpn))
                               if (orX(igpn) < grid(1)%nx.and. igpn >0) then
                                  if (newPhotonPacket%position%x > (grid(1)%xAxis(orX(igpn)) + &
                                       & grid(1)%xAxis(orX(igpn)+1))/2.) then
                                     orX(igpn) = orX(igpn)+1
                                  end if
                               end if
                            elseif (newPhotonPacket%position%x >  grid(1)%xAxis(grid(1)%nx).and. &
                                 & newPhotonPacket%position%x <= &
                                 &(3.*grid(1)%xAxis(grid(1)%nx)-grid(1)%xAxis(grid(1)%nx-1))/2.) then
                               orX(igpn)  = grid(1)%nx
                            elseif (newPhotonPacket%position%x <  grid(1)%xAxis(1) .and. &
                                 & newPhotonPacket%position%x >= &
                                 &(3.*grid(1)%xAxis(1)-grid(1)%xAxis(2))/2.) then
                               orX(igpn)  = 1
                            else
                               orX(igpn) = -1
                            endif

                            if (newPhotonPacket%position%z <=  grid(1)%zAxis(grid(1)%nz)) then
                               call locate( grid(1)%zAxis, newPhotonPacket%position%z,  orZ(igpn))
                               if (orZ(igpn) < grid(1)%nz.and. igpn >0) then
                                  if (newPhotonPacket%position%z > (grid(1)%zAxis(orZ(igpn)) + &
                                    & grid(1)%zAxis(orZ(igpn)+1))/2.) then
                                     orZ(igpn) = orZ(igpn)+1
                                  end if
                               end if
                            elseif (newPhotonPacket%position%z >  grid(1)%zAxis(grid(1)%nz).and. &
                                 & newPhotonPacket%position%z <= &
                                 &(3.*grid(1)%zAxis(grid(1)%nz)-grid(1)%zAxis(grid(1)%nz-1))/2.) then
                               orZ(igpn)  = grid(1)%nz
                            elseif (newPhotonPacket%position%z <  grid(1)%zAxis(1) .and. &
                                 & newPhotonPacket%position%z >= &
                                 &(3.*grid(1)%zAxis(1)-grid(1)%zAxis(2))/2.) then
                               orZ(igpn)  = 1
                            else
                               orZ(igpn) = -1
                            endif

                         end if

                         if (orX(igpn) >= 1 .and. orY(igpn) >= 1 .and. orZ(igpn) >= 1)  exit

                      end do

                   end if

                   if (orX(igpn) < 1 .or. orY(igpn) < 1 .or.orZ(igpn) < 1  ) then
                      nextPacket = 1
                   else

                      newPhotonPacket%xP(igPn) =  orX(igpn)
                      newPhotonPacket%yP(igPn) =  orY(igpn)
                      newPhotonPacket%zP(igPn) =  orZ(igpn)

                      newPhotonPacket%origin(1) = newPhotonPacket%iG
                      newPhotonPacket%origin(2) = &
                           &grid(newPhotonPacket%iG)%active(newPhotonPacket%xP(igpn),&
                           & newPhotonPacket%yP(igpn), newPhotonPacket%zP(igpn))
                   end if

                endif

                if (newPhotonPacket%nu>=1.) then
                   Qphot = Qphot + deltaEUsed/(2.1799153e-11*newPhotonPacket%nu)
                end if

                if (lgFluorescence) then
                   do i = 1, 14
                      if (newPhotonPacket%nu>=fluoThreshArray(i)) then
                         Qfluo(i) = Qfluo(i) + deltaEUsed/(ryd2erg*newPhotonPacket%nu)
                      end if
                      if (newPhotonPacket%nuP == FeKaColdP) &
                           & QfluoBins(1) = QfluoBins(1) + deltaEUsed

                      if (newPhotonPacket%nuP == FeL1P) &
                           & QfluoBins(2) = QfluoBins(2) + deltaEUsed
                      if (newPhotonPacket%nuP == FeL2P) &
                           & QfluoBins(3) = QfluoBins(3) + deltaEUsed
                      if (newPhotonPacket%nuP == CKaP) &
                           & QfluoBins(4) = QfluoBins(4) + deltaEUsed
                      if (newPhotonPacket%nuP == NKaP) &
                           & QfluoBins(5) = QfluoBins(5) + deltaEUsed
                      if (newPhotonPacket%nuP == OKaP) &
                           & QfluoBins(6) = QfluoBins(6) + deltaEUsed
                      if (newPhotonPacket%nuP == NeKaP) &
                           & QfluoBins(7) = QfluoBins(7) + deltaEUsed
                      if (newPhotonPacket%nuP == MgKaP) &
                           & QfluoBins(8) = QfluoBins(8) + deltaEUsed
                      if (newPhotonPacket%nuP == AlKaP) &
                           & QfluoBins(9) = QfluoBins(9) + deltaEUsed
                      if (newPhotonPacket%nuP == SiKaP) &
                           & QfluoBins(10) = QfluoBins(10) + deltaEUsed
                      if (newPhotonPacket%nuP == SKaP) &
                           & QfluoBins(11) = QfluoBins(11) + deltaEUsed
                      if (newPhotonPacket%nuP == ArKaP) &
                           & QfluoBins(12) = QfluoBins(12) + deltaEUsed
                      if (newPhotonPacket%nuP == CaKaP) &
                           & QfluoBins(13) = QfluoBins(13) + deltaEUsed
                      if (newPhotonPacket%nuP == TiKaP) &
                           & QfluoBins(14) = QfluoBins(14) + deltaEUsed
                   end do
                end if

                if (taskid==0 .and. niterateMC==1 .and. lgPlaneIonization) then
                   JnuStart(newPhotonPacket%xP(igpn),newPhotonPacket%zP(igpn),&
                        & newPhotonPacket%nuP) = &
                        &  JnuStart(newPhotonPacket%xP(igpn),newPhotonPacket%zP(igpn),&
                        & newPhotonPacket%nuP)&
                        & + deltaEUsed/widFlx(newPhotonPacket%nuP)
                end if

                ! if the photon is from an extra diffuse source
             case ("diffExt")

                ! check that the grid and cell indeces have been specified
!                if (.not.(present(gp).and.present(difSource))) then
!                    print*, "! newPhotonPacket: grid and cell indeces of the new extra diffuse &
!                         & photon packet have not been specified"
!                    stop
!                end if

                call getNu2(inSpectrumProbDen(0,1:nbins), nuP)

                if (nuP>nbins) then
                   print*, "! newPhotonPacket: insanity occured in extra diffuse photon &
                        & nuP assignment (nuP,gp,activeP)", nuP, gp
                   print*, "difSpectrumProbDen: ", inSpectrumProbDen(0,:)
                   stop
                end if

                if (nuP < 1) then
                   print*, "! newPhotonPacket: insanity occured in extra diffuse photon &
                        & nuP assignment (nuP,gp,activeP)", nuP, gp,grid(gP)%active(xP(igpn),yP(igpn),zP(igpn))
                   print*, "difSpectrumProbDen: ", inSpectrumProbDen(0,:)
                   stop
                end if

                positionLoc%x = grid(gP)%xAxis(difSource(1))
                positionLoc%y = grid(gP)%yAxis(difSource(2))
                positionLoc%z = grid(gP)%zAxis(difSource(3))

                ! initialize the new photon packet
                orX(igPn) = difSource(1)
                orY(igPn) = difSource(2)
                orZ(igPn) = difSource(3)

                ! initialize the new photon packet
                if (grid(gP)%active(orX(igpn), orY(igpn), orZ(igpn)) < 0.) then
                   print*, "! newPhotonPacket: new packet cannot be emitted from re-mapped cell -3-"
                   print*, "chType, nuP, starPosition(iStar), .false., .true., orX,orY,orZ, gp"
                   print*, chType, nuP, starPosition(iStar), .false., .true., orX,orY,orZ, gp
                   stop
                end if
                newPhotonPacket = initPhotonPacket(nuP, positionLoc, .false., .false., orX,&
                     & orY, orZ, gP)

                newPhotonPacket%SEDtype = 0


            ! if the photon is diffuse
            case ("diffuse")

                ! check that gas is present in the grid
                if (.not.lgGas) then
                   print*, "! newPhotonPacket: diffuse packet cannot be created in a no gas grid"
                   stop
                end if

                ! check that the position has been specified
!                if (.not.present(position)) then
!                    print*, "! newPhotonPacket: position of the new diffuse &
!                         & photon packet has not been specified"
!                    stop
!                end if

                ! check that the position indeces have been specified
!                if (.not.(present(xP).and.present(yP).and.present(zP))) then
!                    print*, "! newPhotonPacket: position indeces of the new diffuse &
!                         & photon packet have not been specified"
!                    stop
!                end if
                ! check that the grid indeces have been specified
!                if (.not.(present(gP))) then
!                    print*, "! newPhotonPacket: grid index of the new diffuse &
!                         & photon packet has not been specified"
!                    stop
!                end if

                ! check that the position is not inside the inner region
                if (grid(gP)%active(xP(igPn),yP(igPn),zP(igPn))<= 0) then
                    print*, "! newPhotonPacket: position of the new diffuse &
                         & photon packet is inside the inner region", xP(igPn),yP(igPn),zP(igPn),gP
                    stop
                end if

                ! decide whether continuum or line photon
                call random_number(random)

                random = 1.-random

                if (random <= grid(gP)%totalLines(grid(gP)%active(xP(igPn),yP(igPn),zP(igPn)))) then
                   ! line photon
                   ! line photons escape so don't care which one it is unless debugging
                   if (lgDebug) then
                      call getNu2( grid(gP)%linePDF(grid(gP)%active(xP(igPn),yP(igPn),zP(igPn)),:), nuP )

                      if (nuP < 1) then
                         print*, "! newPhotonPacket: insanity occured in line photon &
                              & nuP assignment"
                         stop
                      end if

                   else

                      nuP = 0

                   end if

                   ! initialize the new photon packet
                   if (grid(gP)%active(xP(igpn), yp(igpn), zp(igpn)) < 0.) then
                      print*, "! newPhotonPacket: new packet cannot be emitted from re-mapped cell -4-"
                      print*, "chType, nuP, starPosition(iStar), .false., .true., xp,yp,zp, gp"
                      print*, chType, nuP, starPosition(iStar), .false., .true., xp,yp,zp, gp
                      stop
                   end if

                   newPhotonPacket = initPhotonPacket(nuP, position, .true., .false., xP, yP, zP, gP)
                   newPhotonPacket%SEDtype = 0

                else
                    ! continuum photon

                    call getNu2(grid(gP)%recPDF(grid(gP)%active(xP(igPn),yP(igPn),zP(igPn)),1:nbins), nuP)

                    if (nuP>nbins) then
                       print*, "! newPhotonPacket: insanity occured in diffuse photon &
                       & nuP assignment (nuP,xP,yP,zP,activeP)", nuP, xP(igPn),yP(igPn),zP(igPn),&
                       &  grid(gP)%active(xP(igPn),yP(igPn),zP(igPn))
                       print*, "recPDF: ", grid(gP)%recPDF(grid(gP)%active(xP(igPn),yP(igPn),zP(igPn)),1:nbins)
                       stop
                    end if

                    if (nuP < 1) then
                        print*, "! newPhotonPacket: insanity occured in diffuse photon &
                             & nuP assignment", nuP, xP(igPn),yP(igPn),zP(igPn), &
                             & grid(gP)%active(xP(igPn),yP(igPn),zP(igPn))
                       print*, "recPDF: ", grid(gP)%recPDF(grid(gP)%active(xP(igPn),yP(igPn),zP(igPn)),:)
                        stop
                    end if

                    ! initialize the new photon packet
                   if (grid(gP)%active(xP(igpn), yp(igpn), zp(igpn)) < 0.) then
                      print*, "! newPhotonPacket: new packet cannot be emitted from re-mapped cell -5-"
                      print*, "chType, nuP, starPosition(iStar), .false., .true., xp,yp,zp, gp"
                      print*, chType, nuP, starPosition(iStar), .false., .true., xp,yp,zp, gp
                      stop
                   end if

                    newPhotonPacket = initPhotonPacket(nuP, position, .false., .false., xP, yP, zP, gP)

                    newPhotonPacket%SEDtype = 0

                end if

            case ("dustEmi")

               ! check dust is present
               if (.not.lgDust) then
                  print*, "! newPhotonPacket: dust emitted packet cannot be created in a &
                       &no dust grid."
                  stop
               end if

               if (lgGas) then
                  print*, "! newPhotonPacket: dustEmi-type packet should be created in a &
                       & grid containing gas."
                  stop
               end if

               ! check that the position has been specified
!               if (.not.present(position)) then
!                  print*, "! newPhotonPacket: position of the new dust emitted &
!                       &photon packet has not been specified"
!                  stop
!               end if

               ! check that the position indeces have been specified
!               if (.not.(present(xP).and.present(yP).and.present(zP))) then
!                  print*, "! newPhotonPacket: position indeces of the new dust emitted &
!                       &photon packet have not been specified"
!                  stop
!               end if
               ! check that the position indeces have been specified
!               if (.not.(present(gP))) then
!                  print*, "! newPhotonPacket: grid index of the new dust emitted &
!                       &photon packet has not been specified"
!                  stop
!               end if

               ! check that the position is not inside the inner region
               if (grid(gP)%active(xP(igPn),yP(igPn),zP(igPn))<= 0) then
                  print*, "! newPhotonPacket: position of the new dust emitted &
                       &photon packet is inside the inner region", xP(igPn),yP(igPn),zP(igPn)
                  stop
               end if

               call getNu2(grid(gP)%dustPDF(grid(gP)%active(xP(igPn),yP(igPn),zP(igPn)),1:nbins), nuP)

               if (nuP>nbins) then
                   print*, "! newPhotonPacket: insanity occured in dust emitted photon &
                       &nuP assignment (iphot, nuP,xP(gP),yP(gP),zP(gP),activeP)", iphot, &
                       & nuP, xP(igPn),yP(igPn),zP(igPn), &
                       & grid(gP)%active(xP(igPn),yP(igPn),zP(igPn))
                   print*, "dustPDF: ", grid(gP)%dustPDF(grid(gP)%active(xP(igPn),yP(igPn),zP(igPn)),1:nbins)
                   print*, "grain T: ", grid(gP)%Tdust(:,0,grid(gP)%active(xP(igPn),yP(igPn),zP(igPn)))
                   stop
                end if

               if (nuP < 1) then
                  print*, "! newPhotonPacket: insanity occured in dust emitted photon &
                       &nuP assignment", nuP,xP(igPn),yP(igPn),zP(igPn),&
                       & grid(gP)%active(xP(igPn),yP(igPn),zP(igPn))
                  print*, "dustPDF: ", grid(gP)%dustPDF(grid(gP)%active(xP(igPn),yP(igPn),zP(igPn)),1:nbins)
                  print*, "grain T: ", grid(gP)%Tdust(:, 0, grid(gP)%active(xP(igPn),yP(igPn),zP(igPn)))
                  stop
               end if

               ! initialize the new photon packet
                   if (grid(gP)%active(xP(igpn), yp(igpn), zp(igpn)) < 0.) then
                      print*, "! newPhotonPacket: new packet cannot be emitted from re-mapped cell -6-"
                      print*, "chType, nuP, starPosition(iStar), .false., .true., xp,yp,zp, gp"
                      print*, chType, nuP, starPosition(iStar), .false., .true., xp,yp,zp, gp
                      stop
                   end if

               newPhotonPacket = initPhotonPacket(nuP, position, .false., .false., xP, yP, zP, gP)
               newPhotonPacket%SEDtype = 0


            ! if the photon packet type is wrong or missing
            case default

                print*, "! newPhotonPacket: wrong photon packet type - 'stellar', 'diffuse' and &
                     & dust emitted types allowed-"
                stop
            end select

        end function newPhotonPacket

        subroutine pathSegment(enPacket)
          implicit none

          ! local variables
          type(vector)                    :: vHat     ! direction vector
          type(vector)                    :: rVec     ! position vector
          type(vector)                    :: rVecloc  ! position vector

          real                            :: costheta !
          real                            :: absTau   ! optical depth
          real                            :: dlLoc    ! local displacement
          real                            :: dx, dy, dz
          real                            :: dSx, dSy, dSz
                                                      ! distances from x,y and z wall
          real                            :: dS       ! distance from nearest wall
          real                            :: dV       ! lume of this cell
          real                            :: passProb ! prob of passing the next segment
          real                            :: probSca  ! prob that the packet scatters
          real                            :: radius   ! radius
          real                            :: random   ! random number
          real                            :: tauCell  ! local tau
          real                            :: deltaEold
          real                            :: fluorescenceE ! energy of fluoresc line [ryd]


          integer                         :: nu1P, highNuP
          integer                         :: idirT,idirP ! direction cosine counters
          integer                         :: i, j, ifluoloc ! counter
          integer                         :: xP,yP,zP ! cartesian axes indeces
          integer                         :: gP       ! grid index
          integer                         :: igpp     ! grid index 1=mother 2=sub
          integer                         :: safeLimit =10000
                                                      ! safe limit for the loop
          integer                         :: ifluopacks
          integer                         :: xploc(2), yploc(2), zploc(2), gploc

          type(photon_packet), intent(inout) :: enPacket ! the energy packet

          character(len=7)                :: packetType ! stellar, diffuse, dustEmitted?

          logical                         :: lgScattered ! is the packet scattering with dust?
          logical                         :: lgReturn




          if (enPacket%iG == 1) then
             igpp = 1
          else if (enPacket%iG>1) then
             igpp = 2
          else
             print*, "! pathSegment: insane grid index"
             stop
          end if

          ! check that the input position is not outside the grid
          if ( (enPacket%iG <= 0).or.(enPacket%iG > nGrids) ) then
             print*, "! pathSegment: starting position not in any defined gridhses",&
                  & enPacket
             stop
          else if ( (enPacket%xP(igpp) <= 0).or.&
               &(enPacket%xP(igpp) > grid(enPacket%iG)%nx) ) then
             print*, "! pathSegment: starting position in x is outside the grid",&
                  & enPacket
             stop
          else if ( (enPacket%yP(igpp) <= 0).or. &
               & (enPacket%yP(igpp) > grid(enPacket%iG)%ny) ) then
             print*, "! pathSegment: starting position in y is outside the grid",&
                  & enPacket
             stop
          else if ( (enPacket%zP(igpp) <= 0).or.&
               & (enPacket%zP(igpp) > grid(enPacket%iG)%nz) ) then
             print*, "! pathSegment: starting position in z is outside the grid",&
                  & enPacket
             stop
          end if

          ! define vHat and rVec
          rVec = enPacket%position
          vHat = enPacket%direction

          ! initialize xP, yP,zP
          xP = enPacket%xP(igpp)
          yP = enPacket%yP(igpp)
          zP = enPacket%zP(igpp)
          gP = enPacket%iG

          ! initialise distance from walls
          dSx = 0.
          dSy = 0.
          dSz = 0.

          if (lg1D) then
             radius = 1.e10*sqrt((rVec%x/1.e10)*(rVec%x/1.e10) + &
                  &                               (rVec%y/1.e10)*(rVec%y/1.e10) + &
                  &                               (rVec%z/1.e10)*(rVec%z/1.e10))
             call locate(grid(1)%xAxis, radius, xP)
             if (nGrids > 1 .or. gP >1) then
                print*, " ! pathSegment: multiple grids are not allowed in a 1D simulation"
                stop
             end if
          end if



          ! initialize optical depth
          absTau = 0.

          ! get a random number
          call random_number(random)

          ! calculate the probability
          passProb = -log(1.-random)

          ! speed up photons that my be trapped
          if (lgPlaneIonization) then
             safeLimit=5000
          else
             safeLimit=500000
!             safeLimit=10000
          end if

          do i = 1, safeLimit

             if (xP <  grid(gP)%nx) then
                if (vHat%x > 0. .and. abs(rVec%x - &
                     & (grid(gP)%xAxis(xP)+grid(gP)%xAxis(xP+1))/2.)/ &
                     & ((grid(gP)%xAxis(xP)+grid(gP)%xAxis(xP+1))/2.) < 2.e-2) &
                     & xP = xP+1
             end if
             if (xP > 1 ) then
                if (vHat%x < 0. .and. abs (rVec%x  -&
                     & (grid(gP)%xAxis(xP)+grid(gP)%xAxis(xP-1))/2.)/ &
                     & ((grid(gP)%xAxis(xP)+grid(gP)%xAxis(xP-1))/2.) < 2.e-2) &
                     & xP = xP-1
             end if

             if (yP <  grid(gP)%ny) then
                if (vHat%y > 0. .and. abs (rVec%y - &
                     & (grid(gP)%yAxis(yP)+grid(gP)%yAxis(yP+1))/2.)/ &
                     & ((grid(gP)%yAxis(yP)+grid(gP)%yAxis(yP+1))/2.) < 2.e-2) &
                     & yP = yP+1
             end if
             if (yP > 1 ) then
                if (vHat%y < 0. .and. abs(rVec%y - &
                     & (grid(gP)%yAxis(yP)+grid(gP)%yAxis(yP-1))/2.)/ &
                     & ((grid(gP)%yAxis(yP)+grid(gP)%yAxis(yP-1))/2.) < 2.e-2) &
                     & yP = yP-1
             end if

             if (zP <  grid(gP)%nz) then
                if (vHat%z > 0. .and. abs(rVec%z - &
                     & (grid(gP)%zAxis(zP)+grid(gP)%zAxis(zP+1))/2.)/ &
                     & ((grid(gP)%zAxis(zP)+grid(gP)%zAxis(zP+1))/2.) < 2.e-2) &
                     & zP = zP+1
             end if
             if (zP > 1 ) then
                if (vHat%z < 0. .and. abs(rVec%z - &
                     & (grid(gP)%zAxis(zP)+grid(gP)%zAxis(zP-1))/2.)/ &
                     & ((grid(gP)%zAxis(zP)+grid(gP)%zAxis(zP-1))/2.) < 2.e-2) &
                     & zP = zP-1
             end if

             do j = 1, safeLimit

                if (xP > grid(gP)%nx .or. xP < 1 .or. &
                     & yP > grid(gP)%ny .or. yP < 1 .or. &
                     & zP > grid(gP)%nz .or. zP < 1 ) then
                   print*, "! pathSegment: insanity [gp,xp,yp,zp,j,i]", &
                        & gp, xp, yp, zp, j, i
                   stop
                end if

                if (grid(gP)%active(xP,yP,zP)<0) then

                   ! packet is entering a subgrid
                   enPacket%xP(1) = xP
                   enPacket%yP(1) = yP
                   enPacket%zP(1) = zP

                   gP = abs(grid(gP)%active(xP,yP,zP))

                   ! where is the packet in the sub-grid?

                   call locate(grid(gP)%xAxis, rVec%x, xP)
                   if (xP==0) xP = xP+1
                   if (xP< grid(gP)%nx) then
                      if (rVec%x >  (grid(gP)%xAxis(xP+1)+grid(gP)%xAxis(xP))/2.) &
                           & xP = xP + 1
                   end if

                   call locate(grid(gP)%yAxis, rVec%y, yP)
                   if (yP==0) yP=yP+1
                   if (yP< grid(gP)%ny) then
                      if (rVec%y >  (grid(gP)%yAxis(yP+1)+grid(gP)%yAxis(yP))/2.) &
                           & yP = yP + 1
                   end if

                   call locate(grid(gP)%zAxis, rVec%z, zP)
                   if (zP==0) zP=zP+1
                   if (zP< grid(gP)%nz) then
                      if (rVec%z >  (grid(gP)%zAxis(zP+1)+grid(gP)%zAxis(zP))/2.) &
                           & zP = zP + 1
                   end if

                end if

                enPacket%iG = gP

                if (gP > 1) then
                   igpp = 2
                else
                   igpp = 1
                end if

                ! find distances from all walls

                if (lgSymmetricXYZ) then
                   if ( rVec%x <= grid(1)%xAxis(1) ) then
                      if (vHat%x<0.) vHat%x = -vHat%x
                      rVec%x = grid(1)%xAxis(1)
                   end if
                   if ( rVec%y <= grid(1)%yAxis(1) ) then
                      if (vHat%y<0.) vHat%y = -vHat%y
                      rVec%y = grid(1)%yAxis(1)
                   end if
                   if ( rVec%z <= grid(1)%zAxis(1) ) then
                      if (vHat%z<0.) vHat%z = -vHat%z
                      rVec%z = grid(1)%zAxis(1)
                   end if
                end if

                if (vHat%x>0.) then
                   if (xP<grid(gP)%nx) then

                      dSx = ( (grid(gP)%xAxis(xP+1)+grid(gP)%xAxis(xP))/2.-rVec%x)/vHat%x

                      if (abs(dSx)<1.e-5) then
                         rVec%x=(grid(gP)%xAxis(xP+1)+grid(gP)%xAxis(xP))/2.
                         xP = xP+1
                      end if
                   else
                      dSx = ( grid(gP)%xAxis(grid(gP)%nx)-rVec%x)/vHat%x
                      if (abs(dSx)<1.e-5) then
                         rVec%x=grid(gP)%xAxis(grid(gP)%nx)
                         if (.not.lgPlaneIonization .and. gP==1) return
                      end if
                   end if
                else if (vHat%x<0.) then
                   if (xP>1) then
                      dSx = ( (grid(gP)%xAxis(xP)+grid(gP)%xAxis(xP-1))/2.-rVec%x)/vHat%x
                      if (abs(dSx)<1.e-5) then
                         rVec%x=(grid(gP)%xAxis(xP)+grid(gP)%xAxis(xP-1))/2.
                         xP = xP-1
                      end if
                   else
                      dSx = (grid(gP)%xAxis(1)-rVec%x)/vHat%x
                      if (abs(dSx)<1.e-5) then
                         rVec%x=grid(gP)%xAxis(1)
                      end if
                   end if
                else if (vHat%x==0.) then
                   dSx = grid(gP)%xAxis(grid(gP)%nx)
                end if

                if (.not.lg1D) then
                   if (vHat%y>0.) then
                      if (yP<grid(gP)%ny) then
                         dSy = ( (grid(gP)%yAxis(yP+1)+grid(gP)%yAxis(yP))/2.-rVec%y)/vHat%y
                         if (abs(dSy)<1.e-5) then
                            rVec%y=(grid(gP)%yAxis(yP+1)+grid(gP)%yAxis(yP))/2.
                            yP = yP+1
                         end if
                      else
                         dSy = (  grid(gP)%yAxis(grid(gP)%ny)-rVec%y)/vHat%y
                         if (abs(dSy)<1.e-5) then
                            rVec%y=grid(gP)%yAxis(grid(gP)%ny)
                            if(gP==1) return
                         end if
                      end if
                   else if (vHat%y<0.) then
                      if (yP>1) then
                         dSy = ( (grid(gP)%yAxis(yP)+grid(gP)%yAxis(yP-1))/2.-rVec%y)/vHat%y
                         if (abs(dSy)<1.e-5) then
                            rVec%y=(grid(gP)%yAxis(yP)+grid(gP)%yAxis(yP-1))/2.
                            yP = yP-1
                         end if
                      else
                         dSy = ( grid(gP)%yAxis(1)-rVec%y)/vHat%y
                         if (abs(dSy)<1.e-5) then
                            rVec%y=grid(gP)%yAxis(1)
                         end if
                      end if
                   else if (vHat%y==0.) then
                      dSy = grid(gP)%yAxis(grid(gP)%ny)
                   end if

                   if (vHat%z>0.) then
                      if (zP<grid(gP)%nz) then
                         dSz = ( (grid(gP)%zAxis(zP+1)+grid(gP)%zAxis(zP))/2.-rVec%z)/vHat%z
                         if (abs(dSz)<1.e-5) then
                            rVec%z=(grid(gP)%zAxis(zP+1)+grid(gP)%zAxis(zP))/2.
                            zP = zP+1
                         end if
                      else
                         dSz = ( grid(gP)%zAxis(grid(gP)%nz)-rVec%z)/vHat%z
                         if (abs(dSz)<1.e-5) then
                            rVec%z=grid(gP)%zAxis(grid(gP)%nz)
                            if (.not.lgPlaneIonization .and. gP==1) return
                         end if
                      end if
                   else if (vHat%z<0.) then
                      if (zP>1) then
                         dSz = ( (grid(gP)%zAxis(zP)+grid(gP)%zAxis(zP-1))/2.-rVec%z)/vHat%z
                         if (abs(dSz)<1.e-5) then
                            rVec%z=(grid(gP)%zAxis(zP)+grid(gP)%zAxis(zP-1))/2.
                            zP = zP-1
                         end if
                      else
                         dSz = ( grid(gP)%zAxis(1)-rVec%z)/vHat%z
                         if (abs(dSz)<1.e-5) then
                            rVec%z=grid(gP)%zAxis(1)
                         end if
                      end if
                   else if (vHat%z==0.) then
                      dSz = grid(gP)%zAxis(grid(gP)%nz)
                   end if

                   if (xP > grid(gP)%nx .or. xP < 1 .or. &
                        & yP > grid(gP)%ny .or. yP < 1 .or. &
                        & zP > grid(gP)%nz .or. zP < 1 ) then
                      print*, "! pathSegment: insanity -2- [gp,xp,yp,zp]", &
                           & gp, xp, yp, zp
                      stop
                   end if

                end if

                if (grid(gP)%active(xP,yP,zP)>=0) exit
             end do

             ! cater for cells on cell wall
             if ( abs(dSx)<1.e-5 ) dSx = grid(gP)%xAxis(grid(gP)%nx)
             if ( abs(dSy)<1.e-5 ) dSy = grid(gP)%yAxis(grid(gP)%ny)
             if ( abs(dSz)<1.e-5 ) dSz = grid(gP)%zAxis(grid(gP)%nz)

             ! find the nearest wall
             dSx = abs(dSx)
             dSy = abs(dSy)
             dSz = abs(dSz)

             if (dSx<=0.) then
                print*, '! pathSegment: [warning] dSx <= 0.',dSx
                print*, 'grid(gP)%xAxis ', grid(gP)%xAxis
                print*, 'gP,xP,grid(gP)%xAxis(xP), rVec%x, vHat%x'
                print*, gP,xP,grid(gP)%xAxis(xP), rVec%x, vHat%x
                dS = min(dSy, dSz)
             else if (dSy<=0.) then
                print*, '! pathSegment: [warning] dSy <= 0.', dSy
                print*, 'grid(gP)%yAxis ', grid(gP)%yAxis
                print*, 'gP,yP,grid(gP)%yAxis(yP), rVec%y, vHat%y'
                print*, gP,yP,grid(gP)%yAxis(yP), rVec%y, vHat%y
                dS = min(dSx, dSz)
             else if (dSz<=0.) then
                print*, '! pathSegment: [warning] dSz <= 0.', dSz
                print*, 'grid(gP)%zAxis ', grid(gP)%zAxis
                print*, 'gP,zP,grid(gP)%zAxis(zP), rVec%z, vHat%z'
                print*, gP,zP,grid(gP)%zAxis(zP), rVec%z, vHat%z
                dS = min(dSx, dSy)
             else
                dS = min(dSx,dSy)
                dS = min(dS, dSz)
             end if

             ! this should now never ever happen
             if (dS <= 0.) then
                print*, 'pathSegment: dS <= 0', dSx, dSy, dSz
                print*, gP, rVec
                stop
             end if

             ! calculate the optical depth to the next cell wall
             tauCell = dS*grid(gP)%opacity(grid(gP)%active(xP,yP,zP), enPacket%nuP)

             if (lg1D) then
                if (nGrids>1) then
                   print*, '! getVolumeLoc: 1D option and multiple grids options are not compatible'
                   stop
                end if

                if (xP == 1) then

                   dV = 4.*Pi* ( (grid(gP)%xAxis(xP+1)/1.e15)**3)/3.

                else if ( xP==grid(gP)%nx) then

                   dV = Pi* ( (3.*(grid(gP)%xAxis(xP)/1.e15)-(grid(gP)%xAxis(xP-1)/1.e15))**3 - &
                        & ((grid(gP)%xAxis(xP)/1.e15)+(grid(gP)%xAxis(xP-1)/1.e15))**3 ) / 6.

                else

                   dV = Pi* ( ((grid(gP)%xAxis(xP+1)/1.e15)+(grid(gP)%xAxis(xP)/1.e15))**3 - &
                        & ((grid(gP)%xAxis(xP-1)/1.e15)+(grid(gP)%xAxis(xP)/1.e15))**3 ) / 6.

                end if

                dV = dV/8.

             else


                if ( (xP>1) .and. (xP<grid(gP)%nx) ) then
                   dx = abs(grid(gP)%xAxis(xP+1)-grid(gP)%xAxis(xP-1))/2.
                else if ( xP==1 ) then
                   if (lgSymmetricXYZ .or. gP>1 .or. lgPlaneIonization) then
                      dx = abs(grid(gP)%xAxis(xP+1)-grid(gP)%xAxis(xP))/2.
                   else
                      dx = abs(grid(gP)%xAxis(xP+1)-grid(gP)%xAxis(xP))
                   end if
                else if ( xP==grid(gP)%nx ) then
                   if (gP > 1 .or. lgPlaneIonization) then
                      dx = abs(grid(gP)%xAxis(xP)  -grid(gP)%xAxis(xP-1))/2.
                   else
                      dx = abs(grid(gP)%xAxis(xP)  -grid(gP)%xAxis(xP-1))
                   end if
                end if

                if ( (yP>1) .and. (yP<grid(gP)%ny) ) then
                   dy = abs(grid(gP)%yAxis(yP+1)-grid(gP)%yAxis(yP-1))/2.
                else if ( yP==1 ) then
                   if (lgSymmetricXYZ .or. gP>1 .or. lgPlaneIonization) then
                      dy = abs(grid(gP)%yAxis(yP+1)-grid(gP)%yAxis(yP))/2.
                   else
                      dy = abs(grid(gP)%yAxis(yP+1)-grid(gP)%yAxis(yP))
                   end if
                else if ( yP==grid(gP)%ny ) then
                   if (gP > 1 .or. lgPlaneIonization) then
                      dy = abs(grid(gP)%yAxis(yP)  -grid(gP)%yAxis(yP-1))/2.
                   else
                      dy = abs(grid(gP)%yAxis(yP)  -grid(gP)%yAxis(yP-1))
                   end if
                end if

                if ( (zP>1) .and. (zP<grid(gP)%nz) ) then
                   dz = abs(grid(gP)%zAxis(zP+1)-grid(gP)%zAxis(zP-1))/2.
                 else if ( zP==1 ) then
                   if (lgSymmetricXYZ .or. gP>1 .or. lgPlaneIonization) then
                      dz = abs(grid(gP)%zAxis(zP+1)-grid(gP)%zAxis(zP))/2.
                   else
                      dz = abs(grid(gP)%zAxis(zP+1)-grid(gP)%zAxis(zP))
                   end if
                else if ( zP==grid(gP)%nz ) then
                   if (gP>1 .or. lgPlaneIonization) then
                      dz = abs(grid(gP)%zAxis(zP)-grid(gP)%zAxis(zP-1))/2.
                   else
                      dz = abs(grid(gP)%zAxis(zP)-grid(gP)%zAxis(zP-1))
                   end if
                end if

                dx = dx/1.e15
                dy = dy/1.e15
                dz = dz/1.e15


                ! calculate the volume
                dV = dx*dy*dz

             end if

             if (rVec%x == 0. .and. rVec%y == 0. .and. rVec%z == 0.) then
                costheta = 1.
             else
                costheta = ((rVec%x*vHat%x)+(rVec%y*vHat%y)+(rVec%z*vHat%z))/&
                     & (1.e10*sqrt( ((rVec%x)/1.e10)*((rVec%x)/1.e10)+&
                     & ((rVec%y)/1.e10)*((rVec%y)/1.e10)+ ((rVec%z)/1.e10)*((rVec%z)/1.e10)))

                if (isnan(costheta) .or. costheta > 1.1) then
                   print*, costheta
                   print*, (1.e10*sqrt( ((rVec%x)/1.e10)*((rVec%x)/1.e10)+&
                        & ((rVec%y)/1.e10)*((rVec%y)/1.e10)+ ((rVec%z)/1.e10)*((rVec%z)/1.e10)))
                   print*, rVec
                   stop
                end if
             end if

             ! check if the paintckets interacts within this cell
             if ((absTau+tauCell > passProb) .and. (grid(gP)%active(xP,yP,zP)>0)) then

                ! packet interacts

                ! calculate where within this cell the packet is absorbed
                dlLoc = (passProb-absTau)/grid(gP)%opacity(grid(gP)%active(xP,yP,zP), enPacket%nuP)

                ! update packet's position
                rVec = rVec + dlLoc*vHat

                if (lgSymmetricXYZ .and. gP==1) then
                   if ( rVec%x <= grid(gP)%xAxis(1) ) then
                      if (vHat%x<0.) vHat%x = -vHat%x
                      rVec%x = grid(gP)%xAxis(1)
                   end if
                   if ( rVec%y <= grid(gP)%yAxis(1) ) then
                      if (vHat%y<0.) vHat%y = -vHat%y
                      rVec%y = grid(gP)%yAxis(1)
                   end if
                   if ( rVec%z <= grid(gP)%zAxis(1) ) then
                      if (vHat%z<0.) vHat%z = -vHat%z
                      rVec%z = grid(gP)%zAxis(1)
                   end if
                end if

                ! add contribution of the packet to the radiation field
                if (enPacket%lgStellar) then
                   grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                        grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dlLoc*deltaEUsed / dV
                   grid(gP)%H(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                        grid(gP)%H(grid(gP)%active(xP,yP,zP),enPacket%nuP) + &
                        &(dlLoc*deltaEUsed)*costheta / dV

                else ! if the energy packet is diffuse
                   if (lgDebug) then
                      grid(gP)%Jdif(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                           & grid(gP)%Jdif(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dlLoc*deltaEUsed / dV
                      grid(gP)%H(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                           grid(gP)%H(grid(gP)%active(xP,yP,zP),enPacket%nuP) + &
                           &(dlLoc*deltaEUsed)*costheta / dV
                   else
                      grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                           & grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dlLoc*deltaEUsed / dV
                      grid(gP)%H(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                           grid(gP)%H(grid(gP)%active(xP,yP,zP),enPacket%nuP) + &
                           & (dlLoc*deltaEUsed)*costheta / dV
                   end if
                end if

                ! check if the position within the cell is still within the outer radius
                if ( sqrt( (rvec%x/1.e10)**2 + (rvec%y/1.e10)**2 + (rvec%z/1.e10)**2)*1.e10 >= R_out &
                     & .and. R_out > 0.) then

                   ! the packet escapes without further interaction

                   idirT = int(acos(enPacket%direction%z)/dTheta)+1
                   if (idirT>totangleBinsTheta) then
                      idirT=totangleBinsTheta
                   end if
                   if (idirT<1 .or. idirT>totAngleBinsTheta) then
                      print*, '! pathSegment: error in theta direction cosine assignment',&
                           &  idirT, enPacket, dTheta, totAngleBinsTheta
                      stop
                   end if

                   if (enPacket%direction%x<1.e-35) then
                      idirP = 0
                   else
                      idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
                   end if
                   if (idirP<0) idirP=totAngleBinsPhi+idirP
                   idirP=idirP+1
                   if (idirP>totangleBinsPhi) then
                      idirP=totangleBinsPhi
                   end if

                   if (idirP<1 .or. idirP>totAngleBinsPhi) then
                      print*, '! pathSegment: error in Phi direction cosine assignment',&
                           &  idirP, enPacket, dPhi, totAngleBinsPhi
                      stop
                   end if


                   if (nAngleBins>0) then
                      if ( (totangleBinsPhi>1 .and. viewPointPtheta(idirT) == &
                           &viewPointPphi(idirP)) .or. &
                           & (totangleBinsPhi==1 .and. viewPointPtheta(idirT)>0) ) then
                         grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,&
                              & viewPointPtheta(idirT)) = &
                              &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                              & enPacket%nuP,viewPointPtheta(idirT)) +  deltaEUsed

                         grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,0)=&
                              &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                              & enPacket%nuP,0) +  deltaEUsed

                         if (lgSeparateSED) then
                            grid(enPacket%origin(1))%escapedPacketsComponents&
                                 &(enPacket%origin(2), enPacket%nuP,&
                              & viewPointPtheta(idirT),enPacket%SEDtype) = &
                              &grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                              & enPacket%nuP,viewPointPtheta(idirT),enPacket%SEDtype) +  deltaEUsed

                            grid(enPacket%origin(1))%escapedPacketsComponents&
                                 &(enPacket%origin(2), enPacket%nuP,0,enPacket%SEDtype)=&
                                 &grid(enPacket%origin(1))%escapedPacketsComponents&
                                 &(enPacket%origin(2), &
                                 & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed

                         endif

                      else
                         grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2),enPacket%nuP,0) = &
                              & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                              & enPacket%nuP,0) +  deltaEUsed
                         if (lgSeparateSED) then
                            grid(enPacket%origin(1))%escapedPacketsComponents&
                                 &(enPacket%origin(2),enPacket%nuP,0,enPacket%SEDtype) = &
                                 & grid(enPacket%origin(1))%escapedPacketsComponents&
                                 &(enPacket%origin(2),  enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                         endif
                      end if
                   else
                      grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                           enPacket%nuP,0) = &
                           & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                           & enPacket%nuP,0) +  deltaEUsed
                      if (lgSeparateSED) then
                         grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                              enPacket%nuP,0,enPacket%SEDtype) = &
                              & grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                              & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                      end if

                   end if


                   return
                end if

                ! check if the packet is absorbed or scattered
                if (lgDust) then

                   probSca = grid(gP)%scaOpac(grid(gP)%active(xP,yP,zP),enPacket%nuP)/&
                        & (grid(gP)%opacity(grid(gP)%active(xP,yP,zP),enPacket%nuP))

                   call random_number(random)

                   random = 1.-random

                   if (random > probSca) then
                      lgScattered = .false.
                   else if (random <= probSca) then
                      lgScattered = .true.
                   else
                      print*, '! pathSegment: insanity occured and scattering/absorption &
                           & decision stage.'
                      stop
                   end if

                   if (.not. lgScattered) then

                      absInt = absInt + 1.

                      if (.not.lgGas) then

                         ! packet is absobed by the dust
                         packetType = "dustEmi"
                         exit

                      else
                         ! HERE YOU MUST ADD COMPTON SCATTERING STUFF!

                         ! packet is absobed by the dust+gas
                         packetType = "diffuse"
                         exit

                      end if

                   else

                      scaInt = scaInt + 1.

!                      if (lgMultiDustChemistry) then
!                         do nS = 1, nSpeciesPart(grid(gP)%dustAbunIndex(grid(gP)%active(xP,yP,zP)))
!                            if (grainabun(grid(gP)%dustAbunIndex(grid(gP)%active(xP,yP,zP)),nS)>0. &
!                                 &.and. grid(gP)%Tdust(nS, 0, &
!                                 & grid(gP)%active(xP,yP,zP))<TdustSublime(dustComPoint(&
!                                 &grid(gP)%dustAbunIndex(grid(gP)%active(xP,yP,zP)))-1+nS)) exit
!                         end do
!                         if (nS>nSpeciesPart(grid(gP)%dustAbunIndex(grid(gP)%active(xP,yP,zP)))) then
!                            print*, "! pathSegment: packet scatters with dust at position where all &
!                                 &grains have sublimed -1-."
!                            print*, xP,yP,zP, grid(gP)%active(xP,yP,zP), tauCell, absTau, passProb
!                            stop
!                         end if
!                      else
!                         do nS = 1, nSpeciesPart(1)
!                            if (grainabun(1,nS)>0. &
!                                 &.and. grid(gP)%Tdust(nS, 0, &
!                                 & grid(gP)%active(xP,yP,zP))<TdustSublime(dustComPoint(&
!                                 &1)-1+nS)) exit
!                         end do
!                         if (nS>nSpeciesPart(1)) then
!                            print*, "! pathSegment: packet scatters with dust at position where all &
!                                 &grains have sublimed -2-."
!                            print*, xP,yP,zP, grid(gP)%active(xP,yP,zP), tauCell, absTau, passProb
!                            stop
!                         end if
!                      end if

!                      do nS = 1, nSpecies
!                         if (grainabun(nS)>0. .and. grid(gP)%Tdust(nS, 0, &
!                              & grid(gP)%active(xP,yP,zP))<TdustSublime(nS)) exit
!                      end do
!                      if (nS>nSpecies) then
!                         print*, "! pathSegment: packet scatters with dust at position where all &
!                              &grains have sublimed."
!                         print*, xP,yP,zP, grid(gP)%active(xP,yP,zP), tauCell, absTau, passProb
!                         stop
!                      end if

                      ! packet is scattered by the grain

                      ! calculate new direction
                      ! for now assume scattering is isotropic, when phase
                      ! function is introduced the following must be changed

                      enPacket%xP(igpp) = xP
                      enPacket%yP(igpp) = yP
                      enPacket%zP(igpp) = zP

                      if (grid(gP)%active(enPacket%xp(igpp), enPacket%yp(igpp), enPacket%zp(igpp)) < 0.) then
                         print*, "! pathSegment: new packet cannot be emitted from re-mapped cell -1-"
                         print*, "nuP, starPosition(iStar), .false., .true., xp,yp,zp, gp"
                         print*, nuP, starPosition(iStar), .false., .true.,  xp,yp,zp, gp
                         stop
                      end if

                      enPacket = initPhotonPacket(enPacket%nuP, rVec, .false., .false., enPacket%xP(1:2), &
                           & enPacket%yP(1:2), enPacket%zP(1:2), gP)


                      vHat%x = enPacket%direction%x
                      vHat%y = enPacket%direction%y
                      vHat%z = enPacket%direction%z

                      ! initialize optical depth
                      absTau = 0.

                      ! get a random number
                      call random_number(random)

                      ! calculate the probability
                      passProb = -log(1.-random)

                   end if

                else


                   ! check if the packet is absorbed or scattered
                   ! NOTE : we are ignoring dust

                   if (lgCompton) then
                      probSca = KNsigmaT(enPacket%nuP)*&
                           &grid(gP)%Ne(grid(gP)%active(xp,yp,zp))/&
                           & (grid(gp)%opacity(grid(gp)%active(xp,yp,zp),&
                           &enPacket%nuP))
                   else
                      probSca = 0.
                   end if

                   call random_number(random)

                   random = 1.-random

                   if (random > probSca) then
                      lgScattered = .false.
                   else if (random <= probSca) then
                      lgScattered = .true.
                   else
                      print*, '! pathSegment: insanity occured and scattering/absorption &
                           & decision stage.'
                      stop
                   end if

                   if (lgScattered) then

                      scaInt = scaInt + 1.
                      ! packet is compton scattered
                      ! calculate new direction
                      ! for now assume scattering is isotropic, when phase
                      ! function is introduced the following must be changed

                      enPacket%xP(igpp) = xP
                      enPacket%yP(igpp) = yP
                      enPacket%zP(igpp) = zP

                      if (grid(gP)%active(enPacket%xp(igpp), enPacket%yp(igpp), enPacket%zp(igpp)) < 0.) then
                         print*, "! pathSegment: new packet cannot be emitted from re-mapped cell -11-"
                         print*, "nuP, .false., .true., xp,yp,zp, gp"
                         print*, nuP, .false., .true.,  xp,yp,zp, gp
                         stop
                      end if

                      call comptonScatter(enPacket)

                      enPacket%origin(1) = gP
                      enPacket%origin(2) = grid(gP)%active(enPacket%xp(igpp), enPacket%yp(igpp), enPacket%zp(igpp))


                      vHat%x = enPacket%direction%x
                      vHat%y = enPacket%direction%y
                      vHat%z = enPacket%direction%z

                      ! initialize optical depth
                      absTau = 0.

                      ! get a random number
                      call random_number(random)

                      ! calculate the probability
                      passProb = -log(1.-random)

                   else

                      absInt = absInt + 1.

                      if (.not.lgGas) then
                         print*, "! pathSegment: Insanity occured - no gas present when no dust interaction"
                         stop
                      end if


                      if (lgFluorescence) then
                      ! transfer fluorescent photons

                         enPacket%xP(igpp) = xp
                         enPacket%yP(igpp) = yp
                         enPacket%zP(igpp) = zp

                         do ifluoloc = 1, nfluo
                            if (fluorescenceLabel(iFluoloc) == 'FeKa') then
                               nu1P    = nu1PFeKa
                               highNuP = highNuPFeKa
                               fluorescenceE = 470.388
                            else if (fluorescenceLabel(iFluoloc) == 'FeL1') then
                               nu1P    = nu1PFeL1
                               highNuP = highNuPFeL1
                               fluorescenceE = 58.2353
                            else if (fluorescenceLabel(iFluoloc) == 'FeL2') then
                               nu1P    = nu1PFeL2
                               highNuP = highNuPFeL2
                               fluorescenceE = 51.8382
                            else if (fluorescenceLabel(iFluoloc) == 'CKa') then

                               nu1P    = nu1PCKa
                               highNuP = highNuPCKa
                               fluorescenceE = 20.3676

                            else if (fluorescenceLabel(iFluoloc) == 'NKa') then
                               nu1P    = nu1PNKa
                               highNuP = highNuPNKa
                               fluorescenceE = 28.8529
                            else if (fluorescenceLabel(iFluoloc) == 'OKa') then
                               nu1P    = nu1POKa
                               highNuP = highNuPOKa
                               fluorescenceE = 38.5956
                            else if (fluorescenceLabel(iFluoloc) == 'NeKa') then
                               nu1P    = nu1PNeKa
                               highNuP = highNuPNeKa
                               fluorescenceE = 62.3971
                            else if (fluorescenceLabel(iFluoloc) == 'MgKa') then
                               nu1P    = nu1PMgKa
                               highNuP = highNuPMgKa
                               fluorescenceE = 92.1765
                            else if (fluorescenceLabel(iFluoloc) == 'AlKa') then
                               nu1P    = nu1PAlKa
                               highNuP = highNuPAlKa
                               fluorescenceE = 109.287
                            else if (fluorescenceLabel(iFluoloc) == 'SiKa') then
                               nu1P    = nu1PSiKa
                               highNuP = highNuPSiKa
                               fluorescenceE = 127.904
                            else if (fluorescenceLabel(iFluoloc) == 'SKa') then
                               nu1P    = nu1PSKa
                               highNuP = highNuPSKa
                               fluorescenceE = 169.632
                            else if (fluorescenceLabel(iFluoloc) == 'ArKa') then
                               nu1P    = nu1PArKa
                               highNuP = highNuPArKa
                               fluorescenceE = 217.397
                            else if (fluorescenceLabel(iFluoloc) == 'CaKa') then
                               nu1P    = nu1PCaKa
                               highNuP = highNuPCaKa
                               fluorescenceE = 271.324
                            else if (fluorescenceLabel(iFluoloc) == 'TiKa') then
                               nu1P    = nu1PTiKa
                               highNuP = highNuPTiKa
                               fluorescenceE = 549.860
                            else
                               print*, '! pathSegment: label not recognised', &
                                    & fluorescenceLabel(iFluoloc)
                               stop
                            end if

                            if (enPacket%nuP >= nu1P .and. enPacket%nuP <= highNuP) then
                               ! calculate the local fluorescence luminosity
                               deltaEOld = deltaEUsed

                               call getFLuorescenceL(deltaEUsed,enPacket%nuP,ifluoloc,&
                                    & gP,xp,yp,zp)

                               fluorescenceTot(ifluoloc) = fluorescenceTot(ifluoloc)+&
                                    &deltaEUsed/(Ryd2erg*fluorescenceE)

                               if (deltaEUsed > deltaEold) then
                                  print*, '! pathSegment:deltaEUsed > deltaEold', &
                                       &deltaEUsed, deltaEold
                                  stop
                               end if

                               deltaEUsed = deltaEUsed/20.
                               if (deltaEUsed< 0.) then
                                  print*, "! pathSegment: deltaEUsed < 0! - 1"
                                  stop
                               end if

                               xploc = enPacket%xP
                               yploc = enPacket%yP
                               zploc = enPacket%zP
                               rvecloc = rVec
                               gploc=gp
                               do ifluopacks = 1, 20
                                  call fluorescencePacketRun(fType=fluorescenceLabel(iFluoloc),&
                                       &  rVecin= rVecloc, &
                                       &xpin=xploc, ypin=yploc, zpin=zploc, gpin=gploc, ifl = ifluoloc)
                               end do
                               deltaEUsed = deltaEUsed*20.

                               deltaEUsed = deltaEOld-deltaEUsed

                               if (deltaEUsed< 0.) then
                                  print*, "! pathSegment: deltaEUsed < 0! - 2", deltaEOld,deltaEUsed
                                  stop
                               end if

                            end if
                         end do

                      end if

                      ! packet interacts with gas
                      packetType = "diffuse"
                      exit

                   end if
                end if

             else

                ! the packet is not absorbed within this cell
                ! add contribution of the packet to the radiation field

                if (enPacket%lgStellar) then
                   grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                        grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dS*deltaEUsed / dV

                   grid(gP)%H(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                        grid(gP)%H(grid(gP)%active(xP,yP,zP),enPacket%nuP) + &
                        &(dS*deltaEUsed)*costheta / dV

                else ! if the energy packet is diffuse
                   if (lgDebug) then
                      grid(gP)%Jdif(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                           & grid(gP)%Jdif(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dS*deltaEUsed / dV

                      grid(gP)%H(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                           grid(gP)%H(grid(gP)%active(xP,yP,zP),enPacket%nuP) + &
                           &(dS*deltaEUsed)*costheta / dV
                   else
                      grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                           & grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dS*deltaEUsed / dV

                      grid(gP)%H(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                           grid(gP)%H(grid(gP)%active(xP,yP,zP),enPacket%nuP) + &
                           & (dS*deltaEUsed)*costheta / dV
                   end if
                end if

                ! update absTau
                absTau = absTau+tauCell

                ! update packet's position
                rVec = rVec+dS*vHat

                ! keep track of where you are on mother grid
                if (gP>1) then

                   if (enPacket%xP(1) <= 0 .or. &
                        & enPacket%yP(1) <= 0 .or. &
                        & enPacket%zP(1) <= 0 ) then

                      ! locate where we are at on the mother grid
                      call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                      if (enPacket%xP(1)<1) enPacket%xP(1)=1
                      if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                         if ( rvec%x >= ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                              & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                              & enPacket%xP(1) = enPacket%xP(1)+1
                      end if

                      call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                      if (enPacket%yP(1)<1) enPacket%yP(1)=1
                      if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                         if ( rvec%y >= ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                              & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                              & enPacket%yP(1) = enPacket%yP(1)+1
                      end if

                      call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                      if (enPacket%zP(1)<1) enPacket%zP(1)=1
                      if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                         if ( rvec%z >= ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                              & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                              & enPacket%zP(1) = enPacket%zP(1)+1
                      end if


                   else

                      if (vHat%x>0.) then
                         if ( enPacket%xP(1) < grid(grid(gP)%motherP)%nx ) then
                            if ( rVec%x > (grid(grid(gP)%motherP)%xAxis(enPacket%xP(1))+&
                                 & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1))/2. ) then
                               enPacket%xP(1) = enPacket%xP(1)+1
                            end if
                         else
                            if ( rVec%x > grid(grid(gP)%motherP)%xAxis(enPacket%xP(1))) then
                            end if
                         end if
                      else
                         if ( enPacket%xP(1) > 1 ) then
                            if ( rVec%x <= (grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)-1)+&
                                 & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)))/2. ) then
                               enPacket%xP(1) = enPacket%xP(1)-1
                            end if
                         else
                            if (rVec%x < grid(grid(gP)%motherP)%xAxis(1)) then
!                            print*, '! pathSegment: insanity occured at mother grid transfer (x axis-)',&
!                                 & rVec%x, gP, grid(gP)%motherP
!                            stop
                            end if
                         end if
                      end if
                      if (vHat%y>0.) then
                         if (  enPacket%yP(1) < grid(grid(gP)%motherP)%ny ) then
                            if ( rVec%y > (grid(grid(gP)%motherP)%yAxis( enPacket%yP(1))+&
                                 & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1))/2. ) then
                               enPacket%yP(1) =  enPacket%yP(1)+1
                            end if
                         else
                            if ( rVec%y > grid(grid(gP)%motherP)%yAxis( enPacket%yP(1))) then
!                            print*, '! pathSegment: insanity occured at mother grid transfer (y axis +)',&
!                                 & rVec%y, gP, grid(gP)%motherP
!                            stop
                            end if
                         end if
                      else
                         if (  enPacket%yP(1) > 1 ) then
                            if ( rVec%y <= (grid(grid(gP)%motherP)%yAxis( enPacket%yP(1)-1)+&
                                 & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)))/2. ) then
                               enPacket%yP(1) =  enPacket%yP(1)-1
                            end if
                         else
                            if (rVec%y < grid(grid(gP)%motherP)%yAxis(1)) then
!                            print*, '! pathSegment: insanity occured at mother grid transfer (y axis -)', &
!                                 & rVec%y, gP, grid(gP)%motherP
!                            stop
                            end if
                         end if
                      end if
                      if (vHat%z>0.) then
                         if (  enPacket%zP(1) < grid(grid(gP)%motherP)%nz ) then
                            if ( rVec%z > (grid(grid(gP)%motherP)%zAxis( enPacket%zP(1))+&
                                 & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1))/2. ) then
                               enPacket%zP(1) =  enPacket%zP(1)+1
                            end if
                         else
                            if ( rVec%z > grid(grid(gP)%motherP)%zAxis( enPacket%zP(1))) then
                            print*, '! pathSegment: insanity occured at mother grid transfer (z axis +)', &
                                 & rVec%z, gP, grid(gP)%motherP
                            stop
                            end if
                         end if
                      else
                         if (  enPacket%zP(1) > 1 ) then
                            if ( rVec%z <= (grid(grid(gP)%motherP)%zAxis( enPacket%zP(1)-1)+&
                                 & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)))/2. ) then
                               enPacket%zP(1) =  enPacket%zP(1)-1
                            end if
                         else
                            if (rVec%z < grid(grid(gP)%motherP)%zAxis(1)) then
!                        print*, '! pathSegment: insanity occured at mother grid transfer (z axis -)', &
!                             & rVec%z, gP, grid(gP)%motherP
!                        stop
                            end if
                         end if
                      end if

                   end if
                end if

                if (.not.lg1D) then
                   if ( (dS == dSx) .and. (vHat%x > 0.)  ) then
                      xP = xP+1
                   else if ( (dS == dSx) .and. (vHat%x < 0.) ) then
                      xP = xP-1
                   else if ( (dS == dSy) .and. (vHat%y > 0.) ) then
                      yP = yP+1
                   else if ( (dS == dSy) .and. (vHat%y < 0.) ) then
                      yP = yP-1
                   else if ( (dS == dSz) .and. (vHat%z > 0.) ) then
                      zP = zP+1
                   else if ( (dS == dSz) .and. (vHat%z < 0.) ) then
                      zP = zP-1
                   else
                      print*, '! pathSegment: insanity occured in dS assignement &
                           & [dS,dSx,dSy,dSz,vHat]', dS,dSx,dSy,dSz,vHat
                   end if
                else
                   radius = 1.e10*sqrt((rVec%x/1.e10)*(rVec%x/1.e10) + &
                        & (rVec%y/1.e10)*(rVec%y/1.e10) + &
                        & (rVec%z/1.e10)*(rVec%z/1.e10))
                   call locate(grid(gP)%xAxis, radius , xP)

                end if

                ! be 6/6/06
                if(.not.lgPlaneIonization.and..not.lgSymmetricXYZ) then
                   lgReturn=.false.

                   if ( rVec%y < grid(gP)%yAxis(1) .or. yP<1) then

                      ! the energy packet escapes this grid
                      if (gP==1) then
                         yP=1
                         lgReturn=.true.
                      else if (gP>1) then

                         ! locate where we are at on the mother grid
                         call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                         if (enPacket%xP(1)<1) enPacket%xP(1)=1
                         if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                            if ( rvec%x >= ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                                 & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                                 & enPacket%xP(1) = enPacket%xP(1)+1
                         end if

                         call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                         if (enPacket%yP(1)<1) enPacket%yP(1)=1
                         if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                            if ( rvec%y >= ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                                 & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                                 & enPacket%yP(1) = enPacket%yP(1)+1
                         end if

                         call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                         if (enPacket%zP(1)<1) enPacket%zP(1)=1
                         if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                            if ( rvec%z >= ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                                 & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                                 & enPacket%zP(1) = enPacket%zP(1)+1
                         end if


                         xP   = enPacket%xP(1)
                         yP   = enPacket%yP(1)
                         zP   = enPacket%zP(1)
                         gP   = 1
                         igpp = 1

                      else
                         print*, '! pathSegment: insanity occured - invalid gP', gP
                         stop
                      end if

                   end if

                   if (rVec%y > grid(gP)%yAxis(grid(gP)%ny) .or. yP>grid(gP)%ny) then

                      if (gP==1) then
                         ! the energy packet escapes
                         yP = grid(gP)%ny
                         lgReturn=.true.
                      else if (gP>1) then

                         ! locate where we are at on the mother grid
                         call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                         if (enPacket%xP(1)<1) enPacket%xP(1)=1
                         if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                            if ( rvec%x >= ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                                 & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                                 & enPacket%xP(1) = enPacket%xP(1)+1
                         end if

                         call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                         if (enPacket%yP(1)<1) enPacket%yP(1)=1
                         if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                            if ( rvec%y >= ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                                 & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                                 & enPacket%yP(1) = enPacket%yP(1)+1
                         end if

                         call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                         if (enPacket%zP(1)<1) enPacket%zP(1)=1
                         if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                            if ( rvec%z >= ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                                 & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                                 & enPacket%zP(1) = enPacket%zP(1)+1
                         end if

                         xP   =  enPacket%xP(1)
                         yP   =  enPacket%yP(1)
                         zP   =  enPacket%zP(1)
                         gP   =  1
                         igpp = 1
                      else
                         print*, '! pathSegment: insanity occured - invalid gP', gP
                         stop
                      end if

                   end if

                   if ( (rVec%x <= grid(gP)%xAxis(1) .or. xP<1) .and. gP==1) then
                      xP=1
                      lgReturn=.true.
                   end if


                   if ( (rVec%x <= grid(gP)%xAxis(1) .or. xP<1) &
                    & .and. gP>1) then

                      ! locate where we are at on the mother grid
                      call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                      if (enPacket%xP(1)<1) enPacket%xP(1)=1
                      if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                         if ( rvec%x >= ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                              & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                              & enPacket%xP(1) = enPacket%xP(1)+1
                      end if

                      call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                      if (enPacket%yP(1)<1) enPacket%yP(1)=1
                      if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                         if ( rvec%y >= ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                              & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                              & enPacket%yP(1) = enPacket%yP(1)+1
                      end if

                      call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                      if (enPacket%zP(1)<1) enPacket%zP(1)=1
                      if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                         if ( rvec%z >= ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                              & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                              & enPacket%zP(1) = enPacket%zP(1)+1
                      end if


                      xP   = enPacket%xP(1)
                      yP   = enPacket%yP(1)
                      zP   = enPacket%zP(1)
                      gP   = 1
                      igpp = 1

                   end if


                   if ( (rVec%x >=  grid(gP)%xAxis(grid(gP)%nx) &
                        & .or. xP>grid(gP)%nx) .and. gP==1 )then
                      xP = grid(gP)%nx
                      lgReturn=.true.

                   end if

                   if ( (rVec%x >=  grid(gP)%xAxis(grid(gP)%nx)&
                        & .or. xP>grid(gP)%nx) .and.  gP>1) then

                      ! locate where we are at on the mother grid
                      call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                      if (enPacket%xP(1)<1) enPacket%xP(1)=1
                      if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                         if ( rvec%x >= ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                              & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                              & enPacket%xP(1) = enPacket%xP(1)+1
                      end if

                      call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                      if (enPacket%yP(1)<1) enPacket%yP(1)=1
                      if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                         if ( rvec%y >= ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                              & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                              & enPacket%yP(1) = enPacket%yP(1)+1
                      end if

                      call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                      if (enPacket%zP(1)<1) enPacket%zP(1)=1
                      if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                         if ( rvec%z >= ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                              & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                              & enPacket%zP(1) = enPacket%zP(1)+1
                      end if

                      xP   = enPacket%xP(1)
                      yP   = enPacket%yP(1)
                      zP   = enPacket%zP(1)
                      gP   = 1
                      igpp = 1
                   end if

                   if ( (rVec%z <= grid(gP)%zAxis(1) .or.zP<1) &
                        & .and. gP==1) then
                      zP=1
                      lgReturn=.true.

                   end if

                   if ( (rVec%z <= grid(gP)%zAxis(1) &
                        & .or.zP<1) .and. gP>1) then

                      ! locate where we are at on the mother grid
                      call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                      if (enPacket%xP(1)<1) enPacket%xP(1)=1
                      if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                         if ( rvec%x >= ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                              & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                              & enPacket%xP(1) = enPacket%xP(1)+1
                      end if

                      call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                      if (enPacket%yP(1)<1) enPacket%yP(1)=1
                      if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                         if ( rvec%y >= ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                              & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                              & enPacket%yP(1) = enPacket%yP(1)+1
                      end if

                      call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                      if (enPacket%zP(1)<1) enPacket%zP(1)=1
                      if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                         if ( rvec%z >= ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                              & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                              & enPacket%zP(1) = enPacket%zP(1)+1
                      end if

                      xP   = enPacket%xP(1)
                      yP   = enPacket%yP(1)
                      zP   = enPacket%zP(1)
                      gP   = 1
                      igpp = 1

                   end if

                   if ( (rVec%z >=  grid(gP)%zAxis(grid(gP)%nz) &
                        & .or. zP>grid(gP)%nz) &
                        & .and. gP==1) then

                      zP = grid(gP)%nz
                      lgReturn=.true.

                   end if

                   if ((rVec%z >=  grid(gP)%zAxis(grid(gP)%nz)&
                        & .or. zP>grid(gP)%nz) .and. gP>1) then

                      ! locate where we are at on the mother grid
                      call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                      if (enPacket%xP(1)<1) enPacket%xP(1)=1
                      if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                         if ( rvec%x >= ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                              & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                              & enPacket%xP(1) = enPacket%xP(1)+1
                      end if

                      call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                      if (enPacket%yP(1)<1) enPacket%yP(1)=1
                      if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                         if ( rvec%y >= ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                              & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                              & enPacket%yP(1) = enPacket%yP(1)+1
                      end if

                      call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                      if (enPacket%zP(1)<1) enPacket%zP(1)=1
                      if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                         if ( rvec%z >= ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                              & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                              & enPacket%zP(1) = enPacket%zP(1)+1
                      end if

                      xP   = enPacket%xP(1)
                      yP   = enPacket%yP(1)
                      zP   = enPacket%zP(1)
                      gP   = 1
                      igpp = 1

                   end if

                   if (lgReturn) then

                      ! the packet escapes without further interaction

                      idirT = int(acos(enPacket%direction%z)/dTheta)+1
                      if (idirT>totangleBinsTheta) then
                         idirT=totangleBinsTheta
                      end if
                      if (idirT<1 .or. idirT>totAngleBinsTheta) then
                         print*, '! pathSegment: error in theta direction cosine assignment',&
                              &  idirT, enPacket, dTheta, totAngleBinsTheta
                         stop
                      end if


                      if (enPacket%direction%x<1.e-35) then
                         idirP = 0
                      else
                         idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
                      end if
                      if (idirP<0) idirP=totAngleBinsPhi+idirP
                      idirP=idirP+1

                      if (idirP>totangleBinsPhi) then
                         idirP=totangleBinsPhi
                      end if
                      if (idirP<1 .or. idirP>totAngleBinsPhi) then
                         print*, '! pathSegment: error in phi direction cosine assignment',&
                              &  idirP, enPacket, dPhi, totAngleBinsPhi
                         stop
                      end if

                      if (nAngleBins>0) then
                         if ( (totangleBinsPhi>1 .and. viewPointPtheta(idirT) == viewPointPphi(idirP)) .or. &
                              & (totangleBinsPhi==1 .and. viewPointPtheta(idirT)>0) ) then
                            grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,&
                                 & viewPointPtheta(idirT)) = &
                                 &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                 & enPacket%nuP,viewPointPtheta(idirT)) +  deltaEUsed
                            grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,0) = &
                                 &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                 & enPacket%nuP,0) +  deltaEUsed

                            if (lgSeparateSED) then
                               grid(enPacket%origin(1))%escapedPacketsComponents&
                                    &(enPacket%origin(2), enPacket%nuP,&
                                    & viewPointPtheta(idirT),enPacket%SEDtype) = &
                                    &grid(enPacket%origin(1))%escapedPacketsComponents&
                                    &(enPacket%origin(2), enPacket%nuP,&
                                    &viewPointPtheta(idirT),enPacket%SEDtype) +  deltaEUsed
                               grid(enPacket%origin(1))%escapedPacketsComponents&
                                    &(enPacket%origin(2), enPacket%nuP,0,enPacket%SEDtype) = &
                                    &grid(enPacket%origin(1))%escapedPacketsComponents&
                                    &(enPacket%origin(2), enPacket%nuP,0,enPacket%SEDtype)+deltaEUsed
                            endif

                         else
                            grid(enPacket%origin(1))%escapedPackets&
                                 &(enPacket%origin(2),enPacket%nuP,0) = &
                                 & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                 & enPacket%nuP,0) +  deltaEUsed
                            if (lgSeparateSED) then
                               grid(enPacket%origin(1))%escapedPacketsComponents&
                                    &(enPacket%origin(2),enPacket%nuP,0,enPacket%SEDtype) = &
                                    & grid(enPacket%origin(1))%escapedPacketsComponents&
                                    &(enPacket%origin(2), &
                                    & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                            endif

                         end if
                      else
                         grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                              enPacket%nuP,0) = &
                              & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                              & enPacket%nuP,0) +  deltaEUsed

                         if (lgSeparateSED) then
                            grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                 enPacket%nuP,0,enPacket%SEDtype) = &
                                 & grid(enPacket%origin(1))%escapedPacketsComponents&
                                 &(enPacket%origin(2), &
                                 & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                         endif

                      end if



                      return
                   end if

                end if

                ! end be 6/6/06


                if(lgPlaneIonization) then
                   lgReturn=.false.

                   if ( rVec%y < grid(gP)%yAxis(1) .or. yP<1) then
                      ! the energy packet escapes this grid
                      if (gP==1) then
                         yP=1
                         lgReturn=.true.
                      else if (gP>1) then

                         ! locate where we are at on the mother grid
                         call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                         if (enPacket%xP(1)<1) enPacket%xP(1)=1
                         if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                            if ( rvec%x >= ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                              & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                              & enPacket%xP(1) = enPacket%xP(1)+1
                         end if

                         call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                         if (enPacket%yP(1)<1) enPacket%yP(1)=1
                         if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                            if ( rvec%y >= ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                                 & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                                 & enPacket%yP(1) = enPacket%yP(1)+1
                         end if

                         call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                         if (enPacket%zP(1)<1) enPacket%zP(1)=1
                         if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                            if ( rvec%z >= ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                                 & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                                 & enPacket%zP(1) = enPacket%zP(1)+1
                         end if

                         xP = enPacket%xP(1)
                         yP = enPacket%yP(1)
                         zP = enPacket%zP(1)
                         gP = 1
                         igpp = 1
                      else
                         print*, '! pathSegment: insanity occured - invalid gP', gP
                         stop
                      end if

                   end if

                   if (rVec%y > grid(gP)%yAxis(grid(gP)%ny) .or. yP>grid(gP)%ny) then

                      if (gP==1) then

                         ! the energy packet escapes
                         yP = grid(gP)%ny
                         lgReturn=.true.

                      else if (gP>1) then

                         ! locate where we are at on the mother grid
                         call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                         if (enPacket%xP(1)<1) enPacket%xP(1)=1
                         if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                            if ( rvec%x >= ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                                 & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                                 & enPacket%xP(1) = enPacket%xP(1)+1
                         end if

                         call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                         if (enPacket%yP(1)<1) enPacket%yP(1)=1
                         if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                            if ( rvec%y >= ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                                 & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                                 & enPacket%yP(1) = enPacket%yP(1)+1
                         end if

                         call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                         if (enPacket%zP(1)<1) enPacket%zP(1)=1
                         if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                            if ( rvec%z >= ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                                 & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                                 & enPacket%zP(1) = enPacket%zP(1)+1
                         end if

                         xP = enPacket%xP(1)
                         yP = enPacket%yP(1)
                         zP = enPacket%zP(1)
                         gP = 1
                         igpp = 1

                      else
                         print*, '! pathSegment: insanity occured - invalid gP', gP
                         stop
                      end if

                   end if

                   if ( (rVec%x <= grid(1)%xAxis(1) .or. xP<1) ) then
                      xP=1
                      rVec%x = grid(gP)%xAxis(1)
                      vHat%x = -vHat%x

                   end if

                   if ( (rVec%x <= grid(gP)%xAxis(1) .or. xP<1) &
                        & .and. gP>1) then

                      ! locate where we are at on the mother grid
                      call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                      if (enPacket%xP(1)<1) enPacket%xP(1)=1
                      if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                         if ( rvec%x >= ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                              & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                              & enPacket%xP(1) = enPacket%xP(1)+1
                      end if

                      call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                      if (enPacket%yP(1)<1) enPacket%yP(1)=1
                      if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                         if ( rvec%y >= ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                              & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                              & enPacket%yP(1) = enPacket%yP(1)+1
                      end if

                      call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                      if (enPacket%zP(1)<1) enPacket%zP(1)=1
                      if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                         if ( rvec%z >= ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                              & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                              & enPacket%zP(1) = enPacket%zP(1)+1
                      end if

                      xP = enPacket%xP(1)
                      yP = enPacket%yP(1)
                      zP = enPacket%zP(1)
                      gP = 1
                      igpp = 1

                   end if

                   if ( (rVec%x >=  grid(1)%xAxis(grid(gP)%nx) &
                        & .or. xP>grid(gP)%nx)  )then
                      xP = grid(gP)%nx
                      rVec%x = grid(gP)%xAxis(grid(gP)%nx)
                      vHat%x = -vHat%x

                   end if

                   if ( (rVec%x >=  grid(gP)%xAxis(grid(gP)%nx)&
                        & .or. xP>grid(gP)%nx) .and.  gP>1) then

                      ! locate where we are at on the mother grid
                      call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                      if (enPacket%xP(1)<1) enPacket%xP(1)=1
                      if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                         if ( rvec%x >= ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                              & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                              & enPacket%xP(1) = enPacket%xP(1)+1
                      end if

                      call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                      if (enPacket%yP(1)<1) enPacket%yP(1)=1
                      if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                         if ( rvec%y >= ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                              & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                              & enPacket%yP(1) = enPacket%yP(1)+1
                      end if

                      call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                      if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                         if (enPacket%zP(1)<1) enPacket%zP(1)=1
                         if ( rvec%z >= ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                              & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                              & enPacket%zP(1) = enPacket%zP(1)+1
                      end if

                      xP = enPacket%xP(1)
                      yP = enPacket%yP(1)
                      zP = enPacket%zP(1)
                      gP = 1
                      igpp = 1
                   end if

                   if ( (rVec%z <= grid(1)%zAxis(1) .or.zP<1) ) then
                      zP=1
                      rVec%z = grid(gP)%yAxis(1)
                      vHat%z = -vHat%z
                   end if

                   if ( (rVec%z <= grid(gP)%zAxis(1) &
                        & .or.zP<1) .and. gP>1) then

                      ! locate where we are at on the mother grid
                      call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                      if (enPacket%xP(1)<1) enPacket%xP(1)=1
                      if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                         if ( rvec%x >= ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                              & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                              & enPacket%xP(1) = enPacket%xP(1)+1
                      end if

                      call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                      if (enPacket%yP(1)<1) enPacket%yP(1)=1
                      if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                         if ( rvec%y >= ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                              & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                              & enPacket%yP(1) = enPacket%yP(1)+1
                      end if

                      call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                      if (enPacket%zP(1)<1) enPacket%zP(1)=1
                      if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                         if ( rvec%z >= ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                              & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                              & enPacket%zP(1) = enPacket%zP(1)+1
                      end if

                      xP = enPacket%xP(1)
                      yP = enPacket%yP(1)
                      zP = enPacket%zP(1)
                      gP = 1
                      igpp = 1

                   end if

                   if ( (rVec%z >=  grid(1)%zAxis(grid(gP)%nz) .or. zP>grid(gP)%nz) &
                        & ) then

                      zP = grid(gP)%nz
                      rVec%z = grid(gP)%zAxis(grid(gP)%nz)
                      vHat%z = -vHat%z

                   end if

                   if ((rVec%z >=  grid(gP)%zAxis(grid(gP)%nz) &
                        & .or. zP>grid(gP)%nz) .and. gP>1) then

                      ! locate where we are at on the mother grid
                      call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                      if (enPacket%xP(1)<1) enPacket%xP(1)=1
                      if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                         if ( rvec%x >= ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                              & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                              & enPacket%xP(1) = enPacket%xP(1)+1
                      end if

                      call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                      if (enPacket%yP(1)<1) enPacket%yP(1)=1
                      if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                         if ( rvec%y >= ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                              & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                              & enPacket%yP(1) = enPacket%yP(1)+1
                      end if

                      call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                      if (enPacket%zP(1)<1) enPacket%zP(1)=1
                      if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                         if ( rvec%z >= ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                              & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                              & enPacket%zP(1) = enPacket%zP(1)+1
                      end if

                      xP = enPacket%xP(1)
                      yP = enPacket%yP(1)
                      zP = enPacket%zP(1)
                      gP = 1
                      igpp = 1

                   end if


                   if (lgReturn) then

                      ! the packet escapes without further interaction

                      idirT = int(acos(enPacket%direction%z)/dTheta)+1
                      if (idirT>totangleBinsTheta) then
                         idirT=totangleBinsTheta
                      end if
                      if (idirT<1 .or. idirT>totAngleBinsTheta) then
                         print*, '! pathSegment: error in theta direction cosine assignment',&
                              &  idirT, enPacket, dTheta, totAngleBinsTheta
                         stop
                      end if


                      if (enPacket%direction%x<1.e-35) then
                         idirP = 0
                      else
                         idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
                      end if
                      if (idirP<0) idirP=totAngleBinsPhi+idirP
                      idirP=idirP+1

                      if (idirP>totangleBinsPhi) then
                         idirP=totangleBinsPhi
                      end if
                      if (idirP<1 .or. idirP>totAngleBinsPhi) then
                         print*, '! pathSegment: error in phi direction cosine assignment',&
                              &  idirP, enPacket, dPhi, totAngleBinsPhi
                         stop
                      end if



                      if (nAngleBins>0) then
                         if ( (totangleBinsPhi>1 .and. viewPointPtheta(idirT) == viewPointPphi(idirP)) .or. &
                              & (totangleBinsPhi==1 .and. viewPointPtheta(idirT)>0) ) then
                            grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,&
                                 & viewPointPtheta(idirT)) = &
                                 &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                 & enPacket%nuP,viewPointPtheta(idirT)) +  deltaEUsed
                            grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,0) = &
                                 &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                 & enPacket%nuP,0) +  deltaEUsed

                            if (lgSeparateSED) then
                               grid(enPacket%origin(1))%escapedPacketsComponents&
                                    &(enPacket%origin(2), enPacket%nuP,&
                                    & viewPointPtheta(idirT),enPacket%SEDtype) = &
                                    & grid(enPacket%origin(1))%escapedPacketsComponents&
                                    &(enPacket%origin(2), &
                                    & enPacket%nuP,viewPointPtheta(idirT),enPacket%SEDtype)+deltaEUsed
                               grid(enPacket%origin(1))%escapedPacketsComponents&
                                    &(enPacket%origin(2), enPacket%nuP,0,enPacket%SEDtype) = &
                                    &grid(enPacket%origin(1))%escapedPacketsComponents&
                                    &(enPacket%origin(2), &
                                    & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                            endif
                         else
                            grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2),enPacket%nuP,0)=&
                                 & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                 & enPacket%nuP,0) +  deltaEUsed
                            if (lgSeparateSED) then
                               grid(enPacket%origin(1))%escapedPacketsComponents&
                                    &(enPacket%origin(2),enPacket%nuP,0,enPacket%SEDtype)=&
                                    & grid(enPacket%origin(1))%escapedPacketsComponents&
                                    &(enPacket%origin(2), &
                                    & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                            endif
                         end if
                      else
                         grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                              enPacket%nuP,0) = &
                              & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                              & enPacket%nuP,0) +  deltaEUsed
                         if (lgSeparateSED) then
                            grid(enPacket%origin(1))%escapedPacketsComponents&
                                 &(enPacket%origin(2), enPacket%nuP,0,enPacket%SEDtype) = &
                                 & grid(enPacket%origin(1))%escapedPacketsComponents&
                                 &(enPacket%origin(2), &
                                 & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                         endif
                      end if


                      return
                   end if

                end if

                ! check if the path is still within the simulation region
                radius = 1.e10*sqrt((rVec%x/1.e10)*(rVec%x/1.e10) + &
                     &                     (rVec%y/1.e10)*(rVec%y/1.e10) + &
                     &                     (rVec%z/1.e10)*(rVec%z/1.e10))

                if (.not.lgPlaneIonization) then

                   if ( (.not.lgSymmetricXYZ .and. (rVec%x<=grid(1)%xAxis(1) .or.&
                        & rVec%y<=grid(1)%yAxis(1) .or. rVec%z<=grid(1)%zAxis(1))) .or.&
                        & (rVec%x >= grid(gP)%xAxis(grid(gP)%nx)) .or.&
                        &(rVec%y >= grid(gP)%yAxis(grid(gP)%ny)) .or.&
                        &(rVec%z >= grid(gP)%zAxis(grid(gP)%nz)) .or. &
                        & xP>grid(gP)%nx .or. yP>grid(gP)%ny .or. zP>grid(gP)%nz ) then

                      if (gP==1) then

                         if (enPacket%xP(1) > grid(1)%nx) xP = grid(1)%nx
                         if (enPacket%yP(1) > grid(1)%ny) yP = grid(1)%ny
                         if (enPacket%zP(1) > grid(1)%nz) zP = grid(1)%nz
                         if (enPacket%xP(1) < 1) xP = 1
                         if (enPacket%yP(1) < 1) yP = 1
                         if (enPacket%zP(1) < 1) zP = 1


                         ! the energy packet escapes

                         ! the packet escapes without further interaction

                         idirT = int(acos(enPacket%direction%z)/dTheta)+1
                         if (idirT>totangleBinsTheta) then
                            idirT=totangleBinsTheta
                         end if
                         if (idirT<1 .or. idirT>totAngleBinsTheta) then
                            print*, '! pathSegment: error in theta direction cosine assignment',&
                                 &  idirT, enPacket, dTheta, totAngleBinsTheta
                            stop
                         end if

                         if (enPacket%direction%x<1.e-35) then
                            idirP = 0
                         else
                            idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
                         end if
                         if (idirP<0) idirP=totAngleBinsPhi+idirP
                         idirP=idirP+1

                         if (idirP>totangleBinsPhi) then
                            idirP=totangleBinsPhi
                         end if

                         if (idirP<1 .or. idirP>totAngleBinsPhi) then
                            print*, '! pathSegment: error in phi direction cosine assignment',&
                                 &  idirP, enPacket, dPhi, totAngleBinsPhi
                            stop
                         end if

                         if (nAngleBins>0) then
                            if ( (totangleBinsPhi>1 .and. viewPointPtheta(idirT) == viewPointPphi(idirP)) .or. &
                                 & (totangleBinsPhi==1 .and. viewPointPtheta(idirT)>0) ) then
                               grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,&
                                    & viewPointPtheta(idirT)) = &
                                    &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                    & enPacket%nuP,viewPointPtheta(idirT)) +  deltaEUsed
                               grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,0) = &
                                    &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                    & enPacket%nuP,0) +  deltaEUsed

                               if (lgSeparateSED) then
                                  grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                       &enPacket%nuP,&
                                   & viewPointPtheta(idirT),enPacket%SEDtype) = &
                                   &grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                   & enPacket%nuP,viewPointPtheta(idirT),enPacket%SEDtype) +  deltaEUsed
                                  grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                       &enPacket%nuP,0,enPacket%SEDtype)= &
                                       &grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                       & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                               endif

                            else
                               grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2),enPacket%nuP,0) = &
                                    & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                    & enPacket%nuP,0) +  deltaEUsed

                               if (lgSeparateSED) then
                                  grid(enPacket%origin(1))%escapedPacketsComponents&
                                       &(enPacket%origin(2),enPacket%nuP,0,enPacket%SEDtype) = &
                                       & grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                       & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                               endif

                            end if
                         else
                            grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                 enPacket%nuP,0) = &
                                 & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                 & enPacket%nuP,0) +  deltaEUsed
                            if (lgSeparateSED) then
                               grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                    enPacket%nuP,0,enPacket%SEDtype) = &
                                    & grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                    & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                            endif
                         end if

                         !b2.005
                         return

                      else if (gP>1) then


                         ! locate where we are at on the mother grid
                         call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                         if (enPacket%xP(1)<1) enPacket%xP(1)=1
                         if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                            if ( rvec%x >= ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                                 & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                                 & enPacket%xP(1) = enPacket%xP(1)+1
                         end if

                         call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                         if (enPacket%yP(1)<1) enPacket%yP(1)=1
                         if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                            if ( rvec%y >= ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                                 & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                                 & enPacket%yP(1) = enPacket%yP(1)+1
                         end if

                         call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                         if (enPacket%zP(1)<1) enPacket%zP(1)=1
                         if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                            if ( rvec%z >= ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                                 & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                                 & enPacket%zP(1) = enPacket%zP(1)+1
                         end if

                         xP = enPacket%xP(1)
                         yP = enPacket%yP(1)
                         zP = enPacket%zP(1)
                         gP = 1
                         igpp = 1


                         if (gP/=1) then
                            print*, '! pathSegment: nested multigrids still not implemented'
                            stop
                         end if

                         if ( (radius >= R_out .and. R_out >= 0.) .or.&
                              & (rVec%x >= grid(1)%xAxis(grid(1)%nx)) .or.&
                              &(rVec%y  >= grid(1)%yAxis(grid(1)%ny)) .or.&
                              &(rVec%z  >= grid(1)%zAxis(grid(1)%nz)) .or. &
                              & (.not.lgSymmetricXYZ .and.  &
                              & (rVec%x <= grid(1)%xAxis(1) .or.&
                              &  rVec%y <= grid(1)%yAxis(1) .or.&
                              &  rVec%z <= grid(1)%zAxis(1)))) then


                            if (xP > grid(gP)%nx) xP = grid(gP)%nx
                            if (yP > grid(gP)%ny) yP = grid(gP)%ny
                            if (zP > grid(gP)%nz) zP = grid(gP)%nz

                            if (xP < 1) xP = 1
                            if (yP < 1) yP = 1
                            if (zP < 1) zP = 1

                            ! the energy packet escapes

                            ! the packet escapes without further interaction

                            idirT = int(acos(enPacket%direction%z)/dTheta)+1
                            if (idirT>totangleBinsTheta) then
                               idirT=totangleBinsTheta
                            end if
                            if (idirT<1 .or. idirT>totAngleBinsTheta) then
                               print*, '! pathSegment: error in theta direction cosine assignment',&
                                    &  idirT, enPacket, dTheta, totAngleBinsTheta
                               stop
                            end if

                            if (enPacket%direction%x<1.e-35) then
                               idirP = 0
                            else
                               idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
                            end if
                            if (idirP<0) idirP=totAngleBinsPhi+idirP
                            idirP=idirP+1

                            if (idirP>totangleBinsPhi) then
                               idirP=totangleBinsPhi
                            end if

                            if (idirP<1 .or. idirP>totAngleBinsPhi) then
                               print*, '! pathSegment: error in phi direction cosine assignment',&
                                    &  idirP, enPacket, dPhi, totAngleBinsPhi
                               stop
                            end if


                            if (nAngleBins>0) then
                               if ( (totangleBinsPhi>1 .and. viewPointPtheta(idirT) == viewPointPphi(idirP)) .or. &
                                    & (totangleBinsPhi==1 .and. viewPointPtheta(idirT)>0) ) then
                              grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,&
                                   & viewPointPtheta(idirT)) = &
                                   &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                   & enPacket%nuP,viewPointPtheta(idirT)) +  deltaEUsed
                              grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,0) = &
                                   &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                   & enPacket%nuP,0) +  deltaEUsed
                              if (lgSeparateSED) then
                                 grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), enPacket%nuP,&
                                      & viewPointPtheta(idirT),enPacket%SEDtype) = &
                                      &grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                      & enPacket%nuP,viewPointPtheta(idirT),enPacket%SEDtype) +  deltaEUsed
                                 grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                      &enPacket%nuP,0,enPacket%SEDtype) = &
                                      &grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                      & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                              endif

                           else
                              grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2),enPacket%nuP,0) = &
                                   & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                   & enPacket%nuP,0) +  deltaEUsed
                              if (lgSeparateSED) then
                                 grid(enPacket%origin(1))%escapedPacketsComponents&
                                      &(enPacket%origin(2),enPacket%nuP,0,enPacket%SEDtype) = &
                                      & grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                      & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                              endif
                           end if
                        else
                           grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                & enPacket%nuP,0) = &
                                & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                & enPacket%nuP,0) +  deltaEUsed
                           if (lgSeparateSED) then
                              grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                   & enPacket%nuP,0,enPacket%SEDtype) = &
                                   & grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                   & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                           endif
                        end if



                            !b2.005
                            return

                         end if
                      else
                         print*, '! pathSegment: insanity occured - invalid gP - ', gP
                         stop
                      end if

                   end if

                   if (lgSymmetricXYZ ) then
                      if (lgPlaneIonization) then
                         print*, '! pathSegment: lgSymmetric and lgPlaneionization flags both raised'
                         stop
                      end if

                      if ( rVec%x <= grid(1)%xAxis(1) .or. (gP==1 .and. xP<1)) then
                         if (vHat%x<0.) vHat%x = -vHat%x
                         enPacket%xP(1) = 1
                         xP = 1
                         rVec%x = grid(gP)%xAxis(1)
                      end if
                      if ( rVec%y <= grid(1)%yAxis(1) .or. (gP==1 .and. yP<1)) then
                         if (vHat%y<0.) vHat%y = -vHat%y
                         enPacket%yP(1)=1
                         yP = 1
                         rVec%y = grid(gP)%yAxis(1)
                      end if
                      if ( rVec%z <= grid(1)%zAxis(1) .or. (gP==1 .and. zP<1)) then
                         if (vHat%z<0.) vHat%z = -vHat%z
                         enPacket%zP(1) = 1
                         zP=1
                         rVec%z = grid(1)%zAxis(1)
                      end if

                   end if

                end if

                if (gP>1) then
!print*, iPhot, rvec, vHat, xp, yp, zp, gp
                   if ( ( (rVec%x <= grid(gP)%xAxis(1) &
                        &.or. xP<1) .and. vHat%x <=0.) .or. &
                        & ( (rVec%y <= grid(gP)%yAxis(1) &
                        & .or. yP<1) .and. vHat%y <=0.) .or. &
                        & ( (rVec%z <= grid(gP)%zAxis(1) &
                        &  .or. zP<1) .and. vHat%z <=0.) .or. &
                        & ( (rVec%x >= grid(gP)%xAxis(grid(gP)%nx) &
                        &.or. xP>grid(gP)%nx) .and. vHat%x >=0.) .or. &
                        & ( (rVec%y >= grid(gP)%yAxis(grid(gP)%ny) &
                        & .or. yP>grid(gP)%ny) .and. vHat%y >=0.) .or. &
                        & ( (rVec%z >= grid(gP)%zAxis(grid(gP)%nz) &
                        &  .or. zP>grid(gP)%nz) .and. vHat%z >=0.) ) then


                      ! locate where we are at on the mother grid
                      call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                      if (enPacket%xP(1)<1) enPacket%xP(1)=1
                      if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                         if ( rvec%x >= ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                              & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                              & enPacket%xP(1) = enPacket%xP(1)+1
                      end if

                      call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                      if (enPacket%yP(1)<1) enPacket%yP(1)=1
                      if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                         if ( rvec%y >= ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                              & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                              & enPacket%yP(1) = enPacket%yP(1)+1
                      end if

                      call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                      if (enPacket%zP(1)<1) enPacket%zP(1)=1
                      if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                         if ( rvec%z >= ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                              & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                              & enPacket%zP(1) = enPacket%zP(1)+1
                      end if

                      gP = 1
                      igpp = 1


                      xP = enPacket%xP(1)
                      yP = enPacket%yP(1)
                      zP = enPacket%zP(1)


                   end if

                end if


             end if

             if (.not. lgPlaneIonization .and. gP==1 .and. (xP > grid(gP)%nx  &
                  & .or. yP > grid(gP)%ny .or. zP > grid(gP)%nz) ) then

                ! the energy packet escapes

                ! the packet escapes without further interaction

                idirT = int(acos(enPacket%direction%z)/dTheta)+1
                if (idirT>totangleBinsTheta) then
                   idirT=totangleBinsTheta
                end if
                if (idirT<1 .or. idirT>totAngleBinsTheta) then
                   print*, '! pathSegment: error in theta direction cosine assignment',&
                        &  idirT, enPacket, dTheta, totAngleBinsTheta
                   stop
                end if

                if (enPacket%direction%x<1.e-35) then
                   idirP = 0
                else
                   idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
                end if
                if (idirP<0) idirP=totAngleBinsPhi+idirP
                idirP=idirP+1

                if (idirP>totangleBinsPhi) then
                   idirP=totangleBinsPhi
                end if

                if (idirP<1 .or. idirP>totAngleBinsPhi) then
                   print*, '! pathSegment: error in phi direction cosine assignment',&
                        &  idirP, enPacket, dPhi, totAngleBinsPhi
                   stop
                end if

            if (nAngleBins>0) then
               if ( (totangleBinsPhi>1 .and. viewPointPtheta(idirT) == viewPointPphi(idirP)) .or. &
                    & (totangleBinsPhi==1 .and. viewPointPtheta(idirT)>0) ) then
                  grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,&
                       & viewPointPtheta(idirT)) = &
                       &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                       & enPacket%nuP,viewPointPtheta(idirT)) +  deltaEUsed
                  grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,0) = &
                       &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                       & enPacket%nuP,0) +  deltaEUsed
                  if (lgSeparateSED) then
                     grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), enPacket%nuP,&
                          & viewPointPtheta(idirT), enPacket%SEDtype) = &
                          &grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                          & enPacket%nuP,viewPointPtheta(idirT), enPacket%SEDtype) +  deltaEUsed
                     grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                          & enPacket%nuP,0, enPacket%SEDtype) = &
                          &grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                          & enPacket%nuP,0, enPacket%SEDtype) +  deltaEUsed
                  endif
               else
                  grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2),enPacket%nuP,0) = &
                       & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                       & enPacket%nuP,0) +  deltaEUsed
                  if (lgSeparateSED) then
                     grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2),&
                          &enPacket%nuP,0, enPacket%SEDtype) = &
                          & grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                          & enPacket%nuP,0, enPacket%SEDtype) +  deltaEUsed
                  endif
               end if
            else
               grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                    enPacket%nuP,0) = &
                    & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                    & enPacket%nuP,0) +  deltaEUsed
               if (lgSeparateSED) then
                  grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                       & enPacket%nuP,0, enPacket%SEDtype) = &
                       & grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                       & enPacket%nuP,0, enPacket%SEDtype) +  deltaEUsed
               endif

            end if


                !b2.005
                return


             end if

          end do ! safelimit loop

          if (i>= safeLimit) then
             if (.not.lgPlaneIonization) then
                print*, '! pathSegment: [warning] packet trajectory has exceeded&
                     &  maximum number of events', safeLimit, gP, xP,yP,zP, grid(gP)%active(xP,yP,zP), &
                     & rvec, vhat, enPacket, iphot
             end if
             return

          end if

          if (gP==1) then
             igpp = 1
          else if (gP>1) then
             igpp = 2
          else
             print*, "! pathSegment: insane grid index "
             stop
          end if

          enPacket%xP(igpp) = xP
          enPacket%yP(igpp) = yP
          enPacket%zP(igpp) = zP

          if (deltaEUsed/(Lstar(istar)/nPhotons(istar)) < 1.e-10) return

!          if (countRecursive > recursionLimit) then
!             trapped = trapped+1
!             return
!          end if



          ! the energy packet has beenid absorbed - reemit a new packet from this position
!          call energyPacketRun(packetType, rVec, enPacket%xP(1:2), enPacket%yP(1:2), &
!               & enPacket%zP(1:2), gP)

          chTypeIn = packetType
          positionIn = rVec
          inX =  enPacket%xP(1:2)
          inY =  enPacket%yP(1:2)
          inZ =  enPacket%zP(1:2)
          gPIn = gP
          reRun = 1
          return

        end subroutine pathSegment

    subroutine getFluorescenceL(deltaE,nuPin,ifluoin,igin,xpin,ypin,zpin)
      implicit none


      real, intent(inout) :: deltaE

      real :: fYieldFeKa(26)                  ! fluorescence yields for FeKa
      real :: fYieldFeL1                      ! fluorescence yields for FeL1
      real :: fYieldFeL2                      ! fluorescence yields for FeL2
      real :: fYieldCKa                       ! fluorescence yields for CKa
      real :: fYieldNKa                       ! fluorescence yields for NKa
      real :: fYieldOKa                       ! fluorescence yields for OKa
      real :: fYieldNeKa                      ! fluorescence yields for NeKa
      real :: fYieldMgKa                      ! fluorescence yields for MgKa
      real :: fYieldAlKa                      ! fluorescence yields for AlKa
      real :: fYieldSiKa                      ! fluorescence yields for SiKa
      real :: fYieldSKa                       ! fluorescence yields for SKa
      real :: fYieldArKa                      ! fluorescence yields for ArK
      real :: fYieldCaKa                      ! fluorescence yields for CaKa
      real :: fYieldTiKa                      ! fluorescence yields for TiKa

      real :: radField
      real :: branch                          ! branching ratio
      real :: phXsec, photoRateFluoFrac

      integer, intent(in) :: ifluoin, nuPin
      integer, intent(in) :: xpin,ypin,zpin, igin
      integer :: ion, icell, elP
      integer :: nu1P, highNuP,xSecP

!      print*, 'in getFluorescenceL ', xpin,ypin,zpin

      ! Krolik&Kallman ApJL 320, 5
      fYieldFeKa = (/0.34, 0.34, 0.35, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, &
           & 0.44, 0.45, 0.46, 0.47, 0.47, 0.48, 0.48, 0.49, 0.49, 0.11, 0.75, 0., 0., 0. /)
      fYieldFeL1 = 1.0e-3
      fYieldFeL2 = 1.26e-2
      fYieldCKa = 2.80e-3
      fYieldNKa = 5.20e-3
      fYieldOKa = 8.30e-3
      fYieldNeKa = 1.39e-2
      fYieldMgKa = 2.77e-2
      fYieldAlKa = 0.039
      fYieldSiKa = 5.0e-2
      fYieldSKa = 7.80e-2
      fYieldArKa = 1.15e-1
      fYieldCaKa = 1.60e-1
      fYieldTiKa = 4.01e-1

      radField = deltaE

      photoRateFluoFrac = 0.
      deltaE = 0.

      icell = grid(igin)%active(xpin,ypin,zpin)
      if (icell <= 0) then
         print*, '! getFluorescenceL: fluorescence lines cannot be produced by an inactive cell.', &
              & igin,xpin,ypin,zpin
         stop
      end if

      select case (fluorescenceLabel(ifluoin))
      case ("FeKa")

         branch = 0.882 ! Bambynek et al., 1972
         do ion = 1, min(nstages, 18)

            nu1P = elementP(26, ion, 1, 1)
            highNuP = elementP(26, ion, 1, 2)
            xSecP = elementP(26, ion, 1, 3)

            phXsec = xSecArray(nuPin+xSecP-nu1P)
            if (phXsec < 1.e-35) then
               phXsec = 0.
            end if

            if (nuPin<highnuP .and. nuPin>nu1P ) then

               photoRateFluoFrac = phXsec*grid(igin)%elemAbun(grid(igin)%&
                    & abFileIndex(xPin,yPin,zPin),26)*grid(igin)%Hden(icell)*&
                    & grid(igin)%ionDen(icell,elementXref(26),ion)/&
                    & ((grid(igin)%opacity(grid(igin)%active(xPin,yPin,zPin),nuPin)-&
                    &KNsigmaT(nuPin)*grid(igin)%Ne(grid(igin)%active(xpin,ypin,zpin)))*&
                    & Ryd2erg*nuArray(nuPin))

!print*, ion, photoRateFluoFrac, grid(igin)%active(xpin,ypin,zpin), xpin,ypin,zpin, nuPin, highNuP, nu1P
! 977           7           1         200
!        1373        1409

!print*, phXsec*grid(igin)%elemAbun(grid(igin)%&
!                 & abFileIndex(xPin,yPin,zPin),26)*grid(igin)%Hden(icell)*&
!                 & grid(igin)%ionDen(icell,elementXref(26),ion), &
!                 & grid(igin)%opacity(grid(igin)%active(xPin,yPin,zPin),nuPin),&
!                 & KNsigmaT(nuPin)*grid(igin)%Ne(grid(igin)%active(xpin,ypin,zpin)), &
!                 & Ryd2erg*nuArray(nuPin)


! phXsec*grid(igin)%elemAbun(grid(igin)%&
!                 & abFileIndex(xPin,yPin,zPin),26)*grid(igin)%Hden(icell)*&
!                 & grid(igin)%ionDen(icell,elementXref(26),ion), &
!                 & grid(igin)%opacity(grid(igin)%active(xPin,yPin,zPin),nuPin)*&
!                 & Ryd2erg*nuArray(nuPin), &
!                 & KNsigmaT(nuPin)*grid(igin)%Ne(grid(igin)%active(xpin,ypin,zpin))*&
!                 & Ryd2erg*nuArray(nuPin)

               photoRateFluoFrac = radField*photoRateFluoFrac
            else
               photoRateFluoFrac = 0.
            end if

            deltaE = deltaE+photoRateFluoFrac*&
                 & fYieldFeKa(ion)*branch

         end do

         deltaE = deltaE*1.03e-8

         ! "hot" iron contribution at 6.9 KeV could be implemented &
         ! here similar fashion as above, just multiply
         ! by correct energy (i.e. 1.11e-8 rather than 1.03e-8 erg)

      case ("FeL1")

         elP = 26

         do ion = 1, min(nstages, elP+1)

            nu1P = elementP(elP, ion, 2, 1)
            highNuP = elementP(elP, ion, 2, 2)
            xSecP = elementP(elP, ion, 2, 3)

            phXsec = xSecArray(nuPin+xSecP-nu1P)
            if (phXsec < 1.e-35) then
               phXsec = 0.
            end if
            photoRateFluoFrac = phXsec*grid(igin)%elemAbun(grid(igin)%&
                 & abFileIndex(xPin,yPin,zPin),elP)*grid(igin)%Hden(icell)*&
                 & grid(igin)%ionDen(icell,elementXref(elP),ion)/&
                 & ((grid(igin)%opacity(grid(igin)%active(xPin,yPin,zPin),nuPin)-&
                 & KNsigmaT(nuPin)*grid(igin)%Ne(grid(igin)%active(xpin,ypin,zpin)))*&
                 & Ryd2erg*nuArray(nuPin))

            photoRateFluoFrac = radField*photoRateFluoFrac

            deltaE = deltaE+photoRateFluoFrac*&
                 & fYieldFeL1
         end do

         deltaE = deltaE*FeL1*Ryd2erg

      case ("FeL2")

         elP = 26

         do ion = 1, min(nstages, elP+1)

            nu1P = elementP(elP, ion, 2, 1)
            highNuP = elementP(elP, ion, 2, 2)
            xSecP = elementP(elP, ion, 2, 3)

            phXsec = xSecArray(nuPin+xSecP-nu1P)
            if (phXsec < 1.e-35) then
               phXsec = 0.
            end if
            photoRateFluoFrac = phXsec*grid(igin)%elemAbun(grid(igin)%&
                 & abFileIndex(xPin,yPin,zPin),elP)*grid(igin)%Hden(icell)*&
                 & grid(igin)%ionDen(icell,elementXref(elP),ion)/&
                 & ((grid(igin)%opacity(grid(igin)%active(xPin,yPin,zPin),nuPin)-&
                 & KNsigmaT(nuPin)*grid(igin)%Ne(grid(igin)%active(xpin,ypin,zpin)))*&
                 & Ryd2erg*nuArray(nuPin))

            photoRateFluoFrac = radField*photoRateFluoFrac

            deltaE = deltaE+photoRateFluoFrac*&
                 & fYieldFeL2
         end do

         deltaE = deltaE*FeL2*Ryd2erg

      case ("CKa")

         elP = 6

         do ion = 1, min(nstages, elP+1)

            nu1P    = elementP(elP, ion, 1, 1)
            highNuP = elementP(elP, ion, 1, 2)
            xSecP   = elementP(elP, ion, 1, 3)

            phXsec = xSecArray(nuPin+xSecP-nu1P)
            if (phXsec < 1.e-35) then
               phXsec = 0.
            end if
            photoRateFluoFrac = phXsec*grid(igin)%elemAbun(grid(igin)%&
                 & abFileIndex(xPin,yPin,zPin),elP)*grid(igin)%Hden(icell)*&
                 & grid(igin)%ionDen(icell,elementXref(elP),ion)/&
                 & ((grid(igin)%opacity(grid(igin)%active(xPin,yPin,zPin),nuPin)-&
                 & KNsigmaT(nuPin)*grid(igin)%Ne(grid(igin)%active(xpin,ypin,zpin)))*&
                 & Ryd2erg*nuArray(nuPin))

            photoRateFluoFrac = radField*photoRateFluoFrac

            deltaE = deltaE+photoRateFluoFrac*&
                 & fYieldCKa
         end do

         deltaE = deltaE*CKa*Ryd2erg

      case ("NKa")

         elP = 7

         do ion = 1, min(nstages, elP+1)

            nu1P    = elementP(elP, ion, 1, 1)
            highNuP = elementP(elP, ion, 1, 2)
            xSecP   = elementP(elP, ion, 1, 3)

            phXsec = xSecArray(nuPin+xSecP-nu1P)
            if (phXsec < 1.e-35) then
               phXsec = 0.
            end if
            photoRateFluoFrac = phXsec*grid(igin)%elemAbun(grid(igin)%&
                 & abFileIndex(xPin,yPin,zPin),elP)*grid(igin)%Hden(icell)*&
                 & grid(igin)%ionDen(icell,elementXref(elP),ion)/&
                 & ((grid(igin)%opacity(grid(igin)%active(xPin,yPin,zPin),nuPin)-&
                 & KNsigmaT(nuPin)*grid(igin)%Ne(grid(igin)%active(xpin,ypin,zpin)))*&
                 & Ryd2erg*nuArray(nuPin))

            photoRateFluoFrac = radField*photoRateFluoFrac

            deltaE = deltaE+photoRateFluoFrac*&
                 & fYieldNKa
         end do

         deltaE = deltaE*NKa*Ryd2erg

      case ("OKa")

         elP = 8

         do ion = 1, min(nstages, elP+1)

            nu1P    = elementP(elP, ion, 1, 1)
            highNuP = elementP(elP, ion, 1, 2)
            xSecP   = elementP(elP, ion, 1, 3)

            phXsec = xSecArray(nuPin+xSecP-nu1P)
            if (phXsec < 1.e-35) then
               phXsec = 0.
            end if
            photoRateFluoFrac = phXsec*grid(igin)%elemAbun(grid(igin)%&
                 & abFileIndex(xPin,yPin,zPin),elP)*grid(igin)%Hden(icell)*&
                 & grid(igin)%ionDen(icell,elementXref(elP),ion)/&
                 & ((grid(igin)%opacity(grid(igin)%active(xPin,yPin,zPin),nuPin)-&
                 & KNsigmaT(nuPin)*grid(igin)%Ne(grid(igin)%active(xpin,ypin,zpin)))*&
                 & Ryd2erg*nuArray(nuPin))

            photoRateFluoFrac = radField*photoRateFluoFrac

            deltaE = deltaE+photoRateFluoFrac*&
                 & fYieldOKa
         end do

         deltaE = deltaE*OKa*Ryd2erg

      case ("NeKa")

         elP = 10

         do ion = 1, min(nstages, elP+1)

            nu1P    = elementP(elP, ion, 1, 1)
            highNuP = elementP(elP, ion, 1, 2)
            xSecP   = elementP(elP, ion, 1, 3)

            phXsec = xSecArray(nuPin+xSecP-nu1P)
            if (phXsec < 1.e-35) then
               phXsec = 0.
            end if
            photoRateFluoFrac = phXsec*grid(igin)%elemAbun(grid(igin)%&
                 & abFileIndex(xPin,yPin,zPin),elP)*grid(igin)%Hden(icell)*&
                 & grid(igin)%ionDen(icell,elementXref(elP),ion)/&
                 & ((grid(igin)%opacity(grid(igin)%active(xPin,yPin,zPin),nuPin)-&
                 & KNsigmaT(nuPin)*grid(igin)%Ne(grid(igin)%active(xpin,ypin,zpin)))*&
                 & Ryd2erg*nuArray(nuPin))

            photoRateFluoFrac = radField*photoRateFluoFrac

            deltaE = deltaE+photoRateFluoFrac*&
                 & fYieldNeKa
         end do

         deltaE = deltaE*NeKa*Ryd2erg

      case ("MgKa")

         elP = 12

         do ion = 1, min(nstages, elP+1)

            nu1P    = elementP(elP, ion, 1, 1)
            highNuP = elementP(elP, ion, 1, 2)
            xSecP   = elementP(elP, ion, 1, 3)

            phXsec = xSecArray(nuPin+xSecP-nu1P)
            if (phXsec < 1.e-35) then
               phXsec = 0.
            end if
            photoRateFluoFrac = phXsec*grid(igin)%elemAbun(grid(igin)%&
                 & abFileIndex(xPin,yPin,zPin),elP)*grid(igin)%Hden(icell)*&
                 & grid(igin)%ionDen(icell,elementXref(elP),ion)/&
                 & ((grid(igin)%opacity(grid(igin)%active(xPin,yPin,zPin),nuPin)-&
                 & KNsigmaT(nuPin)*grid(igin)%Ne(grid(igin)%active(xpin,ypin,zpin)))*&
                 & Ryd2erg*nuArray(nuPin))

            photoRateFluoFrac = radField*photoRateFluoFrac

            deltaE = deltaE+photoRateFluoFrac*&
                 & fYieldMgKa
         end do

         deltaE = deltaE*MgKa*Ryd2erg

      case ("AlKa")

         elP = 13

         do ion = 1, min(nstages, elP+1)

            nu1P    = elementP(elP, ion, 1, 1)
            highNuP = elementP(elP, ion, 1, 2)
            xSecP   = elementP(elP, ion, 1, 3)

            phXsec = xSecArray(nuPin+xSecP-nu1P)
            if (phXsec < 1.e-35) then
               phXsec = 0.
            end if
            photoRateFluoFrac = phXsec*grid(igin)%elemAbun(grid(igin)%&
                 & abFileIndex(xPin,yPin,zPin),elP)*grid(igin)%Hden(icell)*&
                 & grid(igin)%ionDen(icell,elementXref(elP),ion)/&
                 & ((grid(igin)%opacity(grid(igin)%active(xPin,yPin,zPin),nuPin)-&
                 & KNsigmaT(nuPin)*grid(igin)%Ne(grid(igin)%active(xpin,ypin,zpin)))*&
                 & Ryd2erg*nuArray(nuPin))

            photoRateFluoFrac = radField*photoRateFluoFrac

            deltaE = deltaE+photoRateFluoFrac*&
                 & fYieldAlKa
         end do

         deltaE = deltaE*AlKa*Ryd2erg

      case ("SiKa")

         elP = 14

         do ion = 1, min(nstages, elP+1)

            nu1P    = elementP(elP, ion, 1, 1)
            highNuP = elementP(elP, ion, 1, 2)
            xSecP   = elementP(elP, ion, 1, 3)

            phXsec = xSecArray(nuPin+xSecP-nu1P)
            if (phXsec < 1.e-35) then
               phXsec = 0.
            end if
            photoRateFluoFrac = phXsec*grid(igin)%elemAbun(grid(igin)%&
                 & abFileIndex(xPin,yPin,zPin),elP)*grid(igin)%Hden(icell)*&
                 & grid(igin)%ionDen(icell,elementXref(elP),ion)/&
                 & ((grid(igin)%opacity(grid(igin)%active(xPin,yPin,zPin),nuPin)-&
                 & KNsigmaT(nuPin)*grid(igin)%Ne(grid(igin)%active(xpin,ypin,zpin)))*&
                 & Ryd2erg*nuArray(nuPin))

            photoRateFluoFrac = radField*photoRateFluoFrac

            deltaE = deltaE+photoRateFluoFrac*&
                 & fYieldSiKa
         end do

         deltaE = deltaE*SiKa*Ryd2erg

      case ("SKa")

         elP = 16

         do ion = 1, min(nstages, elP+1)

            nu1P    = elementP(elP, ion, 1, 1)
            highNuP = elementP(elP, ion, 1, 2)
            xSecP   = elementP(elP, ion, 1, 3)

            phXsec = xSecArray(nuPin+xSecP-nu1P)
            if (phXsec < 1.e-35) then
               phXsec = 0.
            end if
            photoRateFluoFrac = phXsec*grid(igin)%elemAbun(grid(igin)%&
                 & abFileIndex(xPin,yPin,zPin),elP)*grid(igin)%Hden(icell)*&
                 & grid(igin)%ionDen(icell,elementXref(elP),ion)/&
                 & ((grid(igin)%opacity(grid(igin)%active(xPin,yPin,zPin),nuPin)-&
                 & KNsigmaT(nuPin)*grid(igin)%Ne(grid(igin)%active(xpin,ypin,zpin)))*&
                 & Ryd2erg*nuArray(nuPin))

            photoRateFluoFrac = radField*photoRateFluoFrac

            deltaE = deltaE+photoRateFluoFrac*&
                 & fYieldSKa
         end do

         deltaE = deltaE*SKa*Ryd2erg

      case ("ArKa")

         elP = 18

         do ion = 1, min(nstages, elP+1)

            nu1P    = elementP(elP, ion, 1, 1)
            highNuP = elementP(elP, ion, 1, 2)
            xSecP   = elementP(elP, ion, 1, 3)

            phXsec = xSecArray(nuPin+xSecP-nu1P)
            if (phXsec < 1.e-35) then
               phXsec = 0.
            end if
            photoRateFluoFrac = phXsec*grid(igin)%elemAbun(grid(igin)%&
                 & abFileIndex(xPin,yPin,zPin),elP)*grid(igin)%Hden(icell)*&
                 & grid(igin)%ionDen(icell,elementXref(elP),ion)/&
                 & ((grid(igin)%opacity(grid(igin)%active(xPin,yPin,zPin),nuPin)-&
                 & KNsigmaT(nuPin)*grid(igin)%Ne(grid(igin)%active(xpin,ypin,zpin)))*&
                 & Ryd2erg*nuArray(nuPin))

            photoRateFluoFrac = radField*photoRateFluoFrac

            deltaE = deltaE+photoRateFluoFrac*&
                 & fYieldArKa
         end do

         deltaE = deltaE*ArKa*Ryd2erg

      case ("CaKa")

         elP = 20

         do ion = 1, min(nstages, elP+1)

            nu1P    = elementP(elP, ion, 1, 1)
            highNuP = elementP(elP, ion, 1, 2)
            xSecP   = elementP(elP, ion, 1, 3)

            phXsec = xSecArray(nuPin+xSecP-nu1P)
            if (phXsec < 1.e-35) then
               phXsec = 0.
            end if
            photoRateFluoFrac = phXsec*grid(igin)%elemAbun(grid(igin)%&
                 & abFileIndex(xPin,yPin,zPin),elP)*grid(igin)%Hden(icell)*&
                 & grid(igin)%ionDen(icell,elementXref(elP),ion)/&
                 & ((grid(igin)%opacity(grid(igin)%active(xPin,yPin,zPin),nuPin)-&
                 & KNsigmaT(nuPin)*grid(igin)%Ne(grid(igin)%active(xpin,ypin,zpin)))*&
                 & Ryd2erg*nuArray(nuPin))

            photoRateFluoFrac = radField*photoRateFluoFrac

            deltaE = deltaE+photoRateFluoFrac*&
                 & fYieldCaKa
         end do

         deltaE = deltaE*CaKa*Ryd2erg

      case ("TiKa")

         elP = 22

         do ion = 1, min(nstages, elP+1)

            nu1P    = elementP(elP, ion, 1, 1)
            highNuP = elementP(elP, ion, 1, 2)
            xSecP   = elementP(elP, ion, 1, 3)

            phXsec = xSecArray(nuPin+xSecP-nu1P)
            if (phXsec < 1.e-35) then
               phXsec = 0.
            end if
            photoRateFluoFrac = phXsec*grid(igin)%elemAbun(grid(igin)%&
                 & abFileIndex(xPin,yPin,zPin),elP)*grid(igin)%Hden(icell)*&
                 & grid(igin)%ionDen(icell,elementXref(elP),ion)/&
                 & ((grid(igin)%opacity(grid(igin)%active(xPin,yPin,zPin),nuPin)-&
                 & KNsigmaT(nuPin)*grid(igin)%Ne(grid(igin)%active(xpin,ypin,zpin)))*&
                 & Ryd2erg*nuArray(nuPin))

            photoRateFluoFrac = radField*photoRateFluoFrac

            deltaE = deltaE+photoRateFluoFrac*&
                 & fYieldTiKa
         end do

         deltaE = deltaE*TiKa*Ryd2erg

      case default

         print*, "! setDiffuseFluorescenceL: fluorescence label not recognised", &
              & fluorescenceLabel(ifluoin), ifluoin
         stop

      end select



!print*, icell,  deltaE, ' getFluorescenceL'

!      print*, 'out getFluorescenceL'

    end subroutine getFluorescenceL



      recursive subroutine fluorescencePacketRun(fType, rvecin, xPin, yPin, zPin, gPin, ifl)
        implicit none

        type(vector), intent(in)                  :: rVecin        !
        type(vector)                              :: rVec        !

        character(len=4), intent(inout)                   :: fType ! what resonance line?

        integer, intent(in)                       :: ifl
        integer, dimension(2), intent(in)         :: xPin, yPin, &
             & zPin                                            ! cartesian axes indeces
        integer, dimension(2)                     :: xP, yP, &
             & zP                                            ! cartesian axes indeces
        ! 1= mother; 2=sub

        integer, intent(inout)           :: gPin             ! grid index
        integer                          :: gP               ! grid index
        integer                          :: igpr             ! grid pointer 1= mother 2=sub

        ! local variables

        type(photon_packet)              :: fluoPacket         ! the energu packet

        rvec = rvecin
        xp = xpin
        yp = ypin
        zp = zpin
        gP = gPin

        if (gP==1) then
           igpr = 1
        else if (gP>1) then
           igpr = 2
        else
           print*,  "! fluorescencePacketRun: insane grid index"
           stop
        end if

        fluoPacket = newFluorescencePacket(chType=fType, gP=gP, &
             &orXf=xP,orYf=yP,orZf=zP, rvec=rvec)

        fluoPacket%SEDtype = ifl

        ! compute the next segment of trajectory
        call fluoPathSegment(fluoPacket)

        return

      end subroutine fluorescencePacketRun


        ! this function initializes a photon packet
        function initFluorescencePacket(nuP,  position, lgLine, lgStellar, xP, yP, zP, gP)
          implicit none

          type(photon_packet)      :: initFluorescencePacket  ! the photon packet

          real                     :: random            ! random number

          integer, intent(in)      :: nuP               ! the frequency of the photon packet
          integer, intent(in),dimension(2) :: xP, yP, &
               & zP                                     ! indeces of position on the x, y and z axes
          integer, intent(in)      :: gP                ! grid index
          integer                  :: igpi              ! grid pointer 1=mother, 2=sub


          logical, intent(in)      :: lgLine, lgStellar ! line, stellar packet?

          type(vector), intent(in) :: position          ! the position at which the photon
                                                          ! packet is created
          ! local variables

          integer                  :: irepeat           ! counter
            !print*, 'inin'

!print*, 'init ', gP, xp,yp,zp
            initFluorescencePacket%position = position

            initFluorescencePacket%iG  = gP

            if (gP==1) then
               igpi=1
            else if (gp>1) then
               igpi=2
            else
               print*, "! initFluorescencePacket: insane gridp pointer"
               stop
            end if

            initFluorescencePacket%nuP      = nuP

            initFluorescencePacket%lgStellar = lgStellar

            ! check if photon packen is line or continuum photon
            if ( lgLine ) then
                ! line photon
                initFluorescencePacket%nu       = 0.
                initFluorescencePacket%lgLine   = .true.
            else
                ! continuum photon
                initFluorescencePacket%nu       = nuArray(nuP)
                initFluorescencePacket%lgLine   = .false.
            end if

            initFluorescencePacket%xP  = xP
            initFluorescencePacket%yP  = yP
            initFluorescencePacket%zP  = zP


            ! cater for plane parallel ionization case
            if (initFluorescencePacket%lgStellar .and. lgPlaneIonization) then

               ! get position

               ! x-direction
               call random_number(random)
               random = 1. - random
               initFluorescencePacket%position%x = &
                    & -(grid(gP)%xAxis(2)-grid(gP)%xAxis(1))/2. + random*( &
                    & (grid(gP)%xAxis(2)-grid(gP)%xAxis(1))/2.+&
                    & (grid(gP)%xAxis(grid(gP)%nx)-grid(gP)%xAxis(grid(gP)%nx-1))/2.+&
                    & grid(gP)%xAxis(grid(gP)%nx))
               if (initFluorescencePacket%position%x<grid(gP)%xAxis(1)) &
                    & initFluorescencePacket%position%x=grid(gP)%xAxis(1)
               if (initFluorescencePacket%position%x>grid(gP)%xAxis(grid(gP)%nx)) &
                    initFluorescencePacket%position%x=grid(gP)%xAxis(grid(gP)%nx)

               call locate(grid(gP)%xAxis, initFluorescencePacket%position%x, initFluorescencePacket%xP(igpi))
               if (initFluorescencePacket%xP(igpi) < grid(gP)%nx) then
                  if (initFluorescencePacket%position%x >= (grid(gP)%xAxis(initFluorescencePacket%xP(igpi))+&
                       & grid(gP)%xAxis(initFluorescencePacket%xP(igpi)+1))/2.) &
                       & initFluorescencePacket%xP(igpi) = initFluorescencePacket%xP(igpi)+1
               end if

               ! y-direction
               initFluorescencePacket%position%y = 0.
               initFluorescencePacket%yP(igpi) = 1

               ! z-direction
               call random_number(random)
               random = 1. - random
               initFluorescencePacket%position%z = -(grid(gP)%zAxis(2)-grid(gP)%zAxis(1))/2. + random*( &
                    & (grid(gP)%zAxis(2)-grid(gP)%zAxis(1))/2.+&
                    & (grid(gP)%zAxis(grid(gP)%nz)-grid(gP)%zAxis(grid(gP)%nz-1))/2.+&
                    & grid(gP)%zAxis(grid(gP)%nz))
               if (initFluorescencePacket%position%z<grid(gP)%zAxis(1)) &
                    & initFluorescencePacket%position%z=grid(gP)%zAxis(1)
               if (initFluorescencePacket%position%z>grid(gP)%zAxis(grid(gP)%nz)) &
                    & initFluorescencePacket%position%z=grid(gP)%zAxis(grid(gP)%nz)

              call locate(grid(gP)%zAxis, initFluorescencePacket%position%z, initFluorescencePacket%zP(igpi))
               if (initFluorescencePacket%zP(igpi) < grid(gP)%nz) then
                  if (initFluorescencePacket%position%z >= (grid(gP)%xAxis(initFluorescencePacket%zP(igpi))+&
                       & grid(gP)%zAxis(initFluorescencePacket%zP(igpi)+1))/2.) initFluorescencePacket%zP(igpi) =&
                       & initFluorescencePacket%zP(igpi)+1
               end if

               if (initFluorescencePacket%xP(igpi)<1) initFluorescencePacket%xP(igpi)=1
               if (initFluorescencePacket%zP(igpi)<1) initFluorescencePacket%zP(igpi)=1

               ! direction is parallel to y-axis direction
               initFluorescencePacket%direction%x = 0.
               initFluorescencePacket%direction%y = 1.
               initFluorescencePacket%direction%z = 0.

               if (initFluorescencePacket%xP(igpi) >  grid(gP)%xAxis(grid(gP)%nx) .or. &
                    & initFluorescencePacket%zP(igpi) >  grid(gP)%zAxis(grid(gP)%nz)) then
                  print*, "! initFluorescencePacket: insanity in planeIonisation init"
                  print*, igpi, initFluorescencePacket%xP(igpi),  grid(gP)%xAxis(grid(gP)%nx), &
                       & initFluorescencePacket%zP(igpi), grid(gP)%zAxis(grid(gP)%nz),  &
                       & random, initFluorescencePacket%position%z

                  stop
               end if

                planeIonDistribution(initFluorescencePacket%xP(igpi),initFluorescencePacket%zP(igpi)) = &
                     & planeIonDistribution(initFluorescencePacket%xP(igpi),initFluorescencePacket%zP(igpi)) + 1

             else

                do irepeat = 1, 1000000
                   ! get a random direction
                   initFluorescencePacket%direction = randomUnitVector()
                   if (initFluorescencePacket%direction%x/=0. .and. &
                        & initFluorescencePacket%direction%y/=0. .and. &
                        & initFluorescencePacket%direction%z/=0.) exit
                end do
            end if

            if ((lgSymmetricXYZ) .and. initFluorescencePacket%lgStellar .and. nStars==1) then
                if (initFluorescencePacket%direction%x<0.) &
                     & initFluorescencePacket%direction%x = -initFluorescencePacket%direction%x
                if (initFluorescencePacket%direction%y<0.) &
                     & initFluorescencePacket%direction%y = -initFluorescencePacket%direction%y
                if (initFluorescencePacket%direction%z<0.) &
                     & initFluorescencePacket%direction%z = -initFluorescencePacket%direction%z
            end if

            initFluorescencePacket%origin(1) = gP
            initFluorescencePacket%origin(2) = grid(gP)%active(initFluorescencePacket%xP(igpi),&
                 & initFluorescencePacket%yP(igpi), initFluorescencePacket%zP(igpi))

!print*, 'origin ',  initFluorescencePacket%origin

!print*, 'inout'
          end function initFluorescencePacket

          ! this function creates a new fluorescence packet
          function newFluorescencePacket(chType, gP, orXf,orYf,orZf, rvec)

            type(photon_packet)                :: newFluorescencePacket! the photon packet to be created
            type(vector), intent(in)           :: rVec           !

            character(len=4), intent(in)       :: chType         ! what line?

            ! local variables
            type(vector)                       :: positionLoc    ! the position of the photon
                                                                 ! packet

            integer                            :: nuP            ! the frequency index of the photon packet
            integer, dimension(2),intent(inout):: orXf,orYf,orZf ! dummy
!            integer, intent(in)                :: difSource(3)   ! grid and cell indeces
            integer, intent(inout)             :: gP
            integer                            :: igpn           ! grid pointe 1=motehr, 2=sub

            if (gP==1) then
               igpn = 1
            else if (gp>1) then
               igpn = 2
            else
               print*,  "! newFluorescencePacket: insane grid pointer"
               stop
            end if


            positionLoc%x = rVec%x
            positionLoc%y = rVec%y
            positionLoc%z = rVec%z

            ! initialize the new fluorescence packet
            if (grid(gP)%active(orXf(igpn), orYf(igpn), orZf(igpn)) < 0.) then
               print*, "! newFluorescencePacket: new packet cannot be emitted from re-mapped cell -3-"
               print*, "chType, nuP,  .false., .true., orXf,orYf,orZf, gp"
               print*, chType, nuP,  .false., .true., orXf,orYf,orZf, gp
               stop
            end if

            select case (chType)
            case ("FeKa")

               newFluorescencePacket = initFluorescencePacket(FeKaColdP, positionLoc, .false., .false., orXf,&
                    & orYf, orZf, gP)

            case ("FeL1")

               newFluorescencePacket = initFluorescencePacket(FeL1P, positionLoc, .false., .false., orXf,&
                    & orYf, orZf, gP)

            case ("FeL2")

               newFluorescencePacket = initFluorescencePacket(FeL2P, positionLoc, .false., .false., orXf,&
                    & orYf, orZf, gP)

            case ("CKa")

               newFluorescencePacket = initFluorescencePacket(CKaP, positionLoc, .false., .false., orXf,&
                    & orYf, orZf, gP)

            case ("NKa")

               newFluorescencePacket = initFluorescencePacket(NKaP, positionLoc, .false., .false., orXf,&
                    & orYf, orZf, gP)

            case ("OKa")

               newFluorescencePacket = initFluorescencePacket(OKaP, positionLoc, .false., .false., orXf,&
                    & orYf, orZf, gP)

            case ("NeKa")

               newFluorescencePacket = initFluorescencePacket(NeKaP, positionLoc, .false., .false., orXf,&
                    & orYf, orZf, gP)

            case ("MgKa")

               newFluorescencePacket = initFluorescencePacket(MgKaP, positionLoc, .false., .false., orXf,&
                    & orYf, orZf, gP)

            case ("AlKa")

               newFluorescencePacket = initFluorescencePacket(AlKaP, positionLoc, .false., .false., orXf,&
                    & orYf, orZf, gP)

            case ("SiKa")

               newFluorescencePacket = initFluorescencePacket(SiKaP, positionLoc, .false., .false., orXf,&
                    & orYf, orZf, gP)

            case ("SKa")

               newFluorescencePacket = initFluorescencePacket(SKaP, positionLoc, .false., .false., orXf,&
                    & orYf, orZf, gP)

            case ("ArKa")

               newFluorescencePacket = initFluorescencePacket(ArKaP, positionLoc, .false., .false., orXf,&
                    & orYf, orZf, gP)

            case ("CaKa")

               newFluorescencePacket = initFluorescencePacket(CaKaP, positionLoc, .false., .false., orXf,&
                    & orYf, orZf, gP)

            case ("TiKa")

               newFluorescencePacket = initFluorescencePacket(TiKaP, positionLoc, .false., .false., orXf,&
                    & orYf, orZf, gP)

            ! if the photon packet type is wrong or missing
            case default

                print*, "! newFluorescencePacket: unsupported fluorescent transition.", chType
                stop

             end select

           end function newFluorescencePacket

           subroutine fluoPathSegment(enPacket)
             implicit none

             type(photon_packet), intent(inout) :: enPacket ! the energy packet

             ! local variables
             type(vector)                    :: vHat     ! direction vector
             type(vector)                    :: rVec     ! position vector

             real                            :: absTau   ! optical depth
             real                            :: dlLoc    ! local displacement
             real                            :: dx, dy, dz
             real                            :: dSx, dSy, dSz
                                                      ! distances from x,y and z wall
             real                            :: dS       ! distance from nearest wall
             real                            :: dV       ! lume of this cell
             real                            :: passProb ! prob of passing the next segment
             real                            :: probSca  ! prob that the packet scatters
             real                            :: radius   ! radius
             real                            :: random   ! random number
             real                            :: tauCell  ! local tau

             integer                         :: idirT,idirP ! direction cosine counters
             integer                         :: i, j     ! counter
             integer                         :: xP,yP,zP ! cartesian axes indeces
             integer                         :: gP       ! grid index
             integer                         :: igpp     ! grid index 1=mother 2=sub
             integer                         :: safeLimit =10000 ! safe limit for the loop

             character(len=7)                :: packetType ! what line?

             logical                         :: lgScattered ! is the packet scattering with dust?
             logical                         :: lgReturn
!print*, 'ps1'

             if (enPacket%iG == 1) then
                igpp = 1
             else if (enPacket%iG>1) then
                igpp = 2
             else
                print*, "! fluoPathSegment: insane grid index"
                stop
             end if

             ! check that the input position is not outside the grid
             if ( (enPacket%iG <= 0).or.(enPacket%iG > nGrids) ) then
                print*, "! fluoPathSegment: starting position not in any defined gridhses",&
                     & enPacket
                stop
             else if ( (enPacket%xP(igpp) <= 0).or.&
                  &(enPacket%xP(igpp) > grid(enPacket%iG)%nx) ) then
                print*, "! fluoPathSegment: starting position in x is outside the grid",&
                     & enPacket
                stop
             else if ( (enPacket%yP(igpp) <= 0).or. &
                  & (enPacket%yP(igpp) > grid(enPacket%iG)%ny) ) then
                print*, "! fluoPathSegment: starting position in y is outside the grid",&
                     & enPacket
                stop
             else if ( (enPacket%zP(igpp) <= 0).or.&
                  & (enPacket%zP(igpp) > grid(enPacket%iG)%nz) ) then
                print*, "! fluoPathSegment: starting position in z is outside the grid",&
                     & enPacket
                stop
             end if

             ! define vHat and rVec
             rVec = enPacket%position
             vHat = enPacket%direction

             ! initialize xP, yP,zP
             xP = enPacket%xP(igpp)
             yP = enPacket%yP(igpp)
             zP = enPacket%zP(igpp)
             gP = enPacket%iG

             ! initialise distance from walls
             dSx = 0.
             dSy = 0.
             dSz = 0.

             if (lg1D) then
                radius = 1.e10*sqrt((rVec%x/1.e10)*(rVec%x/1.e10) + &
                     &                               (rVec%y/1.e10)*(rVec%y/1.e10) + &
                     &                               (rVec%z/1.e10)*(rVec%z/1.e10))
                call locate(grid(1)%xAxis, radius, xP)
                if (nGrids > 1 .or. gP >1) then
                   print*, " ! fluorescencePacketRun: multiple grids are not allowed in a 1D simulation"
                   stop
                end if
             end if

             ! initialize optical depth
             absTau = 0.

             ! get a random number
             call random_number(random)

             ! calculate the probability
             passProb = -log(1.-random)

             ! speed up photons that my be trapped
             if (lgPlaneIonization) then
                safeLimit=5000
             else
!             safeLimit=500000
                safeLimit=10000
             end if

             do i = 1, safeLimit
!print*, 'ps2', i

                if (xP <  grid(gP)%nx) then
                   if (vHat%x > 0. .and. abs(rVec%x - &
                        & (grid(gP)%xAxis(xP)+grid(gP)%xAxis(xP+1))/2.)/ &
                        & ((grid(gP)%xAxis(xP)+grid(gP)%xAxis(xP+1))/2.) < 2.e-2) &
                        & xP = xP+1
                end if
                if (xP > 1 ) then
                   if (vHat%x < 0. .and. abs (rVec%x  -&
                        & (grid(gP)%xAxis(xP)+grid(gP)%xAxis(xP-1))/2.)/ &
                        & ((grid(gP)%xAxis(xP)+grid(gP)%xAxis(xP-1))/2.) < 2.e-2) &
                        & xP = xP-1
                end if

                if (yP <  grid(gP)%ny) then
                   if (vHat%y > 0. .and. abs (rVec%y - &
                        & (grid(gP)%yAxis(yP)+grid(gP)%yAxis(yP+1))/2.)/ &
                        & ((grid(gP)%yAxis(yP)+grid(gP)%yAxis(yP+1))/2.) < 2.e-2) &
                     & yP = yP+1
                end if
                if (yP > 1 ) then
                   if (vHat%y < 0. .and. abs(rVec%y - &
                        & (grid(gP)%yAxis(yP)+grid(gP)%yAxis(yP-1))/2.)/ &
                        & ((grid(gP)%yAxis(yP)+grid(gP)%yAxis(yP-1))/2.) < 2.e-2) &
                        & yP = yP-1
                end if

                if (zP <  grid(gP)%nz) then
                   if (vHat%z > 0. .and. abs(rVec%z - &
                        & (grid(gP)%zAxis(zP)+grid(gP)%zAxis(zP+1))/2.)/ &
                        & ((grid(gP)%zAxis(zP)+grid(gP)%zAxis(zP+1))/2.) < 2.e-2) &
                        & zP = zP+1
                end if
                if (zP > 1 ) then
                   if (vHat%z < 0. .and. abs(rVec%z - &
                        & (grid(gP)%zAxis(zP)+grid(gP)%zAxis(zP-1))/2.)/ &
                        & ((grid(gP)%zAxis(zP)+grid(gP)%zAxis(zP-1))/2.) < 2.e-2) &
                        & zP = zP-1
                end if

                do j = 1, safeLimit

                   if (xP > grid(gP)%nx .or. xP < 1 .or. &
                        & yP > grid(gP)%ny .or. yP < 1 .or. &
                        & zP > grid(gP)%nz .or. zP < 1 ) then
                      print*, "! fluoPathSegment: insanity [gp,xp,yp,zp,j,i]", &
                           & gp, xp, yp, zp, j, i
                      stop
                   end if

                   if (grid(gP)%active(xP,yP,zP)<0) then

                      ! packet is entering a subgrid
                      enPacket%xP(1) = xP
                      enPacket%yP(1) = yP
                      enPacket%zP(1) = zP

                      gP = abs(grid(gP)%active(xP,yP,zP))

                      ! where is the packet in the sub-grid?

                      call locate(grid(gP)%xAxis, rVec%x, xP)
                      if (xP==0) xP = xP+1
                      if (xP< grid(gP)%nx) then
                         if (rVec%x >  (grid(gP)%xAxis(xP+1)+grid(gP)%xAxis(xP))/2.) &
                              & xP = xP + 1
                      end if

                      call locate(grid(gP)%yAxis, rVec%y, yP)
                      if (yP==0) yP=yP+1
                      if (yP< grid(gP)%ny) then
                         if (rVec%y >  (grid(gP)%yAxis(yP+1)+grid(gP)%yAxis(yP))/2.) &
                              & yP = yP + 1
                      end if

                      call locate(grid(gP)%zAxis, rVec%z, zP)
                      if (zP==0) zP=zP+1
                      if (zP< grid(gP)%nz) then
                         if (rVec%z >  (grid(gP)%zAxis(zP+1)+grid(gP)%zAxis(zP))/2.) &
                              & zP = zP + 1
                      end if

                   end if

                   enPacket%iG = gP
                   if (gP > 1) then
                      igpp = 2
                   else
                      igpp = 1
                   end if


                   ! find distances from all walls
!print*, 'ps3', xp,yp,zp
                   if (lgSymmetricXYZ) then
                      if ( rVec%x <= grid(1)%xAxis(1) ) then
                         if (vHat%x<0.) vHat%x = -vHat%x
                         rVec%x = grid(1)%xAxis(1)
                      end if
                      if ( rVec%y <= grid(1)%yAxis(1) ) then
                         if (vHat%y<0.) vHat%y = -vHat%y
                         rVec%y = grid(1)%yAxis(1)
                      end if
                      if ( rVec%z <= grid(1)%zAxis(1) ) then
                         if (vHat%z<0.) vHat%z = -vHat%z
                         rVec%z = grid(1)%zAxis(1)
                      end if
                   end if

                   if (vHat%x>0.) then
                      if (xP<grid(gP)%nx) then

                         dSx = ( (grid(gP)%xAxis(xP+1)+grid(gP)%xAxis(xP))/2.-rVec%x)/vHat%x

                         if (abs(dSx)<1.e-10) then
                            rVec%x=(grid(gP)%xAxis(xP+1)+grid(gP)%xAxis(xP))/2.
                            xP = xP+1
                         end if
                      else
                         dSx = ( grid(gP)%xAxis(grid(gP)%nx)-rVec%x)/vHat%x
                         if (abs(dSx)<1.e-10) then
                            rVec%x=grid(gP)%xAxis(grid(gP)%nx)
                            if (.not.lgPlaneIonization .and. gP==1) return
                         end if
                      end if
                   else if (vHat%x<0.) then
                      if (xP>1) then
                         dSx = ( (grid(gP)%xAxis(xP)+grid(gP)%xAxis(xP-1))/2.-rVec%x)/vHat%x
                         if (abs(dSx)<1.e-10) then
                            rVec%x=(grid(gP)%xAxis(xP)+grid(gP)%xAxis(xP-1))/2.
                            xP = xP-1
                         end if
                      else
                         dSx = (grid(gP)%xAxis(1)-rVec%x)/vHat%x
                         if (abs(dSx)<1.e-10) then
                            rVec%x=grid(gP)%xAxis(1)
                         end if
                      end if
                   else if (vHat%x==0.) then
                      dSx = grid(gP)%xAxis(grid(gP)%nx)
                   end if

                   if (.not.lg1D) then
                      if (vHat%y>0.) then
                         if (yP<grid(gP)%ny) then
                            dSy = ( (grid(gP)%yAxis(yP+1)+grid(gP)%yAxis(yP))/2.-rVec%y)/vHat%y
                            if (abs(dSy)<1.e-10) then
                               rVec%y=(grid(gP)%yAxis(yP+1)+grid(gP)%yAxis(yP))/2.
                               yP = yP+1
                            end if
                         else
                            dSy = (  grid(gP)%yAxis(grid(gP)%ny)-rVec%y)/vHat%y
                            if (abs(dSy)<1.e-10) then
                               rVec%y=grid(gP)%yAxis(grid(gP)%ny)
                               if(gP==1) return
                            end if
                         end if
                      else if (vHat%y<0.) then
                         if (yP>1) then
                            dSy = ( (grid(gP)%yAxis(yP)+grid(gP)%yAxis(yP-1))/2.-rVec%y)/vHat%y
                            if (abs(dSy)<1.e-10) then
                               rVec%y=(grid(gP)%yAxis(yP)+grid(gP)%yAxis(yP-1))/2.
                               yP = yP-1
                            end if
                         else
                            dSy = ( grid(gP)%yAxis(1)-rVec%y)/vHat%y
                            if (abs(dSy)<1.e-10) then
                               rVec%y=grid(gP)%yAxis(1)
                            end if
                         end if
                      else if (vHat%y==0.) then
                         dSy = grid(gP)%yAxis(grid(gP)%ny)
                      end if

                      if (vHat%z>0.) then
                         if (zP<grid(gP)%nz) then
                            dSz = ( (grid(gP)%zAxis(zP+1)+grid(gP)%zAxis(zP))/2.-rVec%z)/vHat%z
                            if (abs(dSz)<1.e-10) then
                               rVec%z=(grid(gP)%zAxis(zP+1)+grid(gP)%zAxis(zP))/2.
                               zP = zP+1
                            end if
                         else
                            dSz = ( grid(gP)%zAxis(grid(gP)%nz)-rVec%z)/vHat%z
                            if (abs(dSz)<1.e-10) then
                               rVec%z=grid(gP)%zAxis(grid(gP)%nz)
                               if (.not.lgPlaneIonization .and. gP==1) return
                            end if
                         end if
                      else if (vHat%z<0.) then
                         if (zP>1) then
                            dSz = ( (grid(gP)%zAxis(zP)+grid(gP)%zAxis(zP-1))/2.-rVec%z)/vHat%z
                            if (abs(dSz)<1.e-10) then
                               rVec%z=(grid(gP)%zAxis(zP)+grid(gP)%zAxis(zP-1))/2.
                               zP = zP-1
                            end if
                         else
                            dSz = ( grid(gP)%zAxis(1)-rVec%z)/vHat%z
                            if (abs(dSz)<1.e-10) then
                               rVec%z=grid(gP)%zAxis(1)
                            end if
                         end if
                      else if (vHat%z==0.) then
                         dSz = grid(gP)%zAxis(grid(gP)%nz)
                      end if

                      if (xP > grid(gP)%nx .or. xP < 1 .or. &
                           & yP > grid(gP)%ny .or. yP < 1 .or. &
                           & zP > grid(gP)%nz .or. zP < 1 ) then
                         print*, "! fluoPathSegment: insanity -2- [gp,xp,yp,zp]", &
                              & gp, xp, yp, zp
                         stop
                      end if

                   end if

                   if (grid(gP)%active(xP,yP,zP)>=0) exit
                end do

                ! cater for cells on cell wall
                if ( abs(dSx)<1.e-10 ) dSx = grid(gP)%xAxis(grid(gP)%nx)
                if ( abs(dSy)<1.e-10 ) dSy = grid(gP)%yAxis(grid(gP)%ny)
                if ( abs(dSz)<1.e-10 ) dSz = grid(gP)%zAxis(grid(gP)%nz)

                ! find the nearest wall
                dSx = abs(dSx)
                dSy = abs(dSy)
                dSz = abs(dSz)

                if (dSx<=0.) then
                   print*, '! fluoPathSegment: [warning] dSx <= 0.',dSx
                   print*, 'grid(gP)%xAxis ', grid(gP)%xAxis
                   print*, 'gP,xP,grid(gP)%xAxis(xP), rVec%x, vHat%x'
                   print*, gP,xP,grid(gP)%xAxis(xP), rVec%x, vHat%x
                   dS = min(dSy, dSz)
                else if (dSy<=0.) then
                   print*, '! fluoPathSegment: [warning] dSy <= 0.', dSy
                   print*, 'grid(gP)%yAxis ', grid(gP)%yAxis
                   print*, 'gP,yP,grid(gP)%yAxis(yP), rVec%y, vHat%y'
                   print*, gP,yP,grid(gP)%yAxis(yP), rVec%y, vHat%y
                   dS = min(dSx, dSz)
                else if (dSz<=0.) then
                   print*, '! fluoPathSegment: [warning] dSz <= 0.', dSz
                   print*, 'grid(gP)%zAxis ', grid(gP)%zAxis
                   print*, 'gP,zP,grid(gP)%zAxis(zP), rVec%z, vHat%z'
                   print*, gP,zP,grid(gP)%zAxis(zP), rVec%z, vHat%z
                   dS = min(dSx, dSy)
                else
                   dS = min(dSx,dSy)
                   dS = min(dS, dSz)
                end if

                ! this should now never ever happen
                if (dS <= 0.) then
                   print*, 'fluoPathSegment: dS <= 0', dSx, dSy, dSz
                   print*, gP, rVec
                   stop
                end if
!print*, 'ps4'
                ! calculate the optical depth to the next cell wall
                tauCell = dS*grid(gP)%opacity(grid(gP)%active(xP,yP,zP), enPacket%nuP)

!print*, 'ps5'
                ! find the volume of this cell
!             dV = getVolumeLoc(grid(gP), xP,yP,zP)


                if (lg1D) then
                   if (nGrids>1) then
                      print*, '! getVolumeLoc: 1D option and multiple grids options are not compatible'
                      stop
                   end if

                   if (xP == 1) then

                      dV = 4.*Pi* ( (grid(gP)%xAxis(xP+1)/1.e15)**3)/3.


                   else if ( xP==grid(gP)%nx) then

                      dV = Pi* ( (3.*(grid(gP)%xAxis(xP)/1.e15)-(grid(gP)%xAxis(xP-1)/1.e15))**3 - &
                           & ((grid(gP)%xAxis(xP)/1.e15)+(grid(gP)%xAxis(xP-1)/1.e15))**3 ) / 6.

                   else

                      dV = Pi* ( ((grid(gP)%xAxis(xP+1)/1.e15)+(grid(gP)%xAxis(xP)/1.e15))**3 - &
                           & ((grid(gP)%xAxis(xP-1)/1.e15)+(grid(gP)%xAxis(xP)/1.e15))**3 ) / 6.

                   end if

                   dV = dV/8.

                else

                   if ( (xP>1) .and. (xP<grid(gP)%nx) ) then

                      dx = abs(grid(gP)%xAxis(xP+1)-grid(gP)%xAxis(xP-1))/2.
                   else if ( xP==1 ) then
                      if (lgSymmetricXYZ .or. gP>1 .or. lgPlaneIonization ) then
                         dx = abs(grid(gP)%xAxis(xP+1)-grid(gP)%xAxis(xP))/2.
                      else
                         dx = abs(grid(gP)%xAxis(xP+1)-grid(gP)%xAxis(xP))
                      end if
                   else if ( xP==grid(gP)%nx ) then
                      if (gP > 1 .or. lgPlaneIonization) then
                         dx = abs(grid(gP)%xAxis(xP)  -grid(gP)%xAxis(xP-1))/2.
                      else
                         dx = abs(grid(gP)%xAxis(xP)  -grid(gP)%xAxis(xP-1))
                      end if
                   end if

                   if ( (yP>1) .and. (yP<grid(gP)%ny) ) then
                      dy = abs(grid(gP)%yAxis(yP+1)-grid(gP)%yAxis(yP-1))/2.
                   else if ( yP==1 ) then
                      if (lgSymmetricXYZ .or. gP>1 .or. lgPlaneIonization) then
                         dy = abs(grid(gP)%yAxis(yP+1)-grid(gP)%yAxis(yP))/2.
                      else
                         dy = abs(grid(gP)%yAxis(yP+1)-grid(gP)%yAxis(yP))
                      end if
                   else if ( yP==grid(gP)%ny ) then
                      if (gP > 1 .or. lgPlaneIonization) then
                         dy = abs(grid(gP)%yAxis(yP)  -grid(gP)%yAxis(yP-1))/2.
                      else
                         dy = abs(grid(gP)%yAxis(yP)  -grid(gP)%yAxis(yP-1))
                      end if
                   end if

                   if ( (zP>1) .and. (zP<grid(gP)%nz) ) then
                      dz = abs(grid(gP)%zAxis(zP+1)-grid(gP)%zAxis(zP-1))/2.
                   else if ( zP==1 ) then
                      if (lgSymmetricXYZ .or. gP>1 .or. lgPlaneIonization) then
                         dz = abs(grid(gP)%zAxis(zP+1)-grid(gP)%zAxis(zP))/2.
                      else
                         dz = abs(grid(gP)%zAxis(zP+1)-grid(gP)%zAxis(zP))
                      end if
                   else if ( zP==grid(gP)%nz ) then
                      if (gP>1 .or. lgPlaneIonization) then
                         dz = abs(grid(gP)%zAxis(zP)-grid(gP)%zAxis(zP-1))/2.
                      else
                         dz = abs(grid(gP)%zAxis(zP)-grid(gP)%zAxis(zP-1))
                      end if
                   end if

                   dx = dx/1.e15
                   dy = dy/1.e15
                   dz = dz/1.e15

                   ! calculate the volume
                   dV = dx*dy*dz

                end if
!print*, 'ps4', passProb,   absTau+tauCell
                ! check if the packet interacts within this cell
                if ((absTau+tauCell > passProb) .and. (grid(gP)%active(xP,yP,zP)>0)) then

                   ! packet interacts

                   ! calculate where within this cell the packet interacts
                   dlLoc = (passProb-absTau)/grid(gP)%opacity(grid(gP)%active(xP,yP,zP), enPacket%nuP)

                   ! update packet's position
                   rVec = rVec + dlLoc*vHat

                   if (lgSymmetricXYZ .and. gP==1) then
                      if ( rVec%x <= grid(gP)%xAxis(1) ) then
                         if (vHat%x<0.) vHat%x = -vHat%x
                         rVec%x = grid(gP)%xAxis(1)
                      end if
                      if ( rVec%y <= grid(gP)%yAxis(1) ) then
                         if (vHat%y<0.) vHat%y = -vHat%y
                         rVec%y = grid(gP)%yAxis(1)
                      end if
                      if ( rVec%z <= grid(gP)%zAxis(1) ) then
                         if (vHat%z<0.) vHat%z = -vHat%z
                         rVec%z = grid(gP)%zAxis(1)
                      end if
                   end if

                   ! add contribution of the packet to the radiation field
                   if (enPacket%lgStellar) then
                      grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                           grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dlLoc*deltaEUsed / dV
                   else ! if the energy packet is diffuse
                      if (lgDebug) then
                         grid(gP)%Jdif(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                              & grid(gP)%Jdif(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dlLoc*deltaEUsed / dV
                      else
                         grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                              & grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dlLoc*deltaEUsed / dV
                      end if
                   end if

                   ! check if the position within the cell is still within the outer radius
                   if ( sqrt( (rvec%x/1.e10)**2 + (rvec%y/1.e10)**2 + (rvec%z/1.e10)**2)*1.e10 >= R_out &
                        & .and. R_out > 0.) then

!print*, 'ps5 esc'
                      ! the packet escapes without further interaction
                      idirT = int(acos(enPacket%direction%z)/dTheta)+1
                      if (idirT>totangleBinsTheta) then
                         idirT=totangleBinsTheta
                      end if
                      if (idirT<1 .or. idirT>totAngleBinsTheta) then
                         print*, '! fluorescencePacketRun: error in theta direction cosine assignment',&
                              &  idirT, enPacket, dTheta, totAngleBinsTheta
                         stop
                      end if

                      if (enPacket%direction%x<1.e-35) then
                         idirP = 0
                      else
                         idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
                      end if
                      if (idirP<0) idirP=totAngleBinsPhi+idirP
                      idirP=idirP+1
                      if (idirP>totangleBinsPhi) then
                         idirP=totangleBinsPhi
                      end if

                      if (idirP<1 .or. idirP>totAngleBinsPhi) then
                         print*, '! fluorescencePacketRun: error in Phi direction cosine assignment',&
                              &  idirP, enPacket, dPhi, totAngleBinsPhi
                         stop
                      end if

                      if (nAngleBins>0) then
                         if ( (totangleBinsPhi>1 .and. viewPointPtheta(idirT) == viewPointPphi(idirP)) .or. &
                              & (totangleBinsPhi==1 .and. viewPointPtheta(idirT)>0) ) then
                            grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,&
                                 & viewPointPtheta(idirT)) = &
                                 &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                 & enPacket%nuP,viewPointPtheta(idirT)) +  deltaEUsed
                            grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,0) = &
                                 &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                 & enPacket%nuP,0) +  deltaEUsed

                            if (lgSeparateSED) then
                               grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), enPacket%nuP,&
                                    & viewPointPtheta(idirT),enPacket%SEDtype) = &
                                    &grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                    & enPacket%nuP,viewPointPtheta(idirT),enPacket%SEDtype) +  deltaEUsed
                               grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                    & enPacket%nuP,0,enPacket%SEDtype) = &
                                    &grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                    & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                            endif

                         else
                            grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2),enPacket%nuP,0) = &
                                 & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                 & enPacket%nuP,0) +  deltaEUsed

                            if (lgSeparateSED) then
                               grid(enPacket%origin(1))%escapedPacketsComponents&
                                    &(enPacket%origin(2),enPacket%nuP,0,enPacket%SEDtype) = &
                                    & grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                    & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                            endif

                         end if
                      else
                         grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                              & enPacket%nuP,0) = &
                              & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                              & enPacket%nuP,0) +  deltaEUsed
                         if (lgSeparateSED) then
                            grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                 & enPacket%nuP,0,enPacket%SEDtype) = &
                                 & grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                 & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                         endif

                      end if

!print*, 'ps5 esc out'
                      return
                   end if

!print*, 'ps6 ', probSca, KNsigmaT(enPacket%nuP)*&
!                        &grid(gP)%Ne(grid(gP)%active(xp,yp,zp)), (grid(gp)%opacity(grid(gp)%active(xp,yp,zp),&
!                        &enPacket%nuP))
                   ! check if the packet is absorbed or scattered
                   ! NOTE : we are ignoring dust
                   probSca = KNsigmaT(enPacket%nuP)*&
                        &grid(gP)%Ne(grid(gP)%active(xp,yp,zp))/&
                        & (grid(gp)%opacity(grid(gp)%active(xp,yp,zp),&
                        &enPacket%nuP))

                   call random_number(random)

                   random = 1.-random

                   if (random > probSca) then
                      lgScattered = .false.
                   else if (random <= probSca) then
                      lgScattered = .true.
                   else
                      print*, '! fluoPathSegment: insanity occured and scattering/absorption &
                           & decision stage.'
                      stop
                   end if

!print*, 'ps7 sca', lgscattered
                   if (.not. lgScattered) then

                      absInt = absInt + 1.
                      exit

                   else

                      scaInt = scaInt + 1.
                      ! packet is compton scattered
                      ! calculate new direction
                      ! for now assume scattering is isotropic, when phase
                      ! function is introduced the following must be changed

                      enPacket%xP(igpp) = xP
                      enPacket%yP(igpp) = yP
                      enPacket%zP(igpp) = zP

                      if (grid(gP)%active(enPacket%xp(igpp), enPacket%yp(igpp), enPacket%zp(igpp)) < 0.) then
                         print*, "! fluoPathSegment: new packet cannot be emitted from re-mapped cell -1-"
                         print*, "nuP, .false., .true., xp,yp,zp, gp"
                         print*, nuP, .false., .true.,  xp,yp,zp, gp
                         stop
                      end if

                      call comptonScatter(enPacket)

!                      enPacket = initFluorescencePacket(enPacket%nuP, rVec, .false., .false., enPacket%xP(1:2), &
!                           & enPacket%yP(1:2), enPacket%zP(1:2), gP)

                      enPacket%origin(1) = gP
                      enPacket%origin(2) = grid(gP)%active(enPacket%xp(igpp), enPacket%yp(igpp), enPacket%zp(igpp))


                      vHat%x = enPacket%direction%x
                      vHat%y = enPacket%direction%y
                      vHat%z = enPacket%direction%z

                      ! initialize optical depth
                      absTau = 0.

                      ! get a random number
                      call random_number(random)

                      ! calculate the probability
                      passProb = -log(1.-random)
!print*, 'sca out'
                   end if

                else

                   ! the packet is not absorbed within this cell
                   ! add contribution of the packet to the radiation field

                   if (enPacket%lgStellar) then
                      grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                           grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dS*deltaEUsed / dV
                   else ! if the energy packet is diffuse
                      if (lgDebug) then
                         grid(gP)%Jdif(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                              & grid(gP)%Jdif(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dS*deltaEUsed / dV
                      else
                         grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                              & grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dS*deltaEUsed / dV
                      end if
                   end if

                   ! update absTau
                   absTau = absTau+tauCell

                   ! update packet's position
                   rVec = rVec+dS*vHat

                   ! keep track of where you are on mother grid
                   if (gP>1) then

                      if (enPacket%xP(1) <= 0 .or. &
                           & enPacket%yP(1) <= 0 .or. &
                           & enPacket%zP(1) <= 0 ) then

                         ! locate where we are at on the mother grid
                         call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                         if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                            if ( rvec%x > ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                                 & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                                 & enPacket%xP(1) = enPacket%xP(1)+1
                         end if

                         call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                         if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                            if ( rvec%y > ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                                 & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                                 & enPacket%yP(1) = enPacket%yP(1)+1
                         end if

                         call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                         if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                            if ( rvec%z > ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                                 & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                                 & enPacket%zP(1) = enPacket%zP(1)+1
                         end if


                      else

!print*, 'ps8'
                         if (vHat%x>0.) then
                            if ( enPacket%xP(1) < grid(grid(gP)%motherP)%nx ) then
                               if ( rVec%x > (grid(grid(gP)%motherP)%xAxis(enPacket%xP(1))+&
                                    & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1))/2. ) then
                                  enPacket%xP(1) = enPacket%xP(1)+1
                               end if
                            else
                               if ( rVec%x > grid(grid(gP)%motherP)%xAxis(enPacket%xP(1))) then
!                            print*, '! fluoPathSegment: insanity occured at mother grid transfer (x axis +)', &
!                                 & rVec%x, gP, grid(gP)%motherP
!                            stop
                               end if
                            end if
                         else
                            if ( enPacket%xP(1) > 1 ) then
                               if ( rVec%x <= (grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)-1)+&
                                    & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)))/2. ) then
                                  enPacket%xP(1) = enPacket%xP(1)-1
                               end if
                            else
                               if (rVec%x < grid(grid(gP)%motherP)%xAxis(1)) then
!                            print*, '! fluoPathSegment: insanity occured at mother grid transfer (x axis-)',&
!                                 & rVec%x, gP, grid(gP)%motherP
!                            stop
                               end if
                            end if
                         end if
                         if (vHat%y>0.) then
                            if (  enPacket%yP(1) < grid(grid(gP)%motherP)%ny ) then
                               if ( rVec%y > (grid(grid(gP)%motherP)%yAxis( enPacket%yP(1))+&
                                    & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1))/2. ) then
                                  enPacket%yP(1) =  enPacket%yP(1)+1
                               end if
                            else
                               if ( rVec%y > grid(grid(gP)%motherP)%yAxis( enPacket%yP(1))) then
!                            print*, '! fluoPathSegment: insanity occured at mother grid transfer (y axis +)',&
!                                 & rVec%y, gP, grid(gP)%motherP
!                            stop
                               end if
                            end if
                         else
                            if (  enPacket%yP(1) > 1 ) then
                               if ( rVec%y <= (grid(grid(gP)%motherP)%yAxis( enPacket%yP(1)-1)+&
                                    & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)))/2. ) then
                                  enPacket%yP(1) =  enPacket%yP(1)-1
                               end if
                            else
                               if (rVec%y < grid(grid(gP)%motherP)%yAxis(1)) then
!                            print*, '! fluoPathSegment: insanity occured at mother grid transfer (y axis -)', &
!                                 & rVec%y, gP, grid(gP)%motherP
!                            stop
                               end if
                            end if
                         end if
                         if (vHat%z>0.) then
                            if (  enPacket%zP(1) < grid(grid(gP)%motherP)%nz ) then
                               if ( rVec%z > (grid(grid(gP)%motherP)%zAxis( enPacket%zP(1))+&
                                    & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1))/2. ) then
                                  enPacket%zP(1) =  enPacket%zP(1)+1
                               end if
                            else
                               if ( rVec%z > grid(grid(gP)%motherP)%zAxis( enPacket%zP(1))) then
!                            print*, '! fluoPathSegment: insanity occured at mother grid transfer (z axis +)', &
!                                 & rVec%z, gP, grid(gP)%motherP
!                            stop
                               end if
                            end if
                         else
                            if (  enPacket%zP(1) > 1 ) then
                               if ( rVec%z <= (grid(grid(gP)%motherP)%zAxis( enPacket%zP(1)-1)+&
                                    & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)))/2. ) then
                                  enPacket%zP(1) =  enPacket%zP(1)-1
                               end if
                            else
                               if (rVec%z < grid(grid(gP)%motherP)%zAxis(1)) then
!                        print*, '! fluoPathSegment: insanity occured at mother grid transfer (z axis -)', &
!                             & rVec%z, gP, grid(gP)%motherP
!                        stop
                               end if
                            end if
                         end if
!print*, 'ps9'
                      end if
                   end if

                   if (.not.lg1D) then
                      if ( (dS == dSx) .and. (vHat%x > 0.)  ) then
                         xP = xP+1
                      else if ( (dS == dSx) .and. (vHat%x < 0.) ) then
                         xP = xP-1
                      else if ( (dS == dSy) .and. (vHat%y > 0.) ) then
                         yP = yP+1
                      else if ( (dS == dSy) .and. (vHat%y < 0.) ) then
                         yP = yP-1
                      else if ( (dS == dSz) .and. (vHat%z > 0.) ) then
                         zP = zP+1
                      else if ( (dS == dSz) .and. (vHat%z < 0.) ) then
                         zP = zP-1
                      else
                         print*, '! fluoPathSegment: insanity occured in dS assignement &
                              & [dS,dSx,dSy,dSz,vHat]', dS,dSx,dSy,dSz,vHat
                      end if
                   else
                      radius = 1.e10*sqrt((rVec%x/1.e10)*(rVec%x/1.e10) + &
                           & (rVec%y/1.e10)*(rVec%y/1.e10) + &
                           & (rVec%z/1.e10)*(rVec%z/1.e10))
                      call locate(grid(gP)%xAxis, radius , xP)

                   end if
!print*, 'ps10'
                ! be 6/6/06
                   if(.not.lgPlaneIonization.and..not.lgSymmetricXYZ) then
                      lgReturn=.false.

                      if ( rVec%y <= grid(gP)%yAxis(1)-grid(gP)%geoCorrY .or. yP<1) then

                         ! the energy packet escapes this grid
                         if (gP==1) then
                            yP=1
                            lgReturn=.true.
                         else if (gP>1) then
                            xP = enPacket%xP(grid(gP)%motherP)
                            yP = enPacket%yP(grid(gP)%motherP)
                            zP = enPacket%zP(grid(gP)%motherP)
                            gP = grid(gP)%motherP
                         else
                            print*, '! fluoPathSegment: insanity occured - invalid gP', gP
                            stop
                         end if

                      end if
!print*, 'ps11', lgreturn
                     if (rVec%y > grid(gP)%yAxis(grid(gP)%ny)+grid(gP)%geoCorrY .or. yP>grid(gP)%ny) then

                        if (gP==1) then
                           ! the energy packet escapes
                           yP = grid(gP)%ny
                           lgReturn=.true.
                        else if (gP>1) then
                           xP = enPacket%xP(grid(gP)%motherP)
                           yP =  enPacket%yP(grid(gP)%motherP)
                           zP =  enPacket%zP(grid(gP)%motherP)
                           gP = grid(gP)%motherP
                        else
                           print*, '! fluoPathSegment: insanity occured - invalid gP', gP
                           stop
                        end if
!print*, 'ps12', lgreturn
                     end if

                     if ( (rVec%x <= grid(gP)%xAxis(1)-grid(gP)%geoCorrX .or. xP<1) .and. gP==1) then
                        xP=1
                        lgReturn=.true.
!print*, 'ps12', lgreturn
                     end if


                     if ( (rVec%x <= grid(gP)%xAxis(1)-grid(gP)%geoCorrX .or. xP<1) &
                          & .and. gP>1) then

                        xP = enPacket%xP(grid(gP)%motherP)
                        yP =  enPacket%yP(grid(gP)%motherP)
                        zP =  enPacket%zP(grid(gP)%motherP)
                        gP = grid(gP)%motherP

                     end if


                     if ( (rVec%x >=  grid(gP)%xAxis(grid(gP)%nx)+grid(gP)%geoCorrX &
                          & .or. xP>grid(gP)%nx) .and. gP==1 )then
                        xP = grid(gP)%nx
                        lgReturn=.true.
!print*, 'ps13', lgreturn
                     end if

                     if ( (rVec%x >=  grid(gP)%xAxis(grid(gP)%nx)+grid(gP)%geoCorrX&
                          & .or. xP>grid(gP)%nx) .and.  gP>1) then

                        xP = enPacket%xP(grid(gP)%motherP)
                        yP =  enPacket%yP(grid(gP)%motherP)
                        zP =  enPacket%zP(grid(gP)%motherP)
                        gP = grid(gP)%motherP

                     end if

                     if ( (rVec%z <= grid(gP)%zAxis(1)-grid(gP)%geoCorrZ .or.zP<1) &
                          & .and. gP==1) then
                        zP=1
                        lgReturn=.true.
!print*, 'ps14', lgreturn
                    end if

                    if ( (rVec%z <= grid(gP)%zAxis(1)-grid(gP)%geoCorrZ &
                         & .or.zP<1) .and. gP>1) then

                       xP = enPacket%xP(grid(gP)%motherP)
                       yP =  enPacket%yP(grid(gP)%motherP)
                       zP =  enPacket%zP(grid(gP)%motherP)
                       gP = grid(gP)%motherP

                    end if

                    if ( (rVec%z >=  grid(gP)%zAxis(grid(gP)%nz)+grid(gP)%geoCorrZ &
                         & .or. zP>grid(gP)%nz) &
                         & .and. gP==1) then

                       zP = grid(gP)%nz
                       lgReturn=.true.
!print*, 'ps15', lgreturn
                    end if

                    if ((rVec%z >=  grid(gP)%zAxis(grid(gP)%nz)+grid(gP)%geoCorrZ &
                         & .or. zP>grid(gP)%nz) .and. gP>1) then

                       xP = enPacket%xP(grid(gP)%motherP)
                       yP =  enPacket%yP(grid(gP)%motherP)
                       zP =  enPacket%zP(grid(gP)%motherP)
                       gP = grid(gP)%motherP

                    end if

                    if (lgReturn) then

                       ! the packet escapes without further interaction
!print*, 'ps16', lgreturn
                       idirT = int(acos(enPacket%direction%z)/dTheta)+1
                       if (idirT>totangleBinsTheta) then
                          idirT=totangleBinsTheta
                       end if
                       if (idirT<1 .or. idirT>totAngleBinsTheta) then
                          print*, '! fluorescencePacketRun: error in theta direction cosine assignment',&
                               &  idirT, enPacket, dTheta, totAngleBinsTheta
                          stop
                       end if


                       if (enPacket%direction%x<1.e-35) then
                          idirP = 0
                       else
                          idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
                       end if
                       if (idirP<0) idirP=totAngleBinsPhi+idirP
                       idirP=idirP+1

                       if (idirP>totangleBinsPhi) then
                          idirP=totangleBinsPhi
                       end if
                       if (idirP<1 .or. idirP>totAngleBinsPhi) then
                          print*, '! fluorescencePacketRun: error in phi direction cosine assignment',&
                               &  idirP, enPacket, dPhi, totAngleBinsPhi
                          stop
                       end if



                       if (nAngleBins>0) then
                          if ( (totangleBinsPhi>1 .and. viewPointPtheta(idirT) == viewPointPphi(idirP)) .or. &
                               & (totangleBinsPhi==1 .and. viewPointPtheta(idirT)>0) ) then
                             grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,&
                                  & viewPointPtheta(idirT)) = &
                                  &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                  & enPacket%nuP,viewPointPtheta(idirT)) +  deltaEUsed
                             grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,0) = &
                                  &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                  & enPacket%nuP,0) +  deltaEUsed
                             if (lgSeparateSED) then
                                grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), enPacket%nuP,&
                                     & viewPointPtheta(idirT), enPacket%SEDtype) = &
                                     &grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                     & enPacket%nuP,viewPointPtheta(idirT), enPacket%SEDtype) +  deltaEUsed
                                grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                     &enPacket%nuP,0, enPacket%SEDtype) = &
                                     &grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                     & enPacket%nuP,0, enPacket%SEDtype) +  deltaEUsed
                             endif
                          else
                             grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2),enPacket%nuP,0) = &
                                  & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                  & enPacket%nuP,0) +  deltaEUsed
                             if (lgSeparateSED) then
                                grid(enPacket%origin(1))%escapedPacketsComponents&
                                     &(enPacket%origin(2),enPacket%nuP,0, enPacket%SEDtype) = &
                                     & grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                     & enPacket%nuP,0, enPacket%SEDtype) +  deltaEUsed
                             end if
                          end if
                       else
                          grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                               & enPacket%nuP,0) = &
                               & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                               & enPacket%nuP,0) +  deltaEUsed
                          if (lgSeparateSED) then
                             grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                  & enPacket%nuP,0, enPacket%SEDtype) = &
                                  & grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                  & enPacket%nuP,0, enPacket%SEDtype) +  deltaEUsed
                          endif

                       end if


!print*, 'ps16 out'
                       return
                    end if

                 end if

                 ! end be 6/6/06


                 if(lgPlaneIonization) then
                     lgReturn=.false.

                     if ( rVec%y <= grid(gP)%yAxis(1) .or. yP<1) then

                        ! the energy packet escapes this grid
                        if (gP==1) then
                           yP=1
                           lgReturn=.true.
                        else if (gP>1) then
                           ! locate where we are at on the mother grid
                           call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                           if (enPacket%xP(1)<1) enPacket%xP(1)=1
                           if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                              if ( rvec%x >= ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                                   & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                                   & enPacket%xP(1) = enPacket%xP(1)+1
                           end if

                           call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                           if (enPacket%yP(1)<1) enPacket%yP(1)=1
                           if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                              if ( rvec%y >= ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                                   & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                                   & enPacket%yP(1) = enPacket%yP(1)+1
                           end if

                           call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                           if (enPacket%zP(1)<1) enPacket%zP(1)=1
                           if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                              if ( rvec%z >= ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                                   & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                                   & enPacket%zP(1) = enPacket%zP(1)+1
                           end if

                           xP = enPacket%xP(1)
                           yP = enPacket%yP(1)
                           zP = enPacket%zP(1)
                           gP = 1
                           igpp = 1
                        else
                           print*, '! fluoPathSegment: insanity occured - invalid gP', gP
                           stop
                        end if

                     end if

                     if (rVec%y > grid(gP)%yAxis(grid(gP)%ny) .or. yP>grid(gP)%ny) then

                        if (gP==1) then
                           ! the energy packet escapes
                           yP = grid(gP)%ny
                           lgReturn=.true.
                        else if (gP>1) then

                           ! locate where we are at on the mother grid
                           call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                           if (enPacket%xP(1)<1) enPacket%xP(1)=1
                           if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                              if ( rvec%x >= ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                                   & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                                   & enPacket%xP(1) = enPacket%xP(1)+1
                           end if

                           call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                           if (enPacket%yP(1)<1) enPacket%yP(1)=1
                           if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                              if ( rvec%y >= ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                                   & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                                   & enPacket%yP(1) = enPacket%yP(1)+1
                           end if

                           call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                           if (enPacket%zP(1)<1) enPacket%zP(1)=1
                           if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                              if ( rvec%z >= ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                                 & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                                 & enPacket%zP(1) = enPacket%zP(1)+1
                           end if

                           xP = enPacket%xP(1)
                           yP = enPacket%yP(1)
                           zP = enPacket%zP(1)
                           gP = 1
                           igpp = 1
                        else
                           print*, '! fluoPathSegment: insanity occured - invalid gP', gP
                           stop
                        end if

                     end if

                     if ( (rVec%x <= grid(1)%xAxis(1) .or. xP<1) ) then
                        xP=1
                        rVec%x = grid(gP)%xAxis(1)
                        vHat%x = -vHat%x

                     end if

                     if ( (rVec%x <= grid(gP)%xAxis(1) .or. xP<1) &
                          & .and. gP>1) then


                        ! locate where we are at on the mother grid
                        call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                        if (enPacket%xP(1)<1) enPacket%xP(1)=1
                        if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                           if ( rvec%x >= ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                                & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                                & enPacket%xP(1) = enPacket%xP(1)+1
                        end if

                        call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                        if (enPacket%yP(1)<1) enPacket%yP(1)=1
                        if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                           if ( rvec%y >= ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                              & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                              & enPacket%yP(1) = enPacket%yP(1)+1
                        end if

                        call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                        if (enPacket%zP(1)<1) enPacket%zP(1)=1
                        if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                           if ( rvec%z >= ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                                & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                                & enPacket%zP(1) = enPacket%zP(1)+1
                        end if

                        xP = enPacket%xP(1)
                        yP = enPacket%yP(1)
                        zP = enPacket%zP(1)
                        gP = 1
                        igpp = 1



                     end if

                     if ( (rVec%x >=  grid(1)%xAxis(grid(gP)%nx) &
                          & .or. xP>grid(gP)%nx)  )then
                        xP = grid(gP)%nx
                        rVec%x = grid(gP)%xAxis(grid(gP)%nx)
                        vHat%x = -vHat%x

                     end if

                     if ( (rVec%x >=  grid(gP)%xAxis(grid(gP)%nx)&
                          & .or. xP>grid(gP)%nx) .and.  gP>1) then

                        ! locate where we are at on the mother grid
                        call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                        if (enPacket%xP(1)<1) enPacket%xP(1)=1
                        if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                           if ( rvec%x >= ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                                & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                                & enPacket%xP(1) = enPacket%xP(1)+1
                        end if

                        call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                        if (enPacket%yP(1)<1) enPacket%yP(1)=1
                        if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                           if ( rvec%y >= ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                                & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                                & enPacket%yP(1) = enPacket%yP(1)+1
                        end if

                        call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                        if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                           if (enPacket%zP(1)<1) enPacket%zP(1)=1
                           if ( rvec%z >= ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                                & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                                & enPacket%zP(1) = enPacket%zP(1)+1
                        end if

                        xP = enPacket%xP(1)
                        yP = enPacket%yP(1)
                        zP = enPacket%zP(1)
                        gP = 1
                        igpp = 1

                     end if

                     if ( (rVec%z <= grid(1)%zAxis(1) .or.zP<1) ) then
                        zP=1
                        rVec%z = grid(gP)%yAxis(1)
                        vHat%z = -vHat%z

                     end if

                     if ( (rVec%z <= grid(gP)%zAxis(1) &
                          & .or.zP<1) .and. gP>1) then

                        ! locate where we are at on the mother grid
                        call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                        if (enPacket%xP(1)<1) enPacket%xP(1)=1
                        if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                           if ( rvec%x >= ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                                & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                                & enPacket%xP(1) = enPacket%xP(1)+1
                        end if

                        call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                        if (enPacket%yP(1)<1) enPacket%yP(1)=1
                        if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                           if ( rvec%y >= ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                                & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                                & enPacket%yP(1) = enPacket%yP(1)+1
                        end if

                        call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                        if (enPacket%zP(1)<1) enPacket%zP(1)=1
                        if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                           if ( rvec%z >= ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                                & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                                & enPacket%zP(1) = enPacket%zP(1)+1
                        end if

                        xP = enPacket%xP(1)
                        yP = enPacket%yP(1)
                        zP = enPacket%zP(1)
                        gP = 1
                        igpp = 1

                     end if

                     if ( (rVec%z >=  grid(1)%zAxis(grid(gP)%nz) .or. zP>grid(gP)%nz) &
                          & ) then

                        zP = grid(gP)%nz
                        rVec%z = grid(gP)%zAxis(grid(gP)%nz)
                        vHat%z = -vHat%z

                     end if

                     if ((rVec%z >=  grid(gP)%zAxis(grid(gP)%nz) &
                          & .or. zP>grid(gP)%nz) .and. gP>1) then

                        ! locate where we are at on the mother grid
                        call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                        if (enPacket%xP(1)<1) enPacket%xP(1)=1
                        if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                           if ( rvec%x >= ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                                & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                                & enPacket%xP(1) = enPacket%xP(1)+1
                        end if

                        call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                        if (enPacket%yP(1)<1) enPacket%yP(1)=1
                        if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                           if ( rvec%y >= ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                                & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                              & enPacket%yP(1) = enPacket%yP(1)+1
                        end if

                        call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                        if (enPacket%zP(1)<1) enPacket%zP(1)=1
                        if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                           if ( rvec%z >= ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                                & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                                & enPacket%zP(1) = enPacket%zP(1)+1
                        end if

                        xP = enPacket%xP(1)
                        yP = enPacket%yP(1)
                        zP = enPacket%zP(1)
                        gP = 1
                        igpp = 1

                     end if


                     if (lgReturn) then

                        ! the packet escapes without further interaction

                        idirT = int(acos(enPacket%direction%z)/dTheta)+1
                        if (idirT>totangleBinsTheta) then
                           idirT=totangleBinsTheta
                        end if
                        if (idirT<1 .or. idirT>totAngleBinsTheta) then
                           print*, '! energyPacketRun: error in theta direction cosine assignment',&
                                &  idirT, enPacket, dTheta, totAngleBinsTheta
                           stop
                        end if


                        if (enPacket%direction%x<1.e-35) then
                           idirP = 0
                        else
                           idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
                        end if
                        if (idirP<0) idirP=totAngleBinsPhi+idirP
                        idirP=idirP+1

                        if (idirP>totangleBinsPhi) then
                           idirP=totangleBinsPhi
                        end if
                        if (idirP<1 .or. idirP>totAngleBinsPhi) then
                           print*, '! energyPacketRun: error in phi direction cosine assignment',&
                                &  idirP, enPacket, dPhi, totAngleBinsPhi
                           stop
                        end if



                        if (nAngleBins>0) then
                           if ( (totangleBinsPhi>1 .and. viewPointPtheta(idirT) == viewPointPphi(idirP)) .or. &
                                & (totangleBinsPhi==1 .and. viewPointPtheta(idirT)>0) ) then
                              grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,&
                                   & viewPointPtheta(idirT)) = &
                                   &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                   & enPacket%nuP,viewPointPtheta(idirT)) +  deltaEUsed
                              grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,0) = &
                                   &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                   & enPacket%nuP,0) +  deltaEUsed

                              if (lgSeparateSED) then
                                 grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), enPacket%nuP,&
                                      & viewPointPtheta(idirT),enPacket%SEDtype) = &
                                      &grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                      & enPacket%nuP,viewPointPtheta(idirT),enPacket%SEDtype) +  deltaEUsed
                                 grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                      &enPacket%nuP,0,enPacket%SEDtype) = &
                                      &grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                      & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                              endif

                           else
                              grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2),enPacket%nuP,0) = &
                                   & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                   & enPacket%nuP,0) +  deltaEUsed
                              if (lgSeparateSED) then
                                 grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2),&
                                      &enPacket%nuP,0,enPacket%SEDtype) = &
                                      & grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                      & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                              endif
                           end if

                        else
                           grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                & enPacket%nuP,0) = &
                                & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                & enPacket%nuP,0) +  deltaEUsed
                           if (lgSeparateSED) then
                              grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                   & enPacket%nuP,0,enPacket%SEDtype) = &
                                   & grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                   & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                           endif

                        end if


                        return
                     end if

                  end if

                  ! check if the path is still within the simulation region
                  radius = 1.e10*sqrt((rVec%x/1.e10)*(rVec%x/1.e10) + &
                       &                     (rVec%y/1.e10)*(rVec%y/1.e10) + &
                       &                     (rVec%z/1.e10)*(rVec%z/1.e10))

                  if (.not.lgPlaneIonization) then

                     if ( (.not.lgSymmetricXYZ .and. (rVec%x<=grid(1)%xAxis(1).or.&
                          & rVec%y<=grid(1)%yAxis(1) .or. rVec%z<=grid(1)%zAxis(1)&
                          & )) .or.&
                          & (rVec%x >= grid(gP)%xAxis(grid(gP)%nx)) .or.&
                          &(rVec%y >= grid(gP)%yAxis(grid(gP)%ny)) .or.&
                          &(rVec%z >= grid(gP)%zAxis(grid(gP)%nz)) .or. &
                          & xP>grid(gP)%nx .or. yP>grid(gP)%ny .or. zP>grid(gP)%nz ) then

                        if (gP==1) then

                           if (enPacket%xP(1) > grid(1)%nx) xP = grid(1)%nx
                           if (enPacket%yP(1) > grid(1)%ny) yP = grid(1)%ny
                           if (enPacket%zP(1) > grid(1)%nz) zP = grid(1)%nz
                           if (enPacket%xP(1) < 1) xP = 1
                           if (enPacket%yP(1) < 1) yP = 1
                           if (enPacket%zP(1) < 1) zP = 1


                           ! the energy packet escapes

                           ! the packet escapes without further interaction
!print*, 'ps17'
                           idirT = int(acos(enPacket%direction%z)/dTheta)+1
                           if (idirT>totangleBinsTheta) then
                              idirT=totangleBinsTheta
                           end if
                           if (idirT<1 .or. idirT>totAngleBinsTheta) then
                              print*, '! energyPacketRun: error in theta direction cosine assignment',&
                                   &  idirT, enPacket, dTheta, totAngleBinsTheta
                              stop
                           end if

                           if (enPacket%direction%x<1.e-35) then
                              idirP = 0
                           else
                              idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
                           end if

!print*, idirp, enPacket%direction%x

                           if (idirP<0) idirP=totAngleBinsPhi+idirP
                           idirP=idirP+1
!print*, idirp                          , idirt
                           if (idirP>totangleBinsPhi) then
                              idirP=totangleBinsPhi
                           end if
!print*, idirp                           ,idirt
                           if (idirP<1 .or. idirP>totAngleBinsPhi) then
                              print*, '! energyPacketRun: error in phi direction cosine assignment',&
                                   &  idirP, enPacket, dPhi, totAngleBinsPhi
                              stop
                           end if
!print*,  viewPointPphi(idirp),       viewPointPtheta(idirt)

                           if (nAngleBins>0) then
                              if ( (totangleBinsPhi>1 .and. viewPointPtheta(idirT) == viewPointPphi(idirP)) .or. &
                                   & (totangleBinsPhi==1 .and. viewPointPtheta(idirT)>0) ) then
                                 grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,&
                                      & viewPointPtheta(idirT)) = &
                                      &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                      & enPacket%nuP,viewPointPtheta(idirT)) +  deltaEUsed
                                 grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,0) = &
                                      &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                      & enPacket%nuP,0) +  deltaEUsed
                                 if (lgSeparateSED) then
                                    grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), enPacket%nuP,&
                                         & viewPointPtheta(idirT), enPacket%SEDtype) = &
                                         &grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                         & enPacket%nuP,viewPointPtheta(idirT), enPacket%SEDtype) +  deltaEUsed
                                    grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                         &enPacket%nuP,0, enPacket%SEDtype) = &
                                         &grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                         & enPacket%nuP,0, enPacket%SEDtype) +  deltaEUsed
                                 endif
                              else
                                 grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2),enPacket%nuP,0) = &
                                      & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                      & enPacket%nuP,0) +  deltaEUsed
                                 if (lgSeparateSED) then
                                    grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2),&
                                         &enPacket%nuP,0,enPacket%SEDtype) = &
                                         & grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                         & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                                 endif
                              end if
                           else
                              grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                   & enPacket%nuP,0) = &
                                   & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                   & enPacket%nuP,0) +  deltaEUsed
                              if (lgSeparateSED) then
                                 grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                      & enPacket%nuP,0,enPacket%SEDtype) = &
                                      & grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                      & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                              endif

                           end if

!print*, 'ps17 out'
                           !b2.005
                           return

                        else if (gP>1) then

                           ! locate where we are at on the mother grid
                           call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                           if (enPacket%xP(1)<1) enPacket%xP(1)=1
                           if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                              if ( rvec%x >= ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                                   & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                                   & enPacket%xP(1) = enPacket%xP(1)+1
                           end if

                           call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                           if (enPacket%yP(1)<1) enPacket%yP(1)=1
                           if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                              if ( rvec%y >= ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                                   & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                                   & enPacket%yP(1) = enPacket%yP(1)+1
                           end if

                           call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                           if (enPacket%zP(1)<1) enPacket%zP(1)=1
                           if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                              if ( rvec%z >= ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                                   & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                                   & enPacket%zP(1) = enPacket%zP(1)+1
                           end if

                           xP = enPacket%xP(1)
                           yP = enPacket%yP(1)
                           zP = enPacket%zP(1)
                           gP = 1
                           igpp = 1


                           if (gP/=1) then
                              print*, '! fluoPathSegment: nested multigrids still not implemented'
                              stop
                           end if

                           if ( (radius >= R_out .and. R_out >= 0.) .or.&
                                & (rVec%x >= grid(1)%xAxis(grid(1)%nx)) .or.&
                                &(rVec%y >= grid(1)%yAxis(grid(1)%ny) ) .or.&
                                &(rVec%z >= grid(1)%zAxis(grid(1)%nz) ) .or. &
                                & (.not.lgSymmetricXYZ .and.  (rVec%x<=grid(1)%xAxis(1) .or.&
                                & rVec%y<=grid(1)%yAxis(1) .or. rVec%z<=grid(1)%zAxis(1)&
                                & ))) then

!print*, 'ps18'
                              if (xP > grid(gP)%nx) xP = grid(gP)%nx
                              if (yP > grid(gP)%ny) yP = grid(gP)%ny
                              if (zP > grid(gP)%nz) zP = grid(gP)%nz

                              if (xP < 1) xP = 1
                              if (yP < 1) yP = 1
                              if (zP < 1) zP = 1

                              ! the energy packet escapes

                              ! the packet escapes without further interaction

                              idirT = int(acos(enPacket%direction%z)/dTheta)+1
                              if (idirT>totangleBinsTheta) then
                                 idirT=totangleBinsTheta
                              end if
                              if (idirT<1 .or. idirT>totAngleBinsTheta) then
                                 print*, '! energyPacketRun: error in theta direction cosine assignment',&
                                      &  idirT, enPacket, dTheta, totAngleBinsTheta
                                 stop
                              end if

                              if (enPacket%direction%x<1.e-35) then
                                 idirP = 0
                              else
                                 idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
                              end if
                              if (idirP<0) idirP=totAngleBinsPhi+idirP
                              idirP=idirP+1

                              if (idirP>totangleBinsPhi) then
                                 idirP=totangleBinsPhi
                              end if

                              if (idirP<1 .or. idirP>totAngleBinsPhi) then
                                 print*, '! energyPacketRun: error in phi direction cosine assignment',&
                                      &  idirP, enPacket, dPhi, totAngleBinsPhi
                                 stop
                              end if

                              if (nAngleBins>0) then
                                 if ( (totangleBinsPhi>1 .and. viewPointPtheta(idirT) == viewPointPphi(idirP)) .or. &
                                      & (totangleBinsPhi==1 .and. viewPointPtheta(idirT)>0) ) then
                                    grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,&
                                         & viewPointPtheta(idirT)) = &
                                         &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                         & enPacket%nuP,viewPointPtheta(idirT)) +  deltaEUsed
                                    grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,0) = &
                                         &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                         & enPacket%nuP,0) +  deltaEUsed
                                    if (lgSeparateSED) then
                                       grid(enPacket%origin(1))%escapedPacketsComponents&
                                            &(enPacket%origin(2), enPacket%nuP,&
                                            & viewPointPtheta(idirT),enPacket%SEDtype) = &
                                            &grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                            & enPacket%nuP,viewPointPtheta(idirT),enPacket%SEDtype) +  deltaEUsed
                                       grid(enPacket%origin(1))%escapedPacketsComponents&
                                            &(enPacket%origin(2), enPacket%nuP,0,enPacket%SEDtype) = &
                                            &grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                            & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                                    endif
                                 else
                                    grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2),enPacket%nuP,0) = &
                                         & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                         & enPacket%nuP,0) +  deltaEUsed
                                    if (lgSeparateSED) then
                                       grid(enPacket%origin(1))%escapedPacketsComponents&
                                            &(enPacket%origin(2),enPacket%nuP,0, enPacket%SEDtype) = &
                                            & grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                            & enPacket%nuP,0, enPacket%SEDtype) +  deltaEUsed
                                    endif
                                 end if
                              else
                                 grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                      enPacket%nuP,0) = &
                                      & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                      & enPacket%nuP,0) +  deltaEUsed
                                 if (lgSeparateSED) then
                                    grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                         enPacket%nuP,0,enPacket%SEDtype) = &
                                         & grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                                         & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                                 endif


                              end if


!print*, 'ps 18 out'
                              !b2.005
                              return

                           end if
                        else
                           print*, '! fluoPathSegment: insanity occured - invalid gP - ', gP
                           stop
                        end if

                     end if

!            if (lgSymmetricXYZ .and. gP == 1) then
                      if (lgSymmetricXYZ ) then
                         if (lgPlaneIonization) then
                            print*, '! fluoPathSegment: lgSymmetric and lgPlaneionization flags both raised'
                            stop
                         end if

                         if ( rVec%x <= grid(1)%xAxis(1) .or. (gP==1 .and. xP<1)) then
                            if (vHat%x<0.) vHat%x = -vHat%x
                            enPacket%xP(1) = 1
                            xP = 1
                            rVec%x = grid(gP)%xAxis(1)
                         end if
                         if ( rVec%y <= grid(1)%yAxis(1) .or. (gP==1 .and. yP<1)) then
                            if (vHat%y<0.) vHat%y = -vHat%y
                            enPacket%yP(1)=1
                            yP = 1
                            rVec%y = grid(gP)%yAxis(1)
                         end if
                         if ( rVec%z <= grid(1)%zAxis(1) .or. (gP==1 .and. zP<1)) then
                            if (vHat%z<0.) vHat%z = -vHat%z
                            enPacket%zP(1) = 1
                            zP=1
                            rVec%z = grid(1)%zAxis(1)
                         end if

                      end if

                   end if

                   if (gP>1) then
                      if ( ( (rVec%x <= grid(gP)%xAxis(1) &
                           &.or. xP<1) .and. vHat%x <=0.) .or. &
                           & ( (rVec%y <= grid(gP)%yAxis(1) &
                           & .or. yP<1) .and. vHat%y <=0.) .or. &
                           & ( (rVec%z <= grid(gP)%zAxis(1) &
                           &  .or. zP<1) .and. vHat%z <=0.) .or. &
                           & ( (rVec%x >= grid(gP)%xAxis(grid(gP)%nx) &
                           &.or. xP>grid(gP)%nx) .and. vHat%x >=0.) .or. &
                           & ( (rVec%y >= grid(gP)%yAxis(grid(gP)%ny) &
                           & .or. yP>grid(gP)%ny) .and. vHat%y >=0.) .or. &
                           & ( (rVec%z >= grid(gP)%zAxis(grid(gP)%nz) &
                           &  .or. zP>grid(gP)%nz) .and. vHat%z >=0.) ) then

                         ! locate where we are at on the mother grid
                         call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                         if (enPacket%xP(1)<1) enPacket%xP(1)=1
                         if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                            if ( rvec%x >= ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                                 & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                                 & enPacket%xP(1) = enPacket%xP(1)+1
                         end if

                         call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                         if (enPacket%yP(1)<1) enPacket%yP(1)=1
                         if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                            if ( rvec%y >= ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                                 & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                              & enPacket%yP(1) = enPacket%yP(1)+1
                         end if

                         call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                         if (enPacket%zP(1)<1) enPacket%zP(1)=1
                         if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                            if ( rvec%z >= ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                                 & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                                 & enPacket%zP(1) = enPacket%zP(1)+1
                         end if

                         gP = 1
                         igpp = 1


                         xP = enPacket%xP(1)
                         yP = enPacket%yP(1)
                         zP = enPacket%zP(1)

                      end if

                   end if
!print*, 'ps19'
               end if

               if (.not. lgPlaneIonization .and. gP==1 .and. (xP > grid(gP)%nx  &
                    & .or. yP > grid(gP)%ny .or. zP > grid(gP)%nz) ) then
!print*, 'ps20'
           ! the energy packet escapes

           ! the packet escapes without further interaction

           idirT = int(acos(enPacket%direction%z)/dTheta)+1
           if (idirT>totangleBinsTheta) then
              idirT=totangleBinsTheta
           end if
           if (idirT<1 .or. idirT>totAngleBinsTheta) then
              print*, '! energyPacketRun: error in theta direction cosine assignment',&
                   &  idirT, enPacket, dTheta, totAngleBinsTheta
              stop
           end if

           if (enPacket%direction%x<1.e-35) then
              idirP = 0
           else
              idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
           end if
           if (idirP<0) idirP=totAngleBinsPhi+idirP
           idirP=idirP+1

           if (idirP>totangleBinsPhi) then
              idirP=totangleBinsPhi
           end if

           if (idirP<1 .or. idirP>totAngleBinsPhi) then
              print*, '! energyPacketRun: error in phi direction cosine assignment',&
                   &  idirP, enPacket, dPhi, totAngleBinsPhi
              stop
           end if


           if (nAngleBins>0) then
              if ( (totangleBinsPhi>1 .and. viewPointPtheta(idirT) == viewPointPphi(idirP)) .or. &
                   & (totangleBinsPhi==1 .and. viewPointPtheta(idirT)>0) ) then
                 grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,&
                      & viewPointPtheta(idirT)) = &
                      &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                      & enPacket%nuP,viewPointPtheta(idirT)) +  deltaEUsed
                 grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,0) = &
                      &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                      & enPacket%nuP,0) +  deltaEUsed
                 if (lgSeparateSED) then
                    grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), enPacket%nuP,&
                         & viewPointPtheta(idirT),enPacket%SEDtype) = &
                         &grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                         & enPacket%nuP,viewPointPtheta(idirT),enPacket%SEDtype) +  deltaEUsed
                    grid(enPacket%origin(1))%escapedPacketsComponents&
                         &(enPacket%origin(2), enPacket%nuP,0,enPacket%SEDtype) = &
                         &grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                         & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                 endif
              else
                 grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2),enPacket%nuP,0) = &
                      & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                      & enPacket%nuP,0) +  deltaEUsed
                 if (lgSeparateSED) then
                    grid(enPacket%origin(1))%escapedPacketsComponents&
                         &(enPacket%origin(2),enPacket%nuP,0,enPacket%SEDtype) = &
                         & grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                         & enPacket%nuP,0,enPacket%SEDtype) +  deltaEUsed
                 endif
              end if
           else
              grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                   enPacket%nuP,0) = &
                   & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                   & enPacket%nuP,0) +  deltaEUsed
              if (lgSeparateSED) then
                 grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                      & enPacket%nuP,0, enPacket%SEDtype) = &
                      & grid(enPacket%origin(1))%escapedPacketsComponents(enPacket%origin(2), &
                      & enPacket%nuP,0, enPacket%SEDtype) +  deltaEUsed
              endif
           end if

!print*, 'ps20 out'
           !b2.005
            return


         end if

      end do ! safelimit loop

      if (i>= safeLimit) then
         if (.not.lgPlaneIonization) then
            print*, '! fluoPathSegment: [warning] packet trajectory has exceeded&
                 &  maximum number of events', safeLimit, gP, xP,yP,zP, grid(gP)%active(xP,yP,zP), &
              & rvec, vhat, enPacket
         end if
         return

      end if

      if (gP==1) then
         igpp = 1
      else if (gP>1) then
         igpp = 2
      else
         print*, "! fluoPathSegment: insane grid index "
         stop
      end if

      enPacket%xP(igpp) = xP
      enPacket%yP(igpp) = yP
      enPacket%zP(igpp) = zP


      if (deltaEUsed/(Lstar(istar)*1.e36/nPhotons(istar))<1.e-10) return
      packetType = 'diffuse'

      chTypeIn = packetType
      positionIn = rVec
      inX =  enPacket%xP(1:2)
      inY =  enPacket%yP(1:2)
      inZ =  enPacket%zP(1:2)
      gPIn = gP
      reRun = 1

      ! the energy packet has beenid absorbed - reemit a new packet from this position
!      call energyPacketRun(packetType, rVec, enPacket%xP(1:2), enPacket%yP(1:2), &
!                & enPacket%zP(1:2), gP)
      return

 end subroutine fluoPathSegment

 ! determine energy and direction for compton redistribution of fluorescent lines
 subroutine comptonScatter(fPacket)
   implicit none

   type(photon_packet),intent(inout) :: fPacket            ! the fluorescent packet

   real                              :: newNu,newTheta     ! new energy and theta
   real                              :: random             ! random number
   real                              :: u,v,w,t            ! direction units

   integer                           :: newNuP,newThetaP   ! pointer to new energy and theta


!print*, 'compton :'
   ! calculate new direction
   call getNu2(KNsigmaArray(fPacket%nuP,1:180),newThetaP)
   newTheta = real(newThetaP)*Pi/180.

!print*, newTheta
   ! calculate u,v,w (Harries & Howarth, 1997)
   call random_number(random)
   w = 2.*random - 1.
   t = sqrt(1-w*w)
   u = t*cos(newTheta)
   v = t*sin(newTheta)

   fPacket%direction%x = u
   fPacket%direction%y = v
   fPacket%direction%z = w

   ! calculate new energy
   newNu = fPacket%nu*PcompArray(fPacket%nuP,newThetaP)

   ! find this in the nuArray
   call locate(nuArray,newNu,newNuP)

!print*, fPacket%Nu, fPacket%NuP
   fPacket%Nu  = newNu
   fPacket%NuP = newNuP
!print*, fPacket%Nu, fPacket%NuP
 end subroutine comptonScatter

      end subroutine energyPacketDriver

end module photon_mod
