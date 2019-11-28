module readdata_mod
use common_mod

!subroutines to read in data files to arrays at beginning of simulation

contains

    subroutine readdata()

        implicit none
        integer :: i, j,k, ios               ! counter, io state
        integer :: iup, ilow                 ! H counters
        real    :: reader(8)

    ! initialise all the arrays that we are reading in

        y_dat = 0.
        Ay_dat = 0.
        HIRecLineData = 0.
        HeIILymanData = 0.
        HeIILymanNuData = 0.
        HeIIRecLinedata = 0.
        direc_coeffs%elem = 0
        direc_coeffs%n = 0
        direc_coeffs%a = 0.
        direc_coeffs%b = 0.
        direc_coeffs%c = 0.
        direc_coeffs%d = 0.
        direc_coeffs%f = 0.
        direc_coeffs%g = 0
        aldropequi_coeffs%elem = 0
        aldropequi_coeffs%n = 0
        aldropequi_coeffs%a = 0
        aldropequi_coeffs%b = 0
        aldropequi_coeffs%t0 = 0.
        aldropequi_coeffs%t1 = 0.

    ! read in rates from data/HeI2phot.dat

        open(unit = 93,  action="read", file = PREFIX//"/share/mocassinX/data/HeI2phot.dat", status = "old", position = "rewind", iostat=ios)
        if (ios /= 0) then
            print*, "! readData: can't open file: ",PREFIX,"/share/mocassinX/data/HeI2phot.dat"
            stop
        end if
        do i = 1, 41
            read(unit = 93, fmt = *) y_dat(i), Ay_dat(i)
        end do

        close(93)

        ! dielectronic recombination coefficients

        open (unit=18, file=PREFIX//'/share/mocassinX/data/dielectronic.dat', status='old',position='rewind', iostat = ios, action="read")
        do i = 1, 25
           read(unit=18, fmt=*, iostat=ios) direc_coeffs(i)%elem, direc_coeffs(i)%n, direc_coeffs(i)%a, direc_coeffs(i)%b, direc_coeffs(i)%c, direc_coeffs(i)%d, direc_coeffs(i)%f, direc_coeffs(i)%g
           if (ios < 0) exit ! end of file reached
        enddo
        close(18)

        ! high temperature dielectronic recombination coefficients from
        ! Aldrovandi and Pequignot 1973

        open (unit=17, file=PREFIX//'/share/mocassinX/data/aldrovandi.dat', status='old',position='rewind', iostat = ios, action="read")
        if (ios /= 0) then
           print*, "! readData: can't open file ",PREFIX,"/share/mocassinX/data/alrovandi.dat"
           stop
        end if

        do i = 1, 167
            read(unit=17, fmt=*, iostat=ios) aldropequi_coeffs(i)%elem, aldropequi_coeffs(i)%n, aldropequi_coeffs(i)%a, aldropequi_coeffs(i)%b, aldropequi_coeffs(i)%t0, aldropequi_coeffs(i)%t1
            if (ios < 0) exit ! end of file reached
        end do

        close(17)

        !hydrogenic

        open(file=PREFIX//"/share/mocassinX/data/hydroLinesFiles.dat", unit=19)
        if (ios /= 0) then
           print*, "! readData: can't open file ",PREFIX,"/share/mocassinX/data/hydroLinesFiles.dat"
           stop
        end if

        do i = 1, 9
           do j = 1, 12
              read(19,*) hydroLinesFile(i,j)

              open(unit = 94,  action="read", file = PREFIX//"/share/mocassinX/"//hydroLinesFile(i,j), status = "old", position = "rewind", iostat=ios)
              if (ios /= 0) then
                 print*, "! RecLinesEmission: can't open file: ",PREFIX,"/share/mocassinX/",hydroLinesFile(i,j)
                 stop
              end if

              do k = 1, 8! no more than 8 density points in the data
                 read(unit=94, fmt=*, iostat=ios) hydroLinesData(i,j,k)%dens
                 if (ios<0) exit

                 do iup = 15, 2, -1
                    read(94, fmt=*) (hydroLinesData(i,j,k)%linedata(iup, ilow), ilow = 1, min(8, iup-1))
                 end do
              end do
              close(94)
           end do
        end do

        close(19)

        open(unit = 98,  action="read", file = PREFIX//"/share/mocassinX/data/r2a0100old.dat", status = "old", position = "rewind", iostat=ios)
        if (ios /= 0) then
            print*, "! setDiffusePDF: can't open file: ",PREFIX,"/share/mocassinX/data/r2a0100old.dat"
            stop
        end if
        do i = 1, NHeIILyman
            read(98, fmt=*) HeIILymandata(i), HeIILymanNudata(i)
        end do

        close(98)

        close(19)
        open(file=PREFIX//"/share/mocassinX/data/hydroLinesFiles.dat", unit=19)
        do i = 1, 9
           do j = 1, 12
              read(19,*) hydroLinesFile(i,j)
           end do
        end do
        close(19)

        close(21)
        open(file=PREFIX//"/share/mocassinX/data/gaussfitHb.dat", unit=21)
        do i = 1, 9
           read(21,*) (reader(j), j = 1, 8)
           HbACoeff(i,1:4) = reader(1:4)
           HbBCoeff(i,1:4) = reader(5:8)
        end do
        close(21)

        hydroLinesTemps = (/0.05, 0.10, 0.30, 0.50, 0.75, 1.0, &
             &1.25, 1.50, 2.0, 3.0, 5.0, 10.0/)

        close(22)
        open(file=PREFIX//"/share/mocassinX/data/r1bEdge.dat", unit=22)
        do i = 1, 26
           read(22,*) rbEdge(1,1,i), rbEdge(1,2,i), rbEdge(1,3,i)
        end do
        close(22)

        close(23)
        open(file=PREFIX//"/share/mocassinX/data/r2bEdge.dat", unit=23)
        do i = 1, 26
           read(23,*) rbEdge(2,1,i), rbEdge(2,2,i), rbEdge(2,3,i)
        end do
        close(23)

        close(23)
        open(file=PREFIX//"/share/mocassinX/data/r2aEdge.dat", unit=23)
        do i = 1, 18
           read(23,*) r2aEdge(1,i), r2aEdge(2,i), r2aEdge(3,i)
        end do
        close(23)

        elemLabel = (/'*H','He','Li','Be','*B','*C','*N','*O','*F','Ne','Na','Mg',&
             &'Al','Si','*P','*S','Cl','Ar','*K','Ca','Sc','Ti','*V','Cr','Mn','Fe',&
             &'Co','Ni','Cu','Zn'/)

        ! set up atomic weight array
        aWeight = (/1.0080, 4.0026, 6.941, 9.0122, 10.811, 12.0111, 14.0067, 15.9994, &
            & 18.9984, 20.179, 22.9898, 24.305, 26.9815, 28.086, 30.9738, 32.06, 35.453, &
            & 39.948, 39.102, 40.08, 44.956, 47.90, 50.9414, 51.996, 54.9380, 55.847, 58.9332, &
            & 58.71, 63.546, 65.37 /)

        close(72)
        open (unit= 72,  action="read", file=PREFIX//"/share/mocassinX/dustData/nuDustRyd.dat", status = "old", position = "rewind", &
             & iostat = ios)
        if (ios /= 0) then
           print*, "! initCartesianGrid: can't open dust nu grid file - ",PREFIX,"/share/mocassinX/dustData/nuDustRyd.dat"
           stop
        end if
        close(72)

    end subroutine readdata

end module readdata_mod
