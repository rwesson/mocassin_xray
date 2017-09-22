module readdata_mod
use common_mod

!subroutines to read in data files to arrays at beginning of simulation

contains

    subroutine readdata()

        implicit none
        integer :: i, j,k, ios               ! counter, io state
        integer :: iup, ilow                 ! H counters

    ! initialise all the arrays that we are reading in

        y_dat = 0.
        Ay_dat = 0.
        HIRecLineData = 0.
        HeIILymanData = 0.
        HeIILymanNu = 0.
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

        open(unit = 93,  action="read", file = PREFIX//"/share/mocassin/data/HeI2phot.dat", status = "old", position = "rewind", iostat=ios)
        if (ios /= 0) then
            print*, "! readData: can't open file: ",PREFIX,"/share/mocassin/data/HeI2phot.dat"
            stop
        end if
        do i = 1, 41
            read(unit = 93, fmt = *) y_dat(i), Ay_dat(i)
        end do

        close(93)

        ! dielectronic recombination coefficients

        open (unit=18, file=PREFIX//'/share/mocassin/data/dielectronic.dat', status='old',position='rewind', iostat = ios, action="read")
        do i = 1, 25
           read(unit=18, fmt=*, iostat=ios) direc_coeffs(i)%elem, direc_coeffs(i)%n, direc_coeffs(i)%a, direc_coeffs(i)%b, direc_coeffs(i)%c, direc_coeffs(i)%d, direc_coeffs(i)%f, direc_coeffs(i)%g
           if (ios < 0) exit ! end of file reached
        enddo
        close(18)

        ! high temperature dielectronic recombination coefficients from
        ! Aldrovandi and Pequignot 1973

        open (unit=17, file=PREFIX//'/share/mocassin/data/aldrovandi.dat', status='old',position='rewind', iostat = ios, action="read")
        if (ios /= 0) then
           print*, "! readData: can't open file ",PREFIX,"/share/mocassin/data/alrovandi.dat"
           stop
        end if

        do i = 1, 167
            read(unit=17, fmt=*, iostat=ios) aldropequi_coeffs(i)%elem, aldropequi_coeffs(i)%n, aldropequi_coeffs(i)%a, aldropequi_coeffs(i)%b, aldropequi_coeffs(i)%t0, aldropequi_coeffs(i)%t1
            if (ios < 0) exit ! end of file reached
        end do

        close(17)

        !hydrogenic

        open(file=PREFIX//"/share/mocassin/data/hydroLinesFiles.dat", unit=19)
        if (ios /= 0) then
           print*, "! readData: can't open file ",PREFIX,"/share/mocassin/data/ydroLinesFiles.dat"
           stop
        end if

        do i = 1, 9
           do j = 1, 12
              read(19,*) hydroLinesFile(i,j)

              open(unit = 94,  action="read", file = PREFIX//"/share/mocassin/"//hydroLinesFile(i,j), status = "old", position = "rewind", iostat=ios)
              if (ios /= 0) then
                 print*, "! RecLinesEmission: can't open file: ",PREFIX,"/share/mocassin/",hydroLinesFile(i,j)
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

    end subroutine readdata

end module readdata_mod
