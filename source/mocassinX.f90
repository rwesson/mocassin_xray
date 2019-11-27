! Copyright (C) 2013 Barbara Ercolano
! Department of Physics and Astronomy
! Universitäts Sternwarte München
! Scheinerstr 1
! München, 81679
! Deutschland
!
! ercolano@usm.lmu.de
!
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version. This requires
! that any changes or improvements made to the program should also be
! made freely available.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

!
! MoCaSSiN = MOnte CArlo SImulationS of Nebulae
! Version 3.00
! this is the main driver of the simulation
!
program MoCaSSiN
    use common_mod
    use constants_mod
    use dust_mod
    use voronoi_grid_mod
    use iteration_mod
    use output_mod
    use set_input_mod
    use xSec_mod
    use readdata_mod

    implicit none

    include 'mpif.h'

    type(grid_type) :: grid3D(maxGrids)       ! the 3D Cartesian  grid

    real               :: test                ! test
    integer            :: i, iGrid            ! allocation error status
    character(len=10)  :: time                ! time in text format
    real, dimension(2) :: timing              ! cputimer
    integer            :: nhours, nminutes, nseconds

    call cpu_time(timing(1))
    call date_and_time(TIME=time)
    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, taskid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, numtasks, ierr)

    if (taskid == 0) then
        print*, "MOCASSIN 2007 Version ",VERSION
        print*, "compiled with ",COMPILER
        print *,"data directory: ",PREFIX,"/share/mocassinX"
        if (CO.ne."co") print *,"CO=",CO
        print *,"started running at ",time(1:2),":",time(3:4),":",time(5:6)
        print*, " "
    endif

    ! read the input parameters of the simulation
    call readInput()

    if (taskid == 0) then
        print*, " "
    end if

    if (.not.lgVoronoi) then
      ! initialize the 3D cartesian grid
      do iGrid = 1, nGrids
         call initCartesianGrid(grid3D(iGrid), nxIn(iGrid), nyIn(iGrid), nzIn(iGrid))
      end do
    else
      grid3D(1)%nCellsV = nCellsIn
      call initVoronoiGrid(Grid3D(1))  ! to be included when these routines are ready
    end if

    ! initialize opacities x sections array
    call initXSecArray()

    ! set the ionzing continuum according to the contShape variable
    call setContinuum()

    if (.not. lgVoronoi) then
      call fillGrid(grid3D(1:nGrids))
      call setStarPosition(grid3D(1)%xAxis,grid3D(1)%yAxis,grid3D(1)%zAxis, grid3D(1:nGrids))
    else
      call fillGridV(grid3D(1:nGrids))
      call setStarPositionV(grid3D(1))
    end if

    if (taskid==0) then
      do iGrid = 1, nGrids
        print*, 'Grid : ', iGrid
        print*, 'active cells: ', grid3D(iGrid)%nCells
      end do
    endif

    ! prepare atomica data stuff
    if (lgGas) then
      call readdata()
      call makeElements()
    endif

    if (taskid==0) print*, 'active elements: ', nElementsUsed

    ! if grains are included, calculate the dust opacity
    if (lgDust) then
       print*, '! mocassin: calling dustDriver'
       do iGrid = 1, nGrids
          call dustDriver(grid3D(iGrid))
       end do
       print*, '! mocassin: dustDriver done'
    end if

    call mpi_barrier(mpi_comm_world, ierr)

    ! set local diffuse ionisation field
    if (Ldiffuse>0.) then
       if (lgVoronoi) then
          print*,"! mocassin: diffuse ionisation is not yet implemented with Voronoi grids - contact B Ercolano ercolano@usm.lmu.de"
          stop
       end if

       if (taskid==0) print*, '! mocassin: calling setLdiffuse'
       call setLdiffuse(grid3D(1:nGrids))
       test=0.
       do iGrid = 1, nGrids
          do i = 1, grid3D(igrid)%ncells
             test = test+grid3d(igrid)%LdiffuseLoc(i)
          end do
       end do
       if (taskid==0) then
          print*, '! mocassin: setLdiffuse done, total Ldiffuse: ', test
       end if
    end if

    if (taskid==0) print*, '! mocassin: calling MCIterationDriver'
    ! start the Monte Carlo simulation
    call MCIterationDriver(grid3D(1:nGrids))
    if (taskid==0) print*, '! mocassin: MCIterationDriver done'

    if (taskid ==  0) then
        ! determine final statistics
        if (lgGas) call outputGas(grid3D(1:nGrids))
    end if

    call mpi_barrier(mpi_comm_world, ierr)

    ! free all space allocated to the 3D grid
    if (taskid==0) then
      do iGrid=1, nGrids
         call freeGrid(grid3D(iGrid))
      end do
    endif

    call mpi_finalize(ierr)
    call cpu_time(timing(2))

    timing(1)=timing(1)/60.
    timing(2)=timing(2)/60.
    nhours = int((timing(2)-timing(1))/60.)
    nminutes = int(mod(timing(2)-timing(1),60.))
    nseconds = nint(60.*((timing(2)-timing(1))-real(nhours*60)-real(nminutes)))

    if (taskid==0) then
      print "(A,I3.2,A,I2.2,A,I2.2,A)","total run time per processor ",nhours,":",nminutes,":",nseconds," (HMS)"
      print *,'! MoCaSSin: end simulation reached - clean exit -'
    endif

end program MoCaSSiN
