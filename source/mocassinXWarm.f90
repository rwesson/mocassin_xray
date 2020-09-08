! Copyright (C) 2007 Barbara Ercolano
!
! MoCaSSiNwarm = MOnte CArlo SImulationS of Nebulae
! this is the warm-start driver of the simulation
! (requires grid1.out,grid2.out,grid3.out files)
!
! Version 3.00
program MoCaSSiNwarm
    use common_mod
    use constants_mod
    use dust_mod
    use grid_mod
    use voronoi_grid_mod
    use iteration_mod
    use output_mod
    use set_input_mod
    use xSec_mod
    use readdata_mod

    implicit none

    include 'mpif.h'

    type(grid_type) :: grid3D(maxGrids) ! the 3D Cartesian  grid

    integer         :: err              ! allocation error status
    integer         :: iGrid            !
    character(len=10)  :: time                ! time in text format
    real, dimension(2) :: timing              ! cputimer

    call cpu_time(timing(1))
    call date_and_time(TIME=time)

    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, taskid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, numtasks, ierr)

    lgWarm = .true.

    if (taskid == 0) then
        print*, "mocassinXWarm: version ",VERSION
        print*, "compiled with ",COMPILER
        print *,"data directory: ",PREFIX,"/share/mocassinX"
        if (CO.ne."co") print *,"CO=",CO
        print *,"started running at ",time(1:2),":",time(3:4),":",time(5:6)
        print*, " "
    endif

    if (taskid == 0) then
        print*, " "
    end if

    ! check if voronoi
    close(77)
    open(unit=77, file='output/grid3.out', action="read",position='rewind',  &
         &          status='old', iostat = err)
    if (err /= 0) then
       print*, "! mocassinWarm: error opening file grid3.out"
       stop
    end if

    read(77,*) lgVoronoi

    close(77)

    if (.not. lgVoronoi) then
      ! reset the 3D cartesian grid
      call resetGrid(grid3D)
      call setStarPosition(grid3D(1)%xAxis,grid3D(1)%yAxis,grid3D(1)%zAxis,grid3D)
    else
      call resetGridV(grid3D(1))
      call setStarPositionV(grid3D(1))
    end if

    ! initialize opacities x sections array
    call initXSecArray()

    ! set the ionzing continuum according to the contShape variable
    call setContinuum()

    ! prepare atomica data stuff
    call readdata()
    call makeElements()

    if (.not.lgVoronoi) then
      if (taskid == 0) then
        print*, "Total number of grids (mother+sub-grids): ", nGrids
        print*, "Mother grid used:"
        print*, grid3D(1)%xAxis
        print*, grid3D(1)%yAxis
        print*, grid3D(1)%zAxis
      end if
    end if

    ! if grains are included, calculate the dust opacity
    if (lgDust) then
       do iGrid = 1, nGrids
          call dustDriver(grid3D(iGrid))
       end do
    end if

    call mpi_barrier(mpi_comm_world, ierr)

    ! start the Monte Carlo simulation
    call MCIterationDriver(grid3D)

    if (taskid ==  0 .and. .not.lgOutput) then
        ! determine final statistics
        if (lgGas) call outputGas(grid3D)
        call writeSED(grid3D)
        if (contCube(1)>0. .and. contCube(2)>0. ) &
             & call writeContCube(grid3D, contCube(1),contCube(2))
    end if

    call mpi_barrier(mpi_comm_world, ierr)

    ! free all space allocated to the 3D grid
    do iGrid=1, nGrids
       call freeGrid(grid3D(iGrid))
    end do

!    endtime  = mpi_wtime()
!
!    if (taskid == 0) then
!        print*, "time: ", endtime-starttime
!    end if

    call mpi_finalize(ierr)
    stop 'mpi done'


end program MoCaSSiNwarm
