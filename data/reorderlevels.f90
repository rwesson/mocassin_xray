program reorderLevels
  implicit none

  double precision, dimension(150,150) :: fLineEm     ! forbiddn line emissivity
  
  double precision, dimension(150,150) :: wav         ! wavelength of transition [A]
    

    double precision     :: ax, ex                    ! readers
    double precision     :: constant                  ! calculations constant 
    double precision     :: delTeK                    ! Boltzmann exponent
    double precision     :: expFac                    ! calculations factor
    double precision     :: Eji                       ! energy between levels j and i
    double precision,pointer :: qxar(:,:,:)                        ! reader

    double precision       :: qx                        ! reader
    double precision   :: sumN                      ! normalization factor for populations
    double precision        :: sqrTe                     ! sqrt(Te)
    double precision    :: value                     ! general calculations value
    
    double precision, pointer          :: a(:,:)      ! transition rates array
    double precision, pointer          :: cs(:,:)     ! collisional strengths array
    double precision, pointer          :: e(:)        ! energy levels array
    double precision, pointer          :: qeff(:,:)   ! q eff  array
    double precision, pointer          :: qom(:,:,:)  ! qom array
    double precision, pointer          :: tnij(:,:)   ! tnij array
    double precision, pointer          :: x(:,:)      ! matrix arrays
    double precision, pointer          :: y(:)        !       
    double precision, &
         & pointer :: n(:), n2(:) ! level population arrays

    real                 :: a_r(4),a_d(5),z,br       !
    real,    pointer     :: alphaTotal(:)             ! maximum 100-level ion
    real                 :: a_fit, b_fit              !         
    real                 :: qomInt                    ! interpolated value of qom
    real,    pointer     :: logTemp(:)                ! log10 temperature points array
    real,    pointer     :: qq(:)                     ! qq array
    real,    pointer     :: qq2(:)                    ! 2nd deriv qq2 array

    real, intent(in) &
         & :: Te, &    ! electron temperature [K]
         & Ne, &       ! electron density [cm^-3]
         & ionDenUp    ! ion density of the upper ion stage 

    integer :: il,iu,itt
    
    integer  :: gx               ! stat weight reader  
    integer  :: i, j, k, l, iT   ! counters/indeces
    integer  :: iRats            ! coll strength (iRats=0) or (coll rates)/10**iRats
    integer  :: ios              ! I/O error status
    integer  :: nLev             ! number of levels in atomic data file
    integer  :: numLines         ! number of lines for atomic data refernce
    
    integer  :: nTemp            ! number of temperature points in atomic data file        
    integer  :: err              ! allocation error status
    
    integer, parameter :: safeLim = 10000 ! loop safety limit
    
    integer, pointer :: g(:)              ! statistical weight array
    
    integer, dimension(2) :: ilow, &      ! lower index
         & iup          ! upper index

    character(len = 20), pointer :: &
         & label(:)! labels array

    character(len = 150) :: text  ! lines of text

    character(len = *), &
         & intent(in)  :: file_name   ! ionic data file name


    print*, 'Enter file names -old and new-'
    read*, file_name, file_name_new
    
    ! open file containing atomic data
    close(11)
    open(unit=11,  action="read", file = file_name, status="old", position="rewind", &
         & iostat = ios)
    if (ios /= 0) then
       print*, "! can't open file: ", file_name
       stop
    end if

    close(12)
    open(unit=12,  action="write", file = file_name_new, status="unknown", position="rewind", &
         & iostat = ios)
    if (ios /= 0) then
       print*, "! can't open file: ", file_name_new
       stop
    end if

    ! read reference heading
    read(11, *) numLines
    write(12, *) numLines
    
    do i = 1, numLines
       read(11, '(1A150)') text
       write(12, '(1A150)') text       
    end do
    
    ! read number of levels and temperature points available
    read(11, *) nLev, nTemp
    write(12, *) nLev, nTemp
    
    ! allocate space for labels array
    allocate(label(nLev), stat = err)
    if (err /= 0) then
       print*, "! equilibrium: can't allocate label array memory"
       stop
    end if
    
    ! allocate space for tempertaure points array
    allocate(logTemp(nTemp), stat = err)
    if (err /= 0) then
       print*, "! equilibrium: can't allocate logTemp array memory"
       stop
    end if
    
    ! allocate space for transition probability  array
    allocate(a(nLev, nLev), stat = err)
    if (err /= 0) then 
       print*, "! equilibrium: can't allocate a array memory"
       stop
    end if

    ! allocate space for transition probability  array
    allocate(alphaTotal(nLev), stat = err)
    if (err /= 0) then 
       print*, "! equilibrium: can't allocate a array memory"
       stop
    end if
    
    ! allocate space for energy levels array
    allocate(e(nLev), stat = err)
    if (err /= 0) then 
       print*, "! equilibrium: can't allocate e array memory"
       stop
    end if
    
    ! allocate space for statistical weights array
    allocate(g(nLev), stat = err)
    if (err /= 0) then 
       print*, "! equilibrium: can't allocate g array memory"
       stop
    end if
    
    ! allocate space for qom array
    allocate(qom(nTemp, nLev, nLev), stat = err)
    if (err /= 0) then 
       print*, "! equilibrium: can't allocate qom array memory"
       stop
    end if
    
    ! allocate space for collisional strengths array
    allocate(cs(nLev, nLev), stat = err)
    if (err /= 0) then
       print*, "! equilibrium: can't allocate cm array memory"
       stop
    end if
    
    ! allocate space for q eff array
    allocate(qeff(nLev, nLev), stat = err)
    if (err /= 0) then
       print*, "! equilibrium: can't allocate qeff array memory"
       stop
    end if
        
    ! allocate space for tnij array
    allocate(tnij(nLev, nLev), stat = err)
    if (err /= 0) then
       print*, "! equilibrium: can't allocate tnij array memory"
       stop
    end if

    ! allocate space for x array
    allocate(x(nLev, nLev), stat = err)
    if (err /= 0) then
       print*, "! equilibrium: can't allocate x array memory"
       stop
    end if
         
    ! allocate space for y array
    allocate(y(nLev), stat = err)
    if (err /= 0) then
       print*, "! equilibrium: can't allocate y array memory"
       stop
    end if
    
    ! allocate space for qq array
    allocate(qq(nTemp), stat = err)
    if (err /= 0) then
       print*, "! equilibrium: can't allocate qq array memory"
       stop
    end if
    
    ! allocate space for qq2 array
    allocate(qq2(nTemp), stat = err)
    if (err /= 0) then
       print*, "! equilibrium: can't allocate qq2 array memory"
       stop
    end if
    
    ! allocate space for n array
    allocate(n(nLev), stat = err)
    if (err /= 0) then
       print*, "! equilibrium: can't allocate n array memory"
       stop
    end if

    ! allocate space for n array
    allocate(qxar(nLev,nlev,ntemp), stat = err)
    if (err /= 0) then
       print*, "! equilibrium: can't allocate n array memory"
       stop
    end if
    
    ! allocate space for n2 array
    allocate(n2(nLev), stat = err)
    if (err /= 0) then
       print*, "! equilibrium: can't allocate n2 array memory"
       stop
    end if
    
    ! zero out arrays
    a = 0.
    cs = 0.
    e = 0.
    fLineEm = 0.
    g = 0
    n = 0.
    n2 = 0.
    qom = 0.
    qeff = 0.
    qq = 0.
    qq2 = 0.
    logTemp = 0.
    tnij = 0.
    x = 0.d0
    y = 0.

    ! read labels
    do i = 1, nLev
       read(11, '(A20)') label(i)       
    end do

    ! read temperetaure points
    do i = 1, nTemp
       read(11, *) logTemp(i)
    end do
    
    ! read iRats
    read(11, *) iRats

    qxar=0.
    itt = 1
    do i = 1, safeLim
       ! read in data
       read(11, *) ilow(2), iup(2), qx
 
       if (ilow(2) > 0 ) il = ilow(2)
       if (iup(2) > 0 ) iu = iup(2)
       
       ! check if end of qx dat
       if (qx == 0.d0) exit

       qxar(il,iu,itt) = qx

       itt = itt+1
       if (itt>nTemp) itt = 1
       
    end do


    ! read in transition probabilities
    do k = 1, nLev-1
       do l = k+1, nLev
          read(11, *) i, j, ax
          
          a(j,i) = ax

       end do
    end do
    
    ! read statistical weights, energy levels [1/cm]
    do j = 1, nLev
       read(11, *) i, gx, ex
       
       g(i) = gx
       e(i) = ex
    end do

    ! sort by increasing energy of levels.
    allocate(orderLev(nlev)) 
    orderLev = 0
    
    minlev = 1.e30
    do i = 1, nlev
       if (e(i)<minlev) then
          minlev = e(i)
          orderLev(1) = i
       end if
    end do

    do i = 2, nlev
       minlev = 1.e30
       do j = 1, nlev
          if (e(j)<minlev .and. e(j)>minlev(i-1)) then
             minlev = j
             orderLev(i) = j
          end if
       end do
    end do
    
    do i= 1, nlev 
       print*, i, orderLev(i), e(orderlev(i))
    end do

    


    ! read power law fit coefficients [e-13 cm^3/s]
    ! and calculate total recombination coefficient
    ! (direct + cascades)
    alphaTotal = 0.
    
    do j = 2, nLev
!       read(unit=11,fmt=*,iostat=ios) a_fit, b_fit
       read(unit=11,fmt=*,iostat=ios) br,z,a_r(:),a_d(:) !a_fit, b_fit
       if (ios<0) then
          exit
       else
          alphaTotal(1) = 1.
       end if
                  
!       alphaTotal(j) = a_fit * (TeUsed/1.e4)**(b_fit) * 1.e-13
print*, 'tt ', -a_d(5)*1.e4/Te, a_d(5), Te
       alphaTotal(j)=1.e-13*br*z*a_r(1)*(Te*1.e-4/z**2)**a_r(2)/(1.+a_r(3)*&
            & (Te*1.e-4/Z**2)**a_r(4))+1.0e-12*br*(a_d(1)& 
            & /(Te*1.e-4)+a_d(2)+a_d(3)*Te*1.e-4+a_d(4)*(Te*1.e-4)**2)*&
            & (Te*1.e-4)**(-1.5)*exp(-a_d(5)*1.e4/Te)              

    end do

    ! close atomic data file 
    close(11)
    
    ! form matrices
    ! set up qeff
    do i = 2, nLev
       do j = i, nLev
          do iT = 1, nTemp
             qq(iT) = qom(iT, i-1, j)
          end do
print*, i,j, ntemp
          if (nTemp == 1) then  
             ! collisional strength available only for one temperature - assume constant
             qomInt = qq(1)
          else if (nTemp == 2) then
             ! collisional strength available only for one temperature - linear interpolation
print*, 'aa', (qq(2)-qq(1)),  (logTemp(2)-logTemp(1)), (log10Te - logTemp(1))
             qomInt = qq(1) + &
                  (qq(2)-qq(1)) / (logTemp(2)-logTemp(1)) * (log10Te - logTemp(1))
print*, 'bb'
          else
print*, 'spline'
             ! set up second derivatives for spline interpolation
             call spline(logTemp, qq, 1.e30, 1.e30, qq2)
print*, 'splint'             
             ! get interpolated qom for level
             call splint(logTemp, qq, qq2, log10Te, qomInt)
print*, 'out splint'
          end if
          
          ! set collisional strength
          cs(i-1, j) = qomInt
print*, 'qomint ', qomint          
          ! the calculation constant here is the energy [erg] 
          ! associated to unit wavenumber [1/cm] divided by the 
          ! boltzmann constant k.
          constant = 1.4388463d0

          ! exponential factor 
print*, e(i-1), e(j), constant, file_name
          delTeK = (e(i-1)-e(j))*constant
print*, deltek, ' deltek', Te
          if(expFac < -100.) then
             expFac = 0.
          else
             expFac = exp( delTeK/Te )
          end if
print*, expFac          
          qeff(i-1, j) = 8.63d-6 * cs(i-1, j) * expFac /&
               &(g(i-1)*sqrTe)
print*,qeff(i-1, j) 
          qeff(j, i-1) = 8.63d-6 * cs(i-1, j) / (g(j)*sqrTe)
print*, qeff(j, i-1), 'end'
       end do
    end do

    ! set up x
print*, '    ! set up x'
    do i= 2, nLev
       do j = 1, nLev
          
          x(1,:) = 1.
!          y = 0.
          y(1)   = 1.
          
          if (j /= i) then
             x(i, j) = x(i, j) + Ne*qeff(j, i)
             x(i, i) = x(i, i) - Ne*qeff(i, j)
             if (j > i) then
                x(i, j) = x(i, j) + a(j, i)
             else 
                x(i, i) = x(i, i) - a(i, j)
             end if

          end if
       end do
       
       ! when coefficients are available use the following:
       y(i) = -Ne*ionDenUp*alphaTotal(i)
print*, 'y(i)', y(i)
    end do

    if (nLev < 1) then
       print*, '! equilibrium: nLev < 1!  Check data file', nLev, file_name
    end if

    call luSlv(x, y, nLev)
print*, 'after luslv'
    n = y
    
    sumN = 0.d0
    do i = 1, nLev
print*, 'sumn', sumN, n(i)
       sumN = sumN+n(i)
print*, 'sumn 2', sumN
    end do
    do i = 1, nLev           

       n(i) = n(i)/sumN

    end do
    
    ! now find emissivity and wavelengths
    do i = 1, nLev-1
       do j = i+1, nLev
          if (a(j,i) /= 0.d0) then
             Eji = (e(j)-e(i))

             fLineEm(i,j) = a(j,i) * Eji * n(j)
             if (Eji>0.) then
                if (present(wav)) wav(i,j) = 1.e8/Eji
             else
                if (present(wav)) wav(i,j) = 0.
             end if
          end if
       end do
    end do
print*, 'h'

    ! deallocate arrays
    if( associated(alphaTotal) ) deallocate(alphaTotal)
    if( associated(label) ) deallocate(label)
    if( associated(logTemp) ) deallocate(logTemp)
    if( associated(a) ) deallocate(a)
    if( associated(cs) ) deallocate(cs)
    if( associated(n) ) deallocate(n)
    if( associated(n2) ) deallocate(n2)
    if( associated(qeff) ) deallocate(qeff) 
    if( associated(qq) ) deallocate(qq)
    if( associated(qq2) ) deallocate(qq2)
    if( associated(tnij) ) deallocate(tnij) 
    if( associated(x) ) deallocate(x) 
    if( associated(y) ) deallocate(y) 
    if( associated(e) ) deallocate(e)
    if( associated(g) ) deallocate(g)
    if( associated(qom) ) deallocate(qom)
print*, 'out eq'
  end subroutine 
