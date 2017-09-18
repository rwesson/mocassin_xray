! Copyright (C) 2007 Barbara Ercolano 
!  
! Version 3.00
module common_mod
    use constants_mod
    use vector_mod
    implicit none




    real (kind=8), save,pointer::&
          & forbiddenLines(:,:,:,:)                                 ! emissivity from heavies  rec lines


    real             :: fluoThreshArray(14) ! ID labels of fluorescent transitions

    real, pointer    :: fluorescenceVP(:,:) ! viewing angles for fluorescence cube &
                                            ! (for the first dimension 1=theta;2=phi)
    real, pointer :: Cpr(:,:,:) ! sca cross-section [um^2] for each grain species and size
    real, pointer   :: CTpr(:)              ! total sca cross-section [um^2] for grain mixture
    real, pointer   :: TdustTemp(:,:,:)     ! temporary dust temperature array (species,size,cell)
    real, pointer   :: ionDenTemp(:,:,:)    ! temporary ion density array    
    real, pointer   :: NeTemp(:)            ! temporary electron density array
    real, pointer   :: TeTemp(:)            ! temporary electron Temperature array    
    real            :: HeIrecLineCoeff(34,3,4)
    real, pointer   :: KNsigmaArray(:,:)    ! Klein Nishina calc PDFs (nbins,180)
    real, pointer   :: KNsigmaT(:)          ! Klein Nishina cross-sec integrated over solid angle
    real, pointer   :: PcompArray (:,:) 
    real, pointer   :: TwoDscaleJ(:)

    real            :: starttime      ! start time [sec]
    real            :: endtime        ! end time [sec]    
    real            :: absInt         ! total number of absorption events
    real            :: scaInt         ! total number of scattering events

    real            :: dTheta                  ! 
    real            :: dPhi                    ! 
    real            :: nu0                     ! 
    real            :: nu0Add                  !
    real            :: totalDustMass 
    real            :: convPercent=0.          ! total convergence percentage
    real            :: pwlIndex = 0.           ! power law input spectrum index
    real            :: pwlMin=0.               ! power law lower cutoff [Ryd]
    real            :: pwlMax=0.               ! power law higher cutoff [Ryd]
    real            :: totPercent=0.           ! 
    real            :: minaQHeat = 0.          ! min radius for Q Heat in [um]
    real            :: NeUsed, TeUsed, HdenUsed ! local properties of the gas
    real            :: Tmax                    ! max temp for qheat
    real            :: auger(30,30,10,10)      ! auger yields
    real            :: HbACoeff(9,4)
    real            :: HbBCoeff(9,4)    
    real            :: hydroLinesTemps(12)
    real            :: rbEdge(2,3,26)
    real            :: r2aEdge(3,18)    


    real               :: fluoCubeMineV, fluoCubeMaxeV ! limits of the fluorescence cube band
    real               :: fillingFactor    ! filling factor epsilon
    real               :: contCube(2)      ! continuum cube
    real               :: convIncPercent   ! percentage by  which conv must increase
    real               :: convWriteGrid    ! min conv level before starting to write the grid
    real               :: nPhotIncrease    ! nPhoton increase factor
    real               :: densityLaw(3)    ! density law parameters (R1,n,f,N0)
    real,pointer       :: grainAbun(:)     ! abundance of this species
    real,pointer       :: grainRadius(:)   ! grain radius [um]
    real,pointer       :: grainWeight(:)   ! grain weight normalised to 1
    real,pointer       :: viewPointTheta(:),viewPointPhi(:)     ! viewing angles
    real               :: Hdensity         ! constant H density values (cm-^3)
    real               :: H0Start          ! initial guess at X(H0) for regions I and II
    real               :: Lphot            ! L of ionizing source [e36 phot/sec]
    real               :: Ldiffuse         ! total luminosity of diffuse source [e36 erg/sec]
    real,pointer       :: LStar(:)         ! L of ionizing source [e36 erg/sec]
    real               :: meanField        ! mean ionizing field to be used with plane parallel 
                                           ! geometry [erg/sec/cm^2]
    real               :: meanFieldin
    real               :: minConvergence   ! stop when this level of convergence has been reached  
    real               :: minConvQHeat     ! min convergence level for Qheat routines to run
    real               :: NeStart          ! initial guess at Ne
    real               :: MdMgValue        ! dust-to-gas ratio by mass
    real               :: NdustValue       ! dust number density [cm^-3]
    real               :: nuMax            ! max nu limit [Ryd]
    real               :: nuMin            ! min nu limit [Ryd]
    real               :: nuStepSize       ! 
    real               :: resLinesTransfer ! after what conv percent should the resonant line transfer start?    
    real               :: Rnx,Rny,Rnz      ! edges [cm]
    real               :: R_in             ! inner radius [cm]
    real               :: R_out            ! outer radius [cm]
    real               :: SEDfreq(2)       ! 
    real               :: TeStart          ! initial guess at Te for regions I and II
!    real, pointer      :: deltaE(:)        ! energy carried by a single photon
    real               :: deltaEUsed

    real               :: Tdiffuse         ! energy of diffuse ionising source [K]
    real, pointer      :: Tstellar(:)      ! T of ionizing source [K]    
    real, pointer      :: tStep(:)         ! time step for sb99 inputs [yrs]
    real               :: XHILimit         ! convergence limit on X(HI)

    ! energies temperatures etc..
    real :: Bjump = 0.                        ! Balmer jump


    ! slit stuff
    real :: dxSlit = 0.
    real :: dySlit = 0.

    ! photoionization data 
    real, dimension(5,nElements, nElements) :: CF=0.
    real, dimension(30,30)   :: cionTherm = 0.
    real, dimension(6, 30, 30, 7) :: ph1
    real, dimension(7, 30, 30) :: ph2
    real, dimension(nElements) :: aWeight
    real, dimension(450) :: ionEdge 

    real, pointer :: logGammaHI(:,:),logGammaHeI(:,:),logGammaHeII(:,:), tkGamma(:), nuGammaHI(:),&
         & nuGammaHeI(:),nuGammaHeII(:)

    ! RRDR data Badnell
    real, dimension(1:nElements, 0:nElements) :: RRbA, RRbB, RRbT0, &
         &RRbT1, RRbC, RRbT2 
    real, dimension(9, 1:nElements, 0:nElements) :: DRbC, DRbE


    ! dust arrays
    real, pointer :: TdustSublime(:)        ! sublimation T (nspecies)
    real, pointer :: grainVn(:)             ! potential of neutral grain
    real, pointer :: MsurfAtom(:)           ! mass of surf atom of a grain [amu]
    real, pointer :: rho(:)                 ! intrinsic density (nspecies)    
    real, pointer :: dustHeatingBudget(:,:) ! heating budget of grains (nAbComponents, nResLines+1)       
    real, pointer :: SEDnoExt(:)            ! SED no extinction
    real, pointer :: equivalentTau(:)       ! SED no extinction


    ! linear increments
    real, pointer :: dl(:)

    ! dust arrays
    real, pointer :: dustEmIntegral(:,:,:)   ! dust emission integral (species,size,temperature)
    
    ! allocatable arrays
    real, pointer :: absOpacSpecies(:,:)     ! absOpacSpecies (nspecies, freq)

    real, pointer :: continuum(:)            ! continuum array (freq)

    real, pointer :: gauntFF(:)              ! gaunt factor

    real, pointer :: gauntFFHeII(:)          ! gaunt ff for HeII

    real, pointer :: H2phot(:)               ! H 2-photon emission
 
    real, pointer :: HeI2phot(:)             ! HeI 2-photon emission

    real, pointer :: HeII2phot(:)            ! HeII 2-photon emission

    real, pointer :: ionDenUsed(:,:)         ! local ionDen array

    real, pointer :: nuArray(:)              ! energy array [Ryd]

    real, pointer :: xSecArray(:)            ! x Section array calculated in xSec_mod

    real, pointer :: xSecRecoil(:)           ! compton recoil xSec

    real, pointer :: comXSecH(:)             ! compton Xsec Heat array

    real, pointer :: comXSecC(:)             ! compton Xsec Cool array

    real, pointer :: widFlx(:)               ! widFlx array

    real, pointer :: photoRateFluo(:,:)      ! photoionisation rate (nfluo, nstages) [sec^-1]

    type(vector), pointer :: starPosition(:) ! ionising source(s) position

    integer, pointer  :: HINuEdgeP(:)     ! pointers to the HI, HeI and HeII 
    integer, pointer  :: HeINuEdgeP(:)    ! series edges in nuArray
    integer, pointer  :: HeIINuEdgeP(:)   !

    integer :: nTbins                     ! # of temp bins in qheat
    integer :: nlimGammaHI,nlimGammaHeI, nlimGammaHeII, ntkGamma

    integer :: nIterateMC = 1                  ! current MC interation # 

    integer :: nflux = 1                       ! number of energy bins to the high energy end 
                                               ! of the incident spectrum
    integer :: nLines = 0                      ! total number of lines 

    integer :: iOrigin, jOrigin, kOrigin       ! pointers to the origin of the axes

    ! stack pointers

    integer :: KshellLimitP                    ! KshellLimitP
    integer :: FeKaColdP                       ! pointer to 6.4keV in nuArray
    integer :: FeL1P                           ! pointer tp 0.792 keV in nuArray
    integer :: FeL2P                           ! pointer tp 0.705 keV in nuArray
    integer :: CKaP                            ! pointer tp 0.277 keV in nuArray
    integer :: NKaP                            ! pointer tp 0.3924 keV in nuArray
    integer :: OKaP                            ! pointer tp 0.5924 keV in nuArray
    integer :: NeKaP                           ! pointer tp 0.8486 keV in nuArray
    integer :: MgKaP                           ! pointer tp 1.2536 keV in nuArray
    integer :: AlKaP                           ! pointer tp 1.4863 keV in nuArray
    integer :: SiKaP                           ! pointer tp 1.7395 keV in nuArray
    integer :: SKaP                            ! pointer tp 2.307 keV in nuArray
    integer :: ArKaP                           ! pointer tp 2.9566 keV in nuArray
    integer :: CaKaP                           ! pointer tp 3.690 keV in nuArray
    integer :: TiKaP                           ! pointer tp 7.4781 keV in nuArray
    integer :: secIonP
    integer :: cRecoilP                        ! pointer to compton recoil (194.Ryd) in nuArray
    integer :: xrayP                           ! pointer to carbon k-shell ionization in nuArray
    integer :: cRecoilHP, cRecoilHeP           ! recoil pointers

    integer :: BjumpP                          ! pointer to balmer jump

    integer :: nAuger(30,30,10)                ! max number of auger electrons freed

    integer :: noGpLoc=-1    
    integer, dimension(3) :: noCellLoc=-1
    
    integer, dimension(nElements, nElements, 7, 3) &
         &:: elementP                          ! first dim is atomic number of element,
                                               ! second dimension is ion stage, 1 for atom
                                               ! third is shell, 1 for k shell, up to 7
                                               ! last is 1 for pointer to energy of threshold
                                               ! 2 is highest energy for shell
                                               ! 3 is opacity offset


    integer,pointer :: dustScaXsecP(:,:)       ! pointer to dust scatterring x-sec in xSecArray (species,size)
    integer,pointer :: dustAbsXsecP(:,:)       ! pointer to dust absorption x-sec in xSecArray  (species,size)



    integer, dimension(nHlevel+1) &
         & :: HlevNuP = 0.     ! pointer to the nth H level in nuArray

    integer, dimension(nHeIlevel+1) &
         & :: HeIlevNuP = 0.   ! pointer to the nth HeI level in nuArray


    integer, dimension(nHeIIlevel+1) &
         & :: HeIIlevNuP = 0.  ! pointer to the nth HeII level in nuArray

   

    integer         :: radio4p9GHzP   ! 

    integer         :: ierr           ! MPI error status
    integer         :: numtasks       ! total # of processes
    integer         :: taskid         ! process identification #

    integer, pointer :: lgConvergedTemp(:)  ! temporary converged? flag
    integer, pointer :: lgBlackTemp(:)      ! temporary converged? flag

    integer          :: nFluo ! # of fluorescent lines to be transfered
    integer          :: nVP   ! # of fluorescent viewing angles

    integer, dimension(nElements, nElements)&
         & :: nShells          ! number of shells for each element and ionization stage
                                               ! first dimension is atomic weight, second is ionization stage     


    integer, pointer :: planeIonDistribution(:,:) ! initial distribution of ionising photons 


    integer,pointer    :: viewPointPtheta(:), viewPointPphi(:)       ! viewing angles
    integer            :: nAngleBins=0     ! number of viewing angles for SED
    integer            :: TotAngleBinsTheta=20 ! total # of theta angle bins for SED
    integer            :: TotAngleBinsPhi=20 ! total # of phi angle bins for SED
    integer            :: nGrids           ! total number of grids to be used in the simulation
    integer            :: maxIterateMC     ! limit on number of MC iterations
    integer            :: maxPhotons       ! limit to packets to be used
    integer            :: nAbComponents=1  ! number of abundance components            
    integer            :: nBins            ! number of energy bins 
    integer            :: nElementsUsed    ! actual number of elements used
    integer, pointer   :: nPhotons(:)      ! # of packets to be used in the sim
    integer            :: nPhotonsTot      ! # of packets to be used in the sim
    integer            :: nPhotonsDiffuse  ! # of packets to be used by the diffuse ionisation source
    integer            :: nPhotonsDiffuseLoc! # of packets to be used by the diffuse ionisation source
    integer            :: nPacketsFluo     ! # of packets to be used for each fluorescent transition 
    integer            :: nPacketsFluoLoc  ! # of packets to be used for each fluorescent trans in one cell
    integer, dimension(maxGrids) &
                      &:: nxIn,nyIn, nzIn  ! x, y and z dimensions of the grids
    integer            :: nStars           ! number of ionising sources
    integer            :: elementXref(nElements) ! x reference index array for elements actually used
    integer            :: emittingGrid     ! grid emiting illuminating radiation [0 for all]
    integer, pointer   :: starIndeces(:,:) ! (nstars, 3) 1=x; 2=y; 3=z, 4=gp
    integer            :: nstages          ! # of ionisation stages to be included

    integer             :: lymanP                  ! frequency pointer to Lyman limit
    integer             :: nu1PFeKa, highNuPFeKa
    integer             :: nu1PFeL1, highNuPFeL1
    integer             :: nu1PFeL2, highNuPFeL2
    integer             :: nu1PCKa , highNuPCKa 
    integer             :: nu1PNKa , highNuPNKa 
    integer             :: nu1POKa , highNuPOKa 
    integer             :: nu1PNeKa, highNuPNeKa
    integer             :: nu1PMgKa, highNuPMgKa
    integer             :: nu1PAlKa, highNuPAlKa
    integer             :: nu1PSiKa, highNuPSiKa
    integer             :: nu1PSKa , highNuPSKa 
    integer             :: nu1PArKa, highNuPArKa
    integer             :: nu1PCaKa, highNuPCaKa
    integer             :: nu1PTiKa, highNuPTiKa

    ! dust parameters
    integer            :: nSizes           ! number of grain sizes
    integer            :: nSpecies         ! number of grain species
    integer            :: nResLines=0      ! number of resonant lines to be transfered

    type grid_type             ! derived grid type
        sequence

        character(len=20) :: composition            ! chemical composition

        integer :: nx                               ! size of cartesian grid in x
        integer :: ny                               ! size of cartesian grid in y
        integer :: nz                               ! size of cartesian grid in z
        integer :: nCells                           ! number of active cells
        integer :: motherP                          ! pointer to the motherGrid

        integer, pointer :: resLinePackets(:)       ! extraPackets to be generated cell (cell)
        integer, pointer :: abFileIndex(:,:,:)      ! abundance file index axxay
        integer, pointer :: lgConverged(:)          ! has the model converged at this grid cell?
        integer, pointer :: lgBlack(:)              ! is this to remain a black cell?

        integer, pointer :: active(:,:,:)           ! point to active cell in 1d array; 
                                                    ! returns 0 for inactive cells        

        real    :: geoCorrX,geoCorrY,geoCorrZ       ! geometric correction
        real    :: noHit                            ! cells not sampled by rad field
        real    :: noIonBal                         ! Ion Balance not raeched 
        real    :: noTeBal                          ! Te Balance not reached

        real, pointer :: arad(:)                    ! radiative accelerations - only used if lgRadPress = .t.
        real, pointer :: absOpac(:,:)               ! dust absorption opacity 
        real, pointer :: scaOpac(:,:)               ! dust scattering opacity 
        real, pointer :: dustPDF(:,:)               ! dust emission PDF 
        real, pointer :: fEscapeResPhotons(:,:)     ! frac of res line photons escaping (cell,resline)
        real, pointer :: fluorescenceCube(:,:,:)    ! escaped fluorescence packets (ncell,nFluo,nVP)
        real, pointer :: Hden(:)                    ! density of H [1/cm^3]
        real, pointer :: Ne(:)                      ! electron density [1/cm^3]
        real, pointer :: Te(:)                      ! electron temperature [K]
        real, pointer :: NeOld(:)                   ! electron density [1/cm^3]
        real, pointer :: TeOld(:)                   ! electron temperature [K]
        real, pointer :: Tdust(:,:,:)               ! dust temperature (species,size,cell) [K] 
        real, pointer :: ionDen(:,:,:)              ! fractional ion density (x,y,z,elem,stage)
        real, pointer :: Jste(:,:)                  ! MC estimator of stellar J (cell,nu) 
        real, pointer :: H(:,:)                     ! MC estimator of H (cell,nu) 
        real, pointer :: Jdif(:,:)                  ! MC estimator of diffuse J (cell,nu)
        real, pointer :: JPEots(:,:)                ! OTS line contribution to photoelectric emission
        real, pointer :: escapedPackets(:,:,:)      ! escaped packets (cell,nu, angle)
        real, pointer :: escapedPacketsComponents(:,:,:,:) ! escaped packets (cell,nu, angle)
        real, pointer :: linePackets(:,:)           ! line packets (x,y,z,n)
        real, pointer :: LFluorescence(:,:)         ! local fluorescence luminosity of cell [e36 erg/sec]
        real, pointer :: LdiffuseLoc(:)             ! loc luminosity of diffuse source [e36 erg/sec]
        real, pointer :: Ndust(:)                   ! number density for dust

        real, pointer :: xAxis(:)                   ! x-axis
        real, pointer :: yAxis(:)                   ! y-axiss
        real, pointer :: zAxis(:)                   ! z-axis

        real, pointer :: NeInput(:)                 ! only used if lgNeInput=.t.


        real, pointer :: opacity(:,:)               ! total  opacity (cell,nu) [1/cm]
        real, pointer :: recPDF(:,:)                ! rec prob distribution function
        real, pointer :: linePDF(:,:)               ! line prob distribution function
        real, pointer :: totalLines(:)              ! fraction of non-ionizing line phots
        real, pointer :: totalEmission(:)           !
        real, pointer :: elemAbun(:,:)              ! elemental abundance (specified
                                                    ! by number relative to total 
                                                    ! hydrogen density)
    end type grid_type


    character(len=4), pointer :: fluorescenceLabel(:) ! ID labels of fluorescent transitions
    character(len=4)          :: fluoLabArray(14)     ! ID labels of fluorescent transitions

    character(len=20), dimension(3:nElements, 1:31) :: &
         & dataFile       ! name of the file containing the atomic data for this species

    character(len=20), dimension(1:9, 1:12) :: &
         & hydroLinesFile       ! name of the file containing the emiss for hydrogenic rec lines

    character(len=2)        :: elemLabel(30) 

    character(len=50)  :: gridList ! grid list file name
    character(len=50), pointer  :: abundanceFIle(:) ! abundance file names
    character(len=50)           :: contDiffuse      ! shape of the diffuse ionising spectrum
    character(len=50), pointer  :: contShape(:)     ! continuumShape
    character(len=50), pointer  :: contShapeIn(:)   ! continuumShape
    character(len=50), pointer  :: spID(:)          ! input spectrum generator ID
    character(len=50)  :: densityFile      ! density file
    character(len=50)  :: dustFile(2)      ! dust files
    character(len=50)  :: MdMgFile         ! name of MdMg file
    character(len=50)  :: NdustFile        ! name of Ndust file    
    character(len=50)  :: Qfile            ! name of Qfile 
    character(len=30),pointer       :: grainLabel(:)    ! name of this species

    type resLine_type
       character(len=2)      :: species
       character(len=30)     :: transition, multiplet! label

       integer               :: nuP  ! pointer to freq of line in nuArray
       integer               :: elem ! elemnt
       integer               :: ion  ! ionization stage
       integer               :: glow, ghigh ! lower and upper statistical weights
       integer               :: nmul ! # of components in the multiplet
       integer,pointer       :: moclow(:), mochigh(:) ! mocassin's lower and higher levels code

       real                  :: m_ion ! mass of the ion
       real                  :: wav  ! wavelength [A]
       real                  :: freq ! frequency [Ryd]
       real                  :: Elow,Ehigh ! lower and higher level energies
       real                  :: Aik,fik ! trans prob and oscillator strength
    end type resLine_type

    logical         :: lgRadPress=.false.
    logical         :: lg2D=.false.             ! 2D?
    logical         :: lgMultiStars=.false.
    logical         :: lgEquivalentTau          ! calculate equivalent tau?
    logical         :: lgWarm=.false.           ! warm started?
    logical         :: lgFluorescence=.false.   ! fluorescence run?
    logical         :: lgNeInput=.false.        ! Ne distribution entered
    logical         :: lgResLinesFirst = .true. ! first time the res lines transfer proc is called? 
    logical         :: lgPhotoelectric = .true. ! photoelectric effect on
    logical         :: lgQheat = .false.        ! temperature spiking on? 
    logical         :: lgSeparateSED = .false.  ! write out separate SED components?
    logical         :: lgWritePss=.false.       ! write Pss for qHeat file?
    logical         :: lgTraceHeating = .false. ! trace thermal balance?
    logical         :: lgGrainSpiking = .true.  ! temperature spiking for small grains
    logical            :: lgBadnell        ! badnell rec coeffs?
    logical            :: lgHomoFixed=.false.      ! 
    logical            :: lgCompton        ! evaluate Compton energy exchange?
    logical            :: lgNahar          ! nahar rec coeffs?
    logical            :: lgAutoPackets    ! automatic increase of packets on the fly? 
    logical            :: lgTalk           ! talk on?
    logical            :: lgDfile          ! use an external density file?
    logical            :: lgDebug          ! debigging mode? (memory consuming)
    logical            :: lgDlaw           ! use the density law routine?
    logical            :: lgDust           ! is dust included in this model?
    logical            :: lgDustScattering ! is dust scattering included in this model?
    logical            :: lgDustConstant   ! is the md/Mg constant throughtout the grid?
    logical            :: lgGas            ! is there any gas in teh grid? 
    logical            :: lgHdenConstant   ! use the constant density parameters?
    logical            :: lgMdMg           ! dust to gas mass ratio used>? 
    logical            :: lgMdMh           ! dust to hydrogen mass ratio used>? 
    logical            :: lgOutput         ! output line fluxes,temp & ion struct files?
    logical            :: lg1D             ! 1-D switch
    logical            :: lgMultiChemistry ! do we have a chemically inhomogeneous gas? 
    logical            :: lgNeutral        ! starting from neutral gas?
    logical            :: lgPlaneIonization! plane parallel ionization?    
    logical            :: lgSymmetricXYZ   ! symmetric in x, y, and z?


    logical, pointer :: &
         & lgDataAvailable(:,:)! is the atomic data available for this species?


    ! calculate recombination lines?
    logical :: lgRecombination=.false.

    ! composition switches
    logical, dimension(nElements) :: lgElementOn ! is this element present?


    ! photon packet type
    type photon_packet

        integer       :: nuP       ! pointer
        integer, dimension(2) & 
             & :: xP,yP,zP         ! grids position indeces 1= mother 2=sub
        integer, dimension(2) :: origin ! 1=ig, 2=icell
        integer       :: iG        ! grid index 

        integer       :: SEDtype   ! 0 = direct stellar
                                   ! 1 = compton reflected
                                   ! 2, N = fluorescence

        real          :: nu

        logical       :: lgStellar
        logical       :: lgLine
    
        type(vector)  :: position
        type(vector)  :: direction
  
    end type photon_packet


    ! plot type
    type plot_type
       real, pointer         :: intensity(:,:,:)
       real, pointer         :: xAxis(:), yAxis(:)

       logical               :: lgFilter
       logical               :: lgMap
       logical, pointer      :: lgLine(:)

       integer, pointer      :: lineNumber(:)
       integer               :: nPlots
       integer, pointer      :: nuP(:,:)
       
    end type plot_type

    type(plot_type) :: noPlot
    
  
end module common_mod
   

    

