program makeDustShell
  implicit none 

  real(kind = 8) :: rin, rout, radius, normWeight
  real(kind = 8) :: dx, dy, dz, den, dx1
  real(kind = 8) :: mdmg, sums, Pi = 3.1416
  real(kind = 8), pointer :: xaxis(:), yaxis(:), & 
       & zaxis(:), rho(:), rad(:), a(:), &
       & wa(:), abun(:), rhospec(:), da(:)

  integer :: nx, ny, nz, nrho, iskip, irho
  integer :: i, j, k, nspec, nsizes, ii, ai

  character(len=50) :: sizefile, speciesfile, skip

  open(unit=10, file='dustShellInput.in') 

  read(10, *) rin, rout
  read(10, *) mdmg
  read(10, *) sizefile
  read(10, *) speciesfile
  read(10, *) nx, ny, nz

  close(10)

  open(unit=10, file='rhoClotho.dat')   
  read(10,*) nrho

  allocate(rad(nrho))
  allocate(rho(nrho))

  do i = 1, nrho
     read(10,*) rad(i), rho(i)
  end do

  close(10)

  open(unit=10, file=speciesfile)
  read(10,*) nspec 
  allocate(abun(nspec))
  allocate(rhospec(nspec))
  do i = 1, nspec
     read(10,*) skip, abun(i), rhospec(i)
  end do
  close(10)

  open(unit=10, file=sizefile)
  read(10,*) nsizes
  allocate(a(nsizes))
  allocate(wa(nsizes))
  allocate(da(nsizes))
  do i = 1, nsizes
     read(10,*) iskip, a(i), wa(i)
  end do
  close(10)


  if (nSizes>1) then
     da(1) = a(2)-a(1)                    
     do ai = 2, nSizes-1
        da(ai) = (a(ai+1)-a(ai-1))/2.              
     end do
     da(nSizes) = a(nSizes)-a(nSizes-1)
  end if
  normWeight=  0.
  do ai = 1, nSizes
     normWeight = normWeight+wa(ai)*da(ai) 
  end do
  if (nSizes>1) then
     do ai = 1, nSizes
        wa(ai) = (wa(ai)*da(ai))/normWeight
        if (.not.wa(ai)>=0.) then
           print*, '! makeDustXSec : Invalid grain weight ', wa(ai), ai
           stop
        end if
     end do
  end if
  
  print*, ' index, a [um], da [um], weight '
  do ai = 1, nSizes
     print*, ai, real(a(ai)), real(da(ai)), real(wa(ai))
  end do

  a = a*1.e-4 ! cm
  sums = 0.

  do i = 1, nspec
     do j = 1, nsizes
        sums = sums + rhospec(i)*abun(i)*&
             & (4./3.)*Pi*wa(j)*a(j)**3

     enddo
  enddo

  allocate(xaxis(nx+14))
  allocate(yaxis(ny+14))
  allocate(zaxis(nz+14))

  xaxis = 0.
  yaxis = 0.
  zaxis = 0.


!  dx = (rout-rin)/(nx-1)
!  dy = (rout-rin)/(ny-1)
!  dz = (rout-rin)/(nz-1)  

  dx = (log10(rout)-log10(rin))/(nx-1)
  dy = (log10(rout)-log10(rin))/(ny-1)
  dz = (log10(rout)-log10(rin))/(nz-1)

  dx1 = rin/14.

  do i = 2, 14     
     xaxis(i) = (xaxis(i-1)+dx1)
  end do
  do i = 2, 14
     yaxis(i) = (yaxis(i-1)+dx1)
  end do
  do i = 2, 14  
     zaxis(i) = (zaxis(i-1)+dx1)
  end do

  xaxis(2:14) = log10(xaxis(2:14))
  yaxis(2:14) = log10(yaxis(2:14))
  zaxis(2:14) = log10(zaxis(2:14))

  do i = 15, nx+14     
     xaxis(i) = log10(rad(i-14))
  end do
  do i = 15, ny+14
     yaxis(i) = log10(rad(i-14))
  end do
  do i = 15, nz+14
     zaxis(i) = log10(rad(i-14))
  end do

  xaxis(2:nx+14) = 10.**xaxis(2:nx+14)
  yaxis(2:nx+14) = 10.**yaxis(2:nx+14)
  zaxis(2:nx+14) = 10.**zaxis(2:nx+14)

  open(unit=12, file = 'ndust.out') 

  do i = 1, nx+14
     do j = 1, ny+14
        do k = 1, nz+14

           radius = 1.e10*sqrt( (xaxis(i)/1.e10)*(xaxis(i)/1.e10) + &
                & (yaxis(j)/1.e10)*(yaxis(j)/1.e10) + &
                & (zaxis(k)/1.e10)*(zaxis(k)/1.e10) )
           
           irho = 0

           do ii = 1, nrho

              if(rad(ii) > radius) then
                 irho = ii-1
                 exit
              end if
           end do
           if (irho == 0) then
              den = 0.
           elseif (irho ==nrho) then
              den = rho(nrho)
           elseif (irho > 0 .and. irho < nrho) then
              den = rho(irho+1) + (rho(irho)-rho(irho+1))*&
                   & (rad(irho+1)-radius)/(rad(irho+1)-rad(irho))              
           else
              print*, 'insanity in rho location'
              stop
           end if           

           den = den*0.7/1.67d-24

           den = den*3.32d-10


           write(12, *) real(xaxis(i)), real(yaxis(j)), real(zaxis(k)), real(den)

        end do

     end do
  end do

  close(12)

!  print*, mdmg/sums, 3.32d-10*0.7/1.67d-24, (3.32d-10*0.7/1.67d-24)/(mdmg/sums)

end program makeDustShell
