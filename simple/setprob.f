      subroutine setprob
      implicit double precision (a-h,o-z)
      character*12 fname
      common /comaux/ rho1,amu1,alam1,rho2,amu2,alam2,rho3,amu3,alam3
      common /combc/ t0wall,tperiod,pi2,amplitude,t0force

c
c     # Set the material parameters for the elasticity equations
c
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)
                

c
c     # Piecewise constant medium 
c     # Material parameters

      read(7,*) rho1
      read(7,*) alam1
      read(7,*) amu1

      read(7,*) rho2
      read(7,*) alam2
      read(7,*) amu2

      read(7,*) rho3
      read(7,*) alam3
      read(7,*) amu3

      read(7,*) t0wall
      read(7,*) tperiod
      read(7,*) amplitude
      read(7,*) t0force
      pi2 = 8.d0*datan(1.d0)

      return
      end

