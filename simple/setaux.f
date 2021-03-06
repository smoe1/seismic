c     ============================================
      subroutine setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)
c     ============================================
c
c     # set auxiliary arrays 
c     # variable coefficient acoustics
c     #  aux(1,i,j) = density rho in (i,j) cell
c     #  aux(2,i,j) = lambda in (i,j) cell
c     #  aux(3,i,j) = mu in (i,j) cell
c     #  aux(4,i,j) = cp in (i,j) cell
c     #  aux(5,i,j) = cs in (i,j) cell
c
c     # Piecewise constant medium
c     # Material parameters are set in setprob.f

c
c     
      implicit double precision (a-h,o-z)
      dimension aux(7,1-mbc:mx+mbc,1-mbc:my+mbc)
      common /comaux/ rho1,amu1,alam1,rho2,amu2,alam2,rho3,amu3,alam3


      do 30 j=1-mbc,my+mbc
       do 20 i=1-mbc,mx+mbc
          xl = xlower + (i-1.0d0)*dx
          yl = ylower + (j-1.0d0)*dy
          call cellave(xl,yl,dx,dy,w1)
          w2 = 1.d0 - w1

          if (yl .gt. 0.5d0) then
              aux(1,i,j) = w1*rho1 + w2*rho2
              aux(2,i,j) = w1*alam1 + w2*alam2
              aux(3,i,j) = w1*amu1 + w2*amu2
            else
              aux(1,i,j) = w1*rho3 + w2*rho2
              aux(2,i,j) = w1*alam3 + w2*alam2
              aux(3,i,j) = w1*amu3 + w2*amu2
            endif
          bulk       = aux(2,i,j) + 2.d0*aux(3,i,j)
          aux(4,i,j) = dsqrt(bulk/aux(1,i,j))
          aux(5,i,j) = dsqrt(aux(3,i,j)/aux(1,i,j))
       
          aux(6,i,j)=xl+0.5*dx
          aux(7,i,j)=yl+0.5*dy

   20     continue
   30    continue

       return
       end

