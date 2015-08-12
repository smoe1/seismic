
c
c
c
c     =====================================================
       subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
       implicit double precision (a-h,o-z)
       dimension q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
       dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
       common /comaux/ rho1,amu1,alam1,rho2,amu2,alam2,rho3,amu3,alam3

       double precision Beta,x0,y0,lam,mu,xj,yj,x1,y1,src11,src22,src12
       double precision tbeta,tfun,srcu,srcv


        !Beta=500.d0
        Beta=200.d0
        x0=0.5d0
        y0=0.5d0
        tbeta=50.d0
       cp2 = dsqrt((alam2+2.d0*amu2)/rho2)
       do 20 i=1,mx
          xi = xlower + (i-0.5d0)*dx
          do 20 j=1,my


             lam = aux(2,i,j)
             mu = aux(3,i,j)
             xj = xlower + (i-0.5d0)*dx
             yj = ylower + (j-0.5d0)*dy
             fun1=1.0e-3*exp(-Beta*((xj-x0)**2+(yj-y0)**2))
             if (yj<0.5) then
                fun1=-fun1
                !fun1=fun1
             endif
             tfun=exp(-tbeta*0.d0)

             x1=-(-2.d0*Beta)*(xj-x0)*fun1
             x2=-(-2.d0*Beta)*(yj-y0)*fun1


             src11    = (lam+2.d0*mu)*x1
             src22    = lam*x1
             src12    = mu*x2

             srcu     = fun1*(-tbeta*tfun)
             srcv     = 0.d0


             !q(1,i,j) = src11
             !q(2,i,j) = src22
             !q(3,i,j) = src12
             !q(4,i,j) = srcu
             q(1,i,j) = 0.d0
             q(2,i,j) = 0.d0
             q(3,i,j) = 0.d0
             q(4,i,j) = 0.d0
             q(5,i,j) = 0.d0
             q(6,i,j) = 0.d0
             q(7,i,j) = 0.d0
 20          continue

       return
       end

