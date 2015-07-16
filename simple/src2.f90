subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)

    ! Called to update q by solving source term equation 
    ! $q_t = \psi(q)$ over time dt starting at time t.
    !
    ! This function is attempting to simulate an earthquake by forcing the displacements to change on either side of an interface


    ! Auxiliary variables:
    !       1  density
    !       2  lamda
    !       3  mu
    !       4  cp
    !       5  cs

 
    implicit none
    integer, intent(in) :: mbc,mx,my,meqn,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,dt
    real(kind=8), intent(in) ::  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) ::  q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

    double precision y,v,mu,lam,beta,fun1,sc
    double precision tfun,tbeta,ut,vt
    double precision x,u,x11,x12,x22,x0,y0
    double precision t0wall,tperiod,pi2,amplitude,t0force
    integer j1,j
    integer i1,i
    common /combc/ t0wall,tperiod,pi2,amplitude,t0force

    
        Beta=300.d0
        tbeta=50.d0
        x0=0.5d0
        y0=0.5d0

    do 20 i=1,mx
      do 20 j=1,my
        x=xlower+(i-0.5)*dx
        y=ylower+(j-0.5)*dy

        !if (x.ge.0.3.and.y.ge.0.3.and.x.le.0.8.and.y.le.0.8) then
          lam = aux(2,i,j)
          mu = aux(3,i,j)

          fun1=1.0e-3*exp(-Beta*((x-x0)**2+(y-y0)**2))
          !if (y .le. 0.5) then
          !   fun1=-fun1
          !endif
          tfun=exp(-tbeta*t)

          x11=(-2.d0*Beta+(2.d0*Beta*(x-x0))**2)*fun1*(1.d0-tfun)
          x12=((2.d0*Beta)**2*(x-x0)*(y-y0))*fun1*(1.d0-tfun)
          x22=(-2.d0*Beta+(2.d0*Beta*(y-y0))**2)*fun1*(1.d0-tfun)

          ut=fun1*(-tbeta*tbeta*tfun)
          vt=0.d0
         

          !sc=Max(t,1.0)*1.0e-2
          sc=1.0*1.0e-2


          !u=-1.0e-2*(lam*x11+mu*(2.d0*x11+x22))
          !v=-1.0e-2*(lam*x12+mu*x12)
          u=(lam*x11+mu*(2.d0*x11+x22))
          v=(lam*x12+mu*x12)


          !print *,i1,j1,xlower,ylower," HERE! ",q(4,i1,j1),q(5,i1,j1)
          q(4,i,j)=q(4,i,j)+dt*u-ut*dt
          q(5,i,j)=q(5,i,j)+dt*v
         !print *,i1,j1,xlower,ylower," HERE2! ",q(4,i1,j1),q(5,i1,j1)
        !endif
20      continue    

end subroutine src2
