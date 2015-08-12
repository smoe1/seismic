
subroutine b4step2(mbc,mx,my,meqn,q,xlower,ylower,dx,dy,t,dt,maux,aux)

    ! Called before each call to step2.
    ! Use to set time-dependent aux arrays or perform other tasks.
    !
    ! This default version does nothing. 
 
    implicit none
    integer, intent(in) :: mbc,mx,my,meqn,maux
    integer :: i,j
    double precision xi
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,dt
    real(kind=8), intent(inout) :: q(1:meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

       do 20 i=1,mx
          xi = xlower + (i-0.5d0)*dx
          do 20 j=1,my
               q(6,i,j) = q(6,i,j)+dt*q(4,i,j)
               q(7,i,j) = q(7,i,j)+dt*q(5,i,j)
       20      continue
           
end subroutine b4step2
