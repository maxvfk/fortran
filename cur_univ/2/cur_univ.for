!ѕрограммно-математическое обеспечение, реализующее математический
! алгоритм восстановлени€ пространственного распределени€ тока по магнитному
! полю, измеренному над образцом с током.
        
        module aa
        integer,parameter:: niter=100 !  Iter. number for solution 3-d Laplas eq.
!        integer,parameter:: mmx=10, nx=2**mmx !  nx>2
        integer,parameter:: mmx=9, nx=2**mmx !  nx>2
        integer,parameter:: mmy=9, ny=2**mmy !  ny>2
        integer,parameter:: nz=niter/2 !          
        integer:: nxexp,nyexp, kind
        end module aa
!**********************************
        module aaa
        use aa
        real*8:: dh  ! Thickness of sample
        real*8,parameter:: pi=3.14159265d0
        real*8 :: ax,ay
        real*8 x(nx+1),fx(nx+1),fxnew(nx+1)
        real*8 xk(0:nx),fxk(0:nx)
        real*8 y(ny+1),fy(ny+1),fynew(ny+1)
        real*8 yk(0:ny),fyk(0:ny)
        real*8  hz(nx+1,ny+1)
        real*8  hzexp(nx+1,ny+1)
        real*8  bz3d(nx+1,ny+1,nz+1),dbz3d(nx+1,ny+1,nz+1)
        real*8  xx(nx+1,ny+1),qq(0:nx,0:ny),dxo,dyo
        real*8  qqq(0:nx,0:ny)
        integer  mark(nx+1,ny+1)
        real*8  xxjx(nx+1,ny+1),qqjx(0:nx,0:ny)
        real*8  xxjy(nx+1,ny+1),qqjy(0:nx,0:ny)
        real*8  qqbz(0:nx,0:ny),qqnew(0:nx,0:ny)
        real*8  xxnew(nx+1,ny+1),amj(nx+1,ny+1)
        real*8  bzold(nx+1,ny+1),bznew(nx+1,ny+1),bztes(nx+1,ny+1)
        real*8  b0(nx+1,ny+1)
        real*8  qtest(nx+1,ny+1)
        real*8  a,z,aj
        end module aaa
!************************************************
program subsub
print*,'  '
end program subsub
       subroutine sovsem_main
 use aaa
        implicit real*8(a-h,o-z)
        real*4 ar
        dimension qour(nx+1,ny+1)

        print*,' Nx, Ny ',nx,ny
        call init(dx,dy)
        print*,' dh= ',dh
        print*,'  '

        amax=dsqrt(z**2+dx**2+dy**2) 
        amax=amax/2.

!*********** Creation bznew outside of experem. region due to 
! solution of 3-d Dirichlet problem
        dz=2.d0*z  
        dt=0.5d0/(1.d0/dx**2 +1.d0/dy**2 +1.d0/dz**2)
        do jjj=1,niter ! cycle of Stright iteration
        do iz=1,nz
        do iy=2,ny
        do ix=2,nx
        if(iz==1) then
!********  base lay determination ***************
        d2x=(bz3d(ix-1,iy,iz)-2.d0*bz3d(ix,iy,iz)
     *   +bz3d(ix+1,iy,iz))/dx**2
        d2y=(bz3d(ix,iy-1,iz)-2.d0*bz3d(ix,iy,iz)
     *   +bz3d(ix,iy+1,iz))/dy**2
        d2z=(bz3d(ix,iy,iz)-2.d0*bz3d(ix,iy,iz)
     *   +bz3d(ix,iy,iz+1))/dz**2       
!********  End of base lay determination ***************
        else
        d2x=(bz3d(ix-1,iy,iz)-2.d0*bz3d(ix,iy,iz)
     *   +bz3d(ix+1,iy,iz))/dx**2
        d2y=(bz3d(ix,iy-1,iz)-2.d0*bz3d(ix,iy,iz)
     *   +bz3d(ix,iy+1,iz))/dy**2
        d2z=(bz3d(ix,iy,iz-1)-2.d0*bz3d(ix,iy,iz)
     *   +bz3d(ix,iy,iz+1))/dz**2
        end if
        dbz3d(ix,iy,iz)=dt*(d2x+d2y+d2z)                
        end do
        end do
        end do

        do i=1,nx+1
        do j=1,ny+1
        if(mark(i,j)==0) dbz3d(i,j,1)=0.d0 ! Concervation of exper. Hz
        end do
        end do
        bz3d=bz3d+dbz3d
        end do
        bznew=bz3d(:,:,1)
 
       xx=bznew
         call fstr2
         qq=qqnew
!**************  Filtration  *************************
           qqbz(0,0)=qq(0,0)
         qq(0,0)=0.d0

         do k1=0,nx-1
         do k2=0,ny-1
         ak2=xk(k1)**2+yk(k2)**2
         ak=dsqrt(ak2)
         ak3=ak2*ak
           if(ak>1.d-36) then
         arg=z*ak
         arg1=dh*ak/2.d0
         ex=dexp(arg1)
         sh=0.5d0*(ex-1.d0/ex)
            if(dabs(arg1)<1.d-3 ) then
            ccc=1.d0
            else
            ccc=arg1/sh
            end if
            reg=(dexp(-arg)*ak*ccc)**2     !  Tih filtr
            reg=reg/(reg+(amax*ak2)**2) !


           if(kind==0) then
           correc=1. !  Var 1
!           correc=1.*dsqrt(arg)-1.5*arg !  Var 2
!           correc=1.*(arg)**0.25 !  Var 3           
!           correc=1.*(arg)**0.333 !  Var 4           
!           correc=1.!*(1.+.9*xk(k1)/ak)!*(arg)**0.333 !  Var 5           
!           correc=1.*(arg)**0.333*(1.-.7*yk(k2)/ak) !  Var 6           
  
           qq(k1,k2)=qq(k1,k2)*reg*dexp(arg)*2.d0/(ccc*ak)    ! Mag. mom. determ.
           qq(k1,k2)=qq(k1,k2)*correc
           qqbz(k1,k2)=qq(k1,k2)*dexp(-arg)*ak*ccc/2.d0
           qqbz(k1,k2)=qqbz(k1,k2)/correc
           end if
           if(kind==1) then
           qq(k1,k2)=-qq(k1,k2)*reg*dexp(arg)*2.d0/(ccc*ak2)    ! Mag. mom. determ.
           qqbz(k1,k2)=-qq(k1,k2)*dexp(-arg)*ak2*ccc/2.d0
           end if
           if(kind==2) then
           qq(k1,k2)=qq(k1,k2)*reg*dexp(arg)*2.d0/(ccc*ak3)    ! Mag. mom. determ.
           qqbz(k1,k2)=qq(k1,k2)*dexp(-arg)*ak3*ccc/2.d0
           end if
           if(kind==3) then
           
           correc=1.*(arg)**0.333 !  Var 4           
  
           qq(k1,k2)=qq(k1,k2)*reg*dexp(arg)*2.d0/(ccc*ak)    ! Mag. mom. determ.
           qq(k1,k2)=qq(k1,k2)*correc
           qqbz(k1,k2)=qq(k1,k2)*dexp(-arg)*ak*ccc/2.d0
           qqbz(k1,k2)=qqbz(k1,k2)/correc
           end if
           
           end if
         end do
         end do

         call frev2
         qour=xxnew

         qq=qqbz
         call frev2
         bztes=xxnew
!************** End of filtration  *************************

!***** Jx det. ***********
        xxjx=0.d0
        do i1=2,nx-1
        do i2=2,ny-1
        xxjx(i1,i2)=0.5d0*(qour(i1,i2+1)-qour(i1,i2-1))/dy
        end do
        end do
!***** Jy det. ***********
        xxjy=0.d0
        do i1=2,nx-1
        do i2=2,ny-1
        xxjy(i1,i2)=-0.5d0*(qour(i1+1,i2)-qour(i1-1,i2))/dx
        end do
        end do
!**************** mod J det *********************
        amj=0.d0
        do i1=2,nx-1
        do i2=2,ny-1
        amj(i1,i2)=dsqrt(xxjx(i1,i2)**2+xxjy(i1,i2)**2)
        end do
        end do
!**************** Transport J det *********************
        Ntot=(nx-2)*(ny-2);        Ninner=nx*nyexp
        corrJx=0.d0;        corrJy=0.d0
        nn=0
        do jy=2,ny-1
        do ix=2,nx-1
        if(mark(ix,jy)==0) cycle 
        nn=nn+1
        corrJx=corrJx+xxjx(ix,jy)
        corrJy=corrJy+xxjy(ix,jy)
        end do
        end do
        corrJx=corrJx/dfloat(nn)
        corrJy=corrJy/dfloat(nn)
        
        s=0.;        s1=0.
        do jy=2,ny-1
        do ix=2,nx-1
        s=s+dabs(xxjx(ix,jy))+dabs(xxjy(ix,jy))
        s1=s1+dabs(corrJx)+dabs(corrJy)
        xxjx(ix,jy)=xxjx(ix,jy)-corrJx
        xxjy(ix,jy)=xxjy(ix,jy)-corrJy
        end do
        end do
        cur_norm=s1/s
        tr_cur_X=0.d0;        tr_cur_Y=0.d0
        nn=0
        do jy=2,ny-1
        do ix=2,nx-1
        if(mark(ix,jy).ne.0) cycle 
        nn=nn+1
        tr_cur_X=tr_cur_X+xxjx(ix,jy)
        tr_cur_Y=tr_cur_Y+xxjy(ix,jy)
        end do
        end do
        tr_cur_X=tr_cur_X*dy*dfloat(nyexp)/dfloat(nn)
        tr_cur_Y=tr_cur_Y*dx*dfloat(nx)/dfloat(nn)

        call errdet(err)
          print*,' Err= ',err

!************** Writing ************************
1       format(1x,7e14.5)
        open(1,file='Gmag.dat')
!        do j=1,nyexp
        do j=nyexp/4,nyexp*3/4
!        do j=nyexp/20,nyexp*19/20

        do i=1,nxexp
        i1=i+nx/2-nxexp/2
        i2=j+ny/2-nyexp/2
        write(1,1) x(i1),y(i2),qour(i1,i2)
        end do
        end do
        close(1)


        open(1,file='Jc.dat')
!        do j=1,nyexp
        do j=nyexp/4,nyexp*3/4
!        do j=nyexp/10,nyexp*9/10

        do i=1,nxexp
        i1=i+nx/2-nxexp/2
        i2=j+ny/2-nyexp/2
        write(1,1) x(i1),y(i2),
     *   dsqrt(xxjx(i1,i2)**2+xxjy(i1,i2)**2)
        end do
        end do
        close(1)

        open(1,file='Jx.dat')
!        do j=1,nyexp
        do j=nyexp/4,nyexp*3/4
!        do j=nyexp/10,nyexp*9/10
        do i=1,nxexp
        i1=i+nx/2-nxexp/2
        i2=j+ny/2-nyexp/2
        write(1,1) x(i1),y(i2),xxjx(i1,i2)
        end do
        end do
        close(1)

        open(1,file='Jy.dat')
!        do j=1,nyexp
        do j=nyexp/4,nyexp*3/4
!        do j=nyexp/10,nyexp*9/10

        do i=1,nxexp
        i1=i+nx/2-nxexp/2
        i2=j+ny/2-nyexp/2
        write(1,1) x(i1),y(i2),xxjy(i1,i2)
        end do
        end do
        close(1)


        call cpu_time(ar)
	       print*,' CPU TIME+ ',ar
	       read*
        !end
	end subroutine sovsem_main
***********************************************************
        subroutine init(dx,dy)
        use aaa
        implicit real*8(a-h,o-z)
        dimension hexper(nx,ny)
        dimension bmat(ny,ny+nx)

                open(1,file='input.dat')
        read(1,*)dxo
        read(1,*)nxexp
        read(1,*)dyo
        read(1,*)nyexp
        read(1,*)z
        read(1,*) kind
        read(1,*) dh
        
        close(1)
        
        print*,' dx = ',dxo,' nxexp = ',nxexp
        print*,' dy = ',dyo,' nyexp = ',nyexp
        print*,' Z = ',z
        do
        if(kind==0)then 
        print*,' Hall probe '
        exit
        end if
        if(kind==1) then 
        print*,' AFM ampl probe '
        exit
        end if
        if(kind==2) then
         print*,' AFM phaze probe '
         exit
        end if
        if(kind==3) then
         print*,' AFM empirical phaze probe '
         exit
        end if
        stop ' kind error '
        end do
        

        open(1,file='experem.dat')
        do j=1,nyexp
        do i=1,nxexp
        read (1,*)xxx,yyy, hexper(i,j)
        end do
        end do
        close(1)

!****************  Transformation B(mT) to H(A/mm) ****************
         hexper=hexper*10.d0/(4.d0*pi)
!******************************************************************

        dx=dxo; dy=dyo; dz=2.d0*z
        
        do i=1,nx; x(i)=dx*dfloat(i-nx/2);  end do

        do i=1,ny; y(i)=dy*dfloat(i-ny/2);  end do
  
        bzold=0.d0
!********************   Ribbon period expansion  ***********        
        nblock=nx/nxexp+2   
               
        do iper=0,nblock
        do i=1,nxexp
        ix=2*iper*nxexp+i+nx/2-nxexp/2
        if(ix>nx+1)exit
        do j=1,nyexp
        jy=j+ny/2-nyexp/2
        bzold(ix,jy)=hexper(i,j)
        end do
        end do
        end do
        
        do iper=-1,-nblock,-1
        do i=nxexp,1,-1
        ix=2*iper*nxexp+i+nx/2-nxexp/2
        if(ix<1)exit
        do j=1,nyexp
        jy=j+ny/2-nyexp/2
        bzold(ix,jy)=hexper(i,j)
        end do
        end do
        end do
        
        do iper=0,nblock
        do i=1,nxexp
        ix=2*iper*nxexp+1+2*nxexp-i+nx/2-nxexp/2
        if(ix>nx+1)cycle
        do j=1,nyexp
        jy=j+ny/2-nyexp/2
        bzold(ix,jy)=hexper(i,j)
        end do
        end do
        end do
                
        do iper=-1,-nblock,-1
        do i=nxexp,1,-1
        ix=2*iper*nxexp+1+2*nxexp-i+nx/2-nxexp/2
        if(ix<1)cycle
        do j=1,nyexp
        jy=j+ny/2-nyexp/2
        bzold(ix,jy)=hexper(i,j)
        end do
        end do
        end do
!********************  End of  Ribbon period expansion  ***********
        do iy=1,nyexp
        cy=dy*float(iy-1)
        do iyy=1,nyexp
        cyy=dy*float(iyy-1)
        cr2=z**2+(cy-cyy)**2
        bmat(iy,iyy)=-((cy-cyy)/cr2)*dy/(2.*pi)
        end do
        end do
        
        do i=1,nx
        bmat(1:nyexp,ny+i)=bzold(i,1+ny/2-nyexp/2:nyexp+ny/2-nyexp/2)
        end do

         call sol(bmat)
         
        xxjx=0.d0 
        do j=1,nyexp
        do i=1,nx
        i2=j+ny/2-nyexp/2
        xxjx(i,i2)=bmat(j,ny+i)
        end do
        end do
!************************   Initial volume distr. ********************
        do iz=nz,1,-1
        zh=z+dz*dfloat(iz-1)
        do iy=1,ny
        cy=y(iy)
        do iyy=1,ny
        cyy=y(iyy)
        cr2=zh**2+(cy-cyy)**2
        bmat(iy,iyy)=-((cy-cyy)/cr2)*dy/(2.*pi)
        end do
        end do

         bznew=0.d0 
        do ix=1,nx
        do iy=1,ny
        cy=y(iy)
        do iyy=1,ny
        cyy=y(iyy)
        cr2=z**2+(cy-cyy)**2
        bznew(ix,iy)=bznew(ix,iy)+bmat(iy,iyy)*xxjx(ix,iyy)
        end do
        end do
        end do 
        
        b0=bzold

        mark=1
        mark(:,1+ny/2-nyexp/2:nyexp+ny/2-nyexp/2)=0 

        ax=dx*dfloat(nx);        ay=dy*dfloat(ny)

        xk(0)=0.d0
        do k=0,nx-1; xk(k)=pi*dfloat(k)/ax;  end do

        yk(0)=0.d0
        do k=0,ny-1; yk(k)=pi*dfloat(k)/ay;  end do

22      format(1x,3d14.6)
        bz3d(:,:,iz)=bznew 
         end do ! Z-cycle 
!************************  End of Initial volume distr. ********************        
        end
!***********************************************************
        subroutine sol(bmat)
        use aa
        implicit real*8(a-h,o-z)
        dimension bmat(ny,ny+nx),buf(ny+nx)
        dimension bmat0(ny,ny+nx)
        
        if(mod(nyexp,2)==1) then
        bmat(nyexp,:)=0.d0
        bmat(nyexp,nyexp)=1.d0; bmat(nyexp,nyexp-1)=-1.d0 
        end if 
        
        do i=1,nyexp-1  ! Upper form   creation
          !*****  Main element det.  *****
          smax=0.d0
          jmax=0
          do j=i,nyexp
          s=dabs(bmat(j,i))
          if(s>smax) then
          smax=s
          jmax=j
          end if
          end do
          buf=bmat(i,:)
          bmat(i,:)=bmat(jmax,:)
          bmat(jmax,:)=buf
          !*****  End of Main element det.  *****
          do j=i+1,nyexp
          bmat(j,:)=bmat(j,:)-bmat(i,:)*bmat(j,i)/bmat(i,i)
          end do
        end do     ! End of Upper form creation

          do i=nyexp,2,-1  ! Diag form   creation
          do j=1,i-1
          ccc=bmat(j,i)/bmat(i,i)
          bmat(j,:)=bmat(j,:)-bmat(i,:)*ccc
          end do 
        end do     ! End of Diag form creation
        
        do i=1,nyexp;   buf(i)=bmat(i,i);  end do
          
          do j=1+ny,ny+nx
          bmat(:,j)=bmat(:,j)/buf(1:nyexp)
          end do
          
        end
!************************************************************
        subroutine fstr2  ! Stright 2-D Fourier transform
        use aaa
        implicit real*8(a-h,o-z)
        real*8  xq(nx+1,0:ny)
        ind=1;        xq=0.d0
        do i=1,nx
        fy=xx(i,:)

       call rfty (fy, 1, fyk, 1, 2)

        xq(i,:)=fyk
        end do
        qqnew=0.d0
        do j=0,ny-1
        fx=xq(:,j)
       call rftx (fx, 1, fxk, 1)
        qqnew(:,j)=fxk
        end do
        end
!************************************************************
        subroutine frev2  ! Reverse 2-D Fourier transform
        use aaa
        implicit real*8(a-h,o-z)
        real*8  qx(0:nx,ny+1)
        ind=1;        qx=0.d0
        do j=0,nx-1
        fyk=qq(j,:)
       call rfty (fyk, 1, fynew, 1, 5)
        qx(j,:)=fynew
        end do
        xxnew=0.d0
        do i=1,ny
        fxk=qx(:,i)
       call rftx (fxk, 1, fxnew, 1)
        xxnew(:,i)=fxnew
        end do
        end
********************  RFT X *********************************
      subroutine rftx (Xfft, iX, Yfft, iY)
c     1-D Fast Fourier transForM.
        use aa
        implicit real*8(a-h,o-z)
      dimension xfft(nx+1), yfft(nx+1)
      common /D700Dtx/ n, n2, n4, m, f, rttwo
      common /fWorkx/ W(nx*3+3)

      call rcax (Xfft, iX, Yfft, iY)
      return

      end
****************************************************
      subroutine rcax (Xfft, iX, Yfft, iY)
        use aa
        implicit real*8(a-h,o-z)
      dimension  xfft(nx+1), yfft(nx+1)
      common /d700dtx/ n, n2, n4, m, f, rttWo
      common /fWorkx/ W(nx*3+3)
      data mloc / - 1 /
      m = mmx + 1
      if(m .ne. mloc) call d700sux
      mloc = m;      ni = n4;      no = n4 + n2 + 1
      kX = 1;      nfWa = ni;      nlWa = ni + n2
      do 10 i = nfWa, nlWa
      W(i) = xfft(kX)
   10 kX = kX + iX
      nQ = n2
      do 170 l = 1, m
      nQ2 = nQ / 2;      nQ2m1 = nQ2 - 1;      no1 = no
      no2 = no + n2;      ni1 = ni;      ni2 = ni + nQ
      W(no1) = W(ni1) + W(ni2);      W(no2) = W(ni1) - W(ni2)
      if(nQ2m1) 50, 40, 20
   20 do 30 it = 1, nQ2m1
      nor1 = no1 + it
      nor2 = no1 + nQ2 + it
      nir1 = ni1 + it
      nir2 = ni2 - it
      W(nor1) = W(nir1) + W(nir2)
      W(nor2) = W(nir1) - W(nir2)
   30 continue
   40 nor1 = no1 + nQ2
      nir1 = ni1 + nQ2
      W(nor1) = W(nir1) + W(nir1)
   50 no1 = no1 + nQ
      no2 = no2 - nQ
      ni1 = ni2 + nQ
      ni2 = ni1 + nQ
      if(no1 - no2) 60, 120, 160
   60 nc = 0
      ns = n4
      kc = 0
      ks = n4
   70 continue
      W(no1) = W(ni1) + W(ni2)
      W(no2) = W(ni1) - W(ni2)
      if(nQ2m1) 110, 100, 80
   80 nc = nc + nQ
      ns = ns - nQ
      cc = W(nc)
      ss = W(ns)
      do 90 it = 1, nQ2m1
      nor1 = no1 + it
      nor2 = no2 + it
      nir1 = ni1 + it
      nir2 = ni2 - it
      noi1 = nor1 + nQ2
      noi2 = nor2 + nQ2
      nii1 = nir1 + nQ
      nii2 = nir2 + nQ
      re = cc * W(nir2) - ss * W(nii2)
      ai = ss * W(nir2) + cc * W(nii2)
      W(nor1) = W(nir1) + re
      W(nor2) = W(nir1) - re
      W(noi1) = + W(nii1) - ai
      W(noi2) = - W(nii1) - ai
   90 continue
  100 kc = kc + nQ2
      ks = ks - nQ2
      nor1 = no1 + nQ2
      nor2 = no2 + nQ2
      nir1 = ni1 + nQ2
      nir2 = ni2 + nQ2
      W(nor1) = 2.0 * (W(kc) * W(nir1) - W(ks) * W(nir2))
      W(nor2) = 2.0 * (W(ks) * W(nir1) + W(kc) * W(nir2))
  110 no1 = no1 + nQ
      no2 = no2 - nQ
      ni1 = ni2 + nQ
      ni2 = ni1 + nQ
      if(no1 - no2) 70, 120, 160
  120 W(no1) = W(ni1)
      if(nQ2m1) 160, 150, 130
  130 do 140 it = 1, nQ2m1
      nor1 = no1 + it
      noi1 = nor1 + nQ2
      nir1 = ni + nQ + it
      nir2 = ni + nQ + nQ - it
      W(nor1) = + W(nir1)
      W(noi1) = - W(nir2)
  140 continue
  150 nor1 = no1 + nQ2
      nir1 = ni + nQ + nQ2
      W(nor1) = rttWo * W(nir1)
  160 nt = ni
      ni = no
      no = nt
      nQ = nQ2
  170 continue
      kY = 1;      nfWa = ni;      nlWa = ni + n2
      do 180 i = nfWa, nlWa
      yfft(kY) = W(i) * f
  180 kY = kY + iY
      return
      end
******************************************
********************************************
      subroutine d700sux
        use aa
        implicit real*8(a-h,o-z)
      common /d700dtx/ n, n2, n4, m, f, rttWo
      common /fWorkx/ W(nx*3+3)
      if(m .lt. 2) Go to 90
      n4 = 2 ** (m-2);      n2 = n4 + n4;      n  = n2 + n2
      f  = 1.0d0 / dsQrt(dfloat(n))
      da = 4.0d0 * atan(1.0d0) / dfloat(n2)
      rttWo = sQrt(2.0d0)
      nc = n4 - 1
      if(nc .le. 0) return
      do 10 mc = 1, nc
      W(mc) = cos(da*dfloat(mc))
   10 continue
      return
   90 print 100, m
      stop
  100 format(23H1* error in rft ... m =,i5,20H, sHould be .Ge. 2 *)
      end
********************  rft Y *********************************
      subroutine rfty (Xfft, iX, Yfft, iY, mode)
        use aa
        implicit real*8(a-h,o-z)
      dimension xfft(ny+1), yfft(ny+1)
      common /d700dty/ n, n2, n4, m, f, rttWo
      common /fWorky/ W(ny*3+3)
      Go to (10, 20, 30, 40, 20, 30), mode
   10 call rpay (Xfft, iX, Yfft, iY)
      return
   20 call rcay (Xfft, iX, Yfft, iY)
      return
   30 call rsay (Xfft, iX, Yfft, iY)
      return
   40 call rpsy (Xfft, iX, Yfft, iY)
      return
      end
****************************************************
      subroutine rcay (Xfft, iX, Yfft, iY)
        use aa
        implicit real*8(a-h,o-z)
      dimension  xfft(ny+1), yfft(ny+1)
      common /d700dty/ n, n2, n4, m, f, rttWo
      common /fWorky/ W(ny*3+3)
      data mloc / - 1 /
      m = mmy + 1
      if(m .ne. mloc) call d700suy
      mloc = m;      ni = n4;      no = n4 + n2 + 1
      kX = 1;      nfWa = ni;      nlWa = ni + n2
      do 10 i = nfWa, nlWa
      W(i) = xfft(kX)
   10 kX = kX + iX
      nQ = n2
      do 170 l = 1, m
      nQ2 = nQ / 2;      nQ2m1 = nQ2 - 1;      no1 = no
      no2 = no + n2;      ni1 = ni;      ni2 = ni + nQ
      W(no1) = W(ni1) + W(ni2)
      W(no2) = W(ni1) - W(ni2)
      if(nQ2m1) 50, 40, 20
   20 do 30 it = 1, nQ2m1
      nor1 = no1 + it;      nor2 = no1 + nQ2 + it
      nir1 = ni1 + it;      nir2 = ni2 - it
      W(nor1) = W(nir1) + W(nir2)
      W(nor2) = W(nir1) - W(nir2)
   30 continue
   40 nor1 = no1 + nQ2
      nir1 = ni1 + nQ2
      W(nor1) = W(nir1) + W(nir1)
   50 no1 = no1 + nQ
      no2 = no2 - nQ;      ni1 = ni2 + nQ;      ni2 = ni1 + nQ
      if(no1 - no2) 60, 120, 160
   60 nc = 0
      ns = n4;      kc = 0;      ks = n4
   70 continue
      W(no1) = W(ni1) + W(ni2)
      W(no2) = W(ni1) - W(ni2)
      if(nQ2m1) 110, 100, 80
   80 nc = nc + nQ
      ns = ns - nQ;      cc = W(nc);      ss = W(ns)
      do 90 it = 1, nQ2m1
      nor1 = no1 + it;      nor2 = no2 + it
      nir1 = ni1 + it;      nir2 = ni2 - it
      noi1 = nor1 + nQ2;      noi2 = nor2 + nQ2
      nii1 = nir1 + nQ;      nii2 = nir2 + nQ
      re = cc * W(nir2) - ss * W(nii2)
      ai = ss * W(nir2) + cc * W(nii2)
      W(nor1) = W(nir1) + re;      W(nor2) = W(nir1) - re
      W(noi1) = + W(nii1) - ai;      W(noi2) = - W(nii1) - ai
   90 continue
  100 kc = kc + nQ2
      ks = ks - nQ2;      nor1 = no1 + nQ2;      nor2 = no2 + nQ2
      nir1 = ni1 + nQ2;      nir2 = ni2 + nQ2
      W(nor1) = 2.0 * (W(kc) * W(nir1) - W(ks) * W(nir2))
      W(nor2) = 2.0 * (W(ks) * W(nir1) + W(kc) * W(nir2))
  110 no1 = no1 + nQ
      no2 = no2 - nQ;      ni1 = ni2 + nQ;      ni2 = ni1 + nQ
      if(no1 - no2) 70, 120, 160
  120 W(no1) = W(ni1)
      if(nQ2m1) 160, 150, 130
  130 do 140 it = 1, nQ2m1
      nor1 = no1 + it;      noi1 = nor1 + nQ2
      nir1 = ni + nQ + it;      nir2 = ni + nQ + nQ - it
      W(nor1) = + W(nir1);      W(noi1) = - W(nir2)
  140 continue
  150 nor1 = no1 + nQ2
      nir1 = ni + nQ + nQ2
      W(nor1) = rttWo * W(nir1)
  160 nt = ni
      ni = no;      no = nt;      nQ = nQ2
  170 continue
      kY = 1
      nfWa = ni
      nlWa = ni + n2
      do 180 i = nfWa, nlWa
      yfft(kY) = W(i) * f
  180 kY = kY + iY
      return
      end
******************************************
      subroutine rpay (Xfft, iX, Yfft, iY)
        use aa
        implicit real*8(a-h,o-z)
      dimension  xfft(ny+1), yfft(ny+1)
      common /d700dty/ n, n2, n4, m, f, rttWo
      common /fWorky/ W(ny*3+3)
      data mloc / - 1 /
      m = mmy
      if(m .ne. mloc) call d700suy
      mloc = m;      ni = n4 - 1;      no = ni + n
      kX = 1;      nfWa = ni + 1;      nlWa = ni + n
      do 10 i = nfWa, nlWa
      W(i) = xfft(kX)
   10 kX = kX + iX
      nQ = n2
      do 80 l = 1, m
      no1 = no;      no2 = no + n2
      ni1 = ni;      ni2 = ni + nQ
      do 20 it = 1, nQ
      nor1 = no1 + it;      nor2 = no2 + it
      nir1 = ni1 + it;      nir2 = ni2 + it
      W(nor1) = W(nir1) + W(nir2)
      W(nor2) = W(nir1) - W(nir2)
   20 continue
      nc = 0;      ns = n4;      no1 = no1 + nQ
      no2 = no2 - nQ;      ni1 = ni2 + nQ;      ni2 = ni1 + nQ
      if(no1 - no2) 30, 50, 70
   30 nc = nc + nQ
      ns = ns - nQ;      cc = W(nc);      ss = W(ns)
      do 40 it = 1, nQ
      nor1 = no1 + it;      noi1 = nor1 + n2
      nor2 = no2 + it;      noi2 = nor2 + n2
      nir1 = ni1 + it;      nii1 = nir1 + n2
      nir2 = ni2 + it;      nii2 = nir2 + n2
      re = cc * W(nir2) - ss * W(nii2);
      ai = ss * W(nir2) + cc * W(nii2)
      W(nor1) = W(nir1) + re;      W(nor2) = W(nir1) - re
      W(noi1) = + W(nii1) + ai;    W(noi2) = - W(nii1) + ai
   40 continue
      no1 = no1 + nQ;      no2 = no2 - nQ
      ni1 = ni2 + nQ;      ni2 = ni1 + nQ
      if(no1 - no2) 30, 50, 70
   50 do 60 it = 1, nQ
      nor1 = no1 + it;      noi1 = nor1 + n2
      nir1 = ni1 + it;      nir2 = ni2 + it
      W(nor1) = W(nir1);      W(noi1) = W(nir2)
   60 continue
   70 nt = ni
      ni = no;      no = nt;      nQ = nQ / 2
   80 continue
      kY = 1;      nfWa = ni + 1;      nlWa = ni + n
      do 90 i = nfWa, nlWa
      yfft(kY) = W(i) * f
   90 kY = kY + iY
      yfft(kY) = 0.0
      return
      end
*****************************************************
      subroutine rpsy (Xfft, iX, Yfft, iY)
        use aa
        implicit real*8(a-h,o-z)
      dimension  xfft(ny+1), yfft(ny+1)
      common /d700dty/ n, n2, n4, m, f, rttWo
      common /fWorky/ W(ny*3+3)
      data mloc / - 1 /
      m = mmy
      if(m .ne. mloc) call d700suy
      mloc = m;      ni = n4 - 1;      no = ni + n
      kX = 1;      nfWa = ni + 1;      nlWa = ni + n
      do 10 i = nfWa, nlWa
      W(i) = xfft(kX)
   10 kX = kX + iX
      nQ = 1
      do 80 l = 1, m
      no1 = no;      no2 = no + nQ
      ni1 = ni;      ni2 = ni + n2
      do 20 it = 1, nQ
      nor1 = no1 + it;      nor2 = no2 + it
      nir1 = ni1 + it;      nir2 = ni2 + it
      W(nor1) = W(nir1) + W(nir2)
      W(nor2) = W(nir1) - W(nir2)
   20 continue
      nc = 0;      ns = n4;      no1 = no2 + nQ
      no2 = no1 + nQ;      ni1 = ni1 + nQ;      ni2 = ni2 - nQ
      if(ni1 - ni2) 30, 50, 70
   30 nc = nc + nQ
      ns = ns - nQ;      cc = W(nc);      ss = W(ns)
      do 40 it = 1, nQ
      nor1 = no1 + it;      noi1 = nor1 + n2
      nor2 = no2 + it;      noi2 = nor2 + n2
      nir1 = ni1 + it;      nii1 = nir1 + n2
      nir2 = ni2 + it;      nii2 = nir2 + n2
      W(nor1) = W(nir1) + W(nir2)
      re      = W(nir1) - W(nir2)
      W(noi1) = W(nii1) - W(nii2)
      ai      = W(nii1) + W(nii2)
      W(nor2) = cc * re + ss * ai
      W(noi2) = cc * ai - ss * re
   40 continue
      no1 = no2 + nQ;      no2 = no1 + nQ
      ni1 = ni1 + nQ;      ni2 = ni2 - nQ
      if(ni1 - ni2) 30, 50, 70
   50 do 60 it = 1, nQ
      nor1 = no1 + it;      nor2 = no2 + it
      nir1 = ni1 + it;      nii1 = nir1 + n2
      W(nor1) = W(nir1) + W(nir1)
      W(nor2) = W(nii1) + W(nii1)
   60 continue
   70 nt = ni
      ni = no;      no = nt;      nQ = nQ + nQ
   80 continue
      kY = 1;      nfWa = ni + 1;      nlWa = ni + n
      do 90 i = nfWa, nlWa
      yfft(kY) = W(i) * f
   90 kY = kY + iY
      yfft(kY) = yfft(1)
      return
      end
***************************************************
      subroutine rsay (Xfft, iX, Yfft, iY)
        use aa
        implicit real*8(a-h,o-z)
      dimension  xfft(ny+1), yfft(ny+1)
      common /d700dty/ n, n2, n4, m, f, rttWo
      common /fWorky/ W(ny*3+3)
      data mloc / - 1 /
      m = mmy + 1
      if(m .ne. mloc) call d700suy
      mloc = m;      ni = n4;      no = n4 + n2 + 1
      kX = iX + 1;      nfWa = ni + 1;      nlWa = ni + n2 - 1
      do 10 i = nfWa, nlWa
      W(i) = xfft(kX)
   10 kX = kX + iX
      nQ = n2
      do 170 l = 1, m
      nQ2 = nQ / 2;      nQ2m1 = nQ2 - 1;      no1 = no
      no2 = no + n2;      ni1 = ni;      ni2 = ni + nQ
      if(nQ2m1) 50, 40, 20
   20 do 30 it = 1, nQ2m1
      nor1 = no1 + it;      nor2 = no1 + nQ2 + it
      nir1 = ni1 + it;      nir2 = ni2 - it
      W(nor1) = W(nir1) - W(nir2)
      W(nor2) = W(nir1) + W(nir2)
   30 continue
   40 nor1 = no1 + nQ2
      nir1 = ni1 + nQ2
      W(nor1) = W(nir1) + W(nir1)
   50 no1 = no1 + nQ
      no2 = no2 - nQ;      ni1 = ni2 + nQ;      ni2 = ni1 + nQ
      if(no1 - no2) 60, 120, 160
   60 nc = 0
      ns = n4;      kc = 0;      ks = n4
   70 continue
      W(no1) = + W(ni1) + W(ni2)
      W(no2) = - W(ni1) + W(ni2)
      if(nQ2m1) 110, 100, 80
   80 nc = nc + nQ
      ns = ns - nQ;      cc = W(nc);      ss = W(ns)
      do 90 it = 1, nQ2m1
      noi1 = no1 + it;      noi2 = no2 + it
      nii1 = ni1 + it;      nii2 = ni2 - it
      nor1 = noi1 + nQ2;      nor2 = noi2 + nQ2
      nir1 = nii1 + nQ;      nir2 = nii2 + nQ
      re = cc * W(nir2) - ss * W(nii2)
      ai = ss * W(nir2) + cc * W(nii2)
      W(nor1) = W(nir1) - re;      W(nor2) = W(nir1) + re
      W(noi1) = + W(nii1) + ai;      W(noi2) = - W(nii1) + ai
   90 continue
  100 kc = kc + nQ2
      ks = ks - nQ2
      nor1 = no1 + nQ2;      nor2 = no2 + nQ2
      nir1 = ni1 + nQ2;      nir2 = ni2 + nQ2
      W(nor1) = 2.0 * (+ W(kc) * W(nir1) + W(ks) * W(nir2))
      W(nor2) = 2.0 * (- W(ks) * W(nir1) + W(kc) * W(nir2))
  110 no1 = no1 + nQ
      no2 = no2 - nQ;      ni1 = ni2 + nQ;      ni2 = ni1 + nQ
      if(no1 - no2) 70, 120, 160
  120 nir1 = ni + nQ
      W(no1) = W(nir1)
      if(nQ2m1) 160, 150, 130
  130 do 140 it = 1, nQ2m1
      noi1 = no1 + it
      nor1 = noi1 + nQ2
      nir1 = ni + nQ + it
      nir2 = ni + nQ + nQ - it
      W(nor1) = W(nir1)
      W(noi1) = W(nir2)
  140 continue
  150 nor1 = no1 + nQ2
      nir1 = ni + nQ + nQ2
      W(nor1) = rttWo * W(nir1)
  160 nt = ni
      ni = no;      no = nt;      nQ = nQ2
  170 continue
      yfft(1) = xfft(1)
      kY = iY + 1;      nfWa = ni + 1;      nlWa = ni + n2 - 1
      do 180 i = nfWa, nlWa
      yfft(kY) = W(i) * f
  180 kY = kY + iY
      yfft(kY) = xfft(kX)
      return
      end
********************************************
      subroutine d700suy
        use aa
        implicit real*8(a-h,o-z)
      common /d700dty/ n, n2, n4, m, f, rttWo
      common /fWorky/ W(ny*3+3)
      if(m .lt. 2) Go to 90
      n4 = 2 ** (m-2);      n2 = n4 + n4
      n  = n2 + n2
      f  = 1.0d0 / dsQrt(dfloat(n))
      da = 4.0d0 * atan(1.0d0) / dfloat(n2)
      rttWo = sQrt(2.0d0)
      nc = n4 - 1
      if(nc .le. 0) return
      do 10 mc = 1, nc
      W(mc) = cos(da*dfloat(mc))
   10 continue
      return
   90 print 100, m
      stop
  100 format(23H1* error in rft ... m =,i5,20H, sHould be .Ge. 2 *)
      end
!*******************************************
        subroutine errdet (err)
        use aaa
        implicit real*8 (a-h,o-z)
        s0=0.d0
        s1=0.d0
        do i=1,nxexp
        do j=1,nyexp
        s0=s0+(bztes(i+nx/2-nxexp/2,j+ny/2-nyexp/2)-
     *   bzold(i+nx/2-nxexp/2,j+ny/2-nyexp/2))**2
        s1=s1+(bztes(i+nx/2-nxexp/2,j+ny/2-nyexp/2))**2+
     *   (bzold(i+nx/2-nxexp/2,j+ny/2-nyexp/2))**2
        end do
        end do
        err=(s0/s1)
        end

