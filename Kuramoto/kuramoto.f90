parameter (m=10000,n=m+1)
complex cj,r,r1,r01,r02,r001,r002,realr1,realr2
parameter(cj=(0.,1.))
real min,max,uper,ran  !自然频率w
real::nn,dist ,dis,pppp,pppp1,pppp2
real::k0,k,a
real::dY(n),w(n),phase(n),phase0(n),t,sumw(n),sumww(n),ddy(n),u(n-1)
real::pi,tend
integer::mm0,mm,mmm,jj,pp,ppp,ppp1,ppp2,mmmm
pi=4*atan(1.0)
call random_seed()  
open(1000,file='n_time_phase.dat',status='unknown')
open(1,file='phasetime.dat',status='unknown')
open(11,file='pp_phasetime.dat',status='unknown')
open(111,file='ppp_phasetime.dat',status='unknown')
open(2,file='order0_mod.dat',status='unknown') !记录自然频率大于0的序参量
open(21,file='order1_mod.dat',status='unknown')
open(3,file='order0_fujiao.dat',status='unknown')
open(31,file='order1_fujiao.dat',status='unknown')
open(222,file='realr_mod.dat',status='unknown')
open(333,file='realr_fujiao.dat',status='unknown')
open(4,file='youxiaow.dat',status='unknown')
open(5,file='phaseend.dat',status='unknown')
open(10,file='phase_begin.dat',status='unknown')
open(30,file='W_N_10000.txt') 
open(500,file='ouheqiaodu.dat')
!open(30,file='W3.txt') 
!!!!!!初值!!!!!!!!
!do i=1,m
!call random_number(ran)
!phase(i)=ran
!enddo
!enddo
!do i=1,m
!call random_number(ran)
!phase(i)=ran*2*pi
!write(10,*) phase(i)
!enddo 
do i=1,m
 call random_number(ran)
 if( ran>0) then
  phase(i)=0
  else
  phase(i)=pi
 endif
  write(10,*) phase(i)
enddo
phase(n)=0
!!!!!!自然频率w!!!!!!!!!!!!!!!!!
!双峰分布的随机数
do i=1,m
read(30,*) w(i) !读取服从双峰分布的w
enddo
!!!!!!!!!!!!!耦合强度K!!!!!!!!!!!!!!!!
k0=0.75
k=3
write(500,*) k0,k
!!!!!!!!!!!!!!!! 步长h !!!!!!!!!!
h=0.01
!!!!参数!!!!!
nn=0.0
nnn=0
mmmm=0
!运行的时间
t=2000
tend=1800
!!!!!!计算相位!!!!!!!!!!!!
do while(abs(phase(n))<=t)
 !计算时间的长度
      nnn=nnn+1
     if(mod(nnn,80)==1) write(*,*)  phase(n) ! 写屏
    !  write(*,*)  phase(n),phase(2)  
	  call rungekutta4(h,k0,k,w,phase,ddy,n) 
	         do i=1,m
					do while(phase(i)>2*pi)
						phase(i)=phase(i)-2*pi
					enddo
					do while(phase(i)<=0) !使相位[0,2*pi]
					 phase(i)=phase(i)+2*pi
					enddo
			  enddo 
if(tend<phase(n)) then
	!!!!!!!!!!!!!!!!!!!!!!!!!计算序参量r!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   realr1=(0.,0.)
	   realr2=(0.,0.)
       r=(0.,0.)
	   r1=(0.,0.)
	   r01=(0.,0.)
	   r02=(0.,0.)
	   r001=(0.,0.)
	   r002=(0.,0.)
	   mm=0
	   mmm=0
	   mm0=0
     do j=1,m
        	realr1=realr1+exp(cj*phase(j))
            realr2=realr2+exp(2*cj*phase(j)) !计算自然频率大于0的振子序参量
	   if(w(j)>0) then
        mmm=mmm+1
		r=r+exp(cj*phase(j))
        r1=r1+exp(2*cj*phase(j)) !计算自然频率大于0的振子序参量
	   else
        mm=mm+1
		r01=r01+exp(cj*phase(j))
        r02=r02+exp(2*cj*phase(j)) !计算自然频率小于0的振子序参量
	   endif
   enddo
!!!!!!!!!!!!!!!自然频率记录序参量r,r1的模!!!!!!!!!!!
         call pha(1,realr1/m,theta1)
	     call pha(1,realr2/m,theta)   
         write(222,*) phase(n),abs(realr1/m),abs(realr2/m)   !自然频率大于0记录序参量r,r1的模
         write(333,*) phase(n),theta,theta1          !记录序参量的俯角    
	 !!!!!!!!!!!自然频率大于0记录序参量r,r1的模!!!!!!!!!!!
	   call pha(1,r1/mmm,theta1)
	   call pha(1,r/mmm,theta)
	   write(2,*) phase(n),abs(r/mmm),abs(r1/mmm) !自然频率大于0记录序参量r,r1的模
       write(3,*) phase(n),theta,theta1          !记录序参量的俯角
       !!!!!!!!!!!自然频率小于0记录序参量r,r1的模!!!!!!!!!!!
       call pha(1,r01/mm,theta)
       call pha(1,r02/mm,theta1)
       write(31,*) phase(n),theta,0.5*theta1,0.5*theta1-theta  !记录序参量的俯角
	   write(21,*) phase(n),abs(r01/mm),abs(r02/mm) !记录序参量r,r1的模	
!!!!!!!!!!!!!!!!!!!!!记录每个振子的相位分布随时间演化的情况!!!!!!!!!
      do jj=0,19
            ppp=0
	        pppp=0.
			ppp1=0
	        pppp1=0.
			ppp2=0
			pppp2=0. 
            dis=(2*pi/20.)*jj
            dist=(2*pi/20.)*(jj+1)
	      do i=1,n-1
	         if (phase(i)>=dis.and.phase(i)<dist) then
	          ppp=ppp+1 
	         endif
	      enddo
            if (ppp.eq.0) then
	         write(1,'(F8.3,F8.2,F12.8)')  dis, phase(n),0
	         else
	         pppp=(ppp*1.)/(m*1.0)
	        write(1,'(F8.3,F8.2,F12.8)') dis,phase(n),pppp
	        endif       
        !!!分别统计自然频率大于0；小于0的相振子分布。     
	      do i=1,n-1
            if(phase(i)>=dis.and.phase(i)<dist) then
			   if(w(i)>0) then
	           ppp1=ppp1+1 
			   else
               ppp2=ppp2+1 
			   endif
			endif  
          enddo
		  !! 自然频率大于0的分布
          if (ppp1.eq.0) then
	        write(11,'(F8.3,F8.2,F12.8)')  dis, phase(n),0
	      else
	         pppp1=(ppp1*1.)/(mmm*1.0)
	       write(11,'(F8.3,F8.2,F12.8)')  dis,phase(n),pppp1
	     endif
		 !自然频率小于0的相位分布
         if (ppp2.eq.0) then
	        write(111,'(F8.3,F8.2,F12.8)')  dis, phase(n),0
	     else
	         pppp2=(ppp2*1.)/(mm*1.0)
	        write(111,'(F8.3,F8.2,F12.8)') dis,phase(n),pppp2
	     endif
     enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!计算有效频率!!!!!!!!!!!!!!!!!!
	    summ=summ+1
		sumw=sumw+ddy
!!!!!!!!!!!!!!!!记录每个振子的相位随时间演化情况!!!!!!!!!!!!!!!!!!!!!!!!
  endif	
    	 	  	    	    	  	  	  	  	  	  	  	     
enddo
!!!计算有效频率
do i=1,m
write(5,*) phase(i)  !最后一时刻的相位
write(4,*) w(i),sumw(i)/summ
enddo 
stop  
end


!!!!!!!!!!!!!interation函数定义!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interation(k0,k,w,y,dy,n)
real::k0,k
real::  y(n),dy(n)
complex:: ccc,ouhe,ouhe1
real:: cc,w(n),cc1,a,b
integer:: i,n
ccc=(0,1)
dy(n)=1.
ouhe=0.
ouhe1=0.
        do i=1,n-1
	       ouhe=ouhe+exp(ccc*y(i))
           ouhe1=ouhe1+exp(3.*ccc*y(i))
       enddo
	    do i=1,n-1
           a=(k0/((n-1)*1.))*aimag(ouhe*exp(-ccc*y(i)))
	       b=(k/((n-1)*1.))*aimag(ouhe1*exp(-3*ccc*y(i)))
	       dy(i)=w(i)+a+b
	   enddo
return 
end 

!!!!!!!!!!!!!!!!rungekutta4格式!!!!!!!!!!!!!!!!!!
subroutine  rungekutta4(h,k0,k,w,Y,dyy,n)
	integer:: n !n(方程数)
	real:: k0,k,h
	real::w(n),Y(n),dY(n),k1(n),k2(n),k3(n),k4(n),YY(n),dyy(n)   
	  YY=Y
	 call interation(k0,k,w,YY,dY,n)
      k1=dY*h
	  YY=Y+k1/2
     call interation(k0,k,w,YY,dY,n)
      k2=dY*h
	  YY=Y+k2/2
	 call  interation(k0,k,w,YY,dY,n)
      k3=dY*h
	  YY=Y+k3
	 call interation(k0,k,w,YY,dY,n)
	  k4=dY*h
	  Y=Y+(k1+2*k2+2*k3+k4)/6.
      dyy=(k1+2.*k2+2*k3+k4)/(h*6)
	 return
end 


SUBROUTINE pha(n,y,theta)!计算俯角
	integer n
	real theta(n),PI
	complex y(n)
	PI=atan(1.)*4
	do i=1,n
	  if(real(y(i))<0) then
	      theta(i)=atan(imag(y(i))/real(y(i)))+PI
	  endif
	  if(real(y(i))>0) then
	      if(imag(y(i))>=0) then
	       theta(i)=atan(imag(y(i))/real(y(i)))
	      endif
	      if(imag(y(i))<0) then
	       theta(i)=atan(imag(y(i))/real(y(i)))+2*pi
	      endif		  
	  endif
	  if(real(y(i))==0) then
	      if(imag(y(i))>0)  theta(i)=pi/2
	      if(imag(y(i))<0)  theta(i)=3*pi/2
	  endif
	enddo
	return
	end
