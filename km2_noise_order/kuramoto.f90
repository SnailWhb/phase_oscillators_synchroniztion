parameter (m=10000,n=m+1)
complex cj,r0,r
parameter(cj=(0.,1.))
real min,max,uper,ran  !自然频率w
real::r1,r2 
real::k0,k,a
real::dY(n),w(n),phase(n),phase0(n),t,sumw(n),sumww(n),ddy(n),u(n-1),noise(n-1),D
real::pi,tend
integer nn,nnn,II,summ
CHARACTER(LEN=80) ::FILENAME0, FILENAME1,FILENAME2,FORM
pi=4*atan(1.0)
open(2,file='ordermod.dat',status='unknown')
open(10,file='phase_begin.dat',status='unknown')
open(30,file='W_N_10000.txt') 
!open(10,file='phase0.txt') 
call random_seed()
!open(30,file='W3.txt') 
!!!!!!初值条件!!!!!!!!
do i=1,m
read(10,*) phase(i)
enddo
!enddo
!!!!!!自然频率w!!!!!!!!!!!!!!!!!
!do i=1,m
 !  if(i<m/2) then
!	 call random_number(ran)
  !   w(i)=tan(2*pi*ran)-3.
 !    do while(w(i)<-10.or.w(i)>10) 
!	   call  random_number(ran)
!       w(i)=tan(2*pi*ran)-3.
!	 end do
!  else
!	   call   random_number(ran)
 !      w(i)=tan(2*pi*ran)+3.
  !    do while(w(i)<-10.or.w(i)>10) 
!		  call  random_number(ran)
 !         w(i)=tan(2*pi*ran)+3.
!	   end do
 !  end if
!write(3,*) w(i)
!enddo
!双峰分布的随机数
do i=1,m
read(30,*) w(i) !读取服从单峰分布的w
enddo
!!!!!!!!!!!!!噪声!!!!!!!!!!!!!!!!
D=1.2
!!!!!!!!!!!!!!!! 步长h !!!!!!!!!!
h=0.01
!运行的时间
nn=0
II=0 !文件名
t=2000
tend=1800
k0=1.25
do k=3.75,0.75,-0.15
   II=II+1
  write(*,*) II,k0,k  !写屏
   !!!!每次创建一个文件名!!!!!!!!!
   SELECT CASE (II)  
    CASE (1:9)
      WRITE(FORM,'(I1)') II   
    CASE (10:99)
      WRITE(FORM,'(I2)') II
    CASE (100:999)
      WRITE(FORM,'(I3)') II
    END SELECT
    WRITE(FILENAME1,*) "phase_end_",TRIM(FORM),".dat" !创建记录最后一时刻的相位的文件名
    WRITE(FILENAME2,*) "youxiaopinglv_",TRIM(FORM),".dat"   !创建记录有效频率的文件名
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  !!!!!!!!!!!!!!!!!!!!!初始条件!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   phase(n)=0.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   nnn=0
   !!!!!!!!!!!!!有效频率归零!!!!!!!!!!!
	do i=1,m
	sumw(i)=0.0
	enddo
	!!!!!!!!!!!!!!!!!!!!!!!!序参量!!!!!!!!!!!!!!!!!!!!!!
    r1=0.0
    r2=0.0
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do while(abs(phase(n))<=t) !计算时间的长度
	 nn=nn+1
	 if(mod(nn,4000)==1) write(*,*) phase(n)
	   call rungekutta4(h,k0,k,w,phase,ddy,n,noise,D)
	         do i=1,m
					do while(phase(i)>2*pi)
						phase(i)=phase(i)-2*pi
					enddo
					do while(phase(i)<=0) !使相位[0,2*pi]
					 phase(i)=phase(i)+2*pi
					enddo
			  enddo 
	!!!!!!!!!!!!!抛去暂态!!!!!!!!!!!!!!!!!!!!!
     if(tend<phase(n)) then
	!!!!!!!!!!!!!!!!!!!!!!!!!计算序参量r!!!!!!!!!
	   nnn=nnn+1
       r0=(0.,0.)
	   r=(0.,0.)
       do j=1,m
		r0=r0+exp(cj*phase(j))
        r=r+exp(2*cj*phase(j))
       enddo
	   r1=r1+abs(r0/m)
	   r2=r2+abs(r/m)  
     !!!!!!!!!!!!!!!!!!!!!!!!!!计算有效频率!!!!!!!!!!!!!!!!!!
	   sumw=sumw+ddy     		                  
     endif	 	  	 		 	  	    	    	  	  	  	  	  	  	  	     
enddo

write(2,*) k0,k,r1/nnn,r2/nnn  !记录序参量r1,r2，耦合强度k2

open(99,file=FILENAME1) 
open(100,file=FILENAME2) 
do i=1,m
write(99,*)  phase(i) !记录每个耦合强度下的最后一时刻相位
write(100,*) w(i),sumw(i)/nnn !记录有效频率与自然频率
enddo 
close(99)
close(100)

enddo
stop  
end



subroutine interation(k0,k,w,y,dy,n,noise,D)
real::k0,k
real::  y(n),dy(n),D,noise(n-1)
complex:: ccc,ouhe,ouhe1
real:: cc,w(n),cc1,a,b
integer:: i,n
ccc=(0,1)
ouhe=0.
ouhe1=0.
        do i=1,n-1
	       ouhe=ouhe+exp(ccc*y(i))
           ouhe1=ouhe1+exp(2.*ccc*y(i))
       enddo
	    do i=1,n-1
           a=(k0/((n-1)*1.))*aimag(ouhe*exp(-ccc*y(i)))
	       b=(k/((n-1)*1.))*aimag(ouhe1*exp(-2*ccc*y(i)))
	       dy(i)=w(i)+a+b+sqrt(D)*noise(i)
	   enddo
	   dy(n)=1.
return 
end 
!!!!!!!!!!!!!!加噪声函数!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine niosefun(n,noise)                                         !!子程序：生成高斯白噪声   
implicit none
integer n,i
real u1,u2,pi,noise(n)
call random_seed()
pi=3.141592654
do i=1,n                                                                  !!对于特定时刻，每个振子的噪声不同
  call random_number(u1)
  call random_number(u2)
  noise(i)=sqrt(-2*log(u1))*cos(2*pi*u2)
enddo
return
end




!!!!!!!!!!!!!!!!rungekutta4格式!!!!!!!!!!!!!!!!!!
subroutine  rungekutta4(h,k0,k,w,Y,dyy,n,noise,D)
	integer:: n !n(方程数)
	real:: k0,k,h
	real::w(n),Y(n),dY(n),k1(n),k2(n),k3(n),k4(n),YY(n),dyy(n),noise(n-1),D  
	  YY=Y
	 call  niosefun(n-1,noise) 
	 call interation(k0,k,w,y,dy,n,noise,D)
      k1=dY*h
	  YY=Y+k1/2
     call interation(k0,k,w,y,dy,n,noise,D)
      k2=dY*h
	  YY=Y+k2/2
	 call interation(k0,k,w,y,dy,n,noise,D)
      k3=dY*h
	  YY=Y+k3
	 call interation(k0,k,w,y,dy,n,noise,D)
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
