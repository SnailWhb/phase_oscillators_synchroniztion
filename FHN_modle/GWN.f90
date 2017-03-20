program main                                               !!主程序开始
implicit none
integer i,n,nn,k,m,fk
parameter (n=500)    
character(len=3)cTemp                   
real y(1001),dy(1001),ran,rr,pi,noise(n),D
  nn=2*n+1
  pi=3.141592654
  rr=2.0
  k=0
  m=0  
  do fk=1,20
  write(cTemp,'(i3)')fk
  open(1,file='IC'//trim(adjustl(cTemp))//'.dat', status='unknown')                                             !存初始条件
  open(2,file='Oscillator'//trim(adjustl(cTemp))//'.dat', status='unknown')                                     !存一个点的时间演化图
  open(3,file='1stResult'//trim(adjustl(cTemp))//'.dat', status='unknown')                                      !存指定时刻结果u-i图
  open(4,file='2ndResult'//trim(adjustl(cTemp))//'.dat', status='unknown')                                      !存指定时刻结果u-i图
  open(5,file='3rdResult'//trim(adjustl(cTemp))//'.dat', status='unknown')                                      !存指定时刻结果u-i图
  open(0,file='i_time_phase_'//trim(adjustl(cTemp))//'.dat', status='unknown')                                                      
   call random_seed()
   do i=1,n                                                 !!赋初值    初始条件:u^2+v^2=4  环上随机取值
         call random_number(ran)
         y(i)=rr*sin(pi*ran)
		 y(i+n)=rr*cos(pi*ran)   
   enddo
   y(nn)=0
   write(1,*)y          !存初始条件
    D=0.00005*(fk-1)     !噪声强度的变化
	m=0  
	write(*,*)D,fk                                            
    do while(y(nn)<500) 
	  m=m+1                                         !!截止时间
	     call RK(nn,y,dy,n,noise,D)
	       if(mod(m,10)==1)then
		   if(y(nn)>480)then   
	         do i=1,n                                             !!存储指定时刻的所需值 
	           write(0,*)i,y(nn),y(i)
	         enddo
	        endif
             write(2,*)y(nn),y(1),y(250),y(350),y(500)           !存振子随时间的演化关系图
		    endif
            if(m==47000)then   
	       do i=1,n                                                    !存储指定时刻的所有振子的u值
	         write(3,*)i,y(i)
	       enddo
			 close(3)
		 endif

		  if(m==48000)then
	         do i=1,n                                                 !存储指定时刻的所有振子的u值
	           write(4,*)i,y(i)
	         enddo
			 close(4)
		  endif

		  if(m==49000)then   
	         do i=1,n                                                 !存储指定时刻的所有振子的u值
	           write(5,*)i,y(i)
	         enddo
			 close(5)
	      endif

      enddo   !龙格库塔时间演化结束
    enddo !fk循环结束
end                                                                !!主程序结束
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine  fun(nn,y,dy,n,noise,D)                                 !!子程序：编写演化方程    注意噪声项对于每个时间步长的不同
implicit none
integer i,j,k1,n,nn        
real y(nn),dy(nn),e,thita,rn,buu,bvv,buv,bvu,phase,pi,D,noise(n),a,z(6*n),r
real sum1(nn-1),temp(nn-1),t1(3*n),t2(3*n)
nn=2*n+1
pi=3.141592654
phase=pi/2-0.1
a=1.001
thita=0.4
r=0.12
rn=r*n
!write(*,*)rn
!pause
e=0.05  
buu=cos(phase)
bvv=cos(phase)
bvu=-sin(phase)      
buv=sin(phase)
do i=1,n                                                            !!非局域耦合项的编写
        sum1(i)=0
		sum1(i+n)=0
		temp(i)=0
		temp(i+n)=0                                        
        t1(i)=y(i)
		t2(i)=y(i+n)
		t1(i+n)=t1(i)
		t2(i+n)=t2(i)
		t1(i+2*n)=t1(i)
		t2(i+2*n)=t2(i)
enddo
do i=1,n
		do j=i-rn+n,i+rn+n
		   sum1(i)=sum1(i)+buu*(t1(j)-t1(i))+buv*(t2(j)-t2(i))
	       sum1(i+n)=sum1(i+n)+bvu*(t1(j)-t1(i))+bvv*(t2(j)-t2(i))
		enddo
	temp(i)=sum1(i)*thita/(2*rn)
	temp(i+n)=sum1(i+n)*thita/(2*rn)	
    dy(i)=(y(i)-(y(i)**3)/3-y(i+n)+temp(i))/e
    dy(i+n)=(y(i)+a+temp(i+n)+sqrt(2.*D)*noise(i))
enddo
dy(nn)=1.0d0
return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine  RK(nn,y,dy,n,noise,D)                                         !!子程序：龙格库塔法 
   implicit none
   integer n,nn
   real h
   real Y(nn),k1(nn),k2(nn),k3(nn),k4(nn),YY(nn),dY(nn),noise(n),D
	  YY=Y
	  h=0.01
	  call  niosefun(n,noise) 
	  call fun(nn,YY,dy,n,noise,D)
	    k1=dY*h
	    YY=Y+k1/2.
	  call fun(nn,YY,dy,n,noise,D)
	    k2=dY*h
	    YY=Y+k2/2.
	  call fun(nn,YY,dy,n,noise,D)
	    k3=dY*h
	    YY=Y+k3
	  call fun(nn,YY,dy,n,noise,D)
	    k4=dY*h
		Y=Y+(k1+2*k2+2*k3+k4)/6
		return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!