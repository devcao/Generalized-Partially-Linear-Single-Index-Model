
!Compile: gfortran -O2 Fortran_code.f90 ~/gfortransubs/* -o simu1
module commondata
  implicit none
  save
  integer, parameter :: simus=1000, n=1000, m=15, iseed=99,ng=20
  integer, parameter :: lb1=2, lb2=4, lw1=lb1*(lb1*3+14),lw2=lb2*(lb2*3+14),lw3=2*(2*3+14)
  double precision, dimension(n) :: S, W
  double precision, dimension(4,n) :: Z
  integer, dimension(n) :: Y
  double precision, dimension(m) :: epsdis, p
  double precision, dimension(n) :: weight
  double precision, dimension(ng,2):: Gng
  double precision, dimension(n,2):: G
  double precision, dimension(lb1) :: para 
  double precision, dimension(lb2) :: para2
  double precision,dimension(ng)::Ung
  double precision :: U0
end module commondata

program main
  use commondata
  double precision, dimension(lb1) :: alpha
  double precision, dimension(lb2)::gamma
  double precision, dimension(lb1) :: vfd1,theta1,theta1_tmp,zeta
  double precision, dimension(lb2) :: vfd2,theta2,sec
  double precision, dimension(2) :: ms_eps
  double precision, dimension(lw1) :: wv1
  double precision, dimension(lw2) :: wv2
   double precision, dimension(n,2):: Gs1,Gs2,G0
  integer :: i,k,iter, eflag1,eflag2,ccount
  double precision :: beta, tmp, d, tol
  double precision, dimension(lb2,lb2) :: Bn,An,Cn,Evar,An1
  integer, dimension(simus,lb2):: cover
  double precision :: sec1,sec2,r
  double precision, dimension(n) :: resi 
  external Esti_Eq1
  external Esti_th
  tmp = rand(iseed)
  tol=1.e-5
  alpha = [1.0, 1.0]
  beta = 0.3
  gamma = [1.0, 0.5,1.0,-0.3]
  ms_eps = [0.0, 1.0]
  d=ms_eps(2)*6./(m-1.0)
  tmp=0
  do i=1,m
     epsdis(i)=ms_eps(1)-ms_eps(2)*3+(i-1)*d
     p(i)=exp(-(epsdis(i)-ms_eps(1))**2/2.0/ms_eps(2)**2) 
    !          p(i)=1
     tmp=tmp+p(i)
  end do
   do i=1,m
     p(i)=p(i)/tmp
  end do
  open(1,file='WLSlognor.dat')
  iter = 1
  ccount=0
1000 if (iter .le. simus) then
  call gendata(n,alpha,beta,gamma,S,Z,W,Y,G0)
  print*, 'averagy of Y', dble(sum(Y))/dble(n)
!If OLS use
!  weight = 1.0
!If WLS use
  weight = 1.0/(5./3.+abs(S)**2./2.)
  theta1=alpha
  call hybrd10(Esti_Eq1,lb1,theta1,vfd1,tol,eflag1,wv1,lw1)
  print*, 'iter=',iter,'theta1=',theta1,'vfd1=',vfd1,'eflag1=',eflag1
    r=0.001*rand(0)
  print*,'r',r
    theta2=[0.3+r, 0.5+r,1.0+r,-0.3+r]
  para = theta1
  
 !If alpha known
 !para=alpha
999 call hybrd10(Esti_th,lb2,theta2,vfd2,tol,eflag2,wv2,lw2)
  print*, 'iter=',iter,'theta2=',theta2,'vfd2=',vfd2,'eflag2=',eflag2
  if (vfd2(1).eq.999) then
     eflag2=-1
  end if
  if ((eflag2.ne.1) .or. (sum(abs(theta2)) .gt. 10) ) then
     print*,'iter=',iter,'eflag2=',eflag2,'failed Esti_th','theta2=',theta2
	 ccount=ccount+1
     goto 1000
  end if
 if (minval(abs(theta2-[0.3+r, 0.5+r,1.0+r,-0.3+r]))<1.e-5 ) then
     print*,'iter=',iter,'eflag2=',eflag2,'failed Esti_th to its initial value','theta2=',theta2
     goto 1000
  end if
  write(1,*) iter,theta1,theta2,Gng(:,1),Gng(:,2),Ung
  iter=iter+1
  goto 1000
end if

close(1)
return
end program main



subroutine gendata(n,alpha,beta,gamma,S,Z,W,Y,G)
  integer, intent(in) :: n
  double precision, dimension(2), intent(in) :: alpha
  double precision,dimension(4),intent(in)::gamma
  double precision, intent(in) :: beta
  real, dimension(n) :: t5
  double precision, dimension(n), intent(out) :: S,W
  double precision, dimension(n)::U
  double precision, dimension(4,n), intent(out) :: Z
 integer, dimension(n), intent(out) :: Y
 double precision, dimension(n):: Yt
! double precision, dimension(n), intent(out) :: Y
  double precision, dimension(n,2),intent(out)::G
  integer :: i,kflag,j
  real :: u1,u2,u3,u4,u5,u6
  double precision, dimension(n) :: H_dot, tmp,eps,t,EX,X
  S=0.0; Z=0.0; eps=0.0; X=0.0; W=0.0; Y=0
  tmp=0.0; H_dot=0.0; EX=0.0; G=0.0;
  do i=1,n
     call rnorm(u1,u2)
     call rnorm(u3,u4)
	 call rnorm(u5,u6)
     call rt(t5(i),5)
     S(i) = u1
	 Z(1,i)= u4
     Z(2,i) = u2
	 Z(3,i)=2.0*(rand(0)-0.5)
	 Z(4,i)=u5
	 !Z(4,i)=rand(0)

	 call inprod(4,gamma,Z(:,i),U(i))
	 call trueg(U(i),G(i,1))
	 call truegderi(U(i),G(i,2))
     call mm(alpha,S(i),Z(:,i),EX(i))
     eps(i) = u3*(abs(S(i))*sqrt(1./2.))
     call ownadd(EX(i), eps(i), X(i))
	!If error term follows t distribution 
     !    W(i) = X(i) + dble(t5(i))
     !If error term is normal
     W(i) = X(i) + dble(u6)*sqrt(0.6)
     !If link function is inverse logistic
     call inv_logistic(beta,G(i,1),X(i),H_dot(i))
	 ! If link function is inverse probit
	 !call inv_probit(beta,G(i,1),X(i),H_dot(i))
     tmp(i) = rand(0)

     if (tmp(i) .lt. H_dot(i)) then
        Y(i) = 1
     end if
    end do
  return
end subroutine gendata


subroutine Esti_Eq1(lb1,alpha,vfd1,eflag1)
  use commondata,only : n,S,W,weight
  integer :: i
  integer,intent(in) :: lb1
  integer,intent(out) :: eflag1
  double precision, dimension(lb1),intent(in) :: alpha
  double precision, dimension(lb1),intent(out) :: vfd1
  double precision, dimension(n,lb1) :: vfds
  if (sum(abs(alpha)).gt.100) then
     vfd1(1)=999
     return
  end if
  vfds=0.0; vfd1=0.0
  do i = 1,n
     vfds(i,1) = (W(i)-alpha(1)-alpha(2)*S(i))*weight(i)
     vfds(i,2) = S(i)*(W(i)-alpha(1)-alpha(2)*S(i))*weight(i)
  end do

  vfd1 = sum(vfds,1)

  tmp=0.0
  do i=1,lb1
     tmp=tmp+vfd1(i)
  end do
  if (.not.((tmp.gt.0).or.(tmp.le.0))) then
     vfd1(1)=999 
  end if
  return
end subroutine Esti_Eq1


subroutine Esti_g(df,Gd2,vfd3,eflag3)
 use commondata,only : n,m,S,Z,W,Y,epsdis,p,para,lb2,para2,U0
  integer :: i,j,k
  integer, dimension(m) :: ipvt
  integer,intent(in)::df
  integer,intent(out) :: eflag3
  double precision,dimension(df),intent(in)::Gd2
  double precision, dimension(df),intent(out) :: vfd3
  double precision, dimension(lb2) :: parag
  double precision, dimension(n,df) :: vfds
  double precision, dimension(m) :: wv
  double precision, dimension(n,m,m) :: A
  double precision, dimension(n,m) :: R,alpha0
  double precision, dimension(n,m,2) :: f_Y,deri
  double precision, dimension(n,2) :: int_eps,fYSZ, Sthe,Tthe
  double precision, dimension(n) :: hh1,hh2,U,g_u,Kh !which is about E(ee)
  double precision, dimension(n,m) :: alpha
  double precision :: tmp,cond,h,pi,ker
  pi=3.14159265358979
  A = 0.0; f_Y= 0.0; int_eps = 0.0; hh1 = 0.0; hh2 = 0.0
  fYSZ = 0.0; Sthe=0.0; Sthe2=0.0; alpha0=0.0
  Tthe=0.0; parag=0.0; R=0.0
  vfds=0.0; deri=0.0;g_u=0.0;U=0.0;Kh=0.0
  alpha=0.0;ker=0.0
  h=10./(dble(n)**(1./3.))
  !print*,'h',h
  parag=[1.0d0, para2(2:4)]
   do i=1,n
    call inprod(4,parag,Z(:,i),U(i))
   end do
   !rank U
   do i=1,n
     g_u(i)=Gd2(1)+Gd2(2)*(U(i)-U0)  
     do k=1,2
        do j=1,m
           call Y_condi_other([para,para2],g_u(i),S(i),Z(:,i),(abs(S(i))*sqrt(1./2.))*epsdis(j),k-1,f_Y(i,j,k),deri(i,j,k))
        end do
     end do

     do k=1,2
        do j=1,m
           int_eps(i,k) = int_eps(i,k)+p(j)*(abs(S(i))*sqrt(1./2.))*epsdis(j)*f_Y(i,j,k)
        end do
     end do
	 
     do k=1,2
        do j=1,m
           fYSZ(i,k) = fYSZ(i,k)+p(j)*f_Y(i,j,k)
        end do
     end do
  
     do k=1,2
        do j=1,m
           Sthe(i,k) = Sthe(i,k)+p(j)*deri(i,j,k)
        end do
     end do
	 
     do j=1,m
        hh1(i) = hh1(i)+p(j)*((abs(S(i))*sqrt(1./2.))*epsdis(j))**2.
     end do
     hh2(i)=1./ (hh1(i))
	 
     do j=1,m
        do k=1,m
           A(i,j,k)= p(k)*f_Y(i,k,1)&
           *(f_Y(i,j,1)-hh2(i)*(abs(S(i))*sqrt(1./2.))*epsdis(j)*int_eps(i,1))&
           /(fYSZ(i,1))&
           +p(k)*f_Y(i,k,2)&
           *(f_Y(i,j,2)-hh2(i)*(abs(S(i))*sqrt(1./2.))*epsdis(j)*int_eps(i,2))/(fYSZ(i,2))
        end do
     end do

     do j=1,m
           R(i,j)= Sthe(i,1)*(f_Y(i,j,1)-hh2(i)&
                *(abs(S(i))*sqrt(1./2.))*epsdis(j)*int_eps(i,1))/(fYSZ(i,1))&
                +Sthe(i,2)*(f_Y(i,j,2)-hh2(i)*(abs(S(i))*sqrt(1./2.))*epsdis(j)&
                *int_eps(i,2))/(fYSZ(i,2))
     end do
     
     do j=1,m
        alpha(i,j)=R(i,j)
     end do
    
     do j=1,m
       A(i,j,j)=A(i,j,j)+0.005
     end do

     call decomp(m,m,A(i,:,:),cond,ipvt,wv)
     call solvels(m,m,A(i,:,:),alpha(i,:),ipvt)


     do j=1,m
        alpha0(i,j)=alpha(i,j)
     end do

     do k=1,2
        do j=1,m
           Tthe(i,k) = Tthe(i,k)+p(j)*alpha0(i,j)*f_Y(i,j,k)
        end do
     end do

	call kernelf(h,U(i)-U0,ker)

    Kh(i)=(U(i)-U0)*(1./h)*ker

	vfds(i,1)=ker*(Sthe(i,Y(i)+1)/fYSZ(i,Y(i)+1)-Tthe(i,Y(i)+1)/fYSZ(i,Y(i)+1))

    vfds(i,2) = Kh(i)*(Sthe(i,Y(i)+1)/fYSZ(i,Y(i)+1)-Tthe(i,Y(i)+1)/fYSZ(i,Y(i)+1))

  end do
  vfd3 = sum(vfds,1)
  tmp=0.0
  tmp=tmp+vfd3(1)+vfd3(2)

  if (.not.((tmp.gt.0).or.(tmp.le.0))) then
     vfd3(1)=999 
  end if
return
end subroutine Esti_g 

subroutine Esti_th(lb2,theta2,vfd2,eflag2)
 use commondata,only : n,ng,m,S,Z,W,Y,epsdis,p,para,lw3,para2,U0,Ung,Gng,G
  integer :: i,j,k,eflag3
  integer, dimension(m) :: ipvt
  integer,intent(in) :: lb2
  integer,intent(out) :: eflag2
  double precision, dimension(n,2):: Gt
  double precision, dimension(lb2),intent(in) :: theta2
  double precision, dimension(lb2),intent(out) :: vfd2
  double precision, dimension(n,lb2) :: vfds
  double precision, dimension(m) :: wv
  double precision, dimension(lw3) :: wv3
  double precision, dimension(n,m,m) :: A
  double precision, dimension(n,4,m) :: R,alpha0
  double precision, dimension(n,m,2) :: f_Y,deri
  double precision, dimension(n,2) :: int_eps,fYSZ,Sthe1, Sthe2,Sthe3,Sthe4,Tthe1,Tthe2,Tthe3,Tthe4
  double precision, dimension(2):: Gd2,vfd3
  double precision, dimension(lb2):: parag
  double precision, dimension(n) :: U,hh1,hh2!which is about E(ee)
  double precision, dimension(n,m) :: alpha01,alpha02,alpha03,alpha04
  double precision :: tmp,cond,tol
  double precision, dimension(ng)::U1
  external Esti_g
   if (maxval(abs(theta2-[0.3, 0.5,1.0,-0.3])) .gt. 1.0) then
     vfd2(1)=999
     return
  end if

  A = 0.0; f_Y= 0.0; int_eps = 0.0; hh1 = 0.0; hh2 = 0.0
  fYSZ = 0.0; Sthe1=0.0; Sthe2=0.0; Sthe3=0.0; Sthe4=0.0
  Tthe1=0.0; Tthe2=0.0; Tthe3=0.0; Tthe4=0.0
  alpha01=0.0; alpha02=0.0; alpha03=0.0; alpha04=0.0;alpha0=0.0
  R=0.0; vfds=0.0;  deri=0.0;Gd2=0.0;U1=0.0;Gt=0.0; G=0.0
  tol=1.e-5
  print*,'theta2',theta2
  para2=theta2
  parag=[1.0d0, para2(2:4)]
  do i=1,n
  call inprod(4,parag,Z(:,i),U(i))
  end do

 U0=minval(U)
 do i=1,ng
!initial value of G
 Ung(i)=U0
call trueg(U0,Gd2(1))
call truegderi(U0,Gd2(2))

call hybrd10(Esti_g,2,Gd2,vfd3,tol,eflag3,wv3,lw3)

  if (vfd3(1).eq.999) then
     eflag3=-1
  end if
  if (eflag3.ne.1 ) then
  print*,'failed Esti_g','G',Gd2
  vfd2(1)=999
  eflag2=-1
  return
  end if
  Gng(i,1)=Gd2(1) 
  Gng(i,2)=Gd2(2)
  U1(i)=U0
  U0=U0+(maxval(U)-minval(U))/(ng-1)
  end do

! linear interpolate  
do j=1,n
do i=1,ng-1
 if((U(j) .ge. U1(i)) .and. (U(j) .le. U1(i+1)) )  then
 G(j,1)=Gng(i,1)+Gng(i,2)*(U(j)-U1(i))
G(j,2)=Gng(i+1,2)
 end if
end do
end do


do i=1,n
call trueg(U(i),Gt(i,1))
call truegderi(U(i),Gt(i,2))
end do
print*,'difference',sum(abs(G(:,1)-Gt(:,1)))/1000,sum(abs(G(:,2)-Gt(:,2)))/1000

  do i=1,n 
     do k=1,2
        do j=1,m

           call Y_condi_other([para,theta2],G(i,1),S(i),Z(:,i),(abs(S(i))*sqrt(1./2.))*epsdis(j),k-1,f_Y(i,j,k),deri(i,j,k))
        end do
     end do

     do k=1,2
        do j=1,m
           int_eps(i,k) = int_eps(i,k)+p(j)*(abs(S(i))*sqrt(1./2.))*epsdis(j)*f_Y(i,j,k)
        end do
     end do
	 
	   do k=1,2
        do j=1,m
           fYSZ(i,k) = fYSZ(i,k)+p(j)*f_Y(i,j,k)
        end do
     end do
	 
     do k=1,2
        do j=1,m
           Sthe1(i,k) = Sthe1(i,k)+p(j)*deri(i,j,k)&
                *(para(1)+para(2)*S(i)&
                +(abs(S(i))*sqrt(1./2.))*epsdis(j))
        end do
     end do
	 
     do k=1,2
        do j=1,m
           Sthe2(i,k) = Sthe2(i,k)+p(j)*deri(i,j,k)*G(i,2)*Z(2,i)
        end do
     end do
	 
	do k=1,2
        do j=1,m
           Sthe3(i,k) = Sthe3(i,k)+p(j)*deri(i,j,k)*G(i,2)*Z(3,i)
        end do
     end do
	 
	do k=1,2
        do j=1,m
           Sthe4(i,k) = Sthe4(i,k)+p(j)*deri(i,j,k)*G(i,2)*Z(4,i)
        end do
     end do	 
	 
	 
     do j=1,m
        hh1(i) = hh1(i)+p(j)*((abs(S(i))*sqrt(1./2.))*epsdis(j))**2.
     end do
	 
     hh2(i)=1./ hh1(i)
     do j=1,m
        do k=1,m
           A(i,j,k)= p(k)*f_Y(i,k,1)&
           *(f_Y(i,j,1)-hh2(i)*(abs(S(i))*sqrt(1./2.))*epsdis(j)*int_eps(i,1))&
           /(fYSZ(i,1))&
           +p(k)*f_Y(i,k,2)&
           *(f_Y(i,j,2)-hh2(i)*(abs(S(i))*sqrt(1./2.))*epsdis(j)*int_eps(i,2))/(fYSZ(i,2))
        end do
     end do
     do j=1,m
           R(i,1,j)= Sthe1(i,1)*(f_Y(i,j,1)-hh2(i)&
                *(abs(S(i))*sqrt(1./2.))*epsdis(j)*int_eps(i,1))/(fYSZ(i,1))&
                +Sthe1(i,2)*(f_Y(i,j,2)-hh2(i)*(abs(S(i))*sqrt(1./2.))*epsdis(j)&
                *int_eps(i,2))/(fYSZ(i,2))
     end do
     do j=1,m
           R(i,2,j)= Sthe2(i,1)*(f_Y(i,j,1)-hh2(i)&
                *(abs(S(i))*sqrt(1./2.))*epsdis(j)*int_eps(i,1))/(fYSZ(i,1))&
                +Sthe2(i,2)*(f_Y(i,j,2)-hh2(i)*(abs(S(i))*sqrt(1./2.))*epsdis(j)&
                *int_eps(i,2))/(fYSZ(i,2))
     end do
     do j=1,m
           R(i,3,j)= Sthe3(i,1)*(f_Y(i,j,1)-hh2(i)&
                *(abs(S(i))*sqrt(1./2.))*epsdis(j)*int_eps(i,1))/(fYSZ(i,1))&
                +Sthe3(i,2)*(f_Y(i,j,2)-hh2(i)*(abs(S(i))*sqrt(1./2.))*epsdis(j)&
                *int_eps(i,2))/(fYSZ(i,2))
     end do
	 do j=1,m
           R(i,4,j)= Sthe4(i,1)*(f_Y(i,j,1)-hh2(i)&
                *(abs(S(i))*sqrt(1./2.))*epsdis(j)*int_eps(i,1))/(fYSZ(i,1))&
                +Sthe4(i,2)*(f_Y(i,j,2)-hh2(i)*(abs(S(i))*sqrt(1./2.))*epsdis(j)&
                *int_eps(i,2))/(fYSZ(i,2))
     end do
     do j=1,m
        alpha01(i,j)=R(i,1,j)
        alpha02(i,j)=R(i,2,j)     
		alpha03(i,j)=R(i,3,j)
        alpha04(i,j)=R(i,4,j)
     end do

     do j=1,m
        A(i,j,j)=A(i,j,j)+0.005
     end do
     call decomp(m,m,A(i,:,:),cond,ipvt,wv)

     call solvels(m,m,A(i,:,:),alpha01(i,:),ipvt)
     call solvels(m,m,A(i,:,:),alpha02(i,:),ipvt)
     call solvels(m,m,A(i,:,:),alpha03(i,:),ipvt)
     call solvels(m,m,A(i,:,:),alpha04(i,:),ipvt)
     do j=1,m
        alpha0(i,1,j)=alpha01(i,j)
        alpha0(i,2,j)=alpha02(i,j)
	    alpha0(i,3,j)=alpha03(i,j)
        alpha0(i,4,j)=alpha04(i,j)
     end do

     do k=1,2
        do j=1,m
           Tthe1(i,k) = Tthe1(i,k)+p(j)*alpha0(i,1,j)*f_Y(i,j,k)
        end do
     end do
     do k=1,2
        do j=1,m
           Tthe2(i,k) = Tthe2(i,k)+p(j)*alpha0(i,2,j)*f_Y(i,j,k)
        end do
     end do
	 do k=1,2
        do j=1,m
           Tthe3(i,k) = Tthe3(i,k)+p(j)*alpha0(i,3,j)*f_Y(i,j,k)
        end do
     end do
	 do k=1,2
        do j=1,m
           Tthe4(i,k) = Tthe4(i,k)+p(j)*alpha0(i,4,j)*f_Y(i,j,k)
        end do
     end do
     vfds(i,1) = Sthe1(i,Y(i)+1)/fYSZ(i,Y(i)+1)-Tthe1(i,Y(i)+1)/fYSZ(i,Y(i)+1)
     vfds(i,2) = Sthe2(i,Y(i)+1)/fYSZ(i,Y(i)+1)-Tthe2(i,Y(i)+1)/fYSZ(i,Y(i)+1)
	 vfds(i,3) = Sthe3(i,Y(i)+1)/fYSZ(i,Y(i)+1)-Tthe3(i,Y(i)+1)/fYSZ(i,Y(i)+1)
     vfds(i,4) = Sthe4(i,Y(i)+1)/fYSZ(i,Y(i)+1)-Tthe4(i,Y(i)+1)/fYSZ(i,Y(i)+1)

  end do
  vfd2 = sum(vfds,1)
  tmp=0.0
  do i=1,lb2
     tmp=tmp+vfd2(i)
  end do

  if (.not.((tmp.gt.0).or.(tmp.le.0))) then
     vfd2(1)=999 
  end if
return
end subroutine Esti_th

!If link function is inverse logistic
subroutine Y_condi_other(theta,Gi,S,Z,epsdis,Y,fY,deri)
  double precision, dimension(6), intent(in) :: theta
  double precision, dimension(4), intent(in) :: Z
  double precision, intent(in) :: S,Gi
  double precision, intent(in) :: epsdis
  integer, intent(in) :: Y
  double precision, intent(out) :: fY,deri
  double precision :: tmp1, tmp2,tmp3,tmp4
  tmp1=0.0; tmp2=0.0; tmp3=0.0; tmp4=0.0
  call mm(theta(1:2), S,Z,tmp1)
  call ownadd(tmp1,epsdis,tmp2)
  call inv_logistic(theta(3), Gi,tmp2,tmp3)
  call der_logistic(theta(3), Gi,tmp2,tmp4)
  fY = 1-dble(Y)+(2*dble(Y)-1)*tmp3
  deri = (2*dble(Y)-1)*tmp4
  return
end subroutine Y_condi_other

!If link function is inverse probit
!subroutine Y_condi_other(theta,Gi,S,Z,epsdis,Y,fY,deri)
!  double precision, dimension(6), intent(in) :: theta
!  double precision, dimension(4), intent(in) :: Z
!  double precision, intent(in) :: S,Gi
!  double precision, intent(in) :: epsdis
!  integer, intent(in) :: Y
!  double precision, intent(out) :: fY,deri
!  double precision :: tmp1, tmp2,tmp3,tmp4
!  tmp1=0.0; tmp2=0.0; tmp3=0.0; tmp4=0.0
!  call mm(theta(1:2), S,Z,tmp1)
!  call ownadd(tmp1,epsdis,tmp2)
!  call inv_probit(theta(3), Gi,tmp2,tmp3)
!  call der_probit(theta(3), Gi,tmp2,tmp4)
!  fY = 1-dble(Y)+(2*dble(Y)-1)*tmp3
!  deri = (2*dble(Y)-1)*tmp4
!  return
!end subroutine Y_condi_other

subroutine mm(alpha,S,Z,X)
  double precision, dimension(2), intent(in) :: alpha
  double precision, intent(in):: Z,S
  double precision, intent(out) :: X
     X=alpha(1)+alpha(2)*S
  !   X=alpha(1)+alpha(2)*S+alpha(3)*Z
  return
end subroutine mm

subroutine ownadd(EX,eps,X)
  double precision, intent(in) :: EX, eps
  double precision, intent(out) :: X
  X = EX + eps
  return
end subroutine ownadd

!If link function is inverse logistic
subroutine inv_logistic(beta,Gi, X, H)
  double precision, intent(in) :: beta, X, Gi
  double precision, intent(out) :: H
  H = 1./(1.+exp(-beta*X-Gi))
  return
end subroutine inv_logistic

subroutine der_logistic(beta,G, X,  H)
  double precision, intent(in) :: beta, X,G
!   double precision,dimension(1),intent(in):: gamma
  double precision, intent(out) :: H
  H = exp(beta*X+G)/((1.+exp(beta*X+G))**2.)
  return
end subroutine der_logistic

!If link function is inverse probit
subroutine inv_probit(beta,Gi,X,Phi)
double precision, intent(in) :: beta, X, Gi
  double precision, intent(out) :: Phi
  double precision:: ccum
  call cumnor(beta*X+Gi,Phi,ccum)
  return
end subroutine inv_probit

subroutine der_probit(beta,Gi,X,dPhi)
double precision, intent(in) :: beta, X, Gi
  double precision, intent(out) :: dPhi
  pi=3.14159265358979
  dPhi=1./(sqrt(2.*pi))*exp(-(beta*X+Gi)**2./2.)
  return
end subroutine der_probit

subroutine trueg(u,g)
double precision, intent(in)::u
double precision,intent(out)::g
g=u
end subroutine trueg

subroutine truegderi(u,g)
double precision, intent(in)::u
double precision,intent(out)::g
g=1.0
end subroutine truegderi

subroutine inprod(l,s,t,y)
  integer,intent(in) :: l
  double precision,dimension(l),intent(in) :: s,t
  double precision,intent(out) :: y
  integer :: i
  y=0
  do i=1,l
     y=y+s(i)*t(i)
  end do
  return
end subroutine inprod

subroutine kernelf(h,z,ker)
double precision, intent(in):: h,z
double precision,intent(out):: ker
double precision :: pi
pi=3.14159265358979
ker=(1./h)*(1./(sqrt(2.*pi)))*exp(-(z/h)**2/2.0)
end subroutine kernelf





