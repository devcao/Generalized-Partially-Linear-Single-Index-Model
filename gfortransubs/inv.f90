! invert a matrix
subroutine inv(lb,A,A1)
  integer,intent(in) :: lb
  double precision,dimension(lb,lb),intent(in) :: A 
  double precision,dimension(lb,lb),intent(out) :: A1 
  integer :: i
  integer,dimension(lb) :: ipvt
  double precision,dimension(lb,lb) :: B
  double precision :: cond
  double precision,dimension(lb) :: wv
  B=A
  A1=0
  do i=1,lb
     A1(i,i)=1
  end do
  call decomp(lb,lb,B,cond,ipvt,wv)
  do i=1,lb
     call solvels(lb,lb,B,A1(:,i),ipvt)
  end do
  return
end subroutine inv
