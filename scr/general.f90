subroutine Lowdin(Nbasis,P,S)
implicit none
integer :: i
integer, intent(in) :: Nbasis

complex, dimension(Nbasis,Nbasis) :: P,S,SS
complex, dimension(Nbasis,Nbasis) :: S2
real, dimension(:), allocatable :: eigen

allocate(eigen(Nbasis))

SS=S
S2=(0.0,0.0)
call diag_C(Nbasis,SS,eigen)

do i=1,Nbasis
        S2(i,i)=1.0/sqrt(eigen(i))
end do

S2=matmul(matmul(SS,S2),transpose(conjg(SS)))
P=matmul(S2,matmul(P,S2))

S=S2

deallocate(eigen)

end subroutine Lowdin

subroutine diag_C(NB,H,eigen)
implicit none
integer :: i,NB
real, dimension(NB) :: eigen
complex, dimension(NB,NB) :: H

integer, dimension(:), allocatable :: IW
integer :: LW,LRW,LIW,info
real, dimension(:), allocatable :: RW
complex, dimension(:), allocatable :: W

allocate(W(1))
allocate(RW(1))
allocate(IW(1))
LW=-1
LRW=-1
LIW=-1
W=0.0
RW=0.0d0
IW=0.0d0
info=0 

call cheevd('V','L',NB,H,NB,eigen,W,LW,RW,LRW,IW,LIW,info)

if (info /= 0 ) then
        print*, "cheevd 1 failed"
        stop
end if

LW=int(W(1))
LRW=int(RW(1))
LIW=int(IW(1))

deallocate(W)
deallocate(RW)
deallocate(IW)
allocate(W(LW))
allocate(RW(LRW))
allocate(IW(LIW))

W=0.0
RW=0.0d0
IW=0.0d0

call cheevd('V','U',NB,H,NB,eigen,W,LW,RW,LRW,IW,LIW,info)
if (info /= 0 ) then
        print*, "cheevd 2 failed"
	print*, info
        stop
end if

deallocate(W)
deallocate(RW)
deallocate(IW)
end subroutine diag_C

subroutine cross(a, b, cros) 
  double precision, DIMENSION(3) :: cros
  double precision, DIMENSION(3), INTENT(IN) :: a, b

  cros(1) = a(2) * b(3) - a(3) * b(2)
  cros(2) = a(3) * b(1) - a(1) * b(3)
  cros(3) = a(1) * b(2) - a(2) * b(1)
END subroutine cross

