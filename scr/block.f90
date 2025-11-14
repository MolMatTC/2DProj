subroutine block_proj(NB,Nfrag,NBF,H,S)
implicit none
integer :: i,j,k
integer :: NB,Nfrag
integer :: i0,i1
integer, dimension(0:Nfrag) :: NBF

real, dimension(:), allocatable :: eigen

complex, dimension(NB,NB) :: H,MO,S,SS,HH
complex, dimension(:,:), allocatable :: Htmp,Stmp

HH=H
SS=S
!call lowdin(NB,HH,SS)
MO=0.0d0
do i=1,Nfrag

  allocate(eigen(NBF(i)))
  allocate(Htmp(NBF(i),NBF(i)))
  allocate(Stmp(NBF(i),NBF(i)))
  i0=sum(NBF(0:i-1))+1
  i1=sum(NBF(0:i))
  Htmp=HH(i0:i1,i0:i1)
  Stmp=SS(i0:i1,i0:i1)

  call lowdin(NBF(i),Htmp,Stmp)
  call diag_C(NBF(i),Htmp,eigen)
  MO(i0:i1,i0:i1)=matmul(Stmp,Htmp)
  deallocate(Htmp)
  deallocate(Stmp)
  deallocate(eigen)

end do

allocate(Stmp(NB,NB))

Stmp=matmul(Transpose(conjg(MO)),matmul(S,MO))
H=matmul(Transpose(conjg(MO)),matmul(H,MO))

!call lowdin(NB,H,Stmp)

end subroutine block_proj
