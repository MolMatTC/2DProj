subroutine param_reader(fname,NA,Nbasis,ND,Nelec,mo_o,mo_v,NMat)
implicit none
integer :: i,j,k,ic,l,re
integer :: NA,Nbasis,ND,mo_o,mo_v,Nelec,NMat
integer, dimension(2) :: PBC
character (len=200), intent(in) :: fname
character (len=200) :: line

open(999,file=trim(fname)//'.fchk')

114   read(999,'(A)') line
ic=index(line,'Number of atoms')
if(ic.eq.0)go to 114
backspace(999)
read(999,50) NA

115   read(999,'(A)') line
ic=index(line,'Number of electrons')
if(ic.eq.0)go to 115
backspace(999)
read(999,50) Nelec
read(999,*)
read(999,*)
read(999,50) Nbasis

mo_o=int(Nelec/2)+mod(Nelec,2)
mo_v=Nbasis-mo_o

116   read(999,'(A)') line
ic=index(line,'Integer PBC data')
if(ic.eq.0)go to 116
read(999,*) PBC(:)

ND=PBC(1)
NMat=PBC(2)

rewind(999)
close(999)

50 format(50X,I20)
return 
end subroutine param_reader

subroutine read_PBC(fname,NA,NB,ND,Nmat,Z,VMat,TV)
implicit none
integer :: i,j,k,ic
integer :: NA,NB,ND,NMat

integer, dimension(NA) :: Z
integer, dimension(NMat,3) :: VMat
integer, dimension(NMat*3+3) :: Vaux

real, dimension(NA) :: ZZ

double precision, dimension(ND,3) :: TV
double precision, dimension(3,ND) :: TVP

character(len=200) :: fname
character (len=200) :: line

open(999,file=trim(fname)//'.fchk')

116   read(999,'(A)')line
ic=index(line,'Atomic numbers')
if(ic.eq.0)go to 116
read(999,*) ZZ(1:NA)
Z=int(ZZ)

114   read(999,'(A)') line
ic=index(line,'Translation vectors')
if(ic.eq.0)go to 114
read(999,*) TVP
TV=transpose(TVP)

115   read(999,'(A)') line
ic=index(line,'Integer PBC data')
if(ic.eq.0)go to 115
read(999,*)
read(999,*)
read(999,*)
read(999,*) Vaux

k=4
do i=1,Nmat
  VMat(i,1:3)=Vaux(k:k+2)
  k=k+3 
end do
return
end subroutine read_PBC
