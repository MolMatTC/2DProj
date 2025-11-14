program Diab_PBC
implicit none
integer :: i,j,k
integer :: NA,NB,NMat,NK,ND
integer :: Nelec, mo_o, mo_v, Norb
integer :: NFrag,NAtype,NPath
integer, dimension(:), allocatable :: NAFrag,NBF,AN,AF,Z
integer, dimension(:), allocatable :: Orb_frag, Orb_N
integer, dimension(:,:), allocatable :: VMat

double precision, dimension(:), allocatable :: Path
double precision, dimension(:,:), allocatable :: TV,P,P0
double precision, dimension(:,:,:), allocatable :: F,S

character(len=200) :: fname

namelist /SYS/fname, NFrag, NAtype,NPath,Norb
namelist /FRAG/ NAFrag,AN,AF,Orb_Frag,Orb_N
namelist /BPATH/ NK,Path

Norb=2
NK=20
Nfrag=2
NAtype=0
NPath=1

read(*,NML=SYS)

if (NAtype .eq. 0) then
  print*, "ERROR: Number of atom types (NAtype) = 0"
  stop
end if

allocate(NAFrag(Nfrag))
allocate(NBF(0:Nfrag))
allocate(AN(NAtype))
allocate(AF(NAtype))
allocate(Orb_Frag(Norb))
allocate(Orb_N(Norb))
allocate(Path(6*Npath))
allocate(P(NPath,3))
allocate(P0(NPath,3))

read(*,NML=FRAG)
read(*,NML=BPATH)

call param_reader(fname,NA,NB,ND,Nelec,mo_o,mo_v,NMat)

if (sum(NAfrag) .ne. NA) then
  print*, "ERROR: Sum of NAfrag is different than total number of atoms"
  print*, "NAFrag:",NAfrag 
  print*, "Total number of atoms in calculation:",NA
  stop
end if

allocate(Z(NA))
allocate(VMat(NMat,3))
allocate(TV(ND,3))
allocate(S(NB,NB,NMat))
allocate(F(NB,NB,NMat))

call read_PBC(fname,NA,NB,ND,Nmat,Z,VMat,TV)

call over_reader(fname,NB,Nmat,S)
call fock_reader(fname,NB,Nmat,F)
!print*, "READ OK"
!call flush(6)

call set_NBF(NB,NA,Nfrag,NAtype,NAfrag,AN,AF,Z,NBF)
call calc_path(Npath,Path,P,P0)

call coup_UCell(NB,Nfrag,Norb,NMat,NBF,Orb_Frag,Orb_N,F,S)

call coup_path(NB,ND,NK,Nfrag,Norb,Npath,NMat,Vmat,NBF,Orb_Frag,Orb_N,P0,P,TV,F,S)

print*, " "
print*, "HAPPY END :)"

end program Diab_PBC
