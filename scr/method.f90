subroutine coup_UCell(NB,Nfrag,Norb,NMat,NBF,Orb_Frag,Orb_N,H,S)
implicit none
integer :: i,j,k,t
integer :: NB,Nfrag,NMat,Norb
integer, dimension(0:NFrag) :: NBF
integer, dimension(Norb) :: Orb_Frag, Orb_N

double precision, dimension(NB,NB,NMAT) :: H,S
double precision, parameter :: au2ev=27.2116

complex, dimension(NB,NB) :: Heff,Seff

print*, "**********************************************************************"
print*, "                          Coupling in Gamma                           "
print*, "**********************************************************************"
print*, " "

Heff=0.0d0
Seff=0.0d0
do i=1,NMat
  Heff=Heff+H(:,:,i)
  Seff=Seff+S(:,:,i)
end do
call block_proj(NB,NFrag,NBF,Heff,Seff)
Heff=Heff*au2ev

print*,"Localized Orbitals Matrix elements"
print*, " "
do i=1,Norb
do j=1,Norb
  k=sum(NBF(0:orb_frag(i)-1))+Orb_N(i)
  t=sum(NBF(0:orb_frag(j)-1))+Orb_N(j)

  if ( i .eq. j ) then
    write(*,1000) orb_N(i),orb_frag(i),real(Heff(k,k))

  else if ( orb_frag(i) .ne. orb_frag(j) ) then
    write(*,1200) orb_N(i),orb_N(j),orb_frag(i),orb_frag(j),real(Heff(t,k)),aimag(Heff(t,k))

  else if ( orb_frag(i) .eq. orb_frag(j) ) then
    if ( abs(Heff(t,k))**2 .gt. 1.0d-6 ) then
      write(*,1300) orb_N(i),orb_N(j),orb_frag(i),real(Heff(t,k)),aimag(Heff(t,k))
    end if
  end if

end do
end do

print*, " "
print*, "**********************************************************************"

1000 format("Orbital ",I3," from fragment ",I5," Energy = ", E13.5," eV")
1200 format("Orbitals ",I3," and ",I3," from fragments ",I5," and ", I5," Couplings = ", E13.5," +"E13.5,"i eV")
1300 format("WARNING: coupling between orbitals ",I5," and ",I5," from fragment",I5," too high = ",E13.5," +"E13.5,"i eV")

end subroutine coup_UCell

subroutine coup_path(NB,ND,NK,Nfrag,Norb,Npath,NMat,Vmat,NBF,Orb_Frag,Orb_N,P0,P,TV,H,S)
implicit none
integer :: i,j,k,cc,t,tt
integer :: NB,Nfrag,Norb,Npath,NMat,ND,NK
integer :: r
integer, dimension(0:Nfrag) :: NBF
integer, dimension(Norb) :: Orb_frag, Orb_N
integer, dimension(NMat,3) :: Vmat

double precision, dimension(0:NPath*NK) :: kd 
double precision, dimension(3) :: kvect,kvect0
double precision, dimension(NPath,3) :: P0,P
double precision, dimension(ND,3) :: TV,recv
double precision, dimension(NB,NB,NMat) :: H,S
double precision, parameter :: au2ev=27.2116

complex :: blochf
complex, dimension(NB,NB) :: Heff,Seff
complex, dimension(NK*Npath,Norb,Norb) :: H_print
complex, parameter :: eye= (0.0,-1.0)

character(len=20) :: fout

H_print=0.0d0
cc=0
kvect=0.0d0
kd=0.0d0
recv=0.0d0

print*, "**********************************************************************"
print*, "                      Projection along k-space                        "
print*, "**********************************************************************"
print*, " "

call calc_reciprocal(ND,TV,recv)
print*, "Real Space translation vectors (borh):"
do i=1,ND
print*, TV(i,:)
end do
print*, " "
print*, "Reciprocal space lattice vectors (bohr-1):"
do i=1,ND
print*, recv(i,:)
end do
print*, " "
call flush(6)

kvect=0.0d0
do r=1,Npath

  kvect0=kvect
  write(*,1000) r
  write(*,1010) P0(r,:)
  write(*,1020) P(r,:)+P0(r,:)
  kvect=0.0d0
  do i=1,ND
    kvect=kvect+P0(r,:)*recv(i,:)
    !kvect=kvect+P0(r,:)*recv(r,:)
  end do
  write(*,1030) kvect
  kvect=0.0d0
  do i=1,ND
    kvect=kvect+(P0(r,:)+P(r,:))*recv(i,:)
    !kvect=kvect+(P0(r,i)+P(r,i))*recv(i,:)
  end do
  write(*,1040) kvect
  print*, " "

  do k=1,NK
   
    if (k .ne. 1) kvect0=kvect
    kvect=0.0d0 
    do i=1,ND
      kvect=kvect+(P0(r,i)+P(r,i)*float(k-1)/float(NK-1))*recv(i,:)
    end do

    Heff=(0.0,0.0) 
    Seff=(0.0,0.0)

    do i=1,NMat
      blochf=1.0d0
      do j=1,3
        blochf=blochf*exp(eye*kvect(1)*TV(j,1)*VMat(i,j))
        blochf=blochf*exp(eye*kvect(2)*TV(j,2)*VMat(i,j))
        blochf=blochf*exp(eye*kvect(3)*TV(j,3)*VMat(i,j))
      end do      
      Heff=Heff+H(:,:,i)*blochf
      Seff=Seff+S(:,:,i)*blochf
    end do

do i=2,NB
  do j=1,i-1
    Heff(i,j)=conjg(Heff(j,i))
    Seff(i,j)=conjg(Seff(j,i))
  end do
end do

!!!!!!!!
!print*, "NK:", NB,r,k,NK
!print*, " "
!do i=1,NB
!print*, Heff(i,:)
!end do
!print*, " "
!stop
!!!!!!!!

    call block_proj(NB,NFrag,NBF,Heff,Seff)
!print*, "Diag",r,k
!print*, kvect0
!print*, kvect
!print*, 
call flush(6)
    Heff=Heff*au2ev

    cc=cc+1
    kd(cc)=kd(cc-1)+dsqrt(dot_product(kvect-kvect0,kvect-kvect0))
    kd(cc)=kd(cc)-kd(1)
    do i=1,Norb
    do j=1,Norb
      tt=sum(NBF(0:orb_frag(i)-1))+Orb_N(i)
      t=sum(NBF(0:orb_frag(j)-1))+Orb_N(j)
      H_print(cc,i,j)=Heff(tt,t) 
    end do
    end do

  end do
 
end do


cc=0
do i=1,Norb
	if ( i .lt. 10 ) then
	  write(fout,2000) i
	else if (i .lt. 100 ) then 
	  write(fout,2010) i
	else if (i .lt. 1000 ) then 
	  write(fout,2020) i
	else  
	  write(fout,2030) i
	end if
	
	open(10,file=trim(fout))
	do k=1,NK*Npath
	  write(10,*) kd(k),real(H_print(k,i,i)),aimag(H_print(k,i,i))
	end do
	close(10)

  do j=i+1,Norb
    cc=cc+1
    write(*,4000) cc,orb_frag(i),orb_N(i),orb_frag(j),orb_N(j)
	if ( cc .lt. 10 ) then
	  write(fout,3000) cc
	else if (i .lt. 100 ) then 
	  write(fout,3010) cc
	else if (i .lt. 1000 ) then 
	  write(fout,3020) cc
	else  
	  write(fout,3030) cc
	end if

	open(10,file=trim(fout))
	do k=1,NK*Npath
	  write(10,*) kd(k),real(H_print(k,i,j)),aimag(H_print(k,i,j))
	end do
	close(10)
  end do

end do

print*, " "
print*, "**********************************************************************"

1000 format("Sampling path: ",I5)
1010 format("Initial Path point in fractional coordinates:      ",3(E13.5,2X))
1020 format("Final Path point in fractional coordinates:        ",3(E13.5,2X))
1030 format("Initial Path point in reciprocal lattice (bohr-1): ",3(E13.5,2X))
1040 format("Final Path point in reciprocal lattice (bohr-1):   ",3(E13.5,2X))

2000 format("band_",I1,".dat")
2010 format("band_",I2,".dat")
2020 format("band_",I3,".dat")
2030 format("band_",I4,".dat")

3000 format("coup_",I1,".dat")
3010 format("coup_",I2,".dat")
3020 format("coup_",I3,".dat")
3030 format("coup_",I4,".dat")

4000  format("Coupling ",I5," Frag: ",I5," orb: ",I5, " and Frag: ",I5, " orb: ",I5)
end subroutine coup_path
