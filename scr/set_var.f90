subroutine set_NBF(NB,NA,Nfrag,NAtype,NAfrag,AN,AF,Z,NBF)
implicit none
integer :: i,j,k,a
integer :: NB,NA,Nfrag,NAtype
integer, dimension(Nfrag) :: NAfrag
integer, dimension(0:NFrag) :: NBF
integer, dimension(NAtype) :: AN,AF
integer, dimension(NA) :: Z
NBF=0
k=0
do i=1,Nfrag
	do j=1,NAfrag(i)
		k=k+1
NFB:		do a=1,NAtype
		if ( Z(k) == AN(a) ) then		
			exit NFB
		end if
		end do NFB
		NBF(i)=NBF(i)+AF(a)	
	end do
end do

return 
end subroutine set_NBF

subroutine calc_path(N,Path,V,P0)
implicit none
integer :: i,j,k,N
double precision, dimension(6*N) :: Path
double precision, dimension(N,3) :: V,P0

k=1
do i=1,N
  P0(i,1:3)=Path(k:k+2)
  V(i,1:3)=Path(k+3:k+5)-Path(k:k+2)
  k=k+6
end do

end subroutine calc_path

subroutine calc_reciprocal(ND,TV,recv)
implicit none
integer :: i,j,k
integer :: ND
double precision, dimension(3) :: aa,aux
double precision, dimension(ND,3) :: TV,recv
double precision, parameter :: pi=3.14159265359

aux=0.0d0
aux(3)=1.0d0
recv=0.0d0

select case (ND)
  case(1)
    do i=1,3
      recv(1,i)=1/TV(1,i)
    end do
  case(2)
    do i=1,2
      aa=0.0d0
      !call cross(aux,TV(i,:),aa)
      do j=1,3
        !if ( abs(aa(j)) .ge. 1.0d-8 ) then
        if ( abs(TV(i,j)) .ge. 1.0d-8 ) then
          !recv(i,j)=1/aa(j)
          recv(i,j)=1/TV(i,j)
        end if
      end do
    end do
  case(3)
    print*, "ERROR: 3D periodicity not implemented yet"
    stop
end select

recv=pi*recv

end subroutine calc_reciprocal
