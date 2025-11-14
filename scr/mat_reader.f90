SUBROUTINE over_reader(fname,Nbasis,Nmat,over)
implicit none
double precision, dimension(:), allocatable :: temp      
double precision, dimension(Nbasis,Nbasis,Nmat), intent(out) :: over
integer :: i,j,k,l,error,ic,nci5,re
integer, intent(in) :: Nbasis, Nmat
character (len=200), intent(in) :: fname
character (len=50) :: line

open(1050,file=trim(fname)//'.over')
rewind(1050)

nci5=((Nbasis*Nbasis-Nbasis)/2+Nbasis)*(NMat+1)
allocate(temp(nci5))
!allocate(over(Nbasis,Nbasis))

temp=0.0d0
over=0.0d0

re=mod(nci5,5)
if ( re.eq.0) then
	l=0
else
	l=1
end if

100   read(1050,'(A)') line
   ic=index(line(1:19),' Dump of file   514')
   if(ic.eq.0)go to 100

   do i=1,int(nci5/5)+l
   if (i*5 <= nci5) then
   	read(1050,*) (temp((i-1)*5+j),j=1,5)
   else if (i*5 > nci5) then
   	read(1050,*) (temp((i-1)*5+j),j=1,re)
   end if
   end do

k=((Nbasis*Nbasis-Nbasis)/2+Nbasis)
do l = 1, Nmat
do i = 1, Nbasis
 do j = 1,i
  k=k+1 
  over(i,j,l) = temp(k)
  over(j,i,l) =  over(i,j,l)
 enddo
enddo
enddo


deallocate (temp)

close(1050)
RETURN
End subroutine over_reader


SUBROUTINE fock_reader(fname,Nbasis,NMat,over)
implicit none
double precision, dimension(:), allocatable :: temp      
double precision, dimension(Nbasis,Nbasis,NMat), intent(out) :: over
integer :: i,j,k,l,error,ic,nci5,re
integer, intent(in) :: Nbasis,NMat
character (len=200), intent(in) :: fname
character (len=50) :: line

open(1050,file=trim(fname)//'.fock')
rewind(1050)

nci5=((Nbasis*Nbasis-Nbasis)/2+Nbasis)*(NMat)
allocate(temp(nci5))
!allocate(over(Nbasis,Nbasis))

temp=0.0d0
over=0.0d0

re=mod(nci5,5)
if ( re.eq.0) then
	l=0
else
	l=1
end if

100   read(1050,'(A)') line
   ic=index(line(1:19),' Dump of file   536')
   if(ic.eq.0)go to 100

   do i=1,int(nci5/5)+l
   if (i*5 <= nci5) then
   	read(1050,*) (temp((i-1)*5+j),j=1,5)
   else if (i*5 > nci5) then
   	read(1050,*) (temp((i-1)*5+j),j=1,re)
   end if
   end do
  
k=0
do l = 1, NMat
do i = 1, Nbasis
 do j = 1,i
  k=k+1 
  over(i,j,l) = temp(k)
  over(j,i,l) = temp(k)
 enddo
enddo
enddo

deallocate (temp)

close(1050)
RETURN
End subroutine fock_reader



