program Matmul
integer, parameter :: x=2
integer,dimension(x,x):: a,b,c
do i=1,x
do j=1,x
a(j,i)=i+j
b(j,i)=i*j
enddo
enddo
! Now we intend to do c = a * b 
do i=1,x
do j=1,x
c(j,i)=0
do k=1,x
c(j,i)=c(j,i)+a(j,k)*b(k,i)
enddo
enddo
enddo
print *, "Matrix c is: "
print *, c
print *, "Matrix a is: "
print *, a
print *, "Matrix b is: "
print *, b
end program Matmul
