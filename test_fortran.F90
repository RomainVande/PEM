program test_fortan

  REAL ::  a(100)
  REAL ::  b(100)


   a(:)=10
   b(:)=50
   print *, a(:).LT.0
    
   a(8)=-1
   a(10)=-1
   a(9)=-1

   print *, a(:).LT.0
   print *, "a"
   print *, b(:)*(a(:).LT.0)
   print *, 'b'
   print *, sum(b(:)*(a(:).LT.0))
    print *, "TEST DONE"
end program
