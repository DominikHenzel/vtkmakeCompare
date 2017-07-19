c
c matrix/tridiag.f
c (C) 2004 K. Fukagata 
c
      SUBROUTINE TRIDIAG(A,B,NN)
      IMPLICIT double precision(A-H,O-Z)
      double precision A(NN,3), B(NN)
      DO 20 I=1,NN-1
         A(I,3)=A(I,3)/A(I,2)
         B(I)=B(I)/A(I,2)
         A(I+1,2)=A(I+1,2)-A(I+1,1)*A(I,3)         
         B(I+1)=B(I+1)-A(I+1,1)*B(I)
 20   CONTINUE
      B(NN)=B(NN)/A(NN,2)
      DO 30 I=NN-1,1,-1
         B(I)=B(I)-A(I,3)*B(I+1)
 30   CONTINUE
      RETURN
      END

      SUBROUTINE TRID2(A,B,M,N)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 A(N,3), B(M,N)
      DO 20 I=1,N-1
         A(I,3)=A(I,3)/A(I,2)
         A(I+1,2)=A(I+1,2)-A(I+1,1)*A(I,3)         
         DO 25 J=1,M
            B(J,I)=B(J,I)/A(I,2)
            B(J,I+1)=B(J,I+1)-A(I+1,1)*B(J,I)
 25      CONTINUE
 20   CONTINUE
      DO 27 J=1,M
         B(J,N)=B(J,N)/A(N,2)
 27   CONTINUE
      DO 30 I=N-1,1,-1
         DO 35 J=1,M
            B(J,I)=B(J,I)-A(I,3)*B(J,I+1)
 35      CONTINUE
 30   CONTINUE
      RETURN
      END

      SUBROUTINE CTRIDIAG(A,B,NN)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 A(NN,3)
      COMPLEX*16 B(NN)
      DO 20 I=1,NN-1
         A(I,3)=A(I,3)/A(I,2)
         B(I)=B(I)/A(I,2)
         A(I+1,2)=A(I+1,2)-A(I+1,1)*A(I,3)         
         B(I+1)=B(I+1)-A(I+1,1)*B(I)
 20   CONTINUE
      B(NN)=B(NN)/A(NN,2)
      DO 30 I=NN-1,1,-1
         B(I)=B(I)-A(I,3)*B(I+1)
 30   CONTINUE
      RETURN
      END

