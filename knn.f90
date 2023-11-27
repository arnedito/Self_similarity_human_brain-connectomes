! CALCULO DEL AVERAGE NEAREST NEIGHBOUR DEGREE (ANNd) QUE ES UNA MEDIDA DE LAS CORRELACIONES A PRIMER ORDEN.
       PROGRAM NEAREST
       IMPLICIT NONE
       REAL*8 EXPK,AVERAGE,KK
       INTEGER*4 I,J,K,NODE,SUM,IOSTATUS,SUM2
       REAL*8, ALLOCATABLE :: ANND(:),KNN(:)
       INTEGER*4, ALLOCATABLE :: DEGREE(:), NEIGHBOURS(:)
       INTEGER*4, POINTER :: HEAD(:), COUNT(:)


       NODE = 1014

       ALLOCATE(DEGREE(NODE))
       ALLOCATE(KNN(NODE))
       KNN = 0
       AVERAGE = 0.D0
       KK = 0
       SUM2 = 0

! IMPORTAMOS LOS DEGREES Y CALCULAMOS VALOR ESPERADO
       OPEN(19, FILE = 'degrees0.txt')
       SUM = 0

       DO K = 1,NODE
        READ(19,*) I,J
        DEGREE(K) = J
        SUM = SUM + DEGREE(K)
        SUM2 = SUM2 +DEGREE(K)**2
       ENDDO

       EXPK = DBLE(SUM)/NODE
       KK = DBLE(SUM2)/NODE
       CLOSE(19)

!IMPORTAMOS LOS NEIGHBOURS
       ALLOCATE(NEIGHBOURS(SUM))
       OPEN(37, FILE = 'neighbours0.txt')

       DO K = 1,SUM
        READ(37,*) I
        NEIGHBOURS(K) = I
       ENDDO
       CLOSE(37)
       
! CALCULAMOS EL VECTOR HEAD QUE NOS SERVIRA PARA INDICAR LA POSICIÓN DEL VECINO QUE QUEREMOS ESCOGER
       ALLOCATE(HEAD(NODE+1))

       HEAD(1) = 1
       
       DO I = 2,NODE+1
        HEAD(I) = HEAD(I-1) + DEGREE(I-1)
       ENDDO

! CALCULAMOS NEAREST NEIGHBOUR DEGREE, NORMALIZADO Y DIVIDIDO POR EL VALOR ESPERADO DEL DEGREE
       KNN = 0

       DO I = 1,NODE
        DO J = 1, DEGREE(I)
            KNN(I) = KNN(I) + DBLE(DEGREE(NEIGHBOURS(HEAD(I)+J-1)+1))
        ENDDO
       ENDDO
      
! LO ESCRIBIMOS EN UN ARCHIVO JUNTO CON K/EXPK
       OPEN(68, FILE = 'annd0.txt')

!
       DO I = 1,NODE
        IF (DEGREE(I) /= 0) THEN
            KNN(I) = KNN(I)/DBLE(DEGREE(I))
        ENDIF
        AVERAGE = AVERAGE + KNN(I)
        WRITE(68,*) I, KNN(I)/EXPK
       ENDDO
       CLOSE(68)


       WRITE(*,*) 'layer 0'
       WRITE(*,*) '<k> = ',EXPK
       WRITE(*,*) '<k^2> = ', KK
       WRITE(*,*) 'ANNd =  ',AVERAGE*EXPK/KK/NODE

       ALLOCATE(ANND(NODE))
       ALLOCATE(COUNT(NODE))


! AHORA CALCULAMOS EL AVERAGE NEAREST NEIGHBOUR DEGREE, PARA REPRESENTARLO EN FUNCION DEL DEGREE Y NO DE LOS NODOS
       ANND = 0
       COUNT = 0
      
       DO I = 1, NODE
        ANND(DEGREE(I)) = ANND(DEGREE(I)) + KNN(I)
        COUNT(DEGREE(I)) = COUNT(DEGREE(I))+1
       ENDDO

       DO I = 1, NODE
        IF (COUNT(I) > 0) THEN
            ANND(I) = ANND(I) / COUNT(I)
        ENDIF
       ENDDO

       OPEN(96,FILE = 'annd20.txt')
       DO I = 1, NODE
        WRITE(96,*) I/EXPK, ANND(I)*EXPK/KK
       ENDDO
       CLOSE(96)


       END PROGRAM NEAREST
    
