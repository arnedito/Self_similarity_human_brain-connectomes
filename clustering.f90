! Calculate clustering three-point correlations
       PROGRAM CLUSTERING
       IMPLICIT NONE
       REAL*8 EXPK,EXPC
       INTEGER*4 I,J,K,L,NODE,SUM,IOSTATUS,NA,NB
       REAL*8, ALLOCATABLE :: C(:), CK(:)
       INTEGER*4, ALLOCATABLE :: DEGREE(:), NEIGHBOURS(:),LINKS(:)
       INTEGER*4, POINTER :: HEAD(:),COUNT(:)



       NODE = 1014
       EXPK = 0.D0
       EXPC = 0.D0

       ALLOCATE(DEGREE(NODE))

! IMPORTAMOS LOS DEGREES Y CALCULAMOS VALOR ESPERADO
       OPEN(8, FILE = 'degrees0.txt')
       SUM = 0.D0

       DO K = 1,NODE
        READ(8,*) I,J
        DEGREE(K) = J
        SUM = SUM + DEGREE(K)
       ENDDO

       EXPK = DBLE(SUM)/NODE

       WRITE(*,*) 'LAYER 0'
       WRITE(*,*) '<k> = ',EXPK
       CLOSE(8)

!IMPORTAMOS LOS NEIGHBOURS
       ALLOCATE(NEIGHBOURS(SUM))
       OPEN(26, FILE = 'neighbours0.txt')

       DO K = 1,SUM
        READ(26,*) I
        NEIGHBOURS(K) = I
       ENDDO
       CLOSE(26)

! CALCULAMOS EL VECTOR HEAD QUE NOS SERVIRA PARA INDICAR LA POSICIÓN DEL VECINO QUE QUEREMOS ESCOGER
       ALLOCATE(HEAD(NODE+1))
       ALLOCATE(COUNT(NODE))

       HEAD(1) = 1
       
       DO I = 2,NODE+1
        HEAD(I) = HEAD(I-1) + DEGREE(I-1)
       ENDDO

! CALCULAMOS EL CLUSTERING DE CADA NODO (LOCAL). PRIMERO CALCULAMOS LOS LINKS QUE TIENEN ENTRE ELLOS LOS PRIMEROS VECINOS DEL NODO.
       ALLOCATE(C(NODE))
       ALLOCATE(LINKS(NODE))
       ALLOCATE(CK(NODE))

       C = 0.D0
       LINKS = 0
       CK = 0.D0
       COUNT = 0

! CALCULAMOS EL NUMERO DE TRIANGULOS CERRADOS PARA CADA NODO. SI DOS DE SUS VECINOS ESTAN CONECTADOS ENTONCES SUMAMOS 1 TRIANGULO
       DO I = 1,NODE
        DO J = 1,DEGREE(I)
            DO K = J+1,DEGREE(I)
                NA = NEIGHBOURS(HEAD(I)+J-1) ! VECINO 1 NODO I
                NB = NEIGHBOURS(HEAD(I)+K-1) ! VECINO 2 NODO I
                DO L = 1, DEGREE(NB+1)
                    IF (NA==NEIGHBOURS(HEAD(NB+1)+L-1)) THEN ! SI NA ES VECINO DE NB
                        LINKS(I) = LINKS(I) +1
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
       ENDDO

! CALCULAMOS EL CLUSTERING COEFFICIENT PARA CADA NODO Y EL VALOR ESPERADO
       DO I = 1,NODE
        IF (DEGREE(I) > 1) THEN
            C(I) = 2.*LINKS(I)/(DEGREE(I)*(DEGREE(I)-1))
            EXPC = EXPC + C(I)
        ENDIF
       ENDDO
       

! AHORA LO CALCULAMOS EN FUNCION DEL DEGREE
       DO I = 1,NODE
        CK(DEGREE(I)) = CK(DEGREE(I))+ C(I)
        COUNT(DEGREE(I)) = COUNT(DEGREE(I)) + 1
       ENDDO

       DO I = 1,NODE
        IF (COUNT(I) > 0) THEN
            CK(I) = CK(I)/COUNT(I)
        ENDIF
       ENDDO

! ESCRIBIMOS LOS VALORES EN UN ARCHIVO
       OPEN(94, FILE = 'clustering0.txt')
       DO I = 1,NODE
        WRITE(94,*) I/EXPK,CK(I),C(I)
       ENDDO

       CLOSE(94)
       WRITE(*,*) '<c> = ', DBLE(EXPC)/NODE

       ENDPROGRAM CLUSTERING

