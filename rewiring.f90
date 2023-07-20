! PROGRAMA QUE INTENTA ELIMINAR LAS CORRELACIONES ENTRE LOS NODOS INTERCAMBIANDO CONNEXIONES ALEATORIAS ENTRE DOS LINKS DE LA EDGELIST, MÉTODO CONOCIDO COMO REWIRING
       PROGRAM REWIRING
       IMPLICIT NONE
       REAL*8 RAN2, XNUMBER,EXPK,SUMPROB,ERRORC,EXPC,EXPECTEDK,EXPDEG
       INTEGER*4 IDUM,EDGES,IOSTATUS,LINE,I,J,N,LINK1,LINK2,POS1,POS2,K,AVERAGE,NODE,NEWNODE,SUM
       INTEGER*4, ALLOCATABLE :: EDGELIST(:,:),DEG(:),NEIGHBOURS(:),NEW_DEG(:),NEW_NEIGHBOUR(:)
       REAL*8, ALLOCATABLE :: PROB(:),CPROB(:),CCPROB(:),CK(:)
       
       IDUM = -20
       NODE = 1014
       NEWNODE = 1014

       EDGES = 0
       OPEN(19, FILE = '/Users/arnedojoya/TFG/random/HCP/34/HCP_nobrainstem_34_layer_0_edgelist.txt')
! CALCULAMOS LAS FILAS DE LA EDGE LIST (NUMERO DE CONNEXIONES)
       DO
        READ(19, *, IOSTAT = IOSTATUS) LINE ! LEEMOS LINEA DEL ARCHIVO
        IF (IOSTATUS /= 0) EXIT !SI LLEGA AL FINAL DEL ARCHIVO QUE SALGA
        EDGES = EDGES + 1
       ENDDO
       CLOSE(19)

! REDIMENSION EDGELIST
       ALLOCATE(EDGELIST(EDGES,2))
       EDGELIST = 0
       
! WRITE THE VALUES
       OPEN(23, FILE = '/Users/arnedojoya/TFG/random/HCP/34/HCP_nobrainstem_34_layer_0_edgelist.txt')
       DO N = 1, EDGES
        READ(23,*) I,J
        EDGELIST(N,1) = I
        EDGELIST(N,2) = J
       ENDDO
       CLOSE(23)

!rewiring
       DO N = 1,10*EDGES
        POS1 = INT(RAN2(IDUM)*EDGES+1)
        POS2 = INT(RAN2(IDUM)*EDGES+1)
! REWIRE LINKS
        CALL SWAP(EDGELIST(POS1,2),EDGELIST(POS2,2))
        LINK1 = EDGELIST(POS1,1)
        LINK2 = EDGELIST(POS2,1)

! AVOID SELF-LOOPS
        IF (LINK1 == EDGELIST(POS1,2)) THEN
            POS1 = INT(RAN2(IDUM)*EDGES+1)
! REWIRE LINKS
            CALL SWAP(EDGELIST(POS1,2),EDGELIST(POS2,2))
            LINK1 = EDGELIST(POS1,1)
        ENDIF


        IF (LINK2 == EDGELIST(POS2,2)) THEN
            POS2 = INT(RAN2(IDUM)*EDGES+1)
! REWIRE LINKS
            CALL SWAP(EDGELIST(POS1,2),EDGELIST(POS2,2))
            LINK2 = EDGELIST(POS2,1)
        ENDIF


! AVOID REPEATED LINKS
        DO K = 1,EDGES
            IF ((EDGELIST(POS1,1)==EDGELIST(K,1)).AND.(EDGELIST(POS1,2)==EDGELIST(K,2))) THEN
                POS1 = INT(RAN2(IDUM)*EDGES+1)
! REWIRE LINKS
                CALL SWAP(EDGELIST(POS1,2),EDGELIST(POS2,2))
                LINK1 = EDGELIST(POS1,1)
                LINK2 = EDGELIST(POS2,1)
                ! AVOID SELF-LOOPS AGAIN
                IF (LINK1 == EDGELIST(POS1,2)) THEN
                    POS1 = INT(RAN2(IDUM)*EDGES+1)
! REWIRE LINKS
                    CALL SWAP(EDGELIST(POS1,2),EDGELIST(POS2,2))
                    LINK1 = EDGELIST(POS1,1)
                ENDIF
                IF (LINK2 == EDGELIST(POS2,2)) THEN
                    POS2 = INT(RAN2(IDUM)*EDGES+1)
! REWIRE LINKS
                    CALL SWAP(EDGELIST(POS1,2),EDGELIST(POS2,2))
                    LINK2 = EDGELIST(POS2,1)
                ENDIF
            ENDIF
 
            IF ((EDGELIST(POS1,1)==EDGELIST(K,2)).AND.(EDGELIST(POS1,2)==EDGELIST(K,1))) THEN
                POS1 = INT(RAN2(IDUM)*EDGES+1)
! IF LINK IS REPEATED WE CHANGE THE REPEATED LINK AND APPLY REWIRING AGAIN
                CALL SWAP(EDGELIST(POS1,2),EDGELIST(POS2,2))
                LINK1 = EDGELIST(POS1,1)
                LINK2 = EDGELIST(POS2,1)
            ! AVOID SELF-LOOPS AGAIN
                IF (LINK1 == EDGELIST(POS1,2)) THEN
                    POS1 = INT(RAN2(IDUM)*EDGES+1)
! REWIRE LINKS
                    CALL SWAP(EDGELIST(POS1,2),EDGELIST(POS2,2))
                    LINK1 = EDGELIST(POS1,1)
                ENDIF
                IF (LINK2 == EDGELIST(POS2,2)) THEN
                    POS2 = INT(RAN2(IDUM)*EDGES+1)
! REWIRE LINKS
                    CALL SWAP(EDGELIST(POS1,2),EDGELIST(POS2,2))
                    LINK2 = EDGELIST(POS2,1)
                ENDIF
            ENDIF

            IF ((EDGELIST(POS2,1)==EDGELIST(K,1)).AND.(EDGELIST(POS2,2)==EDGELIST(K,2))) THEN
                POS2 = INT(RAN2(IDUM)*EDGES+1)
! IF LINK IS REPEATED WE CHANGE THE REPEATED LINK AND APPLY REWIRING AGAIN
                CALL SWAP(EDGELIST(POS1,2),EDGELIST(POS2,2))
                ! AVOID SELF-LOOPS AGAIN
                LINK1 = EDGELIST(POS1,1)
                LINK2 = EDGELIST(POS2,1)
                IF (LINK1 == EDGELIST(POS1,2)) THEN
                    POS1 = INT(RAN2(IDUM)*EDGES+1)
! REWIRE LINKS
                    CALL SWAP(EDGELIST(POS1,2),EDGELIST(POS2,2))
                    LINK1 = EDGELIST(POS1,1)
                ENDIF
                IF (LINK2 == EDGELIST(POS2,2)) THEN
                    POS2 = INT(RAN2(IDUM)*EDGES+1)
! REWIRE LINKS
                    CALL SWAP(EDGELIST(POS1,2),EDGELIST(POS2,2))
                    LINK2 = EDGELIST(POS2,1)
                ENDIF
            ENDIF

            IF ((EDGELIST(POS2,1)==EDGELIST(K,2)).AND.(EDGELIST(POS2,2)==EDGELIST(K,1))) THEN
                POS2 = INT(RAN2(IDUM)*EDGES+1)
! IF LINK IS REPEATED WE CHANGE THE REPEATED LINK AND APPLY REWIRING AGAIN
                CALL SWAP(EDGELIST(POS1,2),EDGELIST(POS2,2))
                LINK1 = EDGELIST(POS1,1)
                LINK2 = EDGELIST(POS2,1)
    ! AVOID SELF-LOOPS AGAIN
                IF (LINK1 == EDGELIST(POS1,2)) THEN
                    POS1 = INT(RAN2(IDUM)*EDGES+1)
! REWIRE LINKS
                    CALL SWAP(EDGELIST(POS1,2),EDGELIST(POS2,2))
                    LINK1 = EDGELIST(POS1,1)
                ENDIF

                IF (LINK2 == EDGELIST(POS2,2)) THEN
                    POS2 = INT(RAN2(IDUM)*EDGES+1)
! REWIRE LINKS
                    CALL SWAP(EDGELIST(POS1,2),EDGELIST(POS2,2))
                    LINK2 = EDGELIST(POS2,1)
                ENDIF
            ENDIF
        ENDDO
       ENDDO


! ESCRIBIMOS LA EDGELIST RANDOMIZADA
       OPEN(152, FILE = '/Users/arnedojoya/TFG/random/HCP/34/random_edgelist.txt')
       DO I = 1,EDGES
        WRITE(152,*) EDGELIST(I,1),EDGELIST(I,2)
       ENDDO
       CLOSE(152)

       ALLOCATE(DEG(NODE))
       DEG = 0
       ALLOCATE(CK(NODE))
       CK = 0.D0
       CALL DEGREE(EDGELIST,EDGES,NODE,DEG,EXPK,AVERAGE)
       ALLOCATE(NEIGHBOURS(AVERAGE))
       NEIGHBOURS = 0
       CALL VECINOS(NODE,DEG,EDGES,EDGELIST,NEIGHBOURS,AVERAGE)
       CALL CLUSTERING(SUM,DEG,NODE,NEWNODE,NEIGHBOURS,EXPC,CK,ERRORC)

       WRITE(*,*) EXPK,EXPC,EDGES,NODE
       ENDPROGRAM REWIRING






!--------------SUBRUTINAS Y FUNCIONES-----------------------

!*******************************************************************
!       FUNCTION ran2(idum)
!       INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
!       REAL*8 ran2,AM,EPS,RNMX
!       PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,&
!          IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,&
!          IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
!Long period (> 2 × 1018) random number generator of L’Ecuyer with Bays-Durham shuffle
!and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive
!of the endpoint values). Call with idum a negative integer to initialize; thereafter, do not
!alter idum between successive deviates in a sequence. RNMX should approximate the largest
!floating value that is less than 1.
!       INTEGER idum2,j,k,iv(NTAB),iy
!       SAVE iv,iy,idum2
!       DATA idum2/123456789/, iv/NTAB*0/, iy/0/
!       if (idum.le.0) then !Initialize.
!            idum=max(-idum,1) !Be sure to prevent idum = 0.
!            idum2=idum
!            do j=NTAB+8,1,-1 !Load the shuffle table (after 8 warm-ups).
!                k=idum/IQ1
!                idum=IA1*(idum-k*IQ1)-k*IR1
!                if (idum.lt.0) idum=idum+IM1
 !               if (j.le.NTAB) iv(j)=idum
!            enddo
!            iy = iv(1)
!       endif
!       k=idum/IQ1 !Start here when not initializing.
 !      idum=IA1*(idum-k*IQ1)-k*IR1 !Compute idum=mod(IA1*idum,IM1) without overif (idum.lt.0) idum=idum+IM1 flows by Schrage’s method.
!       k=idum2/IQ2
!       idum2=IA2*(idum2-k*IQ2)-k*IR2 !Compute idum2=mod(IA2*idum2,IM2) likewise.
!       if (idum2.lt.0) idum2=idum2+IM2
!       j=1+iy/NDIV !Will be in the range 1:NTAB.
!       iy=iv(j)-idum2 !Here idum is shuffled, idum and idum2 are comiv(j)=idum bined to generate output.
!       if(iy.lt.1)iy=iy+IMM1
!       ran2=min(AM*iy,RNMX) !Because users don’t expect endpoint values.
!       return
!       END

!******** Uniform Random generator ***********************
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      DOUBLE PRECISION ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,&
        IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,&
        IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
        !Long period (> 2 x 10 18 ) random number generator of L'Ecuyer with Bays-Durham shuffle
        !and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive
        !of the endpoint values). Call with idum a negative integer to initialize; thereafter, do not
        !alter idum between successive deviates in a sequence. RNMX should approximate the largest
        !floating value that is less than 1.
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then               !Initialize.
          idum=max(-idum,1)             !Be sure to prevent idum = 0.
          idum2=idum
          do j=NTAB+8,1,-1           !Load the shuffle table (after 8 warm-ups).
               k=idum/IQ1
               idum=IA1*(idum-k*IQ1)-k*IR1
               if (idum.lt.0) idum=idum+IM1
               if (j.le.NTAB) iv(j)=idum
          enddo
          iy=iv(1)
      endif
      k=idum/IQ1                        !Start here when not initializing.
      idum=IA1*(idum-k*IQ1)-k*IR1       !Compute idum=mod(IA1*idum,IM1) without over-
      if (idum.lt.0) idum=idum+IM1      !flows by Schrage's method.
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2     !Compute idum2=mod(IA2*idum2,IM2) likewise.
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV                       !Will be in the range 1:NTAB.
      iy=iv(j)-idum2                    !Here idum is shuffled, idum and idum2 are com-
      iv(j)=idum                        !bined to generate output.
      if(iy.lt.1)iy=iy+IMM1
      ran2=dmin1(AM*dble(iy),RNMX)              !Because users don't expect endpoint values.
      return
      END





!*******************************************************************
! SUBRUTINA CLUSTERING
       SUBROUTINE CLUSTERING(SUM,DEGREE,NODE,NEWNODE,NEIGHBOURS,EXPC,CK,ERRORC)
       IMPLICIT NONE
       INTEGER*4 NODE,I,J,K,NA,NB,SUM,L,NEWNODE,M,N,O,P,Q
       REAL*8 EXPC,ERRORC,SIGMA
       REAL*8 C(NODE),CK(NODE)
       INTEGER*4 LINKS(NODE),DEGREE(NODE),NEIGHBOURS(SUM)
       INTEGER*4,POINTER:: HEAD(:),COUNT(:)

       ALLOCATE(HEAD(NODE+1))
       ALLOCATE(COUNT(NODE+1))
       C = 0.D0
       LINKS = 0
       CK = 0.D0
       COUNT = 0
       HEAD(1) = 1
       EXPC = 0.D0
       ERRORC = 0.D0

       DO M = 2, NODE+1
        HEAD(M) = HEAD(M-1) + DEGREE(M-1)
       ENDDO


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
       DO N = 1,NODE
        IF (DEGREE(N) > 1) THEN
            C(N) = 2.*LINKS(N)/(DEGREE(N)*(DEGREE(N)-1))
            EXPC = EXPC + C(N)
        ENDIF
       ENDDO
       

! AHORA LO CALCULAMOS EN FUNCION DEL DEGREE
       DO O = 1,NODE
        IF (DEGREE(O) > 0) THEN
            CK(DEGREE(O)) = CK(DEGREE(O))+ C(O)
            COUNT(DEGREE(O)) = COUNT(DEGREE(O)) + 1
        ENDIF
       ENDDO

       DO P = 1,NODE
        IF (COUNT(P) > 0) THEN
            CK(P) = CK(P)/COUNT(P)
        ENDIF
       ENDDO
       EXPC = DBLE(EXPC)/NEWNODE

! ERROR
       ERRORC = 0.D0
       DO Q = 1,NODE
        SIGMA = SIGMA + (C(Q)-EXPC)**2
       ENDDO
       SIGMA = SQRT(SIGMA/DBLE(NEWNODE))
       ERRORC = SIGMA/ SQRT(DBLE(NEWNODE))
       END SUBROUTINE

!*******************************************************************
! SUBRUTINA DEGREES Y VALOR ESPERADO
       SUBROUTINE DEGREE(EDGELIST,EDGES,NODE,DEG,EXPK,AVERAGE)
       IMPLICIT NONE
       INTEGER*4 EDGES,NODE,N,K,I,J,AVERAGE
       REAL*8 EXPK
       INTEGER*4 EDGELIST(EDGES,2), DEG(NODE)
        DO N = 1,EDGES
            I = EDGELIST(N,1)
            J = EDGELIST(N,2)
            DEG(I+1) = DEG(I+1) + 1 !SUMAMOS 1 YA QUE EL PRIMER NODO ES 0
            DEG(J+1) = DEG(J+1) + 1
        ENDDO
! CALCULAR EXP K
       AVERAGE = 0

       DO K = 1,NODE
        AVERAGE = AVERAGE + DEG(K)

       ENDDO
       EXPK = AVERAGE/DBLE(NODE)
    
       END SUBROUTINE

!*******************************************************************
! SUBRUTINA VECTOR NEIGHBOURS
       SUBROUTINE VECINOS(NODE,DEGREE,EDGES,EDGELIST,NEIGHBOURS,SUMDEG)
       IMPLICIT NONE
       INTEGER*4 I,J,K,NODE,SUMDEG,EDGES
       INTEGER*4 HEAD(NODE+1),COUNT(NODE), NEIGHBOURS(SUMDEG),DEGREE(NODE),EDGELIST(EDGES,2)

       HEAD(1) = 1
       COUNT(1) = 1
        
       DO I = 2, NODE+1
        HEAD(I) = HEAD(I-1) + DEGREE(I-1)
        COUNT(I) = 1
       ENDDO


!CALCULAMOS EL VECTOR V RELLENANDO A MEDIDA QUE ENCONTRAMOS LOS VECINOS, SUMAMOS 1 YA QUE EL PRIMER NODO ES EL 0 Y LA POSICIÓN 0 NO EXISTE
       DO K = 1,EDGES
        I = EDGELIST(K,1)
        J = EDGELIST(K,2)
        NEIGHBOURS(HEAD(I+1)+COUNT(I+1)-1) = J
        COUNT(I+1) = COUNT(I+1) + 1
        NEIGHBOURS(HEAD(J+1)+COUNT(J+1)-1) = I
        COUNT(J+1) = COUNT(J+1) + 1
       ENDDO
       END SUBROUTINE

!*******************************************************************
! DISTRIBUCION DE GRADO
       SUBROUTINE DEGREE_DISTRIBUTION(NODE,NEWNODE,DEGREE,PROB,CPROB,CCPROB,SUMPROB)
       IMPLICIT NONE
       INTEGER*4 I,J,K,N,NODE,NEWNODE
       INTEGER*4 DEGREE(NODE)
       REAL*8 SUMPROB
       REAL*8 PROB(NODE),CPROB(NODE),CCPROB(NODE)
!CALCULAMOS LAS PROBABILIDADES (FRECUENCIAS DEL HISTOGRAMA)
       DO I = 1,NODE
        DO J = 1,NODE
            IF (DEGREE(I).EQ.J) THEN
                PROB(J) = PROB(J)+1.d0!DEGREE DISTRIBUTION
            ENDIF
            IF (DEGREE(I).LE.J) THEN
                CPROB(J) = CPROB(J)+1.d0 !CUMULATIVE DEGREE DISTRIBUTION
            ENDIF
            IF (DEGREE(I).GT.J) THEN
                CCPROB(J) = CCPROB(J) +1.d0 ! COMPLEMENTARY CUMULATIVE DEGREE DISTRIBUTION
            ENDIF
        ENDDO
       ENDDO

! ELIMINAMOS VALORES DE LA CCDD SI SON IGUALES PARA QUE SEA MAS NITIDO EL GRAFICO
       DO N = 1,NODE
        IF (CCPROB(N)==CCPROB(N+1))THEN
            CCPROB(N)=0.d0
        ENDIF
       ENDDO
       DO K = 1, NODE
        PROB(K) = PROB(K)/NEWNODE
        CPROB(K) = CPROB(K)/NEWNODE
        CCPROB(K) = CCPROB(K)/NEWNODE
        SUMPROB = SUMPROB + PROB(K)
       ENDDO
       
       END SUBROUTINE

!**********************************
! VALOR ESPERADO DEGREE
       REAL*8 FUNCTION EXPECTEDK(DEGREE,NODE,NEWNODE)
       IMPLICIT NONE
       INTEGER*4 NODE,I,NEWNODE
       INTEGER*4 DEGREE(NODE)
       EXPECTEDK = 0
       DO I = 1,NODE
        EXPECTEDK = EXPECTEDK + DEGREE(I)
       ENDDO
       EXPECTEDK = EXPECTEDK/NEWNODE
       RETURN
       END

!*******************************************************************
! AVERAGE NEAREST NEIGHBOUR DEGREE
       SUBROUTINE AVGNNDEG(NODE,DEGREE,ANND,NEIGHBOURS,EDGES)
       IMPLICIT NONE
       INTEGER*4 I,J,K,N,L,M,O,P,NODE,EDGES
       INTEGER*4 HEAD(NODE+1),DEGREE(NODE),COUNT(NODE),NEIGHBOURS(EDGES)
       REAL*8 SUM,SUM2,EXPK,KK
       REAL*8 KNN(NODE),ANND(NODE)

       HEAD(1) = 1
       
       DO N = 2,NODE+1
        HEAD(I) = HEAD(N-1) + DEGREE(N-1)
       ENDDO

! CALCULAMOS NEAREST NEIGHBOUR DEGREE, NORMALIZADO Y DIVIDIDO POR EL VALOR ESPERADO DEL DEGREE

       DO I = 1,NODE
        DO J = 1, DEGREE(I)
            KNN(I) = KNN(I) + DBLE(DEGREE(NEIGHBOURS(HEAD(I)+J-1)+1))
        ENDDO
       ENDDO

       DO K = 1,NODE
        SUM = SUM + DEGREE(K)
        SUM2 = SUM2 +DEGREE(K)**2
       ENDDO

       EXPK = DBLE(SUM)/NODE
       KK = DBLE(SUM2)/NODE

       DO L = 1,NODE
        IF (DEGREE(L) /= 0) THEN
            KNN(L) = KNN(L)/DBLE(DEGREE(L))
        ENDIF
       ENDDO
    
! AHORA CALCULAMOS EL AVERAGE NEAREST NEIGHBOUR DEGREE, PARA REPRESENTARLO EN FUNCION DEL DEGREE Y NO DE LOS NODOS
       ANND = 0
       COUNT = 0
      
       DO M = 1, NODE
        ANND(DEGREE(M)) = ANND(DEGREE(M)) + KNN(M)
        COUNT(DEGREE(M)) = COUNT(DEGREE(M))+1
       ENDDO

       DO O = 1, NODE
        IF (COUNT(O) > 0) THEN
            ANND(O) = ANND(O) / COUNT(O)
        ENDIF
       ENDDO

       DO P = 1, NODE
        ANND(P) = ANND(P)*EXPK/KK
       ENDDO
       END SUBROUTINE


!************************
! DEGREE THRESHOLDING
       SUBROUTINE DEGREE_THRESHOLDING(EDGELIST,DEGREE,NODE,EDGES,NEW_DEG,THRESHOLD,NEWNODE,SUM)
       IMPLICIT NONE
       INTEGER*4 I,J,N,EDGES,THRESHOLD,K,L,M,NEWNODE,NODE,COUNT,P,SUM
       INTEGER*4 DEGREE(NODE),NEW_DEG(NODE),EDGELIST(EDGES,2)
       INTEGER*4, POINTER :: NULL(:)
       LOGICAL :: EXCLUDE_I, EXCLUDE_J
! UNA VEZ TENEMOS LOS DEGREES DE CADA NODO, APLICAMOS EL THRESHOLD Y GUARDAMOS EN EL VECTOR NULL LOS NODOS QUE ELIMINAMOS
       ALLOCATE(NULL(NODE))
       COUNT = 0
       SUM = 0
       NEWNODE = NODE
       NEW_DEG = DEGREE
       DO K = 1,NODE
        IF (NEW_DEG(K)<=THRESHOLD) THEN
            NEW_DEG(K) = 0
            NULL(COUNT+1) = K !NODOS QUE TIENEN DEGREE NULO
            COUNT = COUNT + 1
            NEWNODE = NEWNODE - 1 !NODOS TRAS APLICAR EL THRESHOLD
        ENDIF
       ENDDO

! RESETEAMOS VECTOR DEGREE
       NEW_DEG = 0
! CALCULAMOS EL VECTOR QUE CONTIENE LOS GRADOS APLICANDO EL THRESHOLD
       DO L = 1, EDGES
        I = EDGELIST(L,1)
        J = EDGELIST(L,2)

! COMPROBAMOS SI LOS NODOS DEL VECTOR NULL ESTAN EN LA FILA QUE LEEMOS
! SI ESTA EN PRIMERA COLUMNA I
        EXCLUDE_I = .FALSE.
        DO M = 1, COUNT
            IF (I == NULL(M)) THEN
                EXCLUDE_I = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI ESTA EN PRIMERA COLUMNA J
        EXCLUDE_J = .FALSE.
        DO N = 1, COUNT
            IF (J == NULL(N)) THEN
                EXCLUDE_J = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI NINGUNO DE LOS DOS NODOS ESTÁ EN EL VECTOR NULL SUMAMOS EL DEGREE
        IF (.NOT.EXCLUDE_I.AND..NOT.EXCLUDE_J) THEN
            NEW_DEG(I+1) = NEW_DEG(I+1) + 1 !SUMAMOS 1 YA QUE EL PRIMER NODO ES 0
            NEW_DEG(J+1) = NEW_DEG(J+1) + 1
        ENDIF
       ENDDO
       DO P = 1,NODE
        SUM = SUM + NEW_DEG(P)
       ENDDO
       DEALLOCATE(NULL)
       END SUBROUTINE

!***************************+
!NEIGHBOUR THREHSOLDING
       SUBROUTINE NEIGHBOUR_THRESHOLDING(EDGELIST,SUM,DEGREE,NODE,EDGES,NEW_NEIGHBOUR,THRESHOLD,NEWNODE)
       IMPLICIT NONE
       INTEGER*4 I,J,N,EDGES,THRESHOLD,K,L,M,NEWNODE,NODE,COUNT,P,SUM
       INTEGER*4 DEGREE(NODE),NEW_DEG(NODE),EDGELIST(EDGES,2),HEAD(NODE+1),INDEX(NODE+1),NEW_NEIGHBOUR(SUM)
       INTEGER*4, POINTER :: NULL(:)
       LOGICAL :: EXCLUDE_I, EXCLUDE_J

       NEW_NEIGHBOUR = 0
       HEAD(1) = 1
       INDEX(1) = 1


       ALLOCATE(NULL(NODE))
       COUNT = 0
       NEWNODE = NODE
       DO K = 1,NODE
        IF (DEGREE(K)<=THRESHOLD) THEN
            DEGREE(K) = 0
            NULL(COUNT+1) = K !NODOS QUE TIENEN DEGREE NULO
            COUNT = COUNT + 1
            NEWNODE = NEWNODE - 1 !NODOS TRAS APLICAR EL THRESHOLD
        ENDIF
       ENDDO

       DEGREE = 0
       
! CALCULAMOS EL VECTOR QUE CONTIENE LOS GRADOS APLICANDO EL THRESHOLD
       DO L = 1, EDGES
        I = EDGELIST(L,1)
        J = EDGELIST(L,2)

! COMPROBAMOS SI LOS NODOS DEL VECTOR NULL ESTAN EN LA FILA QUE LEEMOS
! SI ESTA EN PRIMERA COLUMNA I
        EXCLUDE_I = .FALSE.
        DO M = 1, COUNT
            IF (I == NULL(M)) THEN
                EXCLUDE_I = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI ESTA EN PRIMERA COLUMNA J
        EXCLUDE_J = .FALSE.
        DO N = 1, COUNT
            IF (J == NULL(N)) THEN
                EXCLUDE_J = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI NINGUNO DE LOS DOS NODOS ESTÁ EN EL VECTOR NULL SUMAMOS EL DEGREE
        IF (.NOT.EXCLUDE_I.AND..NOT.EXCLUDE_J) THEN
            DEGREE(I+1) = DEGREE(I+1) + 1 !SUMAMOS 1 YA QUE EL PRIMER NODO ES 0
            DEGREE(J+1) = DEGREE(J+1) + 1
        ENDIF
       ENDDO

       DO N = 2, NODE+1
        HEAD(N) = HEAD(N-1) + DEGREE(N-1)
        INDEX(N) = 1
       ENDDO


! CALCULAMOS EL VECTOR QUE CONTIENE LOS NEIGHBOURS APLICANDO EL THRESHOLD
       DO L = 1, EDGES
        I = EDGELIST(L,1)
        J = EDGELIST(L,2)

! COMPROBAMOS SI LOS NODOS DEL VECTOR NULL ESTAN EN LA FILA QUE LEEMOS
! SI ESTA EN PRIMERA COLUMNA I
        EXCLUDE_I = .FALSE.
        DO M = 1, COUNT
            IF (I == NULL(M)) THEN
                EXCLUDE_I = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI ESTA EN PRIMERA COLUMNA J
        EXCLUDE_J = .FALSE.
        DO N = 1, COUNT
            IF (J == NULL(N)) THEN
                EXCLUDE_J = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI NINGUNO DE LOS DOS NODOS ESTÁ EN EL VECTOR NULL SUMAMOS EL DEGREE
        IF (.NOT.EXCLUDE_I.AND..NOT.EXCLUDE_J) THEN
            NEW_NEIGHBOUR(HEAD(I+1)+INDEX(I+1)-1) = J
            INDEX(I+1) = INDEX(I+1) + 1
            NEW_NEIGHBOUR(HEAD(J+1)+INDEX(J+1)-1) = I
            INDEX(J+1) = INDEX(J+1) + 1
        ENDIF
       ENDDO
       DEALLOCATE(NULL)
       END SUBROUTINE


!**************************************************************
! SUBRUTINA PARA INTERCAMBIAR CONEXIONES
       SUBROUTINE SWAP(A,B)
       IMPLICIT NONE
       INTEGER*4 A,B,C
       C = A
       A = B
       B = C
       END SUBROUTINE
