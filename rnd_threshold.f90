! PROGRAMA EN EL QUE APLICAMOS EL DEGREE THRESHOLDING PARA CONSTRUIR GRAFOS ANIDADOS CON NODOS CON UN DEGREE SUPERIOR AL VALOR UMBRAL, DE LOS CUALES ESTUDIAREMOS SUS PROPIEDADES ESTADÍSTICAS
       PROGRAM DEGREE_THRESHOLDING_RENORMALIZATION
       IMPLICIT NONE
       INTEGER*4 I,J,N,NODE,SUM,K,EDGES,LINE,IOSTATUS,KT,NEWNODE,COUNT
       INTEGER*4, ALLOCATABLE :: DEGREE(:), NEIGHBOURS(:)
       INTEGER*4, POINTER :: HEAD(:), INDEX(:), NULL(:)
       LOGICAL :: EXCLUDE_I, EXCLUDE_J
       
! -----------------------------------------------KT = 5 --------------------------------------------------------
! INICIAMOS LAS VARIABLES
       NODE = 1011
       NEWNODE = 1011 !NUMERO DE NODOS DESPUES DE APLICAR EL UMBRAL
       EDGES = 0
       SUM = 0
       KT = 5
       COUNT = 0

       ALLOCATE(DEGREE(NODE))
       ALLOCATE(NULL(NODE))
    
! LLAMAMOS EL ARCHIVO DE DEGREES
       OPEN(19, FILE = "/Users/arnedojoya/TFG/random/UL/3/degrees.txt")
       DO N = 1,NODE
        READ(19,*) I,J
        DEGREE(N) = J
       ENDDO
       CLOSE(19)

! UNA VEZ TENEMOS LOS DEGREES DE CADA NODO, APLICAMOS EL THRESHOLD Y GUARDAMOS EN EL VECTOR NULL LOS NODOS QUE ELIMINAMOS
       DO I = 1,NODE
        IF (DEGREE(I)<=KT) THEN
            DEGREE(I) = 0
            NULL(COUNT+1) = I-1 !NODOS QUE TIENEN DEGREE NULO
            COUNT = COUNT + 1
            NEWNODE = NEWNODE - 1 !NODOS TRAS APLICAR EL THRESHOLD
        ENDIF
       ENDDO

! CALCULAMOS EL VECTOR DEGREE DESPUÉS DE APLICAR EL THRESHOLD

!FILAS EDGELIST
       OPEN(47, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
       DO
        READ(47, *, IOSTAT = IOSTATUS) LINE ! LEEMOS LINEA DEL ARCHIVO
        IF (IOSTATUS /= 0) EXIT !SI LLEGA AL FINAL DEL ARCHIVO QUE SALGA
        EDGES = EDGES + 1
       ENDDO

       CLOSE(47)
    
! RESETEAMOS DEGREES
       DEGREE = 0

       OPEN(57, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
        
! CALCULAMOS EL VECTOR QUE CONTIENE LOS GRADOS APLICANDO EL THRESHOLD
       DO N = 1, EDGES
        READ(57, *, IOSTAT = IOSTATUS) I,J
        IF (IOSTATUS /= 0) THEN
            WRITE(*,*) 'ERROR AL LEER LA EDGE LIST EN LA FILA ', N
            WRITE(*,*) I, J
            EXIT
        ENDIF
        
! COMPROBAMOS SI LOS NODOS DEL VECTOR NULL ESTAN EN LA FILA QUE LEEMOS
! SI ESTA EN PRIMERA COLUMNA I
        EXCLUDE_I = .FALSE.
        DO K = 1, COUNT
            IF (I == NULL(K)) THEN
                EXCLUDE_I = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI ESTA EN PRIMERA COLUMNA J
        EXCLUDE_J = .FALSE.
        DO K = 1, COUNT
            IF (J == NULL(K)) THEN
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
       CLOSE(57)

       SUM = 0
       DO I = 1,NODE
        SUM = SUM + DEGREE(I)
       ENDDO
     
! AHORA CALCULAREMOS LA LISTA DE ADYACENCIA (VECTOR NEIGHBOURS) SIGUIENDO EL MISMO PROCEDIMIENTO QUE PARA EL DEGREE

       ALLOCATE(NEIGHBOURS(SUM))
       ALLOCATE(HEAD(NODE+1))
       ALLOCATE(INDEX(NODE+1))

       NEIGHBOURS = 0
       HEAD(1) = 1
       INDEX(1) = 1
 
       DO I = 2, NODE+1
        HEAD(I) = HEAD(I-1) + DEGREE(I-1)
        INDEX(I) = 1
       ENDDO

       OPEN(121, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
    
!CALCULAMOS EL VECTOR V RELLENANDO A MEDIDA QUE ENCONTRAMOS LOS VECINOS, SUMAMOS 1 YA QUE EL PRIMER NODO ES EL 0 Y LA POSICIÓN 0 NO EXISTE
       DO N = 1,EDGES
        READ(121,*,IOSTAT = IOSTATUS) I,J
        IF (IOSTATUS /= 0) EXIT
! COMPROBAMOS SI LOS NODOS DEL VECTOR NULL ESTAN EN LA FILA QUE LEEMOS
! SI ESTA EN PRIMERA COLUMNA I
        EXCLUDE_I = .FALSE.
        DO K = 1, COUNT
            IF (I == NULL(K)) THEN
                EXCLUDE_I = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI ESTA EN PRIMERA COLUMNA J
        EXCLUDE_J = .FALSE.
        DO K = 1, COUNT
            IF (J == NULL(K)) THEN
                EXCLUDE_J = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI NINGUNO DE LOS DOS NODOS ESTÁ EN EL VECTOR NULL SUMAMOS EL DEGREE
        IF (.NOT.EXCLUDE_I.AND..NOT.EXCLUDE_J) THEN
            NEIGHBOURS(HEAD(I+1)+INDEX(I+1)-1) = J
            INDEX(I+1) = INDEX(I+1) + 1
            NEIGHBOURS(HEAD(J+1)+INDEX(J+1)-1) = I
            INDEX(J+1) = INDEX(J+1) + 1
        ENDIF
       ENDDO
       CLOSE(121)


! ESCRIBIMOS NUEVA LISTA DE DEGREES
       OPEN(158, FILE = '/Users/arnedojoya/TFG/random/UL/3/degrees5.txt')
       DO N = 1,NODE
        WRITE(158,'(I0,3X,I0)') (N-1),DEGREE(N)
       ENDDO
       CLOSE(158)


! ESCRIBIMOS EN UN ARCHIVO NUEVA LISTA DE NEIGHBOURS
       OPEN(156, FILE = '/Users/arnedojoya/TFG/random/UL/3/neighbours5.txt')
       DO I = 1,SUM
        WRITE(156,'(I0)') NEIGHBOURS(I)
       ENDDO
       CLOSE(156)
    
       WRITE(*,*) 'K_T = 5'
       WRITE(*,*) 'Suma de los degrees: ', SUM
       WRITE(*,*) 'Nodos: ', NEWNODE
       WRITE(*,*)'<k> = ',DBLE(SUM)/NEWNODE

       DEALLOCATE(NEIGHBOURS)
       DEALLOCATE(DEGREE)
       DEALLOCATE(NULL)
     

! -----------------------KT = 10----------------------------------
! INICIAMOS LAS VARIABLES
       NODE = 1011
       NEWNODE = 1011 !NUMERO DE NODOS DESPUES DE APLICAR EL UMBRAL
       EDGES = 0
       SUM = 0
       KT = 10
       COUNT = 0

       ALLOCATE(DEGREE(NODE))
       ALLOCATE(NULL(NODE))
    
! LLAMAMOS EL ARCHIVO DE DEGREES
       OPEN(188, FILE = "/Users/arnedojoya/TFG/random/UL/3/degrees.txt")
       DO N = 1,NODE
        READ(188,*) I,J
        DEGREE(N) = J
       ENDDO
       CLOSE(188)

! UNA VEZ TENEMOS LOS DEGREES DE CADA NODO, APLICAMOS EL THRESHOLD Y GUARDAMOS EN EL VECTOR NULL LOS NODOS QUE ELIMINAMOS
       DO I = 1,NODE
        IF (DEGREE(I)<=KT) THEN
            DEGREE(I) = 0
            NULL(COUNT+1) = I-1 !NODOS QUE TIENEN DEGREE NULO
            COUNT = COUNT + 1
            NEWNODE = NEWNODE - 1 !NODOS TRAS APLICAR EL THRESHOLD
        ENDIF
       ENDDO

! CALCULAMOS EL VECTOR DEGREE DESPUÉS DE APLICAR EL THRESHOLD

!FILAS EDGELIST
       OPEN(208, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
       DO
        READ(208, *, IOSTAT = IOSTATUS) LINE ! LEEMOS LINEA DEL ARCHIVO
        IF (IOSTATUS /= 0) EXIT !SI LLEGA AL FINAL DEL ARCHIVO QUE SALGA
        EDGES = EDGES + 1
       ENDDO

       CLOSE(208)
    
! RESETEAMOS DEGREES
       DEGREE = 0

       OPEN(220, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
        
! CALCULAMOS EL VECTOR QUE CONTIENE LOS GRADOS APLICANDO EL THRESHOLD
       DO N = 1, EDGES
        READ(220, *, IOSTAT = IOSTATUS) I,J
        IF (IOSTATUS /= 0) THEN
            WRITE(*,*) 'ERROR AL LEER LA EDGE LIST EN LA FILA ', N
            WRITE(*,*) I, J
            EXIT
        ENDIF
        
! COMPROBAMOS SI LOS NODOS DEL VECTOR NULL ESTAN EN LA FILA QUE LEEMOS
! SI ESTA EN PRIMERA COLUMNA I
        EXCLUDE_I = .FALSE.
        DO K = 1, COUNT
            IF (I == NULL(K)) THEN
                EXCLUDE_I = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI ESTA EN PRIMERA COLUMNA J
        EXCLUDE_J = .FALSE.
        DO K = 1, COUNT
            IF (J == NULL(K)) THEN
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
       CLOSE(220)

       SUM = 0
       DO I = 1,NODE
        SUM = SUM + DEGREE(I)
       ENDDO
     
! AHORA CALCULAREMOS LA LISTA DE ADYACENCIA (VECTOR NEIGHBOURS) SIGUIENDO EL MISMO PROCEDIMIENTO QUE PARA EL DEGREE

       ALLOCATE(NEIGHBOURS(SUM))
       ALLOCATE(HEAD(NODE+1))
       ALLOCATE(INDEX(NODE+1))

       NEIGHBOURS = 0
       HEAD(1) = 1
       INDEX(1) = 1
 
       DO I = 2, NODE+1
        HEAD(I) = HEAD(I-1) + DEGREE(I-1)
        INDEX(I) = 1
       ENDDO

       OPEN(278, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
    
!CALCULAMOS EL VECTOR V RELLENANDO A MEDIDA QUE ENCONTRAMOS LOS VECINOS, SUMAMOS 1 YA QUE EL PRIMER NODO ES EL 0 Y LA POSICIÓN 0 NO EXISTE
       DO N = 1,EDGES
        READ(278,*,IOSTAT = IOSTATUS) I,J
        IF (IOSTATUS /= 0) EXIT
! COMPROBAMOS SI LOS NODOS DEL VECTOR NULL ESTAN EN LA FILA QUE LEEMOS
! SI ESTA EN PRIMERA COLUMNA I
        EXCLUDE_I = .FALSE.
        DO K = 1, COUNT
            IF (I == NULL(K)) THEN
                EXCLUDE_I = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI ESTA EN PRIMERA COLUMNA J
        EXCLUDE_J = .FALSE.
        DO K = 1, COUNT
            IF (J == NULL(K)) THEN
                EXCLUDE_J = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI NINGUNO DE LOS DOS NODOS ESTÁ EN EL VECTOR NULL SUMAMOS EL DEGREE
        IF (.NOT.EXCLUDE_I.AND..NOT.EXCLUDE_J) THEN
            NEIGHBOURS(HEAD(I+1)+INDEX(I+1)-1) = J
            INDEX(I+1) = INDEX(I+1) + 1
            NEIGHBOURS(HEAD(J+1)+INDEX(J+1)-1) = I
            INDEX(J+1) = INDEX(J+1) + 1
        ENDIF
       ENDDO
       CLOSE(278)


! ESCRIBIMOS NUEVA LISTA DE DEGREES
       OPEN(315, FILE = '/Users/arnedojoya/TFG/random/UL/3/degrees10.txt')
       DO N = 1,NODE
        WRITE(315,'(I0,3X,I0)') (N-1),DEGREE(N)
       ENDDO
       CLOSE(315)


! ESCRIBIMOS EN UN ARCHIVO NUEVA LISTA DE NEIGHBOURS
       OPEN(323, FILE = '/Users/arnedojoya/TFG/random/UL/3/neighbours10.txt')
       DO I = 1,SUM
        WRITE(323,'(I0)') NEIGHBOURS(I)
       ENDDO
       CLOSE(323)
    
       WRITE(*,*) 'K_T = 10'
       WRITE(*,*) 'Suma de los degrees: ', SUM
       WRITE(*,*) 'Nodos: ', NEWNODE
       WRITE(*,*)'<k> = ',DBLE(SUM)/NEWNODE

       DEALLOCATE(NEIGHBOURS)
       DEALLOCATE(DEGREE)
       DEALLOCATE(NULL)

! -----------------------KT = 15----------------------------------
! INICIAMOS LAS VARIABLES
       NODE = 1011
       NEWNODE = 1011 !NUMERO DE NODOS DESPUES DE APLICAR EL UMBRAL
       EDGES = 0
       SUM = 0
       KT = 15
       COUNT = 0

       ALLOCATE(DEGREE(NODE))
       ALLOCATE(NULL(NODE))
    
! LLAMAMOS EL ARCHIVO DE DEGREES
       OPEN(349, FILE = "/Users/arnedojoya/TFG/random/UL/3/degrees.txt")
       DO N = 1,NODE
        READ(349,*) I,J
        DEGREE(N) = J
       ENDDO
       CLOSE(349)

! UNA VEZ TENEMOS LOS DEGREES DE CADA NODO, APLICAMOS EL THRESHOLD Y GUARDAMOS EN EL VECTOR NULL LOS NODOS QUE ELIMINAMOS
       DO I = 1,NODE
        IF (DEGREE(I)<=KT) THEN
            DEGREE(I) = 0
            NULL(COUNT+1) = I-1 !NODOS QUE TIENEN DEGREE NULO
            COUNT = COUNT + 1
            NEWNODE = NEWNODE - 1 !NODOS TRAS APLICAR EL THRESHOLD
        ENDIF
       ENDDO

! CALCULAMOS EL VECTOR DEGREE DESPUÉS DE APLICAR EL THRESHOLD

!FILAS EDGELIST
       OPEN(369, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
       DO
        READ(369, *, IOSTAT = IOSTATUS) LINE ! LEEMOS LINEA DEL ARCHIVO
        IF (IOSTATUS /= 0) EXIT !SI LLEGA AL FINAL DEL ARCHIVO QUE SALGA
        EDGES = EDGES + 1
       ENDDO

       CLOSE(369)
    
! RESETEAMOS DEGREES
       DEGREE = 0

       OPEN(381, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
        
! CALCULAMOS EL VECTOR QUE CONTIENE LOS GRADOS APLICANDO EL THRESHOLD
       DO N = 1, EDGES
        READ(381, *, IOSTAT = IOSTATUS) I,J
        IF (IOSTATUS /= 0) THEN
            WRITE(*,*) 'ERROR AL LEER LA EDGE LIST EN LA FILA ', N
            WRITE(*,*) I, J
            EXIT
        ENDIF
        
! COMPROBAMOS SI LOS NODOS DEL VECTOR NULL ESTAN EN LA FILA QUE LEEMOS
! SI ESTA EN PRIMERA COLUMNA I
        EXCLUDE_I = .FALSE.
        DO K = 1, COUNT
            IF (I == NULL(K)) THEN
                EXCLUDE_I = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI ESTA EN PRIMERA COLUMNA J
        EXCLUDE_J = .FALSE.
        DO K = 1, COUNT
            IF (J == NULL(K)) THEN
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
       CLOSE(381)

       SUM = 0
       DO I = 1,NODE
        SUM = SUM + DEGREE(I)
       ENDDO
     
! AHORA CALCULAREMOS LA LISTA DE ADYACENCIA (VECTOR NEIGHBOURS) SIGUIENDO EL MISMO PROCEDIMIENTO QUE PARA EL DEGREE

       ALLOCATE(NEIGHBOURS(SUM))
       ALLOCATE(HEAD(NODE+1))
       ALLOCATE(INDEX(NODE+1))

       NEIGHBOURS = 0
       HEAD(1) = 1
       INDEX(1) = 1
 
       DO I = 2, NODE+1
        HEAD(I) = HEAD(I-1) + DEGREE(I-1)
        INDEX(I) = 1
       ENDDO

       OPEN(439, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
    
!CALCULAMOS EL VECTOR V RELLENANDO A MEDIDA QUE ENCONTRAMOS LOS VECINOS, SUMAMOS 1 YA QUE EL PRIMER NODO ES EL 0 Y LA POSICIÓN 0 NO EXISTE
       DO N = 1,EDGES
        READ(439,*,IOSTAT = IOSTATUS) I,J
        IF (IOSTATUS /= 0) EXIT
! COMPROBAMOS SI LOS NODOS DEL VECTOR NULL ESTAN EN LA FILA QUE LEEMOS
! SI ESTA EN PRIMERA COLUMNA I
        EXCLUDE_I = .FALSE.
        DO K = 1, COUNT
            IF (I == NULL(K)) THEN
                EXCLUDE_I = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI ESTA EN PRIMERA COLUMNA J
        EXCLUDE_J = .FALSE.
        DO K = 1, COUNT
            IF (J == NULL(K)) THEN
                EXCLUDE_J = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI NINGUNO DE LOS DOS NODOS ESTÁ EN EL VECTOR NULL SUMAMOS EL DEGREE
        IF (.NOT.EXCLUDE_I.AND..NOT.EXCLUDE_J) THEN
            NEIGHBOURS(HEAD(I+1)+INDEX(I+1)-1) = J
            INDEX(I+1) = INDEX(I+1) + 1
            NEIGHBOURS(HEAD(J+1)+INDEX(J+1)-1) = I
            INDEX(J+1) = INDEX(J+1) + 1
        ENDIF
       ENDDO
       CLOSE(439)


! ESCRIBIMOS NUEVA LISTA DE DEGREES
       OPEN(476, FILE = '/Users/arnedojoya/TFG/random/UL/3/degrees15.txt')
       DO N = 1,NODE
        WRITE(476,'(I0,3X,I0)') (N-1),DEGREE(N)
       ENDDO
       CLOSE(476)


! ESCRIBIMOS EN UN ARCHIVO NUEVA LISTA DE NEIGHBOURS
       OPEN(484, FILE = '/Users/arnedojoya/TFG/random/UL/3/neighbours15.txt')
       DO I = 1,SUM
        WRITE(484,'(I0)') NEIGHBOURS(I)
       ENDDO
       CLOSE(484)
    
       WRITE(*,*) 'K_T = 15'
       WRITE(*,*) 'Suma de los degrees: ', SUM
       WRITE(*,*) 'Nodos: ', NEWNODE
       WRITE(*,*)'<k> = ',DBLE(SUM)/NEWNODE

       DEALLOCATE(NEIGHBOURS)
       DEALLOCATE(DEGREE)
       DEALLOCATE(NULL)
! -----------------------KT = 20----------------------------------
! INICIAMOS LAS VARIABLES
       NODE = 1011
       NEWNODE = 1011 !NUMERO DE NODOS DESPUES DE APLICAR EL UMBRAL
       EDGES = 0
       SUM = 0
       KT = 20
       COUNT = 0

       ALLOCATE(DEGREE(NODE))
       ALLOCATE(NULL(NODE))
    
! LLAMAMOS EL ARCHIVO DE DEGREES
       OPEN(511, FILE = "/Users/arnedojoya/TFG/random/UL/3/degrees.txt")
       DO N = 1,NODE
        READ(511,*) I,J
        DEGREE(N) = J
       ENDDO
       CLOSE(511)

! UNA VEZ TENEMOS LOS DEGREES DE CADA NODO, APLICAMOS EL THRESHOLD Y GUARDAMOS EN EL VECTOR NULL LOS NODOS QUE ELIMINAMOS
       DO I = 1,NODE
        IF (DEGREE(I)<=KT) THEN
            DEGREE(I) = 0
            NULL(COUNT+1) = I-1 !NODOS QUE TIENEN DEGREE NULO
            COUNT = COUNT + 1
            NEWNODE = NEWNODE - 1 !NODOS TRAS APLICAR EL THRESHOLD
        ENDIF
       ENDDO

! CALCULAMOS EL VECTOR DEGREE DESPUÉS DE APLICAR EL THRESHOLD

!FILAS EDGELIST
       OPEN(531, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
       DO
        READ(531, *, IOSTAT = IOSTATUS) LINE ! LEEMOS LINEA DEL ARCHIVO
        IF (IOSTATUS /= 0) EXIT !SI LLEGA AL FINAL DEL ARCHIVO QUE SALGA
        EDGES = EDGES + 1
       ENDDO

       CLOSE(531)
    
! RESETEAMOS DEGREES
       DEGREE = 0

       OPEN(543, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
        
! CALCULAMOS EL VECTOR QUE CONTIENE LOS GRADOS APLICANDO EL THRESHOLD
       DO N = 1, EDGES
        READ(543, *, IOSTAT = IOSTATUS) I,J
        IF (IOSTATUS /= 0) THEN
            WRITE(*,*) 'ERROR AL LEER LA EDGE LIST EN LA FILA ', N
            WRITE(*,*) I, J
            EXIT
        ENDIF
        
! COMPROBAMOS SI LOS NODOS DEL VECTOR NULL ESTAN EN LA FILA QUE LEEMOS
! SI ESTA EN PRIMERA COLUMNA I
        EXCLUDE_I = .FALSE.
        DO K = 1, COUNT
            IF (I == NULL(K)) THEN
                EXCLUDE_I = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI ESTA EN PRIMERA COLUMNA J
        EXCLUDE_J = .FALSE.
        DO K = 1, COUNT
            IF (J == NULL(K)) THEN
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
       CLOSE(543)

       SUM = 0
       DO I = 1,NODE
        SUM = SUM + DEGREE(I)
       ENDDO
     
! AHORA CALCULAREMOS LA LISTA DE ADYACENCIA (VECTOR NEIGHBOURS) SIGUIENDO EL MISMO PROCEDIMIENTO QUE PARA EL DEGREE

       ALLOCATE(NEIGHBOURS(SUM))
       ALLOCATE(HEAD(NODE+1))
       ALLOCATE(INDEX(NODE+1))

       NEIGHBOURS = 0
       HEAD(1) = 1
       INDEX(1) = 1
 
       DO I = 2, NODE+1
        HEAD(I) = HEAD(I-1) + DEGREE(I-1)
        INDEX(I) = 1
       ENDDO

       OPEN(601, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
    
!CALCULAMOS EL VECTOR V RELLENANDO A MEDIDA QUE ENCONTRAMOS LOS VECINOS, SUMAMOS 1 YA QUE EL PRIMER NODO ES EL 0 Y LA POSICIÓN 0 NO EXISTE
       DO N = 1,EDGES
        READ(601,*,IOSTAT = IOSTATUS) I,J
        IF (IOSTATUS /= 0) EXIT
! COMPROBAMOS SI LOS NODOS DEL VECTOR NULL ESTAN EN LA FILA QUE LEEMOS
! SI ESTA EN PRIMERA COLUMNA I
        EXCLUDE_I = .FALSE.
        DO K = 1, COUNT
            IF (I == NULL(K)) THEN
                EXCLUDE_I = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI ESTA EN PRIMERA COLUMNA J
        EXCLUDE_J = .FALSE.
        DO K = 1, COUNT
            IF (J == NULL(K)) THEN
                EXCLUDE_J = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI NINGUNO DE LOS DOS NODOS ESTÁ EN EL VECTOR NULL SUMAMOS EL DEGREE
        IF (.NOT.EXCLUDE_I.AND..NOT.EXCLUDE_J) THEN
            NEIGHBOURS(HEAD(I+1)+INDEX(I+1)-1) = J
            INDEX(I+1) = INDEX(I+1) + 1
            NEIGHBOURS(HEAD(J+1)+INDEX(J+1)-1) = I
            INDEX(J+1) = INDEX(J+1) + 1
        ENDIF
       ENDDO
       CLOSE(601)


! ESCRIBIMOS NUEVA LISTA DE DEGREES
       OPEN(638, FILE = '/Users/arnedojoya/TFG/random/UL/3/degrees20.txt')
       DO N = 1,NODE
        WRITE(638,'(I0,3X,I0)') (N-1),DEGREE(N)
       ENDDO
       CLOSE(638)


! ESCRIBIMOS EN UN ARCHIVO NUEVA LISTA DE NEIGHBOURS
       OPEN(646, FILE = '/Users/arnedojoya/TFG/random/UL/3/neighbours20.txt')
       DO I = 1,SUM
        WRITE(646,'(I0)') NEIGHBOURS(I)
       ENDDO
       CLOSE(646)
    
       WRITE(*,*) 'K_T = 20'
       WRITE(*,*) 'Suma de los degrees: ', SUM
       WRITE(*,*) 'Nodos: ', NEWNODE
       WRITE(*,*)'<k> = ',DBLE(SUM)/NEWNODE

       DEALLOCATE(NEIGHBOURS)
       DEALLOCATE(DEGREE)
       DEALLOCATE(NULL)


! -----------------------KT = 25----------------------------------
! INICIAMOS LAS VARIABLES
       NODE = 1011
       NEWNODE = 1011 !NUMERO DE NODOS DESPUES DE APLICAR EL UMBRAL
       EDGES = 0
       SUM = 0
       KT = 25
       COUNT = 0

       ALLOCATE(DEGREE(NODE))
       ALLOCATE(NULL(NODE))
    
! LLAMAMOS EL ARCHIVO DE DEGREES
       OPEN(675, FILE = "/Users/arnedojoya/TFG/random/UL/3/degrees.txt")
       DO N = 1,NODE
        READ(675,*) I,J
        DEGREE(N) = J
       ENDDO
       CLOSE(675)

! UNA VEZ TENEMOS LOS DEGREES DE CADA NODO, APLICAMOS EL THRESHOLD Y GUARDAMOS EN EL VECTOR NULL LOS NODOS QUE ELIMINAMOS
       DO I = 1,NODE
        IF (DEGREE(I)<=KT) THEN
            DEGREE(I) = 0
            NULL(COUNT+1) = I-1 !NODOS QUE TIENEN DEGREE NULO
            COUNT = COUNT + 1
            NEWNODE = NEWNODE - 1 !NODOS TRAS APLICAR EL THRESHOLD
        ENDIF
       ENDDO

! CALCULAMOS EL VECTOR DEGREE DESPUÉS DE APLICAR EL THRESHOLD

!FILAS EDGELIST
       OPEN(695, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
       DO
        READ(695, *, IOSTAT = IOSTATUS) LINE ! LEEMOS LINEA DEL ARCHIVO
        IF (IOSTATUS /= 0) EXIT !SI LLEGA AL FINAL DEL ARCHIVO QUE SALGA
        EDGES = EDGES + 1
       ENDDO

       CLOSE(695)
    
! RESETEAMOS DEGREES
       DEGREE = 0

       OPEN(707, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
        
! CALCULAMOS EL VECTOR QUE CONTIENE LOS GRADOS APLICANDO EL THRESHOLD
       DO N = 1, EDGES
        READ(707, *, IOSTAT = IOSTATUS) I,J
        IF (IOSTATUS /= 0) THEN
            WRITE(*,*) 'ERROR AL LEER LA EDGE LIST EN LA FILA ', N
            WRITE(*,*) I, J
            EXIT
        ENDIF
        
! COMPROBAMOS SI LOS NODOS DEL VECTOR NULL ESTAN EN LA FILA QUE LEEMOS
! SI ESTA EN PRIMERA COLUMNA I
        EXCLUDE_I = .FALSE.
        DO K = 1, COUNT
            IF (I == NULL(K)) THEN
                EXCLUDE_I = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI ESTA EN PRIMERA COLUMNA J
        EXCLUDE_J = .FALSE.
        DO K = 1, COUNT
            IF (J == NULL(K)) THEN
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
       CLOSE(707)

       SUM = 0
       DO I = 1,NODE
        SUM = SUM + DEGREE(I)
       ENDDO
     
! AHORA CALCULAREMOS LA LISTA DE ADYACENCIA (VECTOR NEIGHBOURS) SIGUIENDO EL MISMO PROCEDIMIENTO QUE PARA EL DEGREE

       ALLOCATE(NEIGHBOURS(SUM))
       ALLOCATE(HEAD(NODE+1))
       ALLOCATE(INDEX(NODE+1))

       NEIGHBOURS = 0
       HEAD(1) = 1
       INDEX(1) = 1
 
       DO I = 2, NODE+1
        HEAD(I) = HEAD(I-1) + DEGREE(I-1)
        INDEX(I) = 1
       ENDDO

       OPEN(765, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
    
!CALCULAMOS EL VECTOR V RELLENANDO A MEDIDA QUE ENCONTRAMOS LOS VECINOS, SUMAMOS 1 YA QUE EL PRIMER NODO ES EL 0 Y LA POSICIÓN 0 NO EXISTE
       DO N = 1,EDGES
        READ(765,*,IOSTAT = IOSTATUS) I,J
        IF (IOSTATUS /= 0) EXIT
! COMPROBAMOS SI LOS NODOS DEL VECTOR NULL ESTAN EN LA FILA QUE LEEMOS
! SI ESTA EN PRIMERA COLUMNA I
        EXCLUDE_I = .FALSE.
        DO K = 1, COUNT
            IF (I == NULL(K)) THEN
                EXCLUDE_I = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI ESTA EN PRIMERA COLUMNA J
        EXCLUDE_J = .FALSE.
        DO K = 1, COUNT
            IF (J == NULL(K)) THEN
                EXCLUDE_J = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI NINGUNO DE LOS DOS NODOS ESTÁ EN EL VECTOR NULL SUMAMOS EL DEGREE
        IF (.NOT.EXCLUDE_I.AND..NOT.EXCLUDE_J) THEN
            NEIGHBOURS(HEAD(I+1)+INDEX(I+1)-1) = J
            INDEX(I+1) = INDEX(I+1) + 1
            NEIGHBOURS(HEAD(J+1)+INDEX(J+1)-1) = I
            INDEX(J+1) = INDEX(J+1) + 1
        ENDIF
       ENDDO
       CLOSE(765)


! ESCRIBIMOS NUEVA LISTA DE DEGREES
       OPEN(802, FILE = '/Users/arnedojoya/TFG/random/UL/3/degrees25.txt')
       DO N = 1,NODE
        WRITE(802,'(I0,3X,I0)') (N-1),DEGREE(N)
       ENDDO
       CLOSE(802)


! ESCRIBIMOS EN UN ARCHIVO NUEVA LISTA DE NEIGHBOURS
       OPEN(810, FILE = '/Users/arnedojoya/TFG/random/UL/3/neighbours25.txt')
       DO I = 1,SUM
        WRITE(810,'(I0)') NEIGHBOURS(I)
       ENDDO
       CLOSE(810)
    
       WRITE(*,*) 'K_T = 50'
       WRITE(*,*) 'Suma de los degrees: ', SUM
       WRITE(*,*) 'Nodos: ', NEWNODE
       WRITE(*,*)'<k> = ',DBLE(SUM)/NEWNODE

       DEALLOCATE(NEIGHBOURS)
       DEALLOCATE(DEGREE)
       DEALLOCATE(NULL)

! -----------------------KT = 30----------------------------------
! INICIAMOS LAS VARIABLES
       NODE = 1011
       NEWNODE = 1011 !NUMERO DE NODOS DESPUES DE APLICAR EL UMBRAL
       EDGES = 0
       SUM = 0
       KT = 30
       COUNT = 0

       ALLOCATE(DEGREE(NODE))
       ALLOCATE(NULL(NODE))
    
! LLAMAMOS EL ARCHIVO DE DEGREES
       OPEN(838, FILE = "/Users/arnedojoya/TFG/random/UL/3/degrees.txt")
       DO N = 1,NODE
        READ(838,*) I,J
        DEGREE(N) = J
       ENDDO
       CLOSE(838)

! UNA VEZ TENEMOS LOS DEGREES DE CADA NODO, APLICAMOS EL THRESHOLD Y GUARDAMOS EN EL VECTOR NULL LOS NODOS QUE ELIMINAMOS
       DO I = 1,NODE
        IF (DEGREE(I)<=KT) THEN
            DEGREE(I) = 0
            NULL(COUNT+1) = I-1 !NODOS QUE TIENEN DEGREE NULO
            COUNT = COUNT + 1
            NEWNODE = NEWNODE - 1 !NODOS TRAS APLICAR EL THRESHOLD
        ENDIF
       ENDDO

! CALCULAMOS EL VECTOR DEGREE DESPUÉS DE APLICAR EL THRESHOLD

!FILAS EDGELIST
       OPEN(858, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
       DO
        READ(858, *, IOSTAT = IOSTATUS) LINE ! LEEMOS LINEA DEL ARCHIVO
        IF (IOSTATUS /= 0) EXIT !SI LLEGA AL FINAL DEL ARCHIVO QUE SALGA
        EDGES = EDGES + 1
       ENDDO

       CLOSE(858)
    
! RESETEAMOS DEGREES
       DEGREE = 0

       OPEN(870, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
        
! CALCULAMOS EL VECTOR QUE CONTIENE LOS GRADOS APLICANDO EL THRESHOLD
       DO N = 1, EDGES
        READ(870, *, IOSTAT = IOSTATUS) I,J
        IF (IOSTATUS /= 0) THEN
            WRITE(*,*) 'ERROR AL LEER LA EDGE LIST EN LA FILA ', N
            WRITE(*,*) I, J
            EXIT
        ENDIF
        
! COMPROBAMOS SI LOS NODOS DEL VECTOR NULL ESTAN EN LA FILA QUE LEEMOS
! SI ESTA EN PRIMERA COLUMNA I
        EXCLUDE_I = .FALSE.
        DO K = 1, COUNT
            IF (I == NULL(K)) THEN
                EXCLUDE_I = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI ESTA EN PRIMERA COLUMNA J
        EXCLUDE_J = .FALSE.
        DO K = 1, COUNT
            IF (J == NULL(K)) THEN
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
       CLOSE(870)

       SUM = 0
       DO I = 1,NODE
        SUM = SUM + DEGREE(I)
       ENDDO
     
! AHORA CALCULAREMOS LA LISTA DE ADYACENCIA (VECTOR NEIGHBOURS) SIGUIENDO EL MISMO PROCEDIMIENTO QUE PARA EL DEGREE

       ALLOCATE(NEIGHBOURS(SUM))
       ALLOCATE(HEAD(NODE+1))
       ALLOCATE(INDEX(NODE+1))

       NEIGHBOURS = 0
       HEAD(1) = 1
       INDEX(1) = 1
 
       DO I = 2, NODE+1
        HEAD(I) = HEAD(I-1) + DEGREE(I-1)
        INDEX(I) = 1
       ENDDO

       OPEN(928, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
    
!CALCULAMOS EL VECTOR V RELLENANDO A MEDIDA QUE ENCONTRAMOS LOS VECINOS, SUMAMOS 1 YA QUE EL PRIMER NODO ES EL 0 Y LA POSICIÓN 0 NO EXISTE
       DO N = 1,EDGES
        READ(928,*,IOSTAT = IOSTATUS) I,J
        IF (IOSTATUS /= 0) EXIT
! COMPROBAMOS SI LOS NODOS DEL VECTOR NULL ESTAN EN LA FILA QUE LEEMOS
! SI ESTA EN PRIMERA COLUMNA I
        EXCLUDE_I = .FALSE.
        DO K = 1, COUNT
            IF (I == NULL(K)) THEN
                EXCLUDE_I = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI ESTA EN PRIMERA COLUMNA J
        EXCLUDE_J = .FALSE.
        DO K = 1, COUNT
            IF (J == NULL(K)) THEN
                EXCLUDE_J = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI NINGUNO DE LOS DOS NODOS ESTÁ EN EL VECTOR NULL SUMAMOS EL DEGREE
        IF (.NOT.EXCLUDE_I.AND..NOT.EXCLUDE_J) THEN
            NEIGHBOURS(HEAD(I+1)+INDEX(I+1)-1) = J
            INDEX(I+1) = INDEX(I+1) + 1
            NEIGHBOURS(HEAD(J+1)+INDEX(J+1)-1) = I
            INDEX(J+1) = INDEX(J+1) + 1
        ENDIF
       ENDDO
       CLOSE(928)


! ESCRIBIMOS NUEVA LISTA DE DEGREES
       OPEN(965, FILE = '/Users/arnedojoya/TFG/random/UL/3/degrees30.txt')
       DO N = 1,NODE
        WRITE(965,'(I0,3X,I0)') (N-1),DEGREE(N)
       ENDDO
       CLOSE(965)


! ESCRIBIMOS EN UN ARCHIVO NUEVA LISTA DE NEIGHBOURS
       OPEN(973, FILE = '/Users/arnedojoya/TFG/random/UL/3/neighbours30.txt')
       DO I = 1,SUM
        WRITE(973,'(I0)') NEIGHBOURS(I)
       ENDDO
       CLOSE(973)
    
       WRITE(*,*) 'K_T = 30'
       WRITE(*,*) 'Suma de los degrees: ', SUM
       WRITE(*,*) 'Nodos: ', NEWNODE
       WRITE(*,*)'<k> = ',DBLE(SUM)/NEWNODE

       DEALLOCATE(NEIGHBOURS)
       DEALLOCATE(DEGREE)
       DEALLOCATE(NULL)

! -----------------------KT = 35----------------------------------
! INICIAMOS LAS VARIABLES
       NODE = 1011
       NEWNODE = 1011 !NUMERO DE NODOS DESPUES DE APLICAR EL UMBRAL
       EDGES = 0
       SUM = 0
       KT = 35
       COUNT = 0

       ALLOCATE(DEGREE(NODE))
       ALLOCATE(NULL(NODE))
    
! LLAMAMOS EL ARCHIVO DE DEGREES
       OPEN(1001, FILE = "/Users/arnedojoya/TFG/random/UL/3/degrees.txt")
       DO N = 1,NODE
        READ(1001,*) I,J
        DEGREE(N) = J
       ENDDO
       CLOSE(1001)

! UNA VEZ TENEMOS LOS DEGREES DE CADA NODO, APLICAMOS EL THRESHOLD Y GUARDAMOS EN EL VECTOR NULL LOS NODOS QUE ELIMINAMOS
       DO I = 1,NODE
        IF (DEGREE(I)<=KT) THEN
            DEGREE(I) = 0
            NULL(COUNT+1) = I-1 !NODOS QUE TIENEN DEGREE NULO
            COUNT = COUNT + 1
            NEWNODE = NEWNODE - 1 !NODOS TRAS APLICAR EL THRESHOLD
        ENDIF
       ENDDO

! CALCULAMOS EL VECTOR DEGREE DESPUÉS DE APLICAR EL THRESHOLD

!FILAS EDGELIST
       OPEN(1021, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
       DO
        READ(1021, *, IOSTAT = IOSTATUS) LINE ! LEEMOS LINEA DEL ARCHIVO
        IF (IOSTATUS /= 0) EXIT !SI LLEGA AL FINAL DEL ARCHIVO QUE SALGA
        EDGES = EDGES + 1
       ENDDO

       CLOSE(1021)
    
! RESETEAMOS DEGREES
       DEGREE = 0

       OPEN(1033, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
        
! CALCULAMOS EL VECTOR QUE CONTIENE LOS GRADOS APLICANDO EL THRESHOLD
       DO N = 1, EDGES
        READ(1033, *, IOSTAT = IOSTATUS) I,J
        IF (IOSTATUS /= 0) THEN
            WRITE(*,*) 'ERROR AL LEER LA EDGE LIST EN LA FILA ', N
            WRITE(*,*) I, J
            EXIT
        ENDIF
        
! COMPROBAMOS SI LOS NODOS DEL VECTOR NULL ESTAN EN LA FILA QUE LEEMOS
! SI ESTA EN PRIMERA COLUMNA I
        EXCLUDE_I = .FALSE.
        DO K = 1, COUNT
            IF (I == NULL(K)) THEN
                EXCLUDE_I = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI ESTA EN PRIMERA COLUMNA J
        EXCLUDE_J = .FALSE.
        DO K = 1, COUNT
            IF (J == NULL(K)) THEN
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
       CLOSE(1033)

       SUM = 0
       DO I = 1,NODE
        SUM = SUM + DEGREE(I)
       ENDDO
     
! AHORA CALCULAREMOS LA LISTA DE ADYACENCIA (VECTOR NEIGHBOURS) SIGUIENDO EL MISMO PROCEDIMIENTO QUE PARA EL DEGREE

       ALLOCATE(NEIGHBOURS(SUM))
       ALLOCATE(HEAD(NODE+1))
       ALLOCATE(INDEX(NODE+1))

       NEIGHBOURS = 0
       HEAD(1) = 1
       INDEX(1) = 1
 
       DO I = 2, NODE+1
        HEAD(I) = HEAD(I-1) + DEGREE(I-1)
        INDEX(I) = 1
       ENDDO

       OPEN(1091, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
    
!CALCULAMOS EL VECTOR V RELLENANDO A MEDIDA QUE ENCONTRAMOS LOS VECINOS, SUMAMOS 1 YA QUE EL PRIMER NODO ES EL 0 Y LA POSICIÓN 0 NO EXISTE
       DO N = 1,EDGES
        READ(1091,*,IOSTAT = IOSTATUS) I,J
        IF (IOSTATUS /= 0) EXIT
! COMPROBAMOS SI LOS NODOS DEL VECTOR NULL ESTAN EN LA FILA QUE LEEMOS
! SI ESTA EN PRIMERA COLUMNA I
        EXCLUDE_I = .FALSE.
        DO K = 1, COUNT
            IF (I == NULL(K)) THEN
                EXCLUDE_I = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI ESTA EN PRIMERA COLUMNA J
        EXCLUDE_J = .FALSE.
        DO K = 1, COUNT
            IF (J == NULL(K)) THEN
                EXCLUDE_J = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI NINGUNO DE LOS DOS NODOS ESTÁ EN EL VECTOR NULL SUMAMOS EL DEGREE
        IF (.NOT.EXCLUDE_I.AND..NOT.EXCLUDE_J) THEN
            NEIGHBOURS(HEAD(I+1)+INDEX(I+1)-1) = J
            INDEX(I+1) = INDEX(I+1) + 1
            NEIGHBOURS(HEAD(J+1)+INDEX(J+1)-1) = I
            INDEX(J+1) = INDEX(J+1) + 1
        ENDIF
       ENDDO
       CLOSE(1091)


! ESCRIBIMOS NUEVA LISTA DE DEGREES
       OPEN(1128, FILE = '/Users/arnedojoya/TFG/random/UL/3/degrees35.txt')
       DO N = 1,NODE
        WRITE(1128,'(I0,3X,I0)') (N-1),DEGREE(N)
       ENDDO
       CLOSE(1128)


! ESCRIBIMOS EN UN ARCHIVO NUEVA LISTA DE NEIGHBOURS
       OPEN(1136, FILE = '/Users/arnedojoya/TFG/random/UL/3/neighbours35.txt')
       DO I = 1,SUM
        WRITE(1136,'(I0)') NEIGHBOURS(I)
       ENDDO
       CLOSE(1136)
    
       WRITE(*,*) 'K_T = 35'
       WRITE(*,*) 'Suma de los degrees: ', SUM
       WRITE(*,*) 'Nodos: ', NEWNODE
       WRITE(*,*)'<k> = ',DBLE(SUM)/NEWNODE

       DEALLOCATE(NEIGHBOURS)
       DEALLOCATE(DEGREE)
       DEALLOCATE(NULL)



! -----------------------KT = 40----------------------------------
! INICIAMOS LAS VARIABLES
       NODE = 1011
       NEWNODE = 1011 !NUMERO DE NODOS DESPUES DE APLICAR EL UMBRAL
       EDGES = 0
       SUM = 0
       KT = 40
       COUNT = 0

       ALLOCATE(DEGREE(NODE))
       ALLOCATE(NULL(NODE))
    
! LLAMAMOS EL ARCHIVO DE DEGREES
       OPEN(1165, FILE = "/Users/arnedojoya/TFG/random/UL/3/degrees.txt")
       DO N = 1,NODE
        READ(1165,*) I,J
        DEGREE(N) = J
       ENDDO
       CLOSE(1165)

! UNA VEZ TENEMOS LOS DEGREES DE CADA NODO, APLICAMOS EL THRESHOLD Y GUARDAMOS EN EL VECTOR NULL LOS NODOS QUE ELIMINAMOS
       DO I = 1,NODE
        IF (DEGREE(I)<=KT) THEN
            DEGREE(I) = 0
            NULL(COUNT+1) = I-1!NODOS QUE TIENEN DEGREE NULO
            COUNT = COUNT + 1
            NEWNODE = NEWNODE - 1 !NODOS TRAS APLICAR EL THRESHOLD
        ENDIF
       ENDDO

! CALCULAMOS EL VECTOR DEGREE DESPUÉS DE APLICAR EL THRESHOLD

!FILAS EDGELIST
       OPEN(1185, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
       DO
        READ(1185, *, IOSTAT = IOSTATUS) LINE ! LEEMOS LINEA DEL ARCHIVO
        IF (IOSTATUS /= 0) EXIT !SI LLEGA AL FINAL DEL ARCHIVO QUE SALGA
        EDGES = EDGES + 1
       ENDDO

       CLOSE(1185)
    
! RESETEAMOS DEGREES
       DEGREE = 0

       OPEN(1197, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
        
! CALCULAMOS EL VECTOR QUE CONTIENE LOS GRADOS APLICANDO EL THRESHOLD
       DO N = 1, EDGES
        READ(1197, *, IOSTAT = IOSTATUS) I,J
        IF (IOSTATUS /= 0) THEN
            WRITE(*,*) 'ERROR AL LEER LA EDGE LIST EN LA FILA ', N
            WRITE(*,*) I, J
            EXIT
        ENDIF
        
! COMPROBAMOS SI LOS NODOS DEL VECTOR NULL ESTAN EN LA FILA QUE LEEMOS
! SI ESTA EN PRIMERA COLUMNA I
        EXCLUDE_I = .FALSE.
        DO K = 1, COUNT
            IF (I == NULL(K)) THEN
                EXCLUDE_I = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI ESTA EN PRIMERA COLUMNA J
        EXCLUDE_J = .FALSE.
        DO K = 1, COUNT
            IF (J == NULL(K)) THEN
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
       CLOSE(1197)

       SUM = 0
       DO I = 1,NODE
        SUM = SUM + DEGREE(I)
       ENDDO
     
! AHORA CALCULAREMOS LA LISTA DE ADYACENCIA (VECTOR NEIGHBOURS) SIGUIENDO EL MISMO PROCEDIMIENTO QUE PARA EL DEGREE

       ALLOCATE(NEIGHBOURS(SUM))
       ALLOCATE(HEAD(NODE+1))
       ALLOCATE(INDEX(NODE+1))

       NEIGHBOURS = 0
       HEAD(1) = 1
       INDEX(1) = 1
 
       DO I = 2, NODE+1
        HEAD(I) = HEAD(I-1) + DEGREE(I-1)
        INDEX(I) = 1
       ENDDO

       OPEN(1255, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
    
!CALCULAMOS EL VECTOR V RELLENANDO A MEDIDA QUE ENCONTRAMOS LOS VECINOS, SUMAMOS 1 YA QUE EL PRIMER NODO ES EL 0 Y LA POSICIÓN 0 NO EXISTE
       DO N = 1,EDGES
        READ(1255,*,IOSTAT = IOSTATUS) I,J
        IF (IOSTATUS /= 0) EXIT
! COMPROBAMOS SI LOS NODOS DEL VECTOR NULL ESTAN EN LA FILA QUE LEEMOS
! SI ESTA EN PRIMERA COLUMNA I
        EXCLUDE_I = .FALSE.
        DO K = 1, COUNT
            IF (I == NULL(K)) THEN
                EXCLUDE_I = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI ESTA EN PRIMERA COLUMNA J
        EXCLUDE_J = .FALSE.
        DO K = 1, COUNT
            IF (J == NULL(K)) THEN
                EXCLUDE_J = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI NINGUNO DE LOS DOS NODOS ESTÁ EN EL VECTOR NULL SUMAMOS EL DEGREE
        IF (.NOT.EXCLUDE_I.AND..NOT.EXCLUDE_J) THEN
            NEIGHBOURS(HEAD(I+1)+INDEX(I+1)-1) = J
            INDEX(I+1) = INDEX(I+1) + 1
            NEIGHBOURS(HEAD(J+1)+INDEX(J+1)-1) = I
            INDEX(J+1) = INDEX(J+1) + 1
        ENDIF
       ENDDO
       CLOSE(1255)


! ESCRIBIMOS NUEVA LISTA DE DEGREES
       OPEN(1292, FILE = '/Users/arnedojoya/TFG/random/UL/3/degrees40.txt')
       DO N = 1,NODE
        WRITE(1292,'(I0,3X,I0)') (N-1),DEGREE(N)
       ENDDO
       CLOSE(1292)


! ESCRIBIMOS EN UN ARCHIVO NUEVA LISTA DE NEIGHBOURS
       OPEN(1300, FILE = '/Users/arnedojoya/TFG/random/UL/3/neighbours40.txt')
       DO I = 1,SUM
        WRITE(1300,'(I0)') NEIGHBOURS(I)
       ENDDO
       CLOSE(1300)
    
       WRITE(*,*) 'K_T = 40'
       WRITE(*,*) 'Suma de los degrees: ', SUM
       WRITE(*,*) 'Nodos: ', NEWNODE
       WRITE(*,*)'<k> = ',DBLE(SUM)/NEWNODE

       DEALLOCATE(NEIGHBOURS)
       DEALLOCATE(DEGREE)
       DEALLOCATE(NULL)

! -----------------------------------------------KT = 45 --------------------------------------------------------
! INICIAMOS LAS VARIABLES
       NODE = 1011
       NEWNODE = 1011 !NUMERO DE NODOS DESPUES DE APLICAR EL UMBRAL
       EDGES = 0
       SUM = 0
       KT = 45
       COUNT = 0

       ALLOCATE(DEGREE(NODE))
       ALLOCATE(NULL(NODE))
    
! LLAMAMOS EL ARCHIVO DE DEGREES
       OPEN(1328, FILE = "/Users/arnedojoya/TFG/random/UL/3/degrees.txt")
       DO N = 1,NODE
        READ(1328,*) I,J
        DEGREE(N) = J
       ENDDO
       CLOSE(1328)

! UNA VEZ TENEMOS LOS DEGREES DE CADA NODO, APLICAMOS EL THRESHOLD Y GUARDAMOS EN EL VECTOR NULL LOS NODOS QUE ELIMINAMOS
       DO I = 1,NODE
        IF (DEGREE(I)<=KT) THEN
            DEGREE(I) = 0
            NULL(COUNT+1) = I-1 !NODOS QUE TIENEN DEGREE NULO
            COUNT = COUNT + 1
            NEWNODE = NEWNODE - 1 !NODOS TRAS APLICAR EL THRESHOLD
        ENDIF
       ENDDO

! CALCULAMOS EL VECTOR DEGREE DESPUÉS DE APLICAR EL THRESHOLD

!FILAS EDGELIST
       OPEN(1348, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
       DO
        READ(1348, *, IOSTAT = IOSTATUS) LINE ! LEEMOS LINEA DEL ARCHIVO
        IF (IOSTATUS /= 0) EXIT !SI LLEGA AL FINAL DEL ARCHIVO QUE SALGA
        EDGES = EDGES + 1
       ENDDO

       CLOSE(1348)
    
! RESETEAMOS DEGREES
       DEGREE = 0

       OPEN(1360, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
        
! CALCULAMOS EL VECTOR QUE CONTIENE LOS GRADOS APLICANDO EL THRESHOLD
       DO N = 1, EDGES
        READ(1360, *, IOSTAT = IOSTATUS) I,J
        IF (IOSTATUS /= 0) THEN
            WRITE(*,*) 'ERROR AL LEER LA EDGE LIST EN LA FILA ', N
            WRITE(*,*) I, J
            EXIT
        ENDIF
        
! COMPROBAMOS SI LOS NODOS DEL VECTOR NULL ESTAN EN LA FILA QUE LEEMOS
! SI ESTA EN PRIMERA COLUMNA I
        EXCLUDE_I = .FALSE.
        DO K = 1, COUNT
            IF (I == NULL(K)) THEN
                EXCLUDE_I = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI ESTA EN PRIMERA COLUMNA J
        EXCLUDE_J = .FALSE.
        DO K = 1, COUNT
            IF (J == NULL(K)) THEN
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
       CLOSE(1360)

       SUM = 0
       DO I = 1,NODE
        SUM = SUM + DEGREE(I)
       ENDDO
     
! AHORA CALCULAREMOS LA LISTA DE ADYACENCIA (VECTOR NEIGHBOURS) SIGUIENDO EL MISMO PROCEDIMIENTO QUE PARA EL DEGREE

       ALLOCATE(NEIGHBOURS(SUM))
       ALLOCATE(HEAD(NODE+1))
       ALLOCATE(INDEX(NODE+1))

       NEIGHBOURS = 0
       HEAD(1) = 1
       INDEX(1) = 1
 
       DO I = 2, NODE+1
        HEAD(I) = HEAD(I-1) + DEGREE(I-1)
        INDEX(I) = 1
       ENDDO

       OPEN(1418, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
    
!CALCULAMOS EL VECTOR V RELLENANDO A MEDIDA QUE ENCONTRAMOS LOS VECINOS, SUMAMOS 1 YA QUE EL PRIMER NODO ES EL 0 Y LA POSICIÓN 0 NO EXISTE
       DO N = 1,EDGES
        READ(1418,*,IOSTAT = IOSTATUS) I,J
        IF (IOSTATUS /= 0) EXIT
! COMPROBAMOS SI LOS NODOS DEL VECTOR NULL ESTAN EN LA FILA QUE LEEMOS
! SI ESTA EN PRIMERA COLUMNA I
        EXCLUDE_I = .FALSE.
        DO K = 1, COUNT
            IF (I == NULL(K)) THEN
                EXCLUDE_I = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI ESTA EN PRIMERA COLUMNA J
        EXCLUDE_J = .FALSE.
        DO K = 1, COUNT
            IF (J == NULL(K)) THEN
                EXCLUDE_J = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI NINGUNO DE LOS DOS NODOS ESTÁ EN EL VECTOR NULL SUMAMOS EL DEGREE
        IF (.NOT.EXCLUDE_I.AND..NOT.EXCLUDE_J) THEN
            NEIGHBOURS(HEAD(I+1)+INDEX(I+1)-1) = J
            INDEX(I+1) = INDEX(I+1) + 1
            NEIGHBOURS(HEAD(J+1)+INDEX(J+1)-1) = I
            INDEX(J+1) = INDEX(J+1) + 1
        ENDIF
       ENDDO
       CLOSE(1418)


! ESCRIBIMOS NUEVA LISTA DE DEGREES
       OPEN(1455, FILE = '/Users/arnedojoya/TFG/random/UL/3/degrees45.txt')
       DO N = 1,NODE
        WRITE(1455,'(I0,3X,I0)') (N-1),DEGREE(N)
       ENDDO
       CLOSE(1455)


! ESCRIBIMOS EN UN ARCHIVO NUEVA LISTA DE NEIGHBOURS
       OPEN(1463, FILE = '/Users/arnedojoya/TFG/random/UL/3/neighbours45.txt')
       DO I = 1,SUM
        WRITE(1463,'(I0)') NEIGHBOURS(I)
       ENDDO
       CLOSE(1463)
    
       WRITE(*,*) 'K_T = 45'
       WRITE(*,*) 'Suma de los degrees: ', SUM
       WRITE(*,*) 'Nodos: ', NEWNODE
       WRITE(*,*)'<k> = ',DBLE(SUM)/NEWNODE

       DEALLOCATE(NEIGHBOURS)
       DEALLOCATE(DEGREE)
       DEALLOCATE(NULL)

! -----------------------------------------------KT = 50 --------------------------------------------------------
! INICIAMOS LAS VARIABLES
       NODE = 1011
       NEWNODE = 1011 !NUMERO DE NODOS DESPUES DE APLICAR EL UMBRAL
       EDGES = 0
       SUM = 0
       KT = 50
       COUNT = 0

       ALLOCATE(DEGREE(NODE))
       ALLOCATE(NULL(NODE))
    
! LLAMAMOS EL ARCHIVO DE DEGREES
       OPEN(1491, FILE = "/Users/arnedojoya/TFG/random/UL/3/degrees.txt")
       DO N = 1,NODE
        READ(1491,*) I,J
        DEGREE(N) = J
       ENDDO
       CLOSE(1491)

! UNA VEZ TENEMOS LOS DEGREES DE CADA NODO, APLICAMOS EL THRESHOLD Y GUARDAMOS EN EL VECTOR NULL LOS NODOS QUE ELIMINAMOS
       DO I = 1,NODE
        IF (DEGREE(I)<=KT) THEN
            DEGREE(I) = 0
            NULL(COUNT+1) = I-1 !NODOS QUE TIENEN DEGREE NULO
            COUNT = COUNT + 1
            NEWNODE = NEWNODE - 1 !NODOS TRAS APLICAR EL THRESHOLD
        ENDIF
       ENDDO

! CALCULAMOS EL VECTOR DEGREE DESPUÉS DE APLICAR EL THRESHOLD

!FILAS EDGELIST
       OPEN(1511, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
       DO
        READ(1511, *, IOSTAT = IOSTATUS) LINE ! LEEMOS LINEA DEL ARCHIVO
        IF (IOSTATUS /= 0) EXIT !SI LLEGA AL FINAL DEL ARCHIVO QUE SALGA
        EDGES = EDGES + 1
       ENDDO

       CLOSE(1511)
    
! RESETEAMOS DEGREES
       DEGREE = 0

       OPEN(1523, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
        
! CALCULAMOS EL VECTOR QUE CONTIENE LOS GRADOS APLICANDO EL THRESHOLD
       DO N = 1, EDGES
        READ(1523, *, IOSTAT = IOSTATUS) I,J
        IF (IOSTATUS /= 0) THEN
            WRITE(*,*) 'ERROR AL LEER LA EDGE LIST EN LA FILA ', N
            WRITE(*,*) I, J
            EXIT
        ENDIF
        
! COMPROBAMOS SI LOS NODOS DEL VECTOR NULL ESTAN EN LA FILA QUE LEEMOS
! SI ESTA EN PRIMERA COLUMNA I
        EXCLUDE_I = .FALSE.
        DO K = 1, COUNT
            IF (I == NULL(K)) THEN
                EXCLUDE_I = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI ESTA EN PRIMERA COLUMNA J
        EXCLUDE_J = .FALSE.
        DO K = 1, COUNT
            IF (J == NULL(K)) THEN
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
       CLOSE(1523)

       SUM = 0
       DO I = 1,NODE
        SUM = SUM + DEGREE(I)
       ENDDO
     
! AHORA CALCULAREMOS LA LISTA DE ADYACENCIA (VECTOR NEIGHBOURS) SIGUIENDO EL MISMO PROCEDIMIENTO QUE PARA EL DEGREE

       ALLOCATE(NEIGHBOURS(SUM))
       ALLOCATE(HEAD(NODE+1))
       ALLOCATE(INDEX(NODE+1))

       NEIGHBOURS = 0
       HEAD(1) = 1
       INDEX(1) = 1
 
       DO I = 2, NODE+1
        HEAD(I) = HEAD(I-1) + DEGREE(I-1)
        INDEX(I) = 1
       ENDDO

       OPEN(1581, FILE = "/Users/arnedojoya/TFG/random/UL/3/random_edgelist.txt")
    
!CALCULAMOS EL VECTOR V RELLENANDO A MEDIDA QUE ENCONTRAMOS LOS VECINOS, SUMAMOS 1 YA QUE EL PRIMER NODO ES EL 0 Y LA POSICIÓN 0 NO EXISTE
       DO N = 1,EDGES
        READ(1581,*,IOSTAT = IOSTATUS) I,J
        IF (IOSTATUS /= 0) EXIT
! COMPROBAMOS SI LOS NODOS DEL VECTOR NULL ESTAN EN LA FILA QUE LEEMOS
! SI ESTA EN PRIMERA COLUMNA I
        EXCLUDE_I = .FALSE.
        DO K = 1, COUNT
            IF (I == NULL(K)) THEN
                EXCLUDE_I = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI ESTA EN PRIMERA COLUMNA J
        EXCLUDE_J = .FALSE.
        DO K = 1, COUNT
            IF (J == NULL(K)) THEN
                EXCLUDE_J = .TRUE.
                EXIT
            ENDIF
        ENDDO

! SI NINGUNO DE LOS DOS NODOS ESTÁ EN EL VECTOR NULL SUMAMOS EL DEGREE
        IF (.NOT.EXCLUDE_I.AND..NOT.EXCLUDE_J) THEN
            NEIGHBOURS(HEAD(I+1)+INDEX(I+1)-1) = J
            INDEX(I+1) = INDEX(I+1) + 1
            NEIGHBOURS(HEAD(J+1)+INDEX(J+1)-1) = I
            INDEX(J+1) = INDEX(J+1) + 1
        ENDIF
       ENDDO
       CLOSE(1581)


! ESCRIBIMOS NUEVA LISTA DE DEGREES
       OPEN(1618, FILE = '/Users/arnedojoya/TFG/random/UL/3/degrees50.txt')
       DO N = 1,NODE
        WRITE(1618,'(I0,3X,I0)') (N-1),DEGREE(N)
       ENDDO
       CLOSE(1618)


! ESCRIBIMOS EN UN ARCHIVO NUEVA LISTA DE NEIGHBOURS
       OPEN(1626, FILE = '/Users/arnedojoya/TFG/random/UL/3/neighbours50.txt')
       DO I = 1,SUM
        WRITE(1626,'(I0)') NEIGHBOURS(I)
       ENDDO
       CLOSE(1626)
    
       WRITE(*,*) 'K_T = 50'
       WRITE(*,*) 'Suma de los degrees: ', SUM
       WRITE(*,*) 'Nodos: ', NEWNODE
       WRITE(*,*)'<k> = ',DBLE(SUM)/NEWNODE

       DEALLOCATE(NEIGHBOURS)
       DEALLOCATE(DEGREE)
       DEALLOCATE(NULL)
     
     


       END PROGRAM DEGREE_THRESHOLDING_RENORMALIZATION


