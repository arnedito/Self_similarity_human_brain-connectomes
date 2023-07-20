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
       OPEN(19, FILE = "/Users/arnedojoya/TFG/thresholdren/UL/3/degreesUL.txt")
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
       OPEN(47, FILE = "/Users/arnedojoya/TFG/thresholdren/UL/3/UL_nobrainstem_3_layer_0_edgelist.txt")
       DO
        READ(47, *, IOSTAT = IOSTATUS) LINE ! LEEMOS LINEA DEL ARCHIVO
        IF (IOSTATUS /= 0) EXIT !SI LLEGA AL FINAL DEL ARCHIVO QUE SALGA
        EDGES = EDGES + 1
       ENDDO

       CLOSE(47)
    
! RESETEAMOS DEGREES
       DEGREE = 0

       OPEN(57, FILE = "/Users/arnedojoya/TFG/thresholdren/UL/3/UL_nobrainstem_3_layer_0_edgelist.txt")
        
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
            IF (I/=J) THEN
                DEGREE(I+1) = DEGREE(I+1) + 1 !SUMAMOS 1 YA QUE EL PRIMER NODO ES 0
            ENDIF
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

       OPEN(121, FILE = "/Users/arnedojoya/TFG/thresholdren/UL/3/UL_nobrainstem_3_layer_0_edgelist.txt")
    
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
            IF (I/=J) THEN
                NEIGHBOURS(HEAD(I+1)+INDEX(I+1)-1) = J
                INDEX(I+1) = INDEX(I+1) + 1
            ENDIF
        ENDIF
       ENDDO
       CLOSE(121)


! ESCRIBIMOS NUEVA LISTA DE DEGREES
       OPEN(158, FILE = '/Users/arnedojoya/TFG/thresholdren/UL/3/degrees5UL.txt')
       DO N = 1,NODE
        WRITE(158,'(I0,3X,I0)') (N-1),DEGREE(N)
       ENDDO
       CLOSE(158)


! ESCRIBIMOS EN UN ARCHIVO NUEVA LISTA DE NEIGHBOURS
       OPEN(156, FILE = '/Users/arnedojoya/TFG/thresholdren/UL/3/neighbours5UL.txt')
       DO I = 1,SUM
        WRITE(156,'(I0)') NEIGHBOURS(I)
       ENDDO
       CLOSE(156)
    
       WRITE(*,*) 'K_T = 5'
       WRITE(*,*) 'Suma de los degrees: ', SUM
       WRITE(*,*) 'Nodos: ', NEWNODE
       WRITE(*,*)'<k> = ',DBLE(SUM)/NEWNODE

       END PROGRAM DEGREE_THRESHOLDING_RENORMALIZATION

