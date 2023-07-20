! PROGRAM TO OBTAIN THE DEGREE LIST AND NEIGHBOUR LIST SEPARATELY GIVEN AN UNDIRECTED AND UNWEIGHTED NETWORK'S EDGELIST
       PROGRAM EDGELISTHCP
       IMPLICIT NONE
       INTEGER*8 I,J,N,NODE,SUM,K,EDGES,LINE,IOSTATUS
       INTEGER, ALLOCATABLE :: DEGREE(:), NEIGHBORS(:)
       INTEGER, POINTER :: HEAD(:), COUNT(:)
       
! INTRODUCE THE NUMBER OF NODES OF THE NETWORK
       NODE = 1014


       EDGES = 0 !INITIALIZE TO CALCULATE NUMBER OF LINKS
       SUM = 0 !SUM OF DEGREES

       ALLOCATE(DEGREE(NODE))
       DEGREE = 0
    
! INPUT EDGELIST
       OPEN(19, FILE = "HCP_nobrainstem_22_layer_0_edgelist.txt")

! CALCULATE NUMBER OF LINKS
       DO
        READ(19, *, IOSTAT = IOSTATUS) LINE ! LEEMOS LINEA DEL ARCHIVO
        IF (IOSTATUS /= 0) EXIT !SI LLEGA AL FINAL DEL ARCHIVO QUE SALGA
        EDGES = EDGES + 1
       ENDDO
       CLOSE(19)


! INPUT EDGELIST, AGAIN TO AVOID ERRORS IF NOT CLOSED EARLIER
       OPEN(30, FILE = "HCP_nobrainstem_22_layer_0_edgelist.txt")
        
! FILL DEGREES
       DO N = 1, EDGES
        READ(30, *, IOSTAT = IOSTATUS) I,J
        IF (IOSTATUS /= 0) THEN
            WRITE(*,*) 'ERROR AL LEER LA EDGE LIST EN LA FILA ', N
            WRITE(*,*) I, J
            EXIT
        ENDIF
        DEGREE(I+1) = DEGREE(I+1) + 1 !SUMAMOS 1 YA QUE EL PRIMER NODO ES 0
        DEGREE(J+1) = DEGREE(J+1) + 1
       ENDDO
       CLOSE(30)

       OPEN(48, FILE = "degrees0.txt") !WRITE DEGREES IN A FILE
       DO N = 1, NODE
        SUM = SUM + DEGREE(N)
        WRITE(48,'(I0,3X,I0)') (N-1), DEGREE(N)
       ENDDO

       CLOSE(48)

       WRITE(*,*) 'Suma de los degrees: ', SUM
       WRITE(*,*)'<k> = ',DBLE(SUM)/NODE

!--------------------------------------------------------------------------------------
! NEIGHBOR VECTOR

       ALLOCATE(NEIGHBORS(SUM))
       ALLOCATE(HEAD(NODE+1))
       ALLOCATE(COUNT(NODE+1))
     
       NEIGHBOURS = 0 !INITIALIZE
    
       HEAD(1) = 1 !INITIAL POSITION OF THE FIRST NEIGHBOR OF EACH NODE
       COUNT(1) = 1 !COUNTER TO SET EACH NEIGHBOR IN THE CORRECT PLACE
        
       DO I = 2, NODE+1
        HEAD(I) = HEAD(I-1) + DEGREE(I-1)
        COUNT(I) = 1
       ENDDO

! INPUT EDGELIST
       OPEN(81, FILE = "HCP_nobrainstem_22_layer_0_edgelist.txt")

       DO K = 1,EDGES
        READ(81,*,IOSTAT = IOSTATUS) I,J
        IF (IOSTATUS /= 0) EXIT
        NEIGHBORS(HEAD(I+1)+COUNT(I+1)-1) = J
        COUNT(I+1) = COUNT(I+1) + 1
        NEIGHBORS(HEAD(J+1)+COUNT(J+1)-1) = I
        COUNT(J+1) = COUNT(J+1) + 1
       ENDDO
       CLOSE(81)


       OPEN(101, FILE = 'neighbors0.txt') !WRITE NEIGHBORS IN A FILE
       DO K = 1,SUM
        WRITE(101,'(I0)') NEIGHBORS(K)
       ENDDO
       CLOSE(101)

       END PROGRAM EDGELISTHCP

    

