      SUBROUTINE tmklmf(RO,CO,RC,CC,TLIMIT,A,RBEST,CBEST,VAF,NREPS)
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION TIMEA, TIMEB, A(RO,CO),V(100,100),TLIMIT,VAF,
     1                 C(2000,100),D(2000,100),ZMIN,ZBEST,Z,DMIN,DCOM,
     1                 CMIN,ASUM,ASSE,CENTR(100,2000),CENTC(2000,100)
      DOUBLE PRECISION S1
      INTEGER NR(100),NC(100),CMEM(2000),RMEM(2000),
     1        CBEST(CO),RBEST(RO)
C
C  MULTISTART TWO-MODE KL-MEANS (WCSS CRITERION) BLOCK PLACEMENTS **UNKNOWN**
      call fseedi()
      CALL CPU_TIME(TIMEA)
      NREPS = 0
      NUMOPT = 1
      UNIQUE = 1
      ZBEST = 999999999.0D0
      ASUM = 0.0D0
      ASSE = 0.0D0
      DO I = 1,RO
        DO J = 1,CO
          ASUM = ASUM + A(I,J)
        END DO
      END DO
      ASUM = ASUM/DBLE(RO*CO)
      DO I = 1,RO
        DO J = 1,CO
          ASSE = ASSE + (A(I,J)-ASUM)**2
        END DO
      END DO
 1000 NREPS = NREPS + 1
        DO I = 1,RC
          NR(I) = 0
          DO J = 1,CC
            NC(J) = 0
            V(I,J) = 0.0D0
          END DO
        END DO
        DO I = 1,RO
          RMEM(I) = 0
        END DO
        DO J = 1,CO
          CMEM(J) = 0
        END DO
C
C  ###############################################################
C              GREEDY INITIALIZATION PROCEDURE - PART 1:
C
C    RANDOMLY SELECT K1 OF THE ROW OBJECTS TO PUT ONE OBJECT IN
C    EACH OF THE FIRST K1 ROW CLUSTERS.  THEN, DO THE SAME FOR
C    TO OBTAIN K2 COLUMN OBJECTS TO PUT ONE OBJECT IN EACH OF
C    THE FIRST K2 COLUMN CLUSTERS.  THEN COMPUTE THE OBJECTIVE
C    VALUE.
C  ################################################################
C
        DO K = 1,RC
 98       CALL randr(S1)
          ISEL = S1 * FLOAT(RO) + 1.
          IF(ISEL.GT.RO) ISEL = RO
          IF(RMEM(ISEL).GT.0) GO TO 98
          RMEM(ISEL) = K
          NR(K) = NR(K) + 1
          DO J = 1,CO
            CENTR(K,J) = A(ISEL,J)
          END DO
        END DO
        DO L = 1,CC
 97       CALL randr(S1)
          JSEL = S1 * FLOAT(CO) + 1.
          IF(JSEL.GT.CO) JSEL = CO
          IF(CMEM(JSEL).GT.0) GO TO 97
          CMEM(JSEL) = L
          NC(L) = NC(L) + 1
          DO I = 1,RO
            CENTC(I,L) = A(I,JSEL)
          END DO
        END DO
        DO 300 I = 1,RO
          IF(RMEM(I).NE.0) GO TO 300
          DMIN = 9999999
          DO K = 1,RC
            DCOM = 0
            DO J = 1,CO
              DCOM = DCOM + (A(I,J)-CENTR(K,J))**2
            END DO
            IF(DCOM.LT.DMIN) THEN
              DMIN = DCOM
              KSEL = K
            END IF
          END DO
          RMEM(I) = KSEL
          NR(KSEL) = NR(KSEL) + 1
 300    CONTINUE
        DO 301 J = 1,CO
          IF(CMEM(J).NE.0) GO TO 301
          DMIN = 9999999
          DO K = 1,CC
            DCOM = 0
            DO I = 1,RO
              DCOM = DCOM + (A(I,J)-CENTC(I,K))**2
            END DO
            IF(DCOM.LT.DMIN) THEN
              DMIN = DCOM
              KSEL = K
            END IF
          END DO
          CMEM(J) = KSEL
          NC(KSEL) = NC(KSEL) + 1
 301    CONTINUE
C
                              ! COMPUTE CENTROIDS
        Z = 0.0D0
        DO I = 1,RO
          KR = RMEM(I)
          DO J = 1,CO
            KC = CMEM(J)
            V(KR,KC) = V(KR,KC)+A(I,J)
          END DO
        END DO
        DO K = 1,RC
          DO L = 1,CC
            IF(NR(K).GT.0.AND.NC(L).GT.0) THEN
              V(K,L) = V(K,L)/DBLE(NR(K)*NC(L))
            END IF
          END DO
        END DO
        DO I = 1,RO
          K = RMEM(I)
          DO J = 1,CO
            L = CMEM(J)
            Z = Z + (A(I,J)-V(K,L))**2
          END DO
        END DO
        ZMIN = Z
C
 450    TRIG = 0
        DO I = 1,RO
          DO K = 1,RC
            C(I,K) = 0.0D0
          END DO
        END DO
        DO 451 I = 1,RO
          DO 452 K = 1,RC
            DO 453 J = 1,CO
              L = CMEM(J)
              C(I,K) = C(I,K) + (A(I,J) - V(K,L))**2
 453        CONTINUE
 452      CONTINUE
 451    CONTINUE
        DO K = 1,RC
          NR(K) = 0
        END DO
        DO 463 I = 1,RO
          CMIN = 9.9D+12
          DO 464 K = 1,RC
            IF(C(I,K).LT.CMIN) THEN
              CMIN = C(I,K)
              KSEL = K
            END IF
 464      CONTINUE
          RMEM(I) = KSEL
          NR(KSEL) = NR(KSEL) + 1
 463    CONTINUE
C
        DO K = 1,RC
          DO L = 1,CC
            V(K,L) = 0.0D0
          END DO
        END DO
        Z = 0.0D0
        DO I = 1,RO
          KR = RMEM(I)
          DO J = 1,CO
            KC = CMEM(J)
            V(KR,KC) = V(KR,KC)+A(I,J)
          END DO
        END DO
        DO K = 1,RC
          DO L = 1,CC
            IF(NR(K).GT.0.AND.NC(L).GT.0) THEN
              V(K,L) = V(K,L)/DBLE(NR(K)*NC(L))
            END IF
          END DO
        END DO
        DO I = 1,RO
          K = RMEM(I)
          DO J = 1,CO
            L = CMEM(J)
            Z = Z + (A(I,J)-V(K,L))**2
          END DO
        END DO
C
        DO J = 1,CO
          DO L = 1,CC
            D(J,L) = 0.0D0
          END DO
        END DO
        DO 455 J = 1,CO
          DO 456 L = 1,CC
            DO 457 I= 1,RO
              K = RMEM(I)
              D(J,L) = D(J,L) + (A(I,J) - V(K,L))**2
 457        CONTINUE
 456      CONTINUE
 455    CONTINUE
        DO L = 1,CC
          NC(L) = 0
        END DO
        DO 467 J = 1,CO
          CMIN = 9.9D+12
          DO 468 L = 1,CC
            IF(D(J,L).LT.CMIN) THEN
              CMIN = D(J,L)
              LSEL = L
            END IF
 468      CONTINUE
          CMEM(J) = LSEL
          NC(LSEL) = NC(LSEL) + 1
 467    CONTINUE
C
        DO K = 1,RC
          DO L = 1,CC
            V(K,L) = 0.0D0
          END DO
        END DO
        Z = 0.0D0
        DO I = 1,RO
          KR = RMEM(I)
          DO J = 1,CO
            KC = CMEM(J)
            V(KR,KC) = V(KR,KC)+A(I,J)
          END DO
        END DO
        DO K = 1,RC
          DO L = 1,CC
            IF(NR(K).GT.0.AND.NC(L).GT.0) THEN
              V(K,L) = V(K,L)/DBLE(NR(K)*NC(L))
            END IF
          END DO
        END DO
        DO I = 1,RO
          K = RMEM(I)
          DO J = 1,CO
            L = CMEM(J)
            Z = Z + (A(I,J)-V(K,L))**2
          END DO
        END DO
C
        IF(Z.LT.ZMIN-.000001) THEN
          ZMIN = Z
          GO TO 450
        END IF
C
        IF(ZMIN.LT.ZBEST-.000001) THEN
          DO I = 1,RO
            RBEST(I) = RMEM(I)
          END DO
          DO J = 1,CO
            CBEST(J) = CMEM(J)
          END DO
          ZBEST = ZMIN
          UNIQUE = 1
          NUMOPT = 1
        ELSEIF(ZMIN.GT.ZBEST-.000001.AND.ZMIN.LT.ZBEST+.000001) THEN
          NUMOPT = NUMOPT+1
        END IF
      CALL CPU_TIME(TIMEB)
      IF(TIMEB-TIMEA.LT.TLIMIT) GO TO 1000
      VAF = 1-ZBEST/ASSE
      call fseedo() 
      END
