      SUBROUTINE omkmf(RO,RC,IDIAG,TLIMIT,A,RBEST,ZBEST,VAF,NREPS)
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION TIMEA, TIMEB, A(RO,RO),V(100,100),TLIMIT,VAF,
     1                 C(5000,100),ZMIN,ZBEST,Z,DMIN,DCOM,
     1                 CMIN,ASUM,ASSE,CENTR(100,5000),CENTC(100,5000)
      DOUBLE PRECISION S1
      INTEGER NR(100),RMEM(5000),RBEST(RO)
C
C  MULTISTART ONE-MODE KL-MEANS (WCSS CRITERION) BLOCK PLACEMENTS **UNKNOWN**
      CALL CPU_TIME(TIMEA)
      NREPS = 0
      NUMOPT = 1
      UNIQUE = 1
      ZBEST = 9999999999.0D0
      ASUM = 0.0D0
      ASSE = 0.0D0
      DO I = 1,RO
        DO J = 1,RO
          IF(I.NE.J.OR.IDIAG.NE.0) ASUM = ASUM + A(I,J)
        END DO
      END DO
      IF(IDIAG.NE.0) ASUM = ASUM/DFLOAT(RO*RO)
      IF(IDIAG.EQ.0) ASUM = ASUM/DFLOAT(RO*RO-RO)
      DO I = 1,RO
        DO J = 1,RO
          IF(I.NE.J.OR.IDIAG.NE.0) ASSE = ASSE + (A(I,J)-ASUM)**2
        END DO
      END DO
 1000 NREPS = NREPS + 1
        DO I = 1,RC
          NR(I) = 0
          DO J = 1,RC
            V(I,J) = 0.0D0
          END DO
        END DO
        DO I = 1,RO
          RMEM(I) = 0
        END DO
C
C  ###############################################################
C              GREEDY INITIALIZATION PROCEDURE - PART 1:
C
C    RANDOMLY SELECT K OF THE OBJECTS TO PUT ONE OBJECT IN
C    EACH OF THE ROW CLUSTERS. THEN ASSIGN AND COMPUTE THE OBJECTIVE
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
          DO J = 1,RO
            CENTR(K,J) = A(ISEL,J)
            CENTC(K,J) = A(J,ISEL)
          END DO
        END DO
        DO 300 I = 1,RO
          IF(RMEM(I).NE.0) GO TO 300
          DMIN = 9999999
          DO K = 1,RC
            DCOM = 0
            DO J = 1,RO
              DCOM = DCOM+(A(I,J)-CENTR(K,J))**2+(A(J,I)-CENTC(K,J))**2
            END DO
            IF(DCOM.LT.DMIN) THEN
              DMIN = DCOM
              KSEL = K
            END IF
          END DO
          RMEM(I) = KSEL
          NR(KSEL) = NR(KSEL) + 1
 300    CONTINUE
C
                              ! COMPUTE CENTROIDS
        Z = 0.0D0
        DO I = 1,RO
          KR = RMEM(I)
          DO J = 1,RO
            KC = RMEM(J)
            IF(I.NE.J.OR.IDIAG.NE.0) V(KR,KC) = V(KR,KC)+A(I,J)
          END DO
        END DO
        DO K = 1,RC
          DO L = 1,RC
            IF(NR(K).GT.0.AND.NR(L).GT.0) THEN
              IF(K.NE.L.OR.IDIAG.NE.0) THEN
                V(K,L) = V(K,L)/DFLOAT(NR(K)*NR(L))
              ELSEIF(NR(K).GT.1) THEN
                V(K,L) = V(K,L)/DFLOAT(NR(K)*NR(K)-NR(K))
              END IF  
            END IF
          END DO
        END DO
        DO I = 1,RO
          K = RMEM(I)
          DO J = 1,RO
            L = RMEM(J)
            IF(I.NE.J.OR.IDIAG.NE.0) Z = Z + (A(I,J)-V(K,L))**2
          END DO
        END DO
        ZMIN = Z
C
 450    TRIG = 0
        DO 451 I = 1,RO
          DO 452 K = 1,RC
            C(I,K) = 0
            DO 453 J = 1,RO
              IF(J.EQ.I.AND.IDIAG.EQ.0) GO TO 453
              L = RMEM(J)
              C(I,K) = C(I,K) + (A(I,J)-V(K,L))**2 + (A(J,I)-V(L,K))**2
              IF(J.EQ.I.AND.IDIAG.NE.0) C(I,K)=C(I,K)-(A(J,I)-V(L,K))**2
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
          DO L = 1,RC
            V(K,L) = 0.0D0
          END DO
        END DO
        Z = 0.0D0
        DO I = 1,RO
          KR = RMEM(I)
          DO J = 1,RO
            KC = RMEM(J)
            IF(I.NE.J.OR.IDIAG.NE.0) V(KR,KC) = V(KR,KC)+A(I,J)
          END DO
        END DO
        DO K = 1,RC
          DO L = 1,RC
            IF(NR(K).GT.0.AND.NR(L).GT.0) THEN
              IF(K.NE.L.OR.IDIAG.NE.0) THEN
                V(K,L) = V(K,L)/DFLOAT(NR(K)*NR(L))
              ELSEIF(NR(K).GT.1) THEN
                V(K,L) = V(K,L)/DFLOAT(NR(K)*NR(K)-NR(K))
              END IF  
            END IF
          END DO
        END DO
        DO I = 1,RO
          K = RMEM(I)
          DO J = 1,RO
            L = RMEM(J)
            IF(I.NE.J.OR.IDIAG.NE.0) Z = Z + (A(I,J)-V(K,L))**2
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
          ZBEST = ZMIN
          UNIQUE = 1
          NUMOPT = 1
        ELSEIF(ZMIN.GT.ZBEST-.000001.AND.ZMIN.LT.ZBEST+.000001) THEN
          NUMOPT = NUMOPT+1
        END IF
      CALL CPU_TIME(TIMEB)
      IF(TIMEB-TIMEA.LT.TLIMIT) GO TO 1000
      VAF = 1-ZBEST/ASSE
      END
