      SUBROUTINE tmklmedf(RO,CO,RC,CC,TLIMIT,A,GR,GC,GBEST,NREPS)
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION TIMEA, TIMEB, TLIMIT
      DOUBLE PRECISION S1
      INTEGER NR(100),A(RO,CO),NC(100),S(9300),
     1   CMEM(9300),RMEM(9300),GR(9300),GC(9300),
     1   SUM0(100,100),SUM1(100,100),TRIAL0(100,100),
     1   SUMBEST0(100,100),SUMBEST1(100,100),CENTR(100,9300),
     1   CENTC(9300,20),MED(20,20),TRIAL1(100,100)
C
C  ################################################################
C     LOCAL SEARCH HEURISTIC
C        1. OBJECT TRANSFERS ONLY -- NO EXCHANGES
C        2. RANDOMLY CHOOSE K1/K2 CENTROIDS -- ASSIGN EACH OBJECT TO NEAREST
C  ################################################################
C
      call fseedi()
      CALL CPU_TIME(TIMEA)
      NREPS = 0
      GBEST = 999999
      NBEST = 0
 1000 NREPS = NREPS + 1
        DO I = 1,RC
          NR(I) = 0
        END DO
        DO I = 1,CC
          NC(I) = 0
        END DO
        DO I = 1,RO
          RMEM(I) = 0
        END DO
        DO I = 1,CO
          CMEM(I) = 0
        END DO
        DO I = 1,RC
          DO J = 1,CC
            SUM0(I,J) = 0
            SUM1(I,J) = 0
          END DO
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
        DO I = 1,RC
          DO J = 1,CC
            SUM0(I,J) = 0
            SUM1(I,J) = 0
          END DO
        END DO
        DO I = 1,RO
          K = RMEM(I)
          DO J = 1,CO
            L = CMEM(J)
            IF(A(I,J).LE.0) THEN
              SUM0(K,L) = SUM0(K,L) + 1
            ELSE
              SUM1(K,L) = SUM1(K,L) + 1
            END IF
          END DO
        END DO
        Z = 0
        DO K = 1,RC
          DO L = 1,CC
            IF(SUM1(K,L).LE.SUM0(K,L)) THEN
              Z = Z + SUM1(K,L)
            ELSE
              Z = Z + SUM0(K,L)
            END IF
          END DO
        END DO
C
C       ##################################################
C              BEGIN TWO-MODE K-MEDIAN CLUSTERING
C       ##################################################
C
        ZBEST = Z
 1060   DO K = 1,RC
          DO L = 1,CC
            MED(K,L) = 1
            IF(SUM0(K,L).GT.SUM1(K,L)) MED(K,L) = 0
          END DO
        END DO
        DO K = 1,RC
          NR(K) = 0
        END DO
        KMAXD = 0
        DO I = 1,RO
          KBEST = 9999999
          DO K = 1,RC
            KSET = 0
            DO J = 1,CO
              KSET = KSET + ABS(A(I,J)-MED(K,CMEM(J)))
            END DO
            IF(KSET.LT.KBEST) THEN
              KBEST = KSET
              KSEL = K
            END IF
          END DO
          IF(KBEST.GT.KMAXD) THEN
            KMAXD = KBEST
            ISEL = I
          END IF
          RMEM(I) = KSEL
          NR(KSEL) = NR(KSEL) + 1
        END DO
 1059   IEMPTY = 0
        DO K = 1,RC
          IF(NR(K).EQ.0) THEN
            IEMPTY = 1
            KSEL = K
          END IF
        END DO
        IF(IEMPTY.EQ.0) GO TO 1058
        NR(RMEM(ISEL)) = NR(RMEM(ISEL)) - 1
        RMEM(ISEL) = KSEL
        NR(KSEL) = 1
 1058   DO I = 1,RC
          DO J = 1,CC
            SUM0(I,J) = 0
            SUM1(I,J) = 0
          END DO
        END DO
        DO I = 1,RO
          KR = RMEM(I)
          DO J = 1,CO
            KC = CMEM(J)
            IF(A(I,J).EQ.0) THEN
              SUM0(KR,KC) = SUM0(KR,KC)+1
            ELSE
              SUM1(KR,KC) = SUM1(KR,KC)+1
            END IF
          END DO
        END DO
        DO K = 1,RC
          DO L = 1,CC
            MED(K,L) = 1
            IF(SUM0(K,L).GT.SUM1(K,L)) MED(K,L) = 0
          END DO
        END DO
        DO L = 1,CC
          NC(L) = 0
        END DO
        LMAXD = 0
        DO J = 1,CO
          LBEST = 9999999
          DO L = 1,CC
            LSET = 0
            DO I = 1,RO
              LSET = LSET + ABS(A(I,J)-MED(RMEM(I),L))
            END DO
            IF(LSET.LT.LBEST) THEN
              LBEST = LSET
              LSEL = L
            END IF
          END DO
          IF(LBEST.GT.LMAXD) THEN
            LMAXD = LBEST
            JSEL = J
          END IF
          CMEM(J) = LSEL
          NC(LSEL) = NC(LSEL) + 1
        END DO
 1069   IEMPTY = 0
        DO L = 1,CC
          IF(NC(L).EQ.0) THEN
            IEMPTY = 1
            LSEL = L
          END IF
        END DO
        IF(IEMPTY.EQ.0) GO TO 1068
        NC(CMEM(JSEL)) = NC(CMEM(JSEL)) - 1
        CMEM(JSEL) = LSEL
        NC(LSEL) = 1
 1068   DO I = 1,RC
          DO J = 1,CC
            SUM0(I,J) = 0
            SUM1(I,J) = 0
          END DO
        END DO
        DO I = 1,RO
          KR = RMEM(I)
          DO J = 1,CO
            KC = CMEM(J)
            IF(A(I,J).EQ.0) THEN
              SUM0(KR,KC) = SUM0(KR,KC)+1
            ELSE
              SUM1(KR,KC) = SUM1(KR,KC)+1
            END IF
          END DO
        END DO
        Z = 0
        DO K = 1,RC
          DO L = 1,CC
            IF(SUM1(K,L).LE.SUM0(K,L)) THEN
              Z = Z + SUM1(K,L)
            ELSE
              Z = Z + SUM0(K,L)
            END IF
          END DO
        END DO
        IF(Z.LT.ZBEST) THEN
          ZBEST = Z
          GO TO 1060
        END IF
        go to 1075
C
C       ##################################################
C              BEGIN RELOCATION HEURISTIC
C       ##################################################
C
 1100   DELTA = 0
C
C       ##################################################
C              EVALUATE ALL ROW OBJECT RELOCATIONS
C       ##################################################
C
        DO 600 I = 1,RO
          K = RMEM(I)
          IF(NR(K).EQ.1) GO TO 600
          DO 601 KK = 1,RC
            IF(KK.EQ.K) GO TO 601
            DO K1 = 1,RC
              DO K2 = 1,CC
                TRIAL0(K1,K2) =  SUM0(K1,K2)
                TRIAL1(K1,K2) =  SUM1(K1,K2)
               END DO
            END DO
            DO 602 J = 1,CO
              L = CMEM(J)
              IF(A(I,J).EQ.0) THEN
                TRIAL0(K,L) = TRIAL0(K,L)-1
                TRIAL0(KK,L) = TRIAL0(KK,L)+1
              ELSE
                TRIAL1(K,L) = TRIAL1(K,L)-1
                TRIAL1(KK,L) = TRIAL1(KK,L)+1
              END IF
  602       CONTINUE
            ZTRIAL = 0
            DO K1 = 1,RC
              DO L1 = 1,CC
                IF(TRIAL1(K1,L1).LE.TRIAL0(K1,L1)) THEN
                   ZTRIAL = ZTRIAL + TRIAL1(K1,L1)
                ELSE
                   ZTRIAL = ZTRIAL + TRIAL0(K1,L1)
                END IF
              END DO
            END DO
            IF(ZTRIAL-Z.LT.DELTA) THEN
              BESTMOV = 1
              DELTA = ZTRIAL-Z
              ISEL = I
              KSEL = KK
              DO K1 = 1,RC
                DO L1 = 1,CC
                  SUMBEST0(K1,L1) = TRIAL0(K1,L1)
                  SUMBEST1(K1,L1) = TRIAL1(K1,L1)
                END DO
              END DO
            END IF
 601      CONTINUE
 600    CONTINUE
C
C
C       ##################################################
C              EVALUATE ALL COLUMN OBJECT RELOCATIONS
C       ##################################################
C
        DO 700 J = 1,CO
          L = CMEM(J)
          IF(NC(L).EQ.1) GO TO 700
          DO 701 LL = 1,CC
            IF(LL.EQ.L) GO TO 701
            DO K1 = 1,RC
              DO K2 = 1,CC
                TRIAL0(K1,K2) = SUM0(K1,K2)
                TRIAL1(K1,K2) =  SUM1(K1,K2)
               END DO
            END DO
            DO 702 I = 1,RO
              K = RMEM(I)
              IF(A(I,J).EQ.0) THEN
                TRIAL0(K,L) = TRIAL0(K,L)-1
                TRIAL0(K,LL) = TRIAL0(K,LL)+1
              ELSE
                TRIAL1(K,L) = TRIAL1(K,L)-1
                TRIAL1(K,LL) = TRIAL1(K,LL)+1
              END IF
  702       CONTINUE
            ZTRIAL = 0
            DO K1 = 1,RC
              DO L1 = 1,CC
                IF(TRIAL1(K1,L1).LE.TRIAL0(K1,L1)) THEN
                   ZTRIAL = ZTRIAL + TRIAL1(K1,L1)
                ELSE
                   ZTRIAL = ZTRIAL + TRIAL0(K1,L1)
                END IF
              END DO
            END DO
            IF(ZTRIAL-Z.LT.DELTA) THEN
              BESTMOV = 2
              DELTA = ZTRIAL-Z
              JSEL = J
              LSEL = LL
              DO K1 = 1,RC
                DO L1 = 1,CC
                  SUMBEST0(K1,L1) = TRIAL0(K1,L1)
                  SUMBEST1(K1,L1) = TRIAL1(K1,L1)
                END DO
              END DO
            END IF
 701      CONTINUE
 700    CONTINUE
C
C       ##################################################
C              IDENTIFY THE BEST NEIGHBORHOOD MOVE
C       ##################################################
C
        IF(DELTA.LT.0) THEN
          IF(BESTMOV.EQ.1) THEN
            Z = Z + DELTA
            K=RMEM(ISEL)
            NR(K) = NR(K) - 1
            NR(KSEL) = NR(KSEL)+1
            RMEM(ISEL) = KSEL
            DO K1 = 1,RC
              DO K2 = 1,CC
                SUM0(K1,K2) = SUMBEST0(K1,K2)
                SUM1(K1,K2) = SUMBEST1(K1,K2)
              END DO
            END DO
          ELSE
            Z = Z + DELTA
            L=CMEM(JSEL)
            NC(L) = NC(L) - 1
            NC(LSEL) = NC(LSEL)+1
            CMEM(JSEL) = LSEL
            DO K1 = 1,RC
              DO K2 = 1,CC
                SUM0(K1,K2) = SUMBEST0(K1,K2)
                SUM1(K1,K2) = SUMBEST1(K1,K2)
              END DO
            END DO
          END IF
          GO TO 1100
        END IF
C
 1075   IF(Z.LT.GBEST) THEN
          GBEST = Z
          NBEST = 1
          DO I = 1,RO
            GR(I) = RMEM(I)
          END DO
          DO I = 1,CO
            GC(I) = CMEM(I)
          END DO
        ELSEIF(Z.EQ.GBEST) THEN
          NBEST = NBEST + 1
        END IF
      CALL CPU_TIME(TIMEB)
      IF(TIMEB-TIMEA.LT.TLIMIT) GO TO 1000
C
      DO K = 1,RC
        NR(K) = 0
      END DO
      DO l = 1,CC
        NC(L) = 0
      END DO
      DO I = 1,RO
        NR(GR(I)) = NR(GR(I))+1
      END DO
      DO J = 1,CO
        NC(GC(j)) = NC(GC(j))+1
      END DO
      call fseedo()  
      END
