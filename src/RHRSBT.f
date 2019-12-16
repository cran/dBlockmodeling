      SUBROUTINE rhrsbtf(N, C, TLIMIT, OBJVAL, A, EBEST, NREPS)
      IMPLICIT INTEGER(A-Z)
      REAL TIMEA,TIMEB,TLIMIT,OBJVAL
      DOUBLE PRECISION S1
      INTEGER A(N,N),E(N),M(C),EBEST(N),MBEST(C),
     1  POS(C,C),NEG(C,C),PTRY(C,C),NTRY(C,C),ETRY(N),CENTR(C,N)
C
C *******************************************************************
C  RELOCATION ALGORITHM FOR RELAXED STRUCTURAL BALANCE
C *******************************************************************
C
      call fseedi()
      CALL CPU_TIME(TIMEA)  
      SMAX= 0
      DO I = 1,N
        DO J = 1,N
          SMAX = SMAX + ABS(A(I,J))
        END DO
      END DO
      DO I = 1,N
        ETRY(I) = 0
      END DO
      ZBEST = -999999
      NREPS = 0
 6000 NREPS = NREPS + 1
        DO K = 1,C
          M(K) = 0
        END DO
        DO I = 1,N
          E(I) = 0
        END DO
        DO K = 1,C
 98       CALL randr(S1)
          ISEL = S1 * FLOAT(N) + 1.
          IF(ISEL.GT.N) ISEL = N
          IF(E(ISEL).GT.0) GO TO 98
          DO L = 1,K-1
            IDIFF = 0
            DO J = 1,N
              IDIFF = IDIFF + (A(ISEL,J)-CENTR(L,J))**2
            END DO
            IF(IDIFF.EQ.0) GO TO 98
          END DO
          E(ISEL) = K
          M(K) = M(K) + 1
          DO J = 1,N
            CENTR(K,J) = A(ISEL,J)
          END DO
        END DO
        DO 300 I = 1,N
          IF(E(I).NE.0) GO TO 300
          DMIN = 9999999
          DO K = 1,C
            DCOM = 0
            DO J = 1,N
              DCOM = DCOM + (A(I,J)-CENTR(K,J))**2
            END DO
            IF(DCOM.LT.DMIN) THEN
              DMIN = DCOM
              KSEL = K
            END IF
          END DO
          E(I) = KSEL
          M(KSEL) = M(KSEL) + 1
 300    CONTINUE
        Z = 0
        DO K = 1,C
          DO L = 1,C
            POS(K,L) = 0
            NEG(K,L) = 0
          END DO
        END DO
        DO I = 1,N-1
          DO J = I+1,N
            KI = E(I)
            KJ = E(J)
            IF(A(I,J).GT.0) POS(KI,KJ) = POS(KI,KJ)+A(I,J)
            IF(A(I,J).LT.0) NEG(KI,KJ) = NEG(KI,KJ)-A(I,J)
            IF(A(J,I).GT.0) POS(KJ,KI) = POS(KJ,KI)+A(J,I)
            IF(A(J,I).LT.0) NEG(KJ,KI) = NEG(KJ,KI)-A(J,I)
          END DO
        END DO
        DO K = 1,C
          DO L = 1,C
            IF(POS(K,L).GE.NEG(K,L)) THEN
              Z = Z + POS(K,L)
            ELSE
              Z = Z + NEG(K,L)
            END IF
          END DO
        END DO
 1000   IFLAG = 0
        DO 100 I = 1,N
          KI = E(I)
          IF(M(KI).EQ.1) GO TO 100
          DO L = 1,N
            ETRY(L) = E(L)
          END DO
          DO 120 K = 1,C
            IF(K.EQ.KI) GO TO 120
            ETRY(I) = K
            DO K1 = 1,C
              DO L1 = 1,C
                PTRY(K1,L1) = POS(K1,L1)
                NTRY(K1,L1) = NEG(K1,L1)
              END DO
            END DO
            DO 130 J = 1,N
              IF(J.EQ.I) GO TO 130
              KJ = E(J)
              IF(A(I,J).GT.0) THEN
                PTRY(KI,KJ) = PTRY(KI,KJ) - A(I,J)
                PTRY(K,KJ) = PTRY(K,KJ) + A(I,J)
              ELSE
                NTRY(KI,KJ) = NTRY(KI,KJ) + A(I,J)
                NTRY(K,KJ) = NTRY(K,KJ) - A(I,J)
              END IF
              IF(A(J,I).GT.0) THEN
                PTRY(KJ,KI) = PTRY(KJ,KI) - A(J,I)
                PTRY(KJ,K) = PTRY(KJ,K) + A(J,I)
              ELSE
                NTRY(KJ,KI) = NTRY(KJ,KI) + A(J,I)
                NTRY(KJ,K) = NTRY(KJ,K) - A(J,I)
              END IF
 130        CONTINUE
            ZTRY = 0
            DO K1 = 1,C
              DO L1 = 1,C
                IF(PTRY(K1,L1).GE.NTRY(K1,L1)) THEN
                  ZTRY = ZTRY + PTRY(K1,L1)
                ELSE
                  ZTRY = ZTRY + NTRY(K1,L1)
                END IF
              END DO
            END DO
            IF(ZTRY.GT.Z) THEN
              Z = ZTRY
              IFLAG = 1
              E(I) = ETRY(I)
              DO K1 = 1,C
                DO L1 = 1,C
                  POS(K1,L1) = PTRY(K1,L1)
                  NEG(K1,L1) = NTRY(K1,L1)
                END DO
              END DO
              M(K) = M(K) + 1
              M(KI) = M(KI) - 1
              KI = K
            END IF
 120      CONTINUE
 100    CONTINUE
        IF(IFLAG.EQ.1) GO TO 1000
        IF(Z.GT.ZBEST) THEN
          ZBEST = Z
          DO I = 1,N
            EBEST(I) = E(I)
          END DO
        END IF
      CALL CPU_TIME(TIMEB)
      IF(TIMEB-TIMEA.LT.TLIMIT) GO TO 6000
      OBJVAL = FLOAT(SMAX-ZBEST)
      call fseedo()
      END
