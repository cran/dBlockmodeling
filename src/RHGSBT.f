      SUBROUTINE rhgsbtf(N, C, TLIMIT, OBJVAL, A, EB, NREPS)
      IMPLICIT INTEGER(A-Z)
      REAL TIMEA,TIMEB,TLIMIT,OBJVAL
      DOUBLE PRECISION S1
      INTEGER A(N,N),E(N),M(C),EB(N),CENTR(C,N)
C
C *******************************************************************
C  RELOCATION ALGORITHM FOR GENERALIZED STRUCTURAL BALANCE
C *******************************************************************
C
      call fseedi()
      CALL CPU_TIME(TIMEA)  
      DO I = 1,N
        EB(I) = 0
      END DO
      SMAX = 0
      DO I = 1,N
        DO J = 1,N
          SMAX = SMAX + ABS(A(I,J))
        END DO
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
        DO I = 1,N-1
          DO J = I+1,N
            KI = E(I)
            KJ = E(J)
            IF(KI.EQ.KJ) THEN
              Z = Z + A(I,J) + A(J,I)
            ELSE
              Z = Z - A(I,J) - A(J,I)
            END IF
          END DO
        END DO
 1000   IFLAG = 0
        DO 100 I = 1,N
          KI = E(I)
          IF(M(KI).EQ.1) GO TO 100
          DO 120 K = 1,C
            IF(K.EQ.KI) GO TO 120
            DELTA = 0
            DO 130 J = 1,N
              IF(J.EQ.I) GO TO 130
              KJ = E(J)
              IF(KI.EQ.KJ) THEN
                DELTA = DELTA - A(I,J) - A(J,I)
              ELSE
                DELTA = DELTA + A(I,J) + A(J,I)
              END IF
              IF(K.EQ.KJ) THEN
                DELTA = DELTA + A(I,J) + A(J,I)
              ELSE
                DELTA = DELTA - A(I,J) - A(J,I)
              END IF
 130        CONTINUE
            IF(DELTA.GT.0) THEN
              Z = Z + DELTA
              IFLAG = 1
              E(I) = K
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
            EB(I) = E(I)
          END DO
        END IF
      CALL CPU_TIME(TIMEB)
      IF(TIMEB-TIMEA.LT.TLIMIT) GO TO 6000
      OBJVAL = FLOAT(SMAX-ZBEST)/2
      call fseedo()       
      END
