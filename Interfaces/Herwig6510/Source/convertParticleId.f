      !> Particle ID conversion subroutine
      !! 
      !> Subroutine used to convert hadron IDs used by
      !> MCSTHAR++ into those adopted by Herwig6510
      !!
      !> @author Christopher Bignamini
      !!
      SUBROUTINE CONVERTPARTICLEID
      INCLUDE 'HERWIG65.INC'
      INTEGER I 
      INTEGER IWIG !< Herwig identity code
      CHARACTER*8 NWIG !< Herwig hadron name
C   TODO: omp
      DO I=1,NHEP !< Loop over event record
         IF(IDHW(I).EQ.0) THEN !< Selection of produced hadrons
            CALL HWUIDT(1,IDHEP(I),IWIG,NWIG)
            IDHW(I) = IWIG
         ENDIF
      ENDDO
      RETURN
      END
