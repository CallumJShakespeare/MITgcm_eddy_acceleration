C $Header: /u/gcmpack/MITgcm/pkg/mom_common/mom_calc_diss.F,v 1.2 2013/11/06 00:37:11 jmc Exp $
C $Name: checkpoint65n $

#include "MOM_COMMON_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: MOM_CALC_DISS

C     !INTERFACE:
      SUBROUTINE MOM_CALC_DISS(
     O        hz_diss, vert_diss,
     I        bi, bj, myThid )

C     !DESCRIPTION:
C     Calculates the dissipation of the 3-D flow field

C     !USES:
      IMPLICIT NONE

C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"

C     !INPUT PARAMETERS:
C     bi, bj      :: tile indices
C     myThid      :: my Thread Id number
      INTEGER bi, bj
      INTEGER myThid
      
C     !OUTPUT PARAMETERS:
C     hz_diss   :: hz dissipation terms @ grid-cell center
C     vert_diss :: vert dissipation terms @ grid-cell center
      _RL hz_diss(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL vert_diss(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)

C     !TEMPORARY PARAMETERS:
C     str11       :: u_x @ grid-cell center
C     str22       :: v_y @ grid-cell center
C     str33       :: u_y @ grid-cell center
C     str12       :: v_x @ grid-cell corner
C     str13       :: u_z @ above uVel
C     str23       :: v_z @ above vVel
      _RL str11(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL str22(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL str33(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL str12(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL str13(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr+1)
      _RL str23(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr+1)

C     !LOCAL VARIABLES:
C      i, j, k    :: loop indices
      INTEGER i, j, k
      INTEGER kp1
      _RL maskp1
      LOGICAL freeSlip3d
CEOP
      freeSlip3d = .NOT.( no_slip_sides .AND. no_slip_bottom )

      DO k=1,Nr
       kp1 = MIN(k+1,Nr)
       maskp1 = oneRL
       IF ( k.EQ.Nr ) maskp1 = zeroRL

C-    Fills up array edges:
       i = sNx+OLx
       DO j=1-OLy,sNy+OLy
         str11(i,j,k) = 0. _d 0
       ENDDO
       j = sNy+OLy
       DO i=1-OLx,sNx+OLx
         str22(i,j,k) = 0. _d 0
       ENDDO
       i = 1-OLx
       DO j=1-OLy,sNy+OLy
         str12(i,j,k) = 0. _d 0
         str13(i,j,k) = 0. _d 0
         str33(i,j,k) = 0. _d 0
       ENDDO
       j = 1-OLy
       DO i=1-OLx,sNx+OLx
         str12(i,j,k) = 0. _d 0
         str23(i,j,k) = 0. _d 0
         str33(i,j,k) = 0. _d 0
       ENDDO

C     str11 = u_x
       DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx-1
          str11(i,j,k) = recip_dxF(i,j,bi,bj)
     &           *( uVel(i+1, j , k ,bi,bj)-uVel( i , j , k ,bi,bj) )
        ENDDO
       ENDDO

C     str22 = v_y
       DO j=1-OLy,sNy+OLy-1
        DO i=1-OLx,sNx+OLx
          str22(i,j,k) = recip_dyF(i,j,bi,bj)
     &           *( vVel( i ,j+1, k ,bi,bj)-vVel( i , j , k ,bi,bj) )
        ENDDO
       ENDDO

C     str12 = v_x
       DO j=1-OLy,sNy+OLy
        DO i=2-OLx,sNx+OLx
          str12(i,j,k) = recip_dxV(i,j,bi,bj)
     &           *( vVel( i , j , k ,bi,bj)-vVel(i-1, j , k ,bi,bj) )
        ENDDO
       ENDDO
       
C     str33 = u_y 
       DO j=2-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
          str33(i,j,k) = recip_dyU(i,j,bi,bj)
     &           *( uVel( i , j , k ,bi,bj)-uVel( i ,j-1, k ,bi,bj) )
        ENDDO
       ENDDO

C     str13 & str23 special case: k=1
       IF ( k.EQ.1 .AND. freeSlip3d ) THEN
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          str13(i,j,k) = 0. _d 0
          str23(i,j,k) = 0. _d 0
         ENDDO
        ENDDO
       ELSEIF ( k.EQ.1 ) THEN
C--    should put surface wind-stress if z-coords; but right in p-coords:
        DO j=1-OLy,sNy+OLy
         DO i=2-OLx,sNx+OLx
          str13(i,j,k) = halfRL*
     &      recip_drC(k)*rkSign
     &           *( uVel( i , j , k ,bi,bj)*twoRL )
         ENDDO
        ENDDO
        DO j=2-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          str23(i,j,k) = halfRL*
     &      recip_drC(k)*rkSign
     &           *( vVel( i , j , k ,bi,bj)*twoRL )
         ENDDO
        ENDDO
       ELSE
C     str13 = u_z 
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          str13(i,j,k) = recip_drC(k)*rkSign
     &           *( uVel( i , j , k ,bi,bj)-uVel( i , j ,k-1 ,bi,bj) )
         ENDDO
        ENDDO
C     str23 = v_z
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          str23(i,j,k) = recip_drC(k)*rkSign
     &           *( vVel( i , j , k ,bi,bj)-vVel( i , j ,k-1,bi,bj) )
         ENDDO
        ENDDO
       ENDIF

       IF ( freeSlip3d ) THEN
        DO j=2-OLy,sNy+OLy
         DO i=2-OLx,sNx+OLx
           str12(i,j,k) = str12(i,j,k)
     &                  *maskW(i,j-1,k,bi,bj)*maskW(i,j,k,bi,bj)
         ENDDO
        ENDDO
        IF ( k.GE.2 ) THEN
         DO j=1-OLy,sNy+OLy
          DO i=2-OLx,sNx+OLx
           str13(i,j,k) = str13(i,j,k)
     &                  *maskW(i,j,k-1,bi,bj)*maskW(i,j,k,bi,bj)
          ENDDO
         ENDDO
         DO j=2-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           str23(i,j,k) = str23(i,j,k)
     &                  *maskS(i,j,k-1,bi,bj)*maskS(i,j,k,bi,bj)
          ENDDO
         ENDDO
        ENDIF
       ENDIF

C--   end k loop
      ENDDO

C--   fill-up strain tensor component at the very bottom (k=Nr+1)
      k = Nr+1

      DO j=1-OLy,sNy+OLy
       DO i=1-OLx,sNx+OLx
         str13(i,j,k) = 0. _d 0
         str23(i,j,k) = 0. _d 0
       ENDDO
      ENDDO

      IF ( .NOT.freeSlip3d ) THEN

C     str13 = u_z 
       DO j=1-OLy,sNy+OLy
        DO i=2-OLx,sNx+OLx
          str13(i,j,k) =
     &         recip_drF(Nr)*rkSign
c    &      recip_drC(k)*rkSign
     &            *( 0. _d 0 - uVel( i , j ,k-1 ,bi,bj) )
        ENDDO
       ENDDO

C     str23 = v_z 
       DO j=2-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
        str23(i,j,k) =
     &         recip_drF(Nr)*rkSign
c    &      recip_drC(k)*rkSign
     &            *( 0. _d 0 - vVel( i , j ,k-1,bi,bj) )
        ENDDO
       ENDDO

      ENDIF
      
C AT this point we have all the derivatives needed
C     str11       :: u_x @ grid-cell center
C     str22       :: v_y @ grid-cell center
C     str33       :: u_y @ grid-cell center
C     str12       :: v_x @ grid-cell corner
C     str13       :: u_z @ above uVel
C     str23       :: v_z @ above vVel
C with correct handling of boundaries.

C Now we can compute the terms we're interested in.
      DO k=1,Nr
       DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
C Horizontal dissipation (up to kappa_h)
          hz_diss(i,j,k) = str11(i,j,k)**2.d0
     &                   + str12(i,j,k)**2.d0
     &                   + str22(i,j,k)**2.d0
     &                   + str33(i,j,k)**2.d0
C Vertical profiles: won't work: overlap/single tile
C          hz_diss_z(k) = hz_diss_z(k) 
C     &                 + hz_diss(i,j,k)*rA(i,j,bi,bj)
        ENDDO
       ENDDO
      ENDDO
      
C Now we can compute the terms we're interested in.
      DO k=1,Nr
       DO j=1-OLy,sNy+OLy-1
        DO i=1-OLx,sNx+OLx-1
C Vertical dissipation (up to kappa_v)
C defined at grid edges
          vert_diss(i,j,k) = (0.25d0* (
     &     (str13(i,j,k)+str13(i+1,j,k))
     &   +(str13(i,j,k+1)+str13(i+1,j,k+1))
     &     )  )**2.d0
     &   + (0.25d0* ( 
     &     (str23(i,j,k)+str23(i,j+1,k))
     &   +(str23(i,j,k+1)+str23(i,j+1,k+1))
     &     )  )**2.d0
C Vertical profiles
C          vert_diss_z(k) = vert_diss_z(k) 
C     &               + vert_diss(i,j,k)*rA(i,j,bi,bj)
        ENDDO
       ENDDO
      ENDDO


      RETURN
      END
