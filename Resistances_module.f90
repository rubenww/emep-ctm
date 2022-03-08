! <Resistances_module.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.36>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2020 met.no
!*
!*  Contact information:
!*  Norwegian Meteorological Institute
!*  Box 43 Blindern
!*  0313 OSLO
!*  NORWAY
!*  email: emep.mscw@met.no
!*  http://www.emep.int
!*
!*    This program is free software: you can redistribute it and/or modify
!*    it under the terms of the GNU General Public License as published by
!*    the Free Software Foundation, either version 3 of the License, or
!*    (at your option) any later version.
!*
!*    This program is distributed in the hope that it will be useful,
!*    but WITHOUT ANY WARRANTY; without even the implied warranty of
!*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*    GNU General Public License for more details.
!*
!*    You should have received a copy of the GNU General Public License
!*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*****************************************************************************!
module Resistances_module

  ! NOTATION - follows EMEP/DO3SE  conventions
  ! =======================================================
  !  R, G = canopy or big-leaf values of resistance/conductance
  !  r, g = leaf-level values of resistance/conductance
  !  LC = land-cover
  !  ext = external leaf/veg surface, was 'w' in CEH work
  !  gs = ground surface
  !  igs = in-canopy+gs, e.g. Rigs = Rinc + Rgs = DEPAC Rsoil,eff
  !  ns = non-stomatal, EMEP includes gs+ext, ie Gns=Gext+Gigs
  !  w   = water (do not confuse with CEH w!)
  !  Rsur   = 1/(Gsto + Gigs + Gext)  ! 1/R = 1/Rsto + 1/Rns
  
  implicit none
  private

  ! Main:
  public :: RnsA_ACP2012  ! Rns(NH3),  EMEP ACP2012 method
  public :: RextA_DEPAC   ! Rext(NH3), Rw in CEH/DEPAC notation
  public :: RgsA_DEPAC    ! DEPAC's Rsoil

  ! for testing
  public  :: self_test       ! Used by mk.testX

contains

 !----------------------------------------------------------------------------
 ! From EMEP, ACP2012 (which is based upon Nemitz, Smith etc from CEH)

  function RnsA_ACP2012(degC,fRH,aSN) result(Rns)
      real, intent(in) :: degC, fRH   ! Temp (C), RH, fraction
      real, intent(in) :: aSN    ! acidity ratio, molar [SO4]/[NH3]
      real :: F1, F2
      real :: Rns
     
      F1 = 10.0 * log10(degC+2.0) * exp(100.0*(1.0-fRH)/7.0)
      F2 = 10.0**( 1.6769 - 1.1099 * aSN ) !EMEP's acidity fac

      Rns = min( 1.0/22 * F1 * F2, 200.0 )
      Rns = max( Rns, 10.0 )
  end function RnsA_ACP2012
 !----------------------------------------------------------------------------
 ! Rext as used in DEPAC, from DEPAC Ch.4.
 ! Based upon Sutton & Fowler, 1993, but modified for SAI
 ! DEPAC ch4 notes that this value lower than e.g. Nemitz et al 2001 
 ! due to different pollution climates, and that DEPACE's Xw (Xext) 
 ! accounts for sources in the canopy itself.

  elemental function RextA_DEPAC(SAI,fRH,frozen) result(Rext)
     real, intent(in) :: SAI  ! surface area index, 1-sided, m2/m2
     real, intent(in) :: fRH       ! RH, fraction
     logical, intent(in)  :: frozen
     real, parameter :: SAIHaarweg = 3.5
     real :: Rext
     if ( SAI < 1.0e-6 )  then
       Rext = 99.0e9  ! No uptake since no leaves 
     else if (frozen) then
       Rext = 200.0
     else
       Rext= SAIHaarweg / SAI * 2*exp(100*(1-fRH)/12.0)
     end if
  end function RextA_DEPAC

 !----------------------------------------------------------------------------
 ! 
  subroutine RgsA_DEPAC(water,frozen,RgsDry,RgsWet,dbg)
     !> From DEPAC manual, ch. 5, based upon Erisman et al 1994
     !! Calculate Dry and Wet for later processing. Might make RH-dependent
     !! interpolation rather than simple wet/dry
      logical, intent(in)  :: water
      logical, intent(in)  :: frozen ! to distinguish wet/frozen soil (QUERY ice?)
      real,    intent(out) :: RgsDry, RgsWet
      logical, intent(in), optional  :: dbg

      if ( frozen )  then ! QUERY ICE
        RgsDry = 10000 
        RgsWet = 10000
      else
        RgsDry = 1000
        RgsWet = 10
      end if
      if ( water ) RgsDry = RgsWet ! ;-)
  end subroutine RgsA_DEPAC

 !----------------------------------------------------------------------------
  subroutine self_test()

   real, parameter :: hVeg=20.0, SAI=4.0, ustar=0.5   ! crop example
   real :: Ra=log(45.0/(0.1*hVeg))  / ( 0.41*ustar)   ! z0=0.1 h, neglect d, L
   real :: Rb=6.0/ustar, Rinc=14*hVeg/ustar
   integer :: iRH, iT
   real    :: degC, fRH, aSN
   logical :: water=.false., frozen=.false., dbg=.true.
   real    :: Gext, RgsDry, RgsWet, GnsDry, GnsWet
 
   ! Rgs does not depend on meteo (except wet/dry):

   call  RgsA_DEPAC(water,frozen,RgsDry,RgsWet,dbg)

   aSN = 0.1  ! low SO2/NH3 in eg NL

  ! Print out conductances in cm/s, ie 100/R
   do iRH = 70,100,5
     fRH = 0.01 * iRH
     print *, '====================================================='
     print *, 'RH=', iRH , 'aSN=', aSN
     print '(a,a4,2(2x,4a12))', '> ', 'T', &
        'GnsEMEP', 'GextDEPAC','GnsDEPACdry', 'GnsDEPACwet', &
        'VgEMEP', 'VgDEPACdry', ' VgDEPACwet'

     Gext= 1/RextA_DEPAC(SAI,fRH,frozen)        ! DEPAC
    ! Ggs = 1/(Rinc + Rgs );  Gns = Gext + Ggs
     GnsDry = Gext + 1/(Rinc+RgsDry) ! DEPAC
     GnsWet = Gext + 1/(Rinc+RgsWet) ! DEPAC

     do iT = 0, 30, 10

        degC= real(iT) 

        print '(a,i4,2(2x,4f12.1))', '> ', iT,&
           100/RnsA_ACP2012(degC,fRH,aSN), & ! EMEP
           100*Gext, &                       ! DEPAC, just ext
           100*GnsDry, &                     ! DEPAC dry
           100*GnsWet, &                     ! DEPAC wet
         ! Depositions:
           100/(Ra+Rb+RnsA_ACP2012(degC,fRH,aSN)), & ! EMEP
           100/(Ra+Rb+1/GnsDry), &                   ! DEPAC dry
           100/(Ra+Rb+1/GnsWet)                      ! DEPAC wet
     end do
   end do

  end subroutine self_test
end module Resistances_module
!-----------------------------------------------------------------------------
! And just to test the above....
!TSTEMX program testr
!TSTEMX   use Resistances_module
!TSTEMX   call self_test()
!TSTEMX end program testr
