! <BiDir_module.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.36>
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
module BiDir_module

  ! Bi-directional  exchange module, following methods developed
  ! by Roy Wichink Kruit (RJWK), and DEPAC manual 
  ! Labels : D3.4 means DEPAC chapter 3, eqn 4.
  !          F3   means DEPAC Appendix F, eqn 3.
  !
  ! Refs:
  ! RJWK10: Wichink Kruit, R, van Pul, W, Sauter, F, van den Broek,
  ! M Nemitz, E.  Sutton, M, Krol, M. & Holtslag, A. Modeling the
  ! surface-atmosphere exchange of ammonia Atmos. Environ., 2010, 44,
  ! 945 - 957

  ! RJWK12:  Wichink Kruit, R., J Schaap, M, Sauter, F., J van Zanten,
  ! M. C. & van Pul, W. A. J. Modeling the distribution of ammonia
  ! across Europe including bi-directional surface-atmosphere exchange
  ! Biogeosciences, 2012, 9, 5261-5277

  ! RJWK17: Wichink Kruit, R. J, Aben, J, de Vries, W, Sauter, F, van
  ! der Swaluw, E, van Zanten, M. C. & van Pul, W. A. J., Modelling trends
  ! in ammonia in the Netherlands over the 1990-2014 Atmos. Env., 2017,
  ! 154, 20-30

  ! NOTATION - follows EMEP/DO3SE  conventions
  ! =======================================================
  !    ext = external leaf/veg surface, was 'w' in CEH work
  !    w   = water (do not confuse with CEH w!)
  !    gs = ground surface
  !    R, G = canopy or big-leaf values
  !    r, g = leaf-level values
  !    LC = land-cover
  
  ! ISSUES
  ! =======================================================
  ! 1. Still need to think more about ground-surface Rgs
  ! 2. And water stuff can be added in 2018. Will try to find
  !   ocean concs NH4

  ! FUTURE; Consider Tleaf calculation. EMEP/LOTOS use Tair so far

  use Resistances_module, only : RnsA_ACP2012, RextA_DEPAC, RgsA_DEPAC
  implicit none
  private

  ! Main:
  public  :: BiDirXconcs    ! Sets Gammas, X-terms, Rext, for all LC at once
  public  :: BiDirFluxes     ! Depositions and Emissions from specific LC
  public  :: BiDirResistances ! Rext, Gigs,  Rsur
  ! Sea-stuff Jul 2019
  public :: BiDirXwater

  ! for testing
  public  :: self_test       ! Used by mk.testX


  real, private, save :: Xext,Xgs,Xsto  ! hazelhos deleted Xwater

contains

 !---------------------------------------------------------------------------
  function BiDirXwater(i,j,sstK_nwp,ssNH4,sspH,sstK_monthly,S,water,dbg,method) result(Xwater) !Hazelhos 21-02-2020. Added water boolean and dbg.
   ! cf asman 1994 AE
    integer, intent(in) :: i,j
    real, intent(inout) :: sstK_monthly, ssNH4, sspH, S !Hazelhos 21-02-2020: intent(in) --> intent(inout)
    character(len=*), intent(in) :: method
	
	logical, intent(in) :: water !Hazelhos 21-02-20: added
	logical, intent(in) :: dbg   !Hazelhos 21-02-20: added
	character(len=*), parameter :: dtxt='BiDirXwater:'

    real, intent(in) :: sstK_nwp ! Temps
!     real, intent(in) :: ssNH4    ! umol/L
!     real, intent(in) :: sspH      !
!     real, intent(in) :: S   ! salinity  promille
     real :: Xwater    ! eq conc air

     real :: IonicStrength   ! ionic strength
     real :: gam_nh4, gam_nh3, Knh4, Hnh3, sstK
     real, parameter :: mwNH3 = 17.0, Rgas = 8.2075e-5 ! atm m3 /mol/K

     if ( method ==  'nwpSST') then
       sstK = sstK_nwp
     else
       sstK = sstK_monthly
     end if
	 
	 !if ( dbg ) write(*,'(a,4es12.3,L2)') dtxt//' ssNH4, sspH, sstK, S, water', ssNH4, sspH, sstK, S, water
	 !Hazelhos 21-02-2020+: Values are to be based on RWS waterinfo analysis. These are dummy values for tests.
	 !Hazelhos 20-03-2020:	We should think of a better general way of parameterizing this. Perhaps by reading a global field?
	 !						Also consider splitting the section below up per parameter, not only based on ssNH4.
	 if ( ( water .EQ. .TRUE. ) .AND. ( ssNH4 > 1e35 ) ) then !Invalid flag is 1e36. Is there a more elegant way?
		 ssNH4 = 0.2
		 sspH = 8.2
		 !sstK = 288.0! Get temp from grid? --> sstK_nwp is already from the grid. so do nothing  for sstK is fine here.
		 S = 0.5
		 !if ( dbg ) write(*,'(a,4es12.3,L2)') dtxt//' ssNH4, sspH, sstK, S, water', ssNH4, sspH, sstK, S, water
	 end if
	 !Hazelhos-
	 
	 !Hazelhos 21-02-2020+: Note these functions might not be valid for fresh water. For future work: Reparameterize for fresh water bodies somehow?
     gam_nh4 = 0.883 - 0.0768 * log(S)   !order 0.61 for 35%%

     IonicStrength =  0.00147 + 0.01988 * S + 2.08357e-5 *S*S 

     gam_nh3 = 1 + 0.085 * IonicStrength

     Knh4 = 5.67e-10 * exp( -6286 * (1/sstK - 1/298.15))

     Hnh3 = 56*exp( 4092*(1/sstK - 1/298.15) )  ! M /atm

  !Knh4 * Hnh3 ~ 10 000 from BD

     Xwater =  mwNH3 * ssNH4/ &
           ( Rgas * sstK * Hnh3 * (1/gam_nh3 + ( 10**(-sspH))/ (gam_nh4 * Knh4 )) )
    
  end function BiDirXwater

 !---------------------------------------------------------------------------
  subroutine BiDirXconcs(t2C,sstK,Xwater,aSN,NH3aInst,NH3aLT,dbg)!Hazelhos autumn 2019: removed NHxDep, Xwater_in -> Xwater, !Hazelhos 07-02-2020: added water, Hazelhos 20-03-2020 removed water.

      real, intent(in) :: t2C, sstK ! Temps
      real, intent(inout) :: Xwater ! pre-calculated if monthly data available, otherwise -999   ! Hazelhos, autumn 2020:  Xwater_in -> Xwater, intent(in) -> intent(inout)
      real, intent(in) :: aSN       !  SO2/NH3
      real, intent(in) :: NH3aInst  ! Ammonia at 4 m , ug/m3, instantaneous
      real, intent(in) :: NH3aLT    ! Ammonia at 4 m , ug/m3, monthly or annual
      logical, intent(in) :: dbg
	  !Hazelhos 20-03-2020- logical, intent(in) :: water	! Hazelhos 07-02-2020: added water boolean, 
      real :: Gamma_micromet        ! from micro-met at 4m
      real :: Gamma_sto             ! conc. at stomatal interface
      real :: Gamma_ext             ! conc. at leaf-surface-water interface
      real :: Gamma_water           ! seas, rivers, see RJWK12
      real :: Gamma_gs              ! ground-surface
      character(len=*), parameter :: dtxt='BiDirXconc:'

      real :: nh3i, Tk, fT, fSST, FaSN  !Hazelhos autumn, 2020: Xwater deleted. Is obsolete as it is is now intent(inout)

      Tk   = t2C + 273.15

      nh3i      = max(NH3aInst, 1.0e-6 ) ! Prevents numerical problems

      Gamma_ext = 1.84e3*nh3i  * exp( -0.11*t2C ) - 850 !RJWK10 [13],F.1

      if ( dbg ) write(*,*) 'BiDirGamX:',nh3i, t2C, Gamma_ext
      if ( Gamma_ext > 0.0 ) then  ! co-deposition, RJWK17 [5], aSN <0.83
        FaSN  = 1.1 - 1.32 * aSN 
         if ( dbg ) write(*,*) 'BiDirFaSNA:',nh3i, t2C, Gamma_ext, FaSN
        Gamma_ext =  Gamma_ext * max( 0.0, FaSN)
         if ( dbg ) write(*,*) 'BiDirFaSNB:',nh3i, t2C, Gamma_ext, FaSN
      else
        Gamma_ext = 0.0
      end if

      Gamma_gs =  0.0 ! See DEPAC 5 Set zero as too uncertain
      !Gamma_gs =  Gamma_ext ! PRELIM !!!!!!! CHECK !!!!!
                            ! How do we deal with bare-soil etc??
                            ! and surface area?

      Gamma_micromet = 362 * NH3aLT                        !RJWK10 [15a]

      Gamma_sto = 4.7* Gamma_micromet * exp( -0.071*t2C )  !RJWK10 [16],F.3
	  
      if ( dbg ) write(*,'(a,1es12.3)') dtxt//' Before: Xwater ', Xwater !Hazelhos: added print statement

	  
	 !Hazelhos 20-03-2020:
	 ! The section below should not be in an if statement. There should be an Xwater and a Xext, Xigs and Xsto calculated for every cell.
	 ! In the function BiDirFluxes(), the dependency on LU is used to calculate Xtot and the relative contributions of each LU type.
	 ! Since for LU types that are not water, Xwater is not used, we do not need to calculate it here.
	 ! We changed it, so that  Xsto, Xgs, Xwater and Xext are calculated for every grid cell, independent of LU-type.

	 
    !if( Xwater < 0 ) then		 	!Hazelhos 17-12-2019:  This caused the problem that for land grid cells with values from CMEM, no Xsto and Xgs are calculated, which are needed further on.
	!Hazelhos 20-03-2020- if( water .EQ. .False. ) then	!Hazelhos 07-02-2020:  Fixed the issue described above by having it depend on whether the grid cell is water whether Xsto, Xgs and Xwater are calculated.
	
	  if ( dbg ) write(*,'(a,1es12.3)') dtxt//' Xwater', Xwater
	  !Hazelhos 20-03-2020: we outcommented Gamma_water, fSST and Xwater below, Xwater is already calculated in BiDirXwater().
      !Hazelhos 20-03-2020- Gamma_water = 0.0                                 !TMP for 2017 work

      fT   = 2.75e15/Tk   * exp( -1.04e4/Tk )           ! cf RJWK10 [9]
      !Hazelhos 20-03-2020- fSST = 2.75e15/sstK * exp( -1.04e4/sstK )


      Xext = max(0.0, fT * Gamma_ext)                   !RJWK10 [9],F.2  ug/m3

      if( Xext > 1.0e5 ) then ! TMP checks
          print '(a,f7.2,4es10.2)', 'ERROR? BIGEXT ', tK, fT, aSN, FaSN, nh3i
          Xext = -1 * Xext ! MAKE MINUS TO SIGNAL PROBLEM
      end if

      Xsto   = fT * Gamma_sto                   !RJWK10 [12=9],F.4 ug/m3
	  
	  if (dbg) write(*,'(a,5es12.5)') dtxt//'fT, Gamma_ext, Xext, Gamma_sto, Xsto', fT, Gamma_ext, Xext, Gamma_sto, Xsto !Hazelhos 27-03-2020. added for debugging

      Xgs    = fT * Gamma_gs                    !PRELIM  ug/m3

      !Hazelhos 20-03-2020- Xwater =  fSST * Gamma_water              !RJWK12, adjusted, check!  ug/m3
   !Hazelhos 20-03-2020- else
   !  Xwater = Xwater_in  ! Hazelhos, autumn 2020: this was not needed anymore, since we changed the intent(in) of Xwater to intent(inout)
   
   ! Hazelhos+ 14-02-2020: Two things added here. The first statement was necessary to cover the grid cells with LU water but without a value in CMEMS.
   ! 					   In the second part, we reset Xsto, Xgs and Xext to 0, as they can contain values from the previous Land use loop. 
     if ( Xwater < 0 ) Xwater = 0 

	  !if ( dbg ) write(*,'(a,3es12.3)') dtxt//' Xsto, Xgs, Xext:  ', Xsto, Xgs, Xext
      !Xsto 	= 0 !Hazelhos 20-03-2020: For cells that are partly water, Xsto Xgs and Xext are set to 0. However, in the LU loop, this means that Xtot is only dependent on emissions from water in these cells. Fix by removing this?
	  !Xgs  	= 0
	  !Xext 	= 0
   ! Hazelhos-
   !Hazelhos 20-03-2020- end if ! Cell = Water or Land
   
   if ( dbg ) write(*,'(a,es12.3)') dtxt//' After: Xwater ', Xwater    !Hazelhos: added print statement

      if ( dbg ) then
        write(*,'(/,a,2(1x,a,f5.1),3(1x,a,es8.1),3(1x,a,f5.1))') dtxt, &
          'tC:', t2C,' aSN:', aSN, &
          ' Gam: ext', Gamma_ext, 'sto', Gamma_sto, 'gs', Gamma_gs,  &
          '  Xs: ext', Xext,      'sto', Xsto,      'gs', Xgs
      end if
  end subroutine BiDirXconcs

 !---------------------------------------------------------------------------
 ! Get ground-surface (eg soil), in-canopy and external-leaf resistances,
 !  along with surface canopy resistance
 ! (NB DEPAC has no Rinc for grass etc ... consider later. For now we stick with
 !  EMEP methods which use Rinc whenever SAI significant)

  subroutine BiDirResistances(SAI,fRH,frozen,water,Rinc,Gsto,&
                                 Gext,Gigs,Rsur,debug)
      real, intent(in) :: SAI         ! surface area index (1-sided, m2/m2)
      real, intent(in) :: fRH         ! RH, fraction
      logical, intent(in) :: frozen, water
      real, intent(in) :: Rinc        ! in.canopy resistance  (s/m)   
      real, intent(in) :: Gsto        ! bulk stom. conductance  (m/s)
      real, intent(out) :: Gext       ! bulk external-leaf conductance (m/s)
      real, intent(out) :: Gigs       ! ground-surface (soil, plus Rinc term) (m/s)
      real, intent(out) :: Rsur       ! canopy resistance  (s/m)   
      logical, intent(in), optional :: debug
      logical :: dbg = .false.
      character(len=*), parameter :: dtxt='BiDirRes:'

      real    :: RgsDry, RgsWet, Rgs, Rext, Gns
 
      if ( present(debug) ) dbg = debug
	  
	  if ( dbg ) print '(a, 2es9.2)', 'Bidir resistances Rinc, Gsto: ', Rinc, Gsto !Hazelhos added print statement

      call  RgsA_DEPAC(water,frozen,RgsDry,RgsWet,dbg)

  ! QUERY HERE
  ! TMP! DAVE ASSUMPTION. As EMEP doesn't treat dry/wet separately, we assume
  ! increasing wetness after RH 70%. Revise in future; maybe use exponential
  ! (doesn't seem to matter very much as Rext not affected!?)

      Rgs = RgsDry
      if ( fRH > 0.7 ) Rgs = Rgs + ( fRH-0.7) * (RgsWet-RgsDry)
      Gigs = 1/(Rinc + Rgs )

      Rext  = RextA_DEPAC(SAI,fRH,frozen)

      Gext  = 1/Rext
      Gns   = Gext + Gigs
      Rsur = 1/( Gsto + Gns )

      if ( dbg ) write(*,'(a,2f7.2,99g10.2)') dtxt, SAI, fRH, &
               RgsDry, RgsWet, Gigs, Gext, Gns, Rsur
	  if ( dbg ) write(*,'(a,2L2)') dtxt//' Water, frozen: ',  water, frozen !Hazelhos added print statement

  end subroutine BiDirResistances

 !---------------------------------------------------------------------------
  ! BiDirFluxes - calculates effective surface (???) compensation point  Xtot
  ! and associated emission terms.
  ! The original EMEP NH3 resistances used co-dep terms from Nemitz/Smith
  ! For BiDir NH3 we use different Rext terms since Xtot will take
  ! care of the NH3 resistances (rephrase....)
  ! Need to re-calculate Rsur since Rns changing

  subroutine BiDirFluxes(water,SAI,Ve,Gsto,Gext,Gigs,Xwater,Xtot,Xemis,debug)!Hazelhos: added Xwater as an argument

      logical, intent(in) :: water
      real, intent(in) :: SAI         ! surface area index (1-sided, m2/m2) 	  !Hazelhos: is not used here! Remove?
      real, intent(in) :: Ve          ! Exchange velocity (m/s) =1/(Ra+Rb+Rsur)
      real, intent(in) :: Gsto        ! bulk stom. conductance  (m/s)
      real, intent(in) :: Gext        ! bulk external-leaf conductance (m/s)
      real, intent(in) :: Gigs        ! ground-surface (soil plus Rinc term) (m/s)
      real, intent(in) :: Xwater      ! surface compensation point water (ug/m3)  ! Hazelhos added
      real, intent(out)   :: Xtot     ! effective surface comp. point  (ug/m3)
      real, intent(out)   :: Xemis    !  Emissions from Xtot. We want molec/cm2/s
      logical, intent(in), optional :: debug
      logical :: dbg = .false.

      character(len=20) :: errmsg = 'ok'
      character(len=*),parameter :: dtxt = 'BiDirFlx:'
      real :: Rsur, Gsur
     ! Units F= ug/m3 * m/s = ug/m2/s.  Want molec/cm2/s
      real, parameter :: AVOG = 6.023e23  ! Avogadros
      real, parameter :: toMoleccm2s = 1.0e-6/17.0*AVOG*1.0e-4

      if ( present(debug) ) dbg = debug

      Gsur   = Gsto + Gigs + Gext
      Rsur   = 1/Gsur 
	  if ( dbg ) write(*,'(a,4es9.2)') dtxt//'Gsto, Gigs, Gext, Rsur:', Gsto, Gigs, Gext, Rsur !Hazelhos added print statement

      if ( water ) then
        Xtot = Xwater
        Xemis = Xtot * Ve * toMoleccm2s   ! emis part of D3.6
        !BUG Xemis = Xwater * Gsur * toMoleccm2s
        if ( dbg ) write(*,'(a,3es9.2,L2)') dtxt//'Water:', Xtot, Xemis, Xwater, water
      else
      
      ! Weighted average surface conc./comp. point; Units: ug/m3
      ! Land-based: (water just had Xtot = Xwater)

       !Xtot=( Rsur/Rext * Xext + Rsur/Rgseff * Xgs + Rsur*Gsto * Xsto ) !D3.4
       ! Use Gext and Gsto to simplify zero conductance cases:

        Xtot =  Rsur* ( Gext*Xext + Gigs*Xgs + Gsto*Xsto ) !D3.4 							!Hazelhos 20-03-20: Xext, Xgs and Xsto are apparently 0 for partly water cells. Results in only emissions from the part that is water. Fix!

     ! Emissions :-)

        Xemis = Xtot * Ve * toMoleccm2s   ! emis part of D3.6 								!Hazelhos 10-12-19: Xemis is not used. Why calculate? We have BiDirNHxemissions in output now

        if ( dbg )  then
         write(*,'(a,3(3x,a,f7.4,es8.1),2(1x,a,es9.2))') dtxt//' Gam,X:', &
           'ext G,X', Gext, Xext, 'igs G,X', Gigs, Xgs, &
           'sto G,X', Gsto,Xsto, ' => Xtot:',  Xtot, 'Xemis:', Xemis
        end if
      end if 

      if ( Xtot > 1.0e5 ) then ! had some bug somewhere
           if ( Xtot > 1.0e10 )  errmsg=dtxt//'BIGXTOT:'
           print '(a,f6.2,3(f7.4,es8.1),a,2es9.2)', errmsg, SAI, &
              Gext, Xext, Gigs, Xgs, Gsto,Xsto, ' RX:', Rsur, Xtot
           Xtot = -1 * Xtot ! Will trigger warning
           Xemis = -999.0e9

      end if

  end subroutine BiDirFluxes

  subroutine self_test()

   real, parameter :: hVeg=1.0, SAI=4.0, ustar=0.5   ! crop example
   real    :: Ra=log(45.0/(0.1*hVeg))  / ( 0.41*ustar) ! z0=0.1 h, neglect d,L
   real    :: Rb=6.0/ustar, Rinc=14*hVeg/ustar
   real    :: degC, fRH, aSN, sst = 283.0,Gsto=0.02
   real    :: Ve_BD       ! Exchange velocity (m/s) =1/(Ra+Rb+Rsur)
   real    :: nh3i, nh3lt = 10.0   ! ug/m3  inst, long-term NH3 
   real    :: Xtot, Xemis   ! BiDir terms
   real    :: RsurBD, Rsur_emep, Rns_emep, GextBD, GigsBD, Gsur_emep
   real    :: Xwater ! SEA STUFF
   integer :: iRH, iT, inh3
   logical :: water=.false., frozen=.false., dbg=.true.

  ! 2) RH, T, NH3i ,tests:

   aSN = 0.1
   do iRH = 70,100,5
     fRH = 0.01 * iRH
     write(*,'(9(a6,i5))') "== RH ", iRH , 'h',int(hVeg), &
        'Rinc ', nint(Rinc), ' ==============='

     call BiDirResistances(SAI,fRH,frozen,water,Rinc,Gsto,&
                                 GextBD,GigsBD,RsurBD,debug=.true.)

     Ve_BD = 1/( Ra + Rb + RsurBD )

     do inh3 = 5, 20, 5
       nh3i = real(inh3)
       write(*,'(/,a,9(1x,a,f6.1),/)') '==:', 'NH3inst', nh3i, 'NH3LT', nh3lt,&
          'SAI', SAI, 'h', hVeg, 'Ve(cm/s)', 100*Ve_BD

       do iT = 0, 30, 10
          degC = real(iT)

         ! 3)  Sets Gammas and X values 
         Xwater = 0.1  										! Hazelhos added
          call BiDirXconcs(degC, sst, Xwater, aSN=aSN,  &   ! Hazelhos Xwater_in = 0.1 not permitted anymore
             NH3aInst=nh3i,NH3aLT=nh3lt,  dbg=.true.)		! Hazelhos: removed DepNHx = 5,Hazelhos 20-03-2020 removed water = .true.

!  subroutine BiDirXconcs(t2C,sstK,Xwater,aSN,NH3aInst,NH3aLT,DepNHx, dbg)

          if ( Xext < 0.0 ) stop 'Xext'

         ! 4)  Gets Xtot and new Rsur :
          call BiDirFluxes(water,SAI,Ve_BD,Gsto,GextBD,GigsBD,Xwater,Xtot,Xemis,dbg)  ! Hazelhos: added Xwater


        ! Compare with EMEP, which is temp dependent
          Rns_emep =  RnsA_ACP2012(degC,fRH,aSN)
          Gsur_emep = 1/Rns_emep + Gsto
          Rsur_emep = 1/Gsur_emep 

          write(*,'(a,2f8.2,5x,a,2f7.2)') 'BDcomp (emep,BD)  Rc:', &
           Rsur_emep, RsurBD,'Vg (cm/s):',  100/(Ra+Rb+Rsur_emep), 100*Ve_BD
  
       end do ! iT
     end do ! NH3
   end do ! RH

   ! Jul 2019 SEA STUFF:
    !call BiDirSea(sstK= 285.0,ssNH4=2.6,sspH=8.0, S=35.0,Xwater=Xwater)
    !write(*,*) 'BDsea Xwater max', Xwater
    !call BiDirSea(sstK= 285.0,ssNH4=2.6,sspH=7.0, S=35.0,Xwater=Xwater)
    !write(*,*) 'BDsea Xwater max pH', Xwater
    !call BiDirSea(sstK= 285.0,ssNH4=0.026,sspH=8.0, S=35.0,Xwater=Xwater)
    !write(*,*) 'BDsea Xwater low', Xwater
    !call BiDirSea(sstK= 285.0,ssNH4=0.026,sspH=8.0, S=5.0,Xwater=Xwater)
    !write(*,*) 'BDsea Xwater lowS', Xwater


  end subroutine self_test
end module Bidir_module
!-----------------------------------------------------------------------------
! And just to test the above....
!TSTEMX program testr
!TSTEMX   use Bidir_module
!TSTEMX   call self_test()
!TSTEMX end program testr
