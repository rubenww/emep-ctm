! <BiDir_emep.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.36>
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
module BiDir_emep
  use AllocInits,            only: AllocInit
!  use BiDir_module, only: BiDirConcs
  use CheckStop_mod
  use ChemSpecs_mod,         only: NH3
  use Config_module,         only: BiDirInputFile, BiDirInputDir, &
                                    KMAX_MID,MasterProc, USES, &
                                    step_main, iyr_trend ! QUERY yrtrend or meteo yr?
  use Debug_module,          only: DEBUG
  use GridValues_mod,        only: debug_proc, debug_li,debug_lj, glon, glat
  use LocalVariables_mod,    only: Grid     ! Meteorology as e.g. Grid%t2C
  use NetCDF_mod,            only: ReadField_CDF, printCDF
  use PhysicalConstants_mod, only: AVOG
  use ZchemData_mod,         only: xn_2d
  use PAR_mod,               only: LIMAX, LJMAX, me
  use SmallUtils_mod,        only: key2str
  use TimeDate_mod,          only: date, current_date

  use Config_module, only : USES, MasterProc
  implicit none
  private

  public :: Init_BiDir ! Allocates and sets long-term NH3, updates monthly fields
  public :: load_input_data_BiDir ! Allocates and sets long-term NH3, updates monthly fields
  public :: set_BiDirSea
  private :: fixdates
  integer, public :: old_month
  
  !long term NH3 (monthly or annual), ug/m3
  real, public, allocatable, dimension(:,:), save :: BiDir_NH3aLT
  real, public, allocatable, dimension(:,:), save :: BiDir_SO2aLT !Hazelhos 15-11-2019: added, used to calculate aSN
  real, public, allocatable, dimension(:,:), save :: BiDir_NH3inst ! NH3 @ 3-4m
! real, public, allocatable, dimension(:,:), save :: BiDir_NHxDep !Hazelhos 15-11-2019: deleted BiDir_NHxDep, is not used anymore
  real, public, allocatable, dimension(:,:), save :: BiDir_NHxEmis
  real, public, allocatable, dimension(:,:), save :: BiDir_aSN
  real, public, allocatable, dimension(:,:), save :: BiDir_Xtot ! Save results
  real, public, allocatable, dimension(:,:), save :: BiDir_NH3_3m ! NH3 @ 3-4m
  real, public, allocatable, dimension(:,:), save :: BiDir_XH3_3m ! NH3 @ 3-4m

  !JUL2019 SEA STUFF
  real, public, allocatable, dimension(:,:), save :: BiDir_sea_nh4 ! want umole/L
  real, public, allocatable, dimension(:,:), save :: BiDir_sea_pH !
  real, public, allocatable, dimension(:,:), save :: BiDir_sea_gridT !
  real, public, allocatable, dimension(:,:), save :: BiDir_sea_gridS !
! TMP until netcdf read sorted out
!  real, public, allocatable, dimension(:,:), save :: BiDir_NH3acc
!  real, public, allocatable, dimension(:,:), save :: BiDir_SO2acc

contains
!----------------------------------------------------------------------------
  ! quick function to replace YYYY MM etc with dates
  function fixdates(fname) result (f)
    character(len=*), intent(in) :: fname
    character(len=200) :: f
     f=key2str(fname,'YYYY',iyr_trend)
     f=key2str(f,'MM',current_date%month )
     if ( current_date%month < 12 ) then 
       f=key2str(f,'ZZZZ',iyr_trend)
       f=key2str(f,'NN',current_date%month+1 )
     else
       f=key2str(f,'ZZZZ',iyr_trend+1)
       f=key2str(f,'NN', 1 ) !For Jan.
     end if
  end function  fixdates
!----------------------------------------------------------------------------
  subroutine Init_BiDir()
     character(len=*), parameter :: dtxt='BIDIR:Init:'
     character(len=200) :: tmpdir !for SEASTUFF
     character(len=200) :: ifile  !for SEASTUFF
     integer :: i,j
     
     if ( .not. USES%BIDIR ) return
     if ( MasterProc) write(*,*) dtxt//'called DATE',step_main, iyr_trend, current_date
     
     allocate (BiDir_NH3aLT(LIMAX,LJMAX))
	 allocate (BiDir_SO2aLT(LIMAX,LJMAX))!Hazelhos 15-11-2019: added this line
   ! allocate (BiDir_NHxDep(LIMAX,LJMAX))
     allocate (BiDir_aSN(LIMAX,LJMAX))
     allocate (BiDir_Xtot(LIMAX,LJMAX))
     allocate (BiDir_NH3inst(LIMAX,LJMAX))
     allocate (BiDir_NHxEmis(LIMAX,LJMAX))
     allocate (BiDir_NH3_3m(LIMAX,LJMAX))
     allocate (BiDir_XH3_3m(LIMAX,LJMAX))
     
     allocate (BiDir_sea_nh4(LIMAX,LJMAX))
     allocate (BiDir_sea_ph(LIMAX,LJMAX))
     allocate (BiDir_sea_gridS(LIMAX,LJMAX))
     allocate (BiDir_sea_gridT(LIMAX,LJMAX))
	 
     BiDir_NH3inst = -999.  ! Initial value as flag  0.0
     old_month = -999
	 
    !Hazelhos 15-11-2019:
    !A standard simulation, specified by BiDirInputFile in the namelist file, is used for reading SO2 and NH3 
	!values. This is a temporary solution, as we want to obtain SURF_ug_SO2 and SURF_ug_NH3 from the global
    !model for the first domain, and use for the other domains the values from domain 1
	
	
     ! 1) Annual NH3, ug/m3
     call ReadField_CDF(BiDirInputFile, 'SURF_ug_NH3',&
        BiDir_NH3aLT, nstart=1,interpol='zero_order',needed=.true.,&
        UnDef=0.0, debug_flag=.false.)

     call printCDF('BIDIR_NH3',BiDir_NH3aLT,'ugm3')

	 ! 	2) Annual SO2, ug/m3 !Hazelhos 15-11-2019: added
     call ReadField_CDF(BiDirInputFile, 'SURF_ug_SO2',&
        BiDir_SO2aLT, nstart=1,interpol='zero_order',needed=.true.,&
        UnDef=0.0, debug_flag=.false.)

     call printCDF('BIDIR_SO2',BiDir_SO2aLT,'ugm3')
	 
	 
	 !Hazelhos 15-11-2019: Outdated. is not used. Can be removed?
	 
     !  2) Annual DDEP NHx, mgN/m2
     ! WHY? call ReadField_CDF(BiDirInputFile,'DDEP_RDN_m2Water_D', BiDir_NHxDep,&
  !     call ReadField_CDF(BiDirInputFile,'TDEP_RDN', BiDir_NHxDep,& 
  !         nstart=1,interpol='zero_order',needed=.true.,UnDef=0.0, &
  !           debug_flag=.false.)
 

       ! from mgN/m2 to kg(NH4???)/ha:

  !   BiDir_NHxDep = BiDir_NHxDep * 1.0e-6*1.0e4*18.0/14  
  

  

	 
    ! 3) Annual aSN = [SO2]/[NH3]
     !call ReadField_CDF(BiDirInputFile,'aSN', BiDir_aSN, nstart=1,&
     !   interpol='zero_order',needed=.true.,UnDef=0.0, debug_flag=.true.)
     
	 !Hazelhos 20191115: calculate aSN based on SO2 and NH3 (convert ug m-3 to mol m-3), rather than reading it directly from a global file
	 BiDir_aSN = (BiDir_SO2aLT / 64.066) / (BiDir_NH3aLT / 17.031)
     call printCDF('BIDIR_aSN',BiDir_aSN,'ratio')
     
     if( debug_proc ) then
       i=debug_li
       j=debug_lj
       write(*, "(a,3i4,2f8.2)") dtxt//' coords/lonlat', &
          me, i,j, glon(i,j), glat(i,j)
       write(*, "(a,3a12)") dtxt//"    ",'ijVal','min','max'
       write(*, "(a,3f12.3)") dtxt//"NH3LT ",BiDir_NH3aLT(i, j),&
            minval(BiDir_NH3aLT(:,:)), maxval(BiDir_NH3aLT(:,:))
	   write(*, "(a,3f12.3)") dtxt//"SO2LT ",BiDir_SO2aLT(i, j),& !Hazelhos 15-11-2019: added
            minval(BiDir_SO2aLT(:,:)), maxval(BiDir_SO2aLT(:,:))
       ! write(*, "(a,3f12.3)") dtxt//"Dep ",BiDir_NHxDep(i,j),&
       !      minval(BiDir_NHxDep(:,:)), maxval(BiDir_NHxDep(:,:))
       write(*, "(a,3f12.3)") dtxt//"aSN ",BiDir_aSN(i,j),&
            minval(BiDir_aSN(:,:)), maxval(BiDir_aSN(:,:))
      end if

  end subroutine Init_BiDir

  subroutine load_input_data_BiDir()
     character(len=*), parameter :: dtxt='BIDIR:load monthly input data:'
     character(len=200) :: tmpdir !for SEASTUFF
     character(len=200) :: ifile  !for SEASTUFF
	   integer :: i,j

     if ( .not. USES%BIDIR ) return
     if ( current_date%month == old_month ) return
	 if ( MasterProc) write(*,*) dtxt//'called DATE',step_main, iyr_trend, current_date 
 
     ! ============== MONTHLY UPDATES ========================
     !JUL2019 SEA STUFF
     ! 4) 
     ! NH4 is in mmole/m3 which is umole/L
     ! dates: use YYYY, MM and ZZZZ, NN (for end date if needed)

     ifile=trim(BiDirInputDir)//'BiDirSea/FREEBIORYS2V4/YYYY/'// &
              'emep-ext-FREEBIORYS2V4_1mAV_YYYYMM01_ZZZZNN01_nh4.nc'
     ifile= fixdates(ifile)  ! replaces YYYY etc
       
     call ReadField_CDF(ifile, 'nh4', BiDir_sea_nh4, nstart=1,&
       needed=.true.,debug_flag=.true.)
       !TMP interpol='zero_order',needed=.true.,UnDef=0.0, debug_flag=.true.)
       
     ifile=trim(BiDirInputDir)//'BiDirSea/FREEBIORYS2V4/YYYY/'// &
              'emep-ext-FREEBIORYS2V4_1mAV_YYYYMM01_ZZZZNN01_ph.nc'
     ifile= fixdates(ifile)

     call ReadField_CDF(ifile, 'ph', BiDir_sea_ph, nstart=1,&
       needed=.true.,debug_flag=.true.)

if( debug_proc ) then
    
   print *, 'TEST STR bidir',  trim(ifile), debug_li, debug_lj !
   print *, 'TEST STR nh4', maxval(BiDir_sea_nh4), minval(BiDir_sea_nh4)
   print *, 'TEST STR ph',  maxval(BiDir_sea_ph), minval(BiDir_sea_ph)
end if

     ifile=trim(BiDirInputDir)//'BiDirSea/FREEGLORYS2V4/YYYY/'// &
        'emep-ext-GLORYS2V4_ORCA025_YYYYMM_gridT.nc'
     ifile= fixdates(ifile)

     call ReadField_CDF(ifile,'votemper', BiDir_sea_gridT, nstart=1,&
       needed=.true.,debug_flag=.true.)

     !ifile=fix_dates( trim(BiDirInputDir)//'BiDirSea/FREEGLORYS2V4/YYYY/'// &
     !   'emep-ext-GLORYS2V4_ORCA025_YYYYMM_gridS.nc' )
     ifile= trim(BiDirInputDir)//'BiDirSea/FREEGLORYS2V4/YYYY/'// &
        'emep-ext-GLORYS2V4_ORCA025_YYYYMM_gridS.nc' 
     ifile= fixdates(ifile)
     call ReadField_CDF(ifile,'vosaline', BiDir_sea_gridS, nstart=1,&
       needed=.true.,debug_flag=.true.)

     if ( debug_proc ) write(*,*) minval(BiDir_sea_nh4), minval(BiDir_sea_gridT)
     call printCDF('BIDIR_SEA_nh4',BiDir_sea_nh4,'seaNH4units')
     call printCDF('BIDIR_SEA_ph',BiDir_sea_ph,'phUnits')
     call printCDF('BIDIR_SEA_gridT',BiDir_sea_gridT,'seagridT')
     call printCDF('BIDIR_SEA_gridS',BiDir_sea_gridS,'seagridS')

     old_month = current_date%month

  end subroutine load_input_data_BiDir
  
  
  subroutine set_BiDirSea(i,j,ssNH4,sspH,sstK,S)
    integer, intent(in) :: i,j
    real, intent(out) :: sstK, ssNH4, sspH, S
     ssNH4  = BiDir_sea_nh4(i,j)
     sspH   = BiDir_sea_ph(i,j)
     sstK = BiDir_sea_gridT(i,j) + 273.15  ! or sstK_nwp  TO DO
     S    = BiDir_sea_gridS(i,j)
  end subroutine set_BiDirSea

end module BiDir_emep
