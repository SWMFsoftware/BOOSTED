!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission For more information, see
!  http://csem.engin.umich.edu/tools/swmf
module ModUser

  ! General routines to update primitive variables in boosted frame in which
  ! the stacetime is tramsformed: the time offset is introduced:
  ! Time = Time_boosted + tOffsetPerR*(R - R_0) for R > R_0
  ! Follows description and notation in Sokolov, I.V. and Gombosi, T.I. (2024)
  ! Physics-Based Forecasting of Tomorrow's Solar Wind at 1 AU
  ! This module compiles and works if the folder THREAD is copied from
  ! repository to the root directory of the SWMF and the boosted frame
  ! is installed as follows:
  !
  ! cd THREAD
  ! make installboost
  !
  ! Revision history:
  ! 12/25/2024 - Igor Sokolov, implemented vesion as described in the
  !              2024 arxiv paper

  use BATL_lib, ONLY: &
       test_start, test_stop, iProc, lVerbose
  use ModIO,        ONLY: write_prefix, write_myname, iUnitOut
  use ModMain, ONLY: nI, nJ, nK
  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_init_session,               &
       IMPLEMENTED3 => user_set_resistivity,            &
       IMPLEMENTED4 => user_update_states,              &
       IMPLEMENTED5 => user_calc_timestep

  include 'user_module.h' ! list of public methods

  character (len=*), parameter :: NameUserFile = "ModUserBoosted.f90"
  character (len=*), parameter :: NameUserModule = 'Boosted frame model'

  ! Input parameters for two-temperature effects
  real    :: TeFraction, TiFraction
  real    :: EtaPerpSi
  real    :: tOffsetSi = -1.0, rMinTimeOffset, rTimeOffset
  real    :: BoostFactor = 1.0

contains
  !============================================================================
  subroutine user_read_inputs

    use ModReadParam, ONLY: read_line, read_command, read_var
    use ModMain,      ONLY: UseRotatingFrame

    character (len=100) :: NameCommand
    integer :: iDir
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_read_inputs'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    if(iProc == 0 .and. lVerbose > 0)then
       call write_prefix;
       write(iUnitOut,*)'User read_input for boosted frame module starts'
    endif

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE

       select case(NameCommand)

       case('#ROTATINGFRAME')
          call read_var('UseRotatingFrame',UseRotatingFrame)

       case('#BOOSTEDFRAME')
          call read_var('tOffsetSi', tOffsetSi)
          call read_var('rMinTimeOffset', rMinTimeOffset)
          call read_var('rTimeOffset',rTimeOffset)

       case('#USERINPUTEND')
          if(iProc == 0 .and. lVerbose > 0)then
             call write_prefix;
             write(iUnitOut,*)'User read_input SOLAR CORONA ends'
          endif
          EXIT

       case default
          if(iProc == 0) then
             call write_myname; write(*,*) &
                  'ERROR: Invalid user defined #COMMAND in user_read_inputs. '
             write(*,*) '--Check user_read_inputs for errors'
             write(*,*) '--Check to make sure a #USERINPUTEND command was used'
             write(*,*) '  *Unrecognized command was: '//NameCommand
             call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
          end if
       end select
    end do

    call test_stop(NameSub, DoTest)
  end subroutine user_read_inputs
  !============================================================================
  subroutine user_init_session

    use ModBoostedFrame, ONLY: tOffsetPerR, set_param
    use ModIO,         ONLY: write_prefix, iUnitOut
    use ModAdvance,    ONLY: UseElectronPressure
    use ModMultiFluid, ONLY: MassIon_I
    use ModConst,      ONLY: cElectronCharge, cLightSpeed, cBoltzmann, cEps, &
         cElectronMass
    use ModNumConst,   ONLY: cTwoPi, cDegToRad
    use ModPhysics,    ONLY: ElectronTemperatureRatio, AverageIonCharge, &
         Si2No_V, No2Io_V, UnitU_, UnitT_, Gamma, GammaElectron

    real, parameter :: CoulombLog = 20.0
    real :: tOffsetPerRin
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_init_session'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if(iProc == 0)then
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'user_init_session:'
       call write_prefix; write(iUnitOut,*) ''
    end if

        ! TeFraction is used for ideal EOS:
    if(UseElectronPressure)then
       ! Pe = ne*Te (dimensionless) and n=rho/ionmass
       ! so that Pe = ne/n *n*Te = (ne/n)*(rho/ionmass)*Te
       ! TeFraction is defined such that Te = Pe/rho * TeFraction
       TiFraction = MassIon_I(1)
       TeFraction = MassIon_I(1)/AverageIonCharge
    else
       ! p = n*T + ne*Te (dimensionless) and n=rho/ionmass
       ! so that p=rho/massion *T*(1+ne/n Te/T)
       ! TeFraction is defined such that Te = p/rho * TeFraction
       TiFraction = MassIon_I(1) &
            /(1 + AverageIonCharge*ElectronTemperatureRatio)
       TeFraction = TiFraction*ElectronTemperatureRatio
    end if

    ! perpendicular resistivity, used for temperature relaxation
    ! Note EtaPerpSi is divided by cMu.
    EtaPerpSi = sqrt(cElectronMass)*CoulombLog &
         *(cElectronCharge*cLightSpeed)**2/(3*(cTwoPi*cBoltzmann)**1.5*cEps)
    if(tOffsetSi > 0.0)then
       tOffsetPerRin = tOffsetSi*Si2No_V(UnitT_)/&
            (rTimeOffset - rMinTimeOffset)
       if(iProc == 0)then
          call write_prefix; write(iUnitOut,*) 'Boosted frame is used'
          call write_prefix; write(iUnitOut,'(a,f7.0),a') 'Lambda_max=',&
               No2Io_V(UnitU_)/tOffsetPerRin,' km/s'
       end if
       call set_param(tOffsetPerRin=tOffsetPerRin, GammaIn=Gamma, &
            GammaElectronIn=GammaElectron)
    end if
    if(iProc == 0)then
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'user_init_session finished'
       call write_prefix; write(iUnitOut,*) ''
    end if

    call test_stop(NameSub, DoTest)
  end subroutine user_init_session
  !============================================================================
  subroutine user_set_resistivity(iBlock, Eta_G)

    use ModAdvance,    ONLY: State_VGB
    use ModPhysics,    ONLY: No2Si_V, Si2No_V, UnitTemperature_, UnitX_, UnitT_
    use ModVarIndexes, ONLY: Rho_, Pe_

    integer, intent(in) :: iBlock
    real,    intent(out):: Eta_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

    integer :: i, j, k
    real :: Te, TeSi

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_set_resistivity'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    do k = MinK,MaxK; do j = MinJ,MaxJ; do i = MinI,MaxI
       Te = TeFraction*State_VGB(Pe_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
       TeSi = Te*No2Si_V(UnitTemperature_)

       Eta_G(i,j,k) = EtaPerpSi/TeSi**1.5 *Si2No_V(UnitX_)**2/Si2No_V(UnitT_)
    end do; end do; end do

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine user_set_resistivity
  !============================================================================
  subroutine user_update_states(iBlock)

    use ModUpdateState, ONLY: update_state_normal
    use ModAdvance,    ONLY: Source_VC, StateOld_VGB, State_VGB, DtMax_CB, &
         UseElectronPressure, Source_VC, Flux_VXI, Flux_VYI, Flux_VZI,     &
         nVar, nVarUpdate, iVarUpdate_I
    use ModBoostedFrame, ONLY: tOffsetPerR, update_mhd_var, update_mhd_pe_var
    use ModVarIndexes, ONLY:  Rho_, RhoUx_, RhoUy_, RhoUz_, Ux_, Uy_, Uz_, &
         Bx_, By_, Bz_, Pe_, p_, nFluid, Energy_, ScalarFirst_, ScalarLast_
    use ModMain,        ONLY: IsTimeAccurate, iStage, nStage, Dt, Cfl, &
         UseBufferGrid, UseHalfStep, UseFlic, UseUserSourceImpl
    use ModPhysics,     ONLY: InvGammaMinus1, InvGammaElectronMinus1
    use ModEnergy,      ONLY: limit_pressure
    use ModWaves,       ONLY: nWave, WaveFirst_,WaveLast_, &
         UseWavePressure, UseWavePressureLtd, UseAlfvenWaves, DoAdvectWaves, &
         update_wave_group_advection
    use BATL_lib,       ONLY: Xyz_DGB, CellVolume_GB
    use ModGeometry,    ONLY: r_GB, Used_GB
    use ModBuffer,      ONLY: fix_buffer_grid
    use ModIonElectron, ONLY: ion_electron_source_impl, HypEDecay

    integer,intent(in):: iBlock

    real    :: DtLocal, DtFactor, ScalarFactor, Buffer8_V(8), Buffer9_V(9)
    ! Old and new primitive variables
    real    :: Prim_VC(nVar,nI,nJ,nK), PrimOld_VC(nVar,nI,nJ,nK)
    ! Unit vector of radial direction
    real    :: Radial_DC(3,nI,nJ,nK)
    integer :: i, j, k, iVar
    ! true if the update is State = StateOld + Source
    logical:: DoAddToStateOld

    logical, parameter :: UseScalar = ScalarLast_ >= ScalarFirst_

    logical :: DoTest

    character(len=*), parameter:: NameSub = 'user_update_states'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    if(tOffsetPerR <= 0.0)then
       call update_state_normal(iBlock)
       RETURN
    end if

    ! Nothing to do if time step is zero
    if(IsTimeAccurate .and. Dt == 0.0) RETURN

    DoAddToStateOld = UseHalfStep .or. nStage == 1 .or. nStage == 4

    ! Calculate partial step size compared to largest stable time step
    if(nStage==4.or.UseFlic)then
       ! Classical 4th order Runge-Kutta scheme
       select case(iStage)
       case(1)
          DtFactor = Cfl/2
       case(2)
          DtFactor = Cfl/2
       case(3)
          DtFactor = Cfl
       case(4)
          DtFactor = Cfl/6
       end select
    elseif(UseHalfStep)then
       DtFactor = (Cfl*iStage)/nStage
    else
       DtFactor = Cfl
    end if

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       DtLocal = DtFactor*DtMax_CB(i,j,k,iBlock)
       do iVar = 1, nVar+nFluid
          Source_VC(iVar,i,j,k) = &
               DtLocal* (Source_VC(iVar,i,j,k) + &
               ( Flux_VXI(iVar,i,j,k,1)  - Flux_VXI(iVar,i+1,j,k,1)  &
               + Flux_VYI(iVar,i,j,k,1)  - Flux_VYI(iVar,i,j+1,k,1)  &
               + Flux_VZI(iVar,i,j,k,1)  - Flux_VZI(iVar,i,j,k+1,1)  ) &
               /CellVolume_GB(i,j,k,iBlock) )
       end do
    end do; end do; end do
    ! Source of internal energy
    Source_VC(p_,:,:,:) = Source_VC(p_,:,:,:)*InvGammaMinus1
    if(UseElectronPressure)then
       ! Source of electron internal energy
       Source_VC(Pe_,:,:,:) = Source_VC(Pe_,:,:,:)*InvGammaElectronMinus1
       ! Add electron internal energy to the total energy ???
       ! Source_VC(Energy_,:,:,:) = &
       !     Source_VC(Energy_,:,:,:) + Source_VC(Pe_,:,:,:)
    end if

    ! Vector of radial direction
    Radial_DC = Xyz_DGB(:,1:nI,1:nJ,1:nK,iBlock)
    Prim_VC = State_VGB(:,1:nI,1:nJ,1:nK,iBlock)
    if(DoAddToStateOld)then
       PrimOld_VC = StateOld_VGB(:,1:nI,1:nJ,1:nK,iBlock)
    else
       PrimOld_VC = Prim_VC
    end if

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not. Used_GB(i,j,k,iBlock)) CYCLE
       ! Express velocity in terms of momentum density
       Prim_VC(Ux_:Uz_,i,j,k) = Prim_VC(RhoUx_:RhoUz_,i,j,k)/&
            Prim_VC(Rho_,i,j,k)
       PrimOld_VC(Ux_:Uz_,i,j,k) = PrimOld_VC(RhoUx_:RhoUz_,i,j,k)/&
            PrimOld_VC(Rho_,i,j,k)
       if(UseScalar)then
          ScalarFactor = 1 - tOffsetPerR*&
               sum(Radial_DC(:,i,j,k)*PrimOld_VC(Ux_:Uz_,i,j,k))
          Prim_VC(ScalarFirst_:ScalarLast_,i,j,k) = &
               PrimOld_VC(ScalarFirst_:ScalarLast_,i,j,k)*ScalarFactor + &
               Source_VC(ScalarFirst_:ScalarLast_,i,j,k)
       end if
       ! Normalize vector of radial direction
       Radial_DC(:,i,j,k) = Radial_DC(:,i,j,k)/r_GB(i,j,k,iBlock)
       if(UseElectronPressure)then
          Buffer9_V = Prim_VC(&
               [Rho_,&
               Ux_,&
               Uy_,&
               Uz_,&
               Bx_,&
               By_,&
               Bz_,&
               Pe_,&
               p_],i,j,k)
          call update_mhd_pe_var(Old_V=PrimOld_VC([Rho_,&
               Ux_,&
               Uy_,&
               Uz_,&
               Bx_,&
               By_,&
               Bz_,&
               Pe_,&
               p_],i,j,k),&
               Source_V=Source_VC([Rho_,&
               RhoUx_,&
               RhoUy_,&
               RhoUz_,&
               Bx_,&
               By_,&
               Bz_,&
               Pe_,&
               Energy_],i,j,k), Radial_D=Radial_DC(:,i,j,k), New_V=Buffer9_V)
          Prim_VC(&
               [Rho_,&
               Ux_,&
               Uy_,&
               Uz_,&
               Bx_,&
               By_,&
               Bz_,&
               Pe_,&
               p_],i,j,k) = Buffer9_V
       else
          Buffer8_V = Prim_VC(&
               [Rho_,&
               Ux_,&
               Uy_,&
               Uz_,&
               Bx_,&
               By_,&
               Bz_,&
               p_],i,j,k)
          call update_mhd_var(Old_V=PrimOld_VC([Rho_,&
               Ux_,&
               Uy_,&
               Uz_,&
               Bx_,&
               By_,&
               Bz_,&
               p_],i,j,k),&
               Source_V=Source_VC([Rho_,&
               RhoUx_,&
               RhoUy_,&
               RhoUz_,&
               Bx_,&
               By_,&
               Bz_,&
               Energy_],i,j,k), Radial_D=Radial_DC(:,i,j,k), New_V=Buffer8_V)
          Prim_VC(&
               [Rho_,&
               Ux_,&
               Uy_,&
               Uz_,&
               Bx_,&
               By_,&
               Bz_,&
               p_],i,j,k) = Buffer8_V
       end if
    end do; end do; end do
    if(nVarUpdate == nVar)then
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(UseScalar)then
             ScalarFactor = 1/(1 - tOffsetPerR*&
                  sum(Radial_DC(:,i,j,k)*Prim_VC(Ux_:Uz_,i,j,k)))
             Prim_VC(ScalarFirst_:ScalarLast_,i,j,k) = &
                  Prim_VC(ScalarFirst_:ScalarLast_,i,j,k)*ScalarFactor
          end if
          State_VGB(:,i,j,k,iBlock) = Prim_VC(:,i,j,k)
          State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = &
               State_VGB(Rho_,i,j,k,iBlock)*&
               State_VGB(Ux_:Uz_,i,j,k,iBlock)
       end do; end do; end do
    else
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          State_VGB(iVarUpdate_I,i,j,k,iBlock) = Prim_VC(iVarUpdate_I,i,j,k)
          if(UseScalar)then
             ScalarFactor = 1/(1 - tOffsetPerR*&
                  sum(Radial_DC(:,i,j,k)*State_VGB(Ux_:Uz_,i,j,k,iBlock)))
             State_VGB(ScalarFirst_:ScalarLast_,i,j,k,iBlock) = ScalarFactor*&
                  State_VGB(ScalarFirst_:ScalarLast_,i,j,k,iBlock)
          end if
          State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = &
               State_VGB(Rho_,i,j,k,iBlock)*&
               State_VGB(Ux_:Uz_,i,j,k,iBlock)
       end do; end do; end do
    end if
    if(UseBufferGrid) call fix_buffer_grid(iBlock)
    call limit_pressure(1, nI, 1, nJ, 1, nK, iBlock, 1, nFluid, State_VGB)
    call test_stop(NameSub, DoTest, iBlock)
  end subroutine user_update_states
  !============================================================================
  subroutine user_calc_timestep(iBlock)

    use ModBoostedFrame, ONLY: tOffsetPerR, fast_mhd_speed, fast_mhd_pe_speed
    use ModAdvance,  ONLY: State_VGB, DtMax_CB, UseElectronPressure
    use ModMain,     ONLY: DtMax_B, Unused_B, nBlock
    use ModVarIndexes, ONLY: nVar, Rho_, RhoUx_, RhoUz_, Ux_, Uy_, Uz_, &
         Bx_, By_, Bz_, Pe_, p_
    use BATL_size,   ONLY: nI, nJ, nK
    use BATL_lib,    ONLY: Xyz_DGB, nProc, iComm
    use ModGeometry, ONLY: r_GB, Used_GB
    use ModMpi

    integer, intent(in) :: iBlock
    real :: Radial_DC(3,nI,nJ,nK), Primitive_VC(nVar,nI,nJ,nK)
    integer :: i, j, k, iBlockLoop, iError
    logical :: IsFirstBlock
    logical :: IsLastBlock
    real    :: BoostFactorPe = 1.0
    !--------------------------------------------------------------------------

    if(iBlock==1)then
       IsFirstBlock = .true.
    else
       IsFirstBlock = all(Unused_B(1:iBlock-1))
    end if
    if(IsFirstBlock)BoostFactorPe = 1.0

    Primitive_VC = State_VGB(:,1:nI,1:nJ,1:nK,iBlock)
    Radial_DC = Xyz_DGB(:,1:nI,1:nJ,1:nK,iBlock)
    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not. Used_GB(i,j,k,iBlock)) CYCLE
       ! Radialize radial vector
       Radial_DC(:,i,j,k) = Radial_DC(:,i,j,k)/r_GB(i,j,k,iBlock)
       ! Express velocity in terms of momentum density
       Primitive_VC(Ux_:Uz_,i,j,k) = Primitive_VC(RhoUx_:RhoUz_,i,j,k)/&
            Primitive_VC(Rho_,i,j,k)
       ! Calculate boost factor
       if(UseElectronPressure)then
          BoostFactorPe = min(BoostFactorPe, &
               (1 - tOffsetPerR*(&
               sum(Radial_DC(:,i,j,k)*Primitive_VC(Ux_:Uz_,i,j,k)) +&
               fast_mhd_pe_speed(&
               Primitive_VC([Rho_,Ux_,Uy_,Uz_,Bx_,By_,Bz_,Pe_,p_],i,j,k),&
               Radial_DC(:,i,j,k) ) ) ) )
       else
          BoostFactorPe = min(BoostFactorPe, &
               (1 - tOffsetPerR*(&
               sum(Radial_DC(:,i,j,k)*Primitive_VC(Ux_:Uz_,i,j,k)) +&
               fast_mhd_speed(&
               Primitive_VC([Rho_,Ux_,Uy_,Uz_,Bx_,By_,Bz_,p_],i,j,k),&
               Radial_DC(:,i,j,k) ) ) ) )
       end if
    end do; end do; end do
    if(BoostFactorPe <= 0.0)call STOP_mpi(&
         'Negative time step, reduce time offset')
    if(iBlock==nBlock)then
       IsLastBlock = .true.
    else
       IsLastBlock = all(Unused_B(iBlock+1:nBlock))
    end if
    if(IsLastBlock)then
       if(nProc > 1)then
          call MPI_allreduce(&
               BoostFactorPe, BoostFactor, 1, MPI_REAL, MPI_MIN, iComm, iError)
       else
          BoostFactor = BoostFactorPe
       end if
       do iBlockLoop = 1, iBlock-1
          if(Unused_B(iBlockLoop))CYCLE
          DtMax_CB(1:nI,1:nJ,1:nK,iBlockLoop) = &
               DtMax_CB(1:nI,1:nJ,1:nK,iBlockLoop)*BoostFactor
          DtMax_B(iBlockLoop) = DtMax_B(iBlockLoop)*BoostFactor
       end do
       DtMax_CB(1:nI,1:nJ,1:nK,iBlock) = &
            DtMax_CB(1:nI,1:nJ,1:nK,iBlock)*BoostFactor
    end if
  end subroutine user_calc_timestep
  !============================================================================
end module ModUser
!==============================================================================
