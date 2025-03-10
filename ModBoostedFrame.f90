!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModBoostedFrame

  ! General routines to update primitive variables in boosted frame in which
  ! the stacetime is tramsformed: the time offset is introduced:
  ! Time = Time_boosted + tOffsetPerR*(R - R_0) for R > R_0
  ! Follows description and notation in Sokolov, I.V. and Gombosi, T.I. (2024)
  ! Physics-Based Forecasting of Tomorrow's Solar Wind at 1 AU
  ! Revision history:
  ! 12/18/2024 - Igor Sokolov, implemented vesion as described in the
  !              2024 arxiv paper

  use ModUtilities, ONLY: CON_stop
  implicit none
  PRIVATE ! Except
  real, public :: tOffsetPerR = 0
  real :: Gamma = 1.6666666666666666
  real :: GammaMinus1 = 0.6666666666666666
  real :: InvGammaMinus1 = 1.50
  real :: GammaElectron = 1.6666666666666666
  real :: GammaElectronMinus1 = 0.6666666666666666
  real :: InvGammaElectronMinus1 = 1.50
  real :: Tolerance = 1.0e-6
  real :: Mu0 = 1.0
  integer :: nIter = 40
  ! Public members
  public :: set_param
  public :: update_hydro_var
  public :: update_mhd_var
  public :: update_mhd_pe_var
  public :: fast_mhd_speed
  public :: fast_mhd_pe_speed
  logical, parameter :: DoTest = .false.
  character(LEN=*), parameter :: NameModule = 'ModBoostedFrame'
contains
  !============================================================================
  subroutine set_param(tOffsetPerRin, GammaIn, GammaElectronIn, &
       ToleranceIn, UseSi, nIterMax)
    use ModConst, ONLY: cMu
    ! Prepares the module to update primitive variables in the boosted frame
    real,           intent(in) :: tOffsetPerRin
    real, OPTIONAL, intent(in) :: GammaIn
    real, OPTIONAL, intent(in) :: GammaElectronIn
    real, OPTIONAL, intent(in) :: ToleranceIn
    integer, OPTIONAL, intent(in) :: nIterMax
    logical, OPTIONAL, intent(in) :: UseSi
    character(len=*), parameter:: NameSub = 'set_param'
    !--------------------------------------------------------------------------
    tOffsetPerR = tOffsetPerRin
    if(present(GammaIn))then
       Gamma = GammaIn
       GammaMinus1 = Gamma - 1
       InvGammaMinus1 = 1/GammaMinus1
    end if
    if(present(GammaElectronIn))then
       GammaElectron = GammaElectronIn
       GammaElectronMinus1 = GammaElectron - 1
       InvGammaElectronMinus1 = 1/GammaElectronMinus1
    end if
    if(present(ToleranceIn)) Tolerance = ToleranceIn
    if(present(nIterMax)) nIter = nIterMax
    if(present(UseSi))then
       if(.not.UseSi)call CON_stop(NameModule//':'//NameSub//&
            ': UseSi either should not be present, or should be True')
       Mu0 = cMu
    end if
  end subroutine set_param
  !============================================================================
  subroutine update_hydro_var(Old_V, Source_V, Radial_D, New_V)
    ! Components of primitive variable vectors
    integer, parameter :: Rho_ = 1, Ux_ = 2, Uz_ = 4, P_ = 5
    ! Components of the source vector
    integer, parameter :: RhoUx_ = 2, RhoUz_ = 4, Energy_ = 5
    ! INPUTS:
    ! Primitive variables at the time level t^n
    real, intent(in) :: Old_V(Rho_:P_)
    ! Sum of numerical fluxes and sources, for full or half time step
    real, intent(in) :: Source_V(Rho_:Energy_)
    ! Unit vector of radial direction
    real, intent(in) :: Radial_D(3)
    ! OUTPUTS:
    ! Primitive variables at the time level t^{n+1/2} or t^{n+1}
    real, intent(out) :: New_V(Rho_:P_)
    ! Misc:
    real :: RhoBar, UrOld, UrNew, UperpOld_D(3)
    real :: SourceRho, SourceRhoUr, SourceRhoUperp_D(3), SourceEnergy
    real :: DeltaP
    ! Coefficients of quadratic equation to solve delta p
    real :: Alpha, Beta, Epsilon, Discr

    ! For testing hydro as MHD with zero magnetic field
    real :: OldMhd_V(9), SourceMhd_V(9), NewMhd_V(9)
    character(len=*), parameter:: NameSub = 'update_hydro_var'
    !--------------------------------------------------------------------------

    if(DoTest)then
       OldMhd_V = 0.0; SourceMhd_V = 0.0; NewMhd_V = 0.0
       OldMhd_V(1:4) = Old_V(1:4)
       OldMhd_V(8)   = 0.5*Old_V(5); OldMhd_V(9)   = 0.5*Old_V(5)
       NewMhd_V(1:4) = Old_V(1:4)
       NewMhd_V(8)   = 0.5*Old_V(5); NewMhd_V(9) = 0.5*Old_V(5)
       SourceMhd_V(1:4) = Source_V(1:4)
       SourceMhd_V(8) = 0.2*Source_V(5); SourceMhd_V(9) = 0.8*Source_V(5)
       call update_mhd_pe_var(OldMhd_V, SourceMhd_V, Radial_D, NewMhd_V)
       New_V(1:4) = NewMhd_V(1:4); New_V(5) = NewMhd_V(8) + NewMhd_V(9)
       RETURN
    end if
    SourceRho = Source_V(Rho_)

    ! Separate radial and perpendicular velocity
    UrOld = sum(Radial_D*Old_V(Ux_:Uz_))
    UperpOld_D = Old_V(Ux_:Uz_) - UrOld*Radial_D

    ! Separate radial and perpendicular momentum source
    SourceRhoUr = sum(Radial_D*Source_V(RhoUx_:RhoUz_))
    SourceRhoUperp_D = Source_V(RhoUx_:RhoUz_) - SourceRhoUr*Radial_D

    ! \bar{\rho} = \rho^{n+1}*(1 - u_R^{n+1)*tOffsetPerR), may be found
    ! explicitly, althouth neither \rho^{n+1} nor u_R^{n+1) are not known yet
    RhoBar = SourceRho + Old_V(Rho_)*(1 - UrOld*tOffsetPerR)

    ! As long as both \rho^{n+1} and \rho^{n} are expressed in terms of
    ! \bar{\rho), the sources need to be corrected to account for this
    ! substitution. Radial momentum source corrected and divided by RhoBar:
    SourceRhoUr = (SourceRhoUr - UrOld*SourceRho)/RhoBar
    ! Perpendicular momentum source corrected and divided by RhoBar:
    SourceRhoUperp_D = (SourceRhoUperp_D - UperpOld_D*SourceRho)/RhoBar
    ! Perpendicular velocity may be solved explicitly and put in place:
    New_V(Ux_:Uz_) = SourceRhoUperp_D + UperpOld_D

    ! Radial velocity can be expressed in terms of \delta P = P^{n+1} - P^{n}
    ! u_r^{n+1} = u_r^n + (SourceRhoUr + tOffsetPerR*DeltaP)/RhoBar
    ! On eliminating the velocities, u_R^{n+1} and u_\perp^{n+1) from
    ! the energy equation, we arrive at a quadratic equation for \delta P:
    ! \alpha (\delpa P)**2 - \beta \delta P + \varepsilon = 0
    ! The coefficients in this equation are:
    Alpha = (Gamma + 1)*InvGammaMinus1*tOffsetPerR**2/(2*RhoBar)
    Beta  = (1 - tOffsetPerR*(UrOld + SourceRhoUr) - &
         Gamma*tOffsetPerR**2*Old_V(P_)/RhoBar)*InvGammaMinus1
    if(Beta <= 0.0) call CON_stop(NameModule//':'//NameSub//&
         ': Negative \beta, reduce tOffsetPerR')
    Epsilon = Source_V(Energy_) + 0.5*sum(Old_V(Ux_:Uz_)**2)*SourceRho - &
         sum(Old_V(Ux_:Uz_)*Source_V(RhoUx_:RhoUz_)) + Gamma*&
         InvGammaMinus1*Old_V(P_)*SourceRhoUr*tOffsetPerR - &
         sum((Source_V(RhoUx_:RhoUz_) - Old_V(Ux_:Uz_)*SourceRho)**2)/&
         (2*RhoBar)
    Discr = Beta**2 - 4*Alpha*Epsilon
    if(Discr < 0.0) call CON_stop(NameModule//':'//NameSub//&
         ': Negative determinant, reduce time step or tOffsetPerR')

    ! The only meaningful solution is:
    DeltaP = 2*Epsilon/(Beta + sqrt(Discr))

    ! Finalize:
    ! Solve p^{n+1}:
    New_V(P_) = Old_V(P_) + DeltaP
    ! Solve and add radial velocity:
    UrNew = UrOld + SourceRhoUr + tOffsetPerR*DeltaP/RhoBar
    New_V(Ux_:Uz_) = New_V(Ux_:Uz_) + UrNew*Radial_D
    ! Solve density:
    New_V(Rho_) = RhoBar/(1 - tOffsetPerR*UrNew)
  end subroutine update_hydro_var
  !============================================================================
  subroutine update_mhd_var(Old_V, Source_V, Radial_D, New_V)
    ! Components of primitive variable vectors
    integer, parameter :: Rho_ = 1, Ux_ = 2, Uz_ = 4, Bx_= 5, Bz_ = 7, P_ = 8
    ! Components of the source vector
    integer, parameter :: RhoUx_ = 2, RhoUz_ = 4, Energy_ = 8
    ! Components of modified conserved variables, reduced primitive variables
    ! and conserved variable defects:
    integer, parameter :: Ur_ = 1, BperpX_ = 2, BperpZ_ = 4, IterateP_ = 5, &
         IterateE_ = IterateP_, RhoUr_=Ur_, Bperp_ = BperpZ_
    ! INPUTS:
    ! Primitive variables at the time level t^n
    real, intent(in) :: Old_V(Rho_:P_)
    ! Sum of numerical fluxes and sources, for full or half time step
    real, intent(in) :: Source_V(Rho_:Energy_)
    ! Unit vector of radial direction
    real, intent(in) :: Radial_D(3)
    ! INPUT/OUTPUT:
    ! In: primitive variables used to calculate div B source terms
    ! Out: Primitive variables at the time level t^{n+1/2} or t^{n+1}
    real, intent(inout) :: New_V(Rho_:P_)
    ! Primitive vars at n th time level:
    real :: UrOld, UperpOld_D(3), BrOld, BperpOld_D(3)
    ! Sources:
    real :: SourceRhoUr, SourceRhoUperp_D(3), SourceBr, SourceBperp_D(3)

    ! Iterated variables:
    ! Primitive variables at h th iteration:
    real :: ReducedPrim_V(Ur_:IterateP_)
    ! Increments in primitive variables solved from linear equations
    real :: DeltaPrim_V(Ur_:IterateP_)
    ! Modified conserved variables at h th iteration
    real :: ConsOld_V(Ur_:IterateP_)
    ! Modified conserved variables at h+1 th iteration
    real :: ConsNew_V(Ur_:IterateP_)
    ! Defects in the conserved variables:
    real :: Defect_V(Ur_:IterateP_)
    ! Used in the iterative processes
    real :: DeltaUr, UrIter, InvD_V(Bperp_:IterateP_), &
         AurPrime_V(BperpX_:IterateP_), InvDdotDefect_V(BperpX_:IterateP_), &
         InvDdotAprimeUr_V(BperpX_:IterateP_)

    ! Misc:
    real :: RhoBar, SourceRho, UrNew, BrNew, BrBar, UrBar, vAlfven2R, Discr
    real :: EnergyDens, CrMax, BoostFactor
    integer :: iIter
    character(len=*), parameter:: NameSub = 'update_mhd_var'
    !--------------------------------------------------------------------------

    SourceRho = Source_V(Rho_)

    ! Separate radial and perpendicular velocity
    UrOld = sum(Radial_D*Old_V(Ux_:Uz_))
    UperpOld_D = Old_V(Ux_:Uz_) - UrOld*Radial_D

    ! Separate radial and perpendicular momentum source
    SourceRhoUr = sum(Radial_D*Source_V(RhoUx_:RhoUz_))
    SourceRhoUperp_D = Source_V(RhoUx_:RhoUz_) - SourceRhoUr*Radial_D

     ! Separate radial and perpendicular field
    BrOld = sum(Radial_D*Old_V(Bx_:Bz_))
    BperpOld_D = Old_V(Bx_:Bz_) - BrOld*Radial_D

    ! Separate radial and perpendicular field source
    SourceBr = sum(Radial_D*Source_V(Bx_:Bz_))
    SourceBperp_D = Source_V(Bx_:Bz_) - SourceBr*Radial_D

    ! Solve radial field with Ur used in div B source term (now in New_V!)
    UrBar = sum(New_V(Ux_:Uz_)*Radial_D)
    BrNew = BrOld + SourceBr/(1 - tOffsetPerR*UrBar)
    ! BrBar is an arithmetic average of the old and new Br:
    BrBar = 0.5*(BrOld + BrNew)

    CrMax = UrBar + fast_mhd_speed(New_V, Radial_D)
    BoostFactor = 1/(1 -  tOffsetPerR*CrMax)

    ! \bar{\rho} = \rho^{n+1}*(1 - u_R^{n+1)*tOffsetPerR), may be found
    ! explicitly, althouth neither \rho^{n+1} nor u_R^{n+1) are not known yet
    RhoBar = SourceRho + Old_V(Rho_)*(1 - UrOld*tOffsetPerR)

    ! As long as both \rho^{n+1} and \rho^{n} are expressed in terms of
    ! \bar{\rho), the sources need to be corrected to account for this
    ! substitution. Radial momentum source corrected:
    SourceRhoUr = SourceRhoUr - UrOld*SourceRho
    ! Perpendicular momentum source corrected and divided by RhoBar:
    SourceRhoUperp_D = (SourceRhoUperp_D - UperpOld_D*SourceRho)/RhoBar

    ! Misc: a square of radial Alfven wave speed multiplied by tOffsetPerR**2:
    vAlfven2R = tOffsetPerR**2*BrBar**2/(Mu0*RhoBar)
    ! Prepare iterations:
    ! At 0th iteration the primitive variables are taken from nth time level
    ReducedPrim_V(Ur_) = UrOld
    ReducedPrim_V(BperpX_:BperpZ_) = BperpOld_D
    ReducedPrim_V(IterateP_) = Old_V(P_)

    ! Get modified conserved variables:
    call get_modcons(ReducedPrim_V, ConsOld_V)
    ! Get initial defects in the conserved variables:
    Defect_V(RhoUr_) = SourceRhoUr
    Defect_V(BperpX_:BperpZ_) = SourceBperp_D - &
         BrBar*tOffsetPerR*SourceRhoUperp_D
    Defect_V(IterateE_) = Source_V(Energy_)                   + &
         0.5*SourceRho*sum(Old_V(Ux_:Uz_)**2)                 - &
         sum(Old_V(Ux_:Uz_)*Source_V(RhoUx_:RhoUz_))          - &
         (BrBar*SourceBr + sum(SourceBperp_D*BperpOld_D))/Mu0 - &
         0.5*RhoBar*sum(SourceRhoUperp_D**2)
    ! write(*,*)'Internal energy defect=',Defect_V(IterateE_)
    ! Recover true energy defect from modified one
    Defect_V(IterateE_) = Defect_V(IterateE_) + &
         ReducedPrim_V(Ur_)*Defect_V(RhoUr_)  + &
         sum(ReducedPrim_V(BperpX_:BperpZ_)*Defect_V(BperpX_:BperpZ_))/Mu0

    ! write(*,*)'Defects=', Defect_V
    EnergyDens = RhoBar*0.50*ReducedPrim_V(Ur_)**2      + &
         sum(ReducedPrim_V(BperpX_:BperpZ_)**2)/(2*Mu0) + &
         InvGammaMinus1*ReducedPrim_V(IterateP_)        + &
         Defect_V(IterateE_)
    ReducedPrim_V(Ur_) = ReducedPrim_V(Ur_) + Defect_V(RhoUr_)/RhoBar
    ReducedPrim_V(BperpX_:BperpZ_) = ReducedPrim_V(BperpX_:BperpZ_) + &
         Defect_V(BperpX_:BperpZ_)
    ReducedPrim_V(IterateP_) = GammaMinus1*(EnergyDens  - &
         RhoBar*0.50*ReducedPrim_V(Ur_)**2              - &
         sum(ReducedPrim_V(BperpX_:BperpZ_)**2)/(2*Mu0))
    ! Modify the energy defect
    Defect_V = Defect_V + ConsOld_V
    ! Get modified conserved variables:
    call get_modcons(ReducedPrim_V, ConsOld_V)
    ! Modify the energy defect
    Defect_V = Defect_V - ConsOld_V
    ! write(*,*)'Defects=', Defect_V
    ! Recover modified energy defect from true one
    Defect_V(IterateE_) = Defect_V(IterateE_) - &
         ReducedPrim_V(Ur_)*Defect_V(RhoUr_)  - &
         sum(ReducedPrim_V(BperpX_:BperpZ_)*Defect_V(BperpX_:BperpZ_))/Mu0
    ! write(*,*)'Internal energy defect=',Defect_V(IterateE_)
    iIter = 0
    ITERATIONS: do
       ! Solve a system of linear equations
       ! RhoBar*\delta u_R + A_{u_R,Prim}\cdot\delta Prime = Defect_{\rho u_R}
       ! A_{Prim,u_R)\delta u_R + D\cdot \delta Prime = Defect_{Prim},
       ! where D is a diagonal 4x4 matrix

       ! Get reusable u_R*tOffsetPerR:
       UrNew = tOffsetPerR*ReducedPrim_V(Ur_)

       ! Calculate elements of the invere diagonal D
       InvD_V(Bperp_) = 1/(1 - UrNew - vAlfven2R)
       InvD_V(IterateP_) = GammaMinus1/(1 - UrNew)

       ! Calculate vectors present in the formal solution of the last 4 eqs,
       ! A_{Prim,u_R)\delta u_R \delta Prime = D^{-1}\cdot Defect_{Prim}
       InvDdotDefect_V(BperpX_:BperpZ_) = InvD_V(Bperp_)*&
            Defect_V(BperpX_:BperpZ_)
       InvDdotDefect_V(IterateP_) = InvD_V(IterateP_)*Defect_V(IterateP_)

       InvDdotAprimeUr_V(BperpX_:BperpZ_) = InvD_V(Bperp_)*&
            (-tOffsetPerR*ReducedPrim_V(BperpX_:BperpZ_))
       InvDdotAprimeUr_V(IterateP_) = InvD_V(IterateP_)*(-Gamma*&
            tOffsetPerR*ReducedPrim_V(IterateP_)*InvGammaMinus1)

       ! Vector of derivatives of total pressure over primitives:
       AurPrime_V(BperpX_:BperpZ_) = -tOffsetPerR*&
            ReducedPrim_V(BperpX_:BperpZ_)/Mu0
       AurPrime_V(IterateP_) = -tOffsetPerR

       Discr = RhoBar - sum(AurPrime_V*InvDdotAprimeUr_V)
       ! write(*,*)'iIter=', iIter,' Discr/RhoBar', Discr/RhoBar
       if (Discr <= 0.0) then
          write(*,'(a,i6)')'At iter=',iIter
          call CON_stop(NameModule//':'//NameSub//&
               ': Negative discriminator, reduce tOffsetPerR')
       end if

       ! Solve \delta u_R
       DeltaUr = (Defect_V(RhoUr_) - sum(AurPrime_V*InvDdotDefect_V))/Discr
       ! Solve all other \delta Prim:
       DeltaPrim_V(BperpX_:) = InvDdotDefect_V - InvDdotAprimeUr_V*DeltaUr
       DeltaPrim_V(Ur_) = DeltaUr
       ! write(*,*)'DeltaPrim_V=', DeltaPrim_V
       if(tOffsetPerR*abs(DeltaUr) <= Tolerance) then
          ReducedPrim_V = ReducedPrim_V + DeltaPrim_V
          EXIT ITERATIONS
       else
          ! Prepare next iteration
          iIter = iIter + 1
          if(iIter > nIter) call CON_stop(NameModule//':'//NameSub//&
               ': No iteration convergence')
          ! Recover true energy defect from modified one
          Defect_V(IterateE_) = Defect_V(IterateE_) + &
               ReducedPrim_V(Ur_)*Defect_V(RhoUr_)  + &
               sum(ReducedPrim_V(BperpX_:BperpZ_)*Defect_V(BperpX_:BperpZ_))&
               /Mu0
          ! Update reduced primitive variables:
          ReducedPrim_V = ReducedPrim_V + DeltaPrim_V
          ! Update modified conserved variables:
          call get_modcons(ReducedPrim_V, ConsNew_V)
          ! Update defects in the modified conserved variables
          Defect_V = Defect_V + (ConsOld_V - ConsNew_V)
          ! write(*,*)'Defects=', Defect_V
          ! Modify the energy equation and, accordingly, the energy defect,
          ! to facilitate linearization:
          Defect_V(IterateE_) = Defect_V(IterateE_) - &
               ReducedPrim_V(Ur_)*Defect_V(RhoUr_)  - &
               sum(ReducedPrim_V(BperpX_:BperpZ_)*Defect_V(BperpX_:BperpZ_))&
               /Mu0
          ! write(*,*)'Internal energy defect=',Defect_V(IterateE_)
          ! Store the modified conserved variables:
          ConsOld_V = ConsNew_V
       end if
    end do ITERATIONS
    ! Finalize:
    ! Put radial velocity:
    UrNew = ReducedPrim_V(Ur_)
    New_V(Ux_:Uz_) = Radial_D*UrNew
    ! Put perpendicular field:
    New_V(Bx_:Bz_) = ReducedPrim_V(BperpX_:BperpZ_)
    ! Add radial field:
    New_V(Bx_:Bz_) = New_V(Bx_:Bz_) + Radial_D*BrNew
    ! Put pressure:
    New_V(P_) = ReducedPrim_V(IterateP_)
    ! Add perpendicular velocity
    New_V(Ux_:Uz_) = New_V(Ux_:Uz_) + UperpOld_D + SourceRhoUperp_D &
         - tOffsetPerR*BrBar/(Mu0*RhoBar)*&
         (ReducedPrim_V(BperpX_:BperpZ_) - BperpOld_D)
    ! Solve density:
    New_V(Rho_) = RhoBar/(1 - tOffsetPerR*UrNew)
  contains
    !==========================================================================
    subroutine get_modcons(ReducedPrim_V, ModifiedCons_V)
      ! Solves iteration for the modified conservative variables
      ! in terms of the reduced set of primitive variables:
      ! INPUTS:
      real, intent(in) :: ReducedPrim_V(Ur_:IterateP_)
      ! OUTPUTS:
      real, intent(out) :: ModifiedCons_V(RhoUr_:IterateE_)
      ! Misc:
      real :: Ur
      !------------------------------------------------------------------------

      Ur = ReducedPrim_V(Ur_)
      ModifiedCons_V(RhoUr_) = RhoBar*Ur - tOffsetPerR*(&
           sum(ReducedPrim_V(BperpX_:BperpZ_)**2)/(2*Mu0) + &
           ReducedPrim_V(IterateP_))
      ! Kinetic energy
      ModifiedCons_V(IterateE_) = 0.5*RhoBar*Ur**2
      ! To simplify computations:
      Ur = Ur*tOffsetPerR
      ModifiedCons_V(BperpX_:BperpZ_) = ReducedPrim_V(BperpX_:BperpZ_)*&
           (1 - Ur - vAlfven2R)
      ! Add magnetic and internal energy
      ModifiedCons_V(IterateE_) =  ModifiedCons_V(IterateE_) &
           + (1 - 2*Ur - vAlfven2R)*&
           sum(ReducedPrim_V(BperpX_:BperpZ_)**2)/(2*Mu0) &
           + (1 - Gamma*Ur)*ReducedPrim_V(IterateP_)*InvGammaMinus1
    end subroutine get_modcons
    !==========================================================================
  end subroutine update_mhd_var
  !============================================================================
  subroutine update_mhd_pe_var(Old_V, Source_V, Radial_D, New_V)
    ! Components of primitive variable vectors
    integer, parameter :: Rho_ = 1, Ux_ = 2, Uz_ = 4, Bx_= 5, Bz_ = 7, &
         Pe_ = 8, P_ = 9
    ! Components of the source vector
    integer, parameter :: RhoUx_ = 2, RhoUz_ = 4, Energy_ = 9
    ! Components of modified conserved variables, reduced primitive variables
    ! and conserved variable defects:
    integer, parameter :: Ur_ = 1, BperpX_ = 2, BperpZ_ = 4, IteratePe_ = 5, &
         IterateP_ = 6, IterateE_ = IterateP_, RhoUr_=Ur_, Bperp_ = BperpZ_
    ! INPUTS:
    ! Primitive variables at the time level t^n
    real, intent(in) :: Old_V(Rho_:P_)
    ! Sum of numerical fluxes and sources, for full or half time step
    real, intent(in) :: Source_V(Rho_:Energy_)
    ! Unit vector of radial direction
    real, intent(in) :: Radial_D(3)
    ! INPUT/OUTPUT:
    ! In: primitive variables used to calculate div B source terms
    ! Out: Primitive variables at the time level t^{n+1/2} or t^{n+1}
    real, intent(inout) :: New_V(Rho_:P_)
    ! Primitive vars at n th time level:
    real :: UrOld, UperpOld_D(3), BrOld, BperpOld_D(3)
    ! Sources:
    real :: SourceRhoUr, SourceRhoUperp_D(3), SourceBr, SourceBperp_D(3)

    ! Iterated variables:
    ! Primitive variables at h th iteration:
    real :: ReducedPrim_V(Ur_:IterateP_)
    ! Increments in primitive variables solved from linear equations
    real :: DeltaPrim_V(Ur_:IterateP_)
    ! Modified conserved variables at h th iteration
    real :: ConsOld_V(Ur_:IterateP_)
    ! Modified conserved variables at h+1 th iteration
    real :: ConsNew_V(Ur_:IterateP_)
    ! Defects in the conserved variables:
    real :: Defect_V(Ur_:IterateP_)
    ! Used in the iterative processes
    real :: DeltaUr, UrIter, InvD_V(Bperp_:IterateP_), &
         AurPrime_V(BperpX_:IterateP_), InvDdotDefect_V(BperpX_:IterateP_), &
         InvDdotAprimeUr_V(BperpX_:IterateP_)

    ! Misc:
    real :: RhoBar, PeBar, SourceRho, UrNew, BrNew, BrBar, UrBar, vAlfven2R, &
         Discr, EnergyDens, CrMax, BoostFactor
    integer :: iIter

    character(len=*), parameter:: NameSub = 'update_mhd_pe_var'
    !--------------------------------------------------------------------------
    SourceRho = Source_V(Rho_)
    PeBar = New_V(Pe_)
    ! Separate radial and perpendicular velocity
    UrOld = sum(Radial_D*Old_V(Ux_:Uz_))
    UperpOld_D = Old_V(Ux_:Uz_) - UrOld*Radial_D

    ! Separate radial and perpendicular momentum source
    SourceRhoUr = sum(Radial_D*Source_V(RhoUx_:RhoUz_))
    SourceRhoUperp_D = Source_V(RhoUx_:RhoUz_) - SourceRhoUr*Radial_D

     ! Separate radial and perpendicular field
    BrOld = sum(Radial_D*Old_V(Bx_:Bz_))
    BperpOld_D = Old_V(Bx_:Bz_) - BrOld*Radial_D

    ! Separate radial and perpendicular field source
    SourceBr = sum(Radial_D*Source_V(Bx_:Bz_))
    SourceBperp_D = Source_V(Bx_:Bz_) - SourceBr*Radial_D

    ! Solve radial field with Ur used in div B source term (now in New_V!)
    UrBar = sum(Radial_D*New_V(Ux_:Uz_))
    BrNew = BrOld + SourceBr/(1 - tOffsetPerR*UrBar)
    ! BrBar is an arithmetic average of the old and new Br:
    BrBar = 0.5*(BrOld + BrNew)

    CrMax = UrBar + fast_mhd_pe_speed(New_V, Radial_D)
    BoostFactor = 1!/(1 -  tOffsetPerR*CrMax)

    ! \bar{\rho} = \rho^{n+1}*(1 - u_R^{n+1)*tOffsetPerR), may be found
    ! explicitly, althouth neither \rho^{n+1} nor u_R^{n+1) are not known yet
    RhoBar = SourceRho + Old_V(Rho_)*(1 - UrOld*tOffsetPerR)

    ! As long as both \rho^{n+1} and \rho^{n} are expressed in terms of
    ! \bar{\rho), the sources need to be corrected to account for this
    ! substitution. Radial momentum source corrected:
    SourceRhoUr = SourceRhoUr - UrOld*SourceRho
    ! Perpendicular momentum source corrected and divided by RhoBar:
    SourceRhoUperp_D = (SourceRhoUperp_D - UperpOld_D*SourceRho)/RhoBar

    ! Misc: a square of radial Alfven wave speed multiplied by tOffsetPerR**2:
    vAlfven2R = tOffsetPerR**2*BrBar**2/(Mu0*RhoBar)
    ! Prepare iterations:
    ! At 0th iteration the primitive variables are taken from nth time level
    ReducedPrim_V(Ur_) = UrOld
    ReducedPrim_V(BperpX_:BperpZ_) = BperpOld_D
    ReducedPrim_V(IteratePe_) = Old_V(Pe_)
    ReducedPrim_V(IterateP_) = Old_V(P_)

    ! Get modified conserved variables:
    call get_modcons(ReducedPrim_V, ConsOld_V)
    ! Get initial defects in the conserved variables:
    Defect_V(RhoUr_) = SourceRhoUr
    Defect_V(BperpX_:BperpZ_) = SourceBperp_D - &
         BrBar*tOffsetPerR*SourceRhoUperp_D
    Defect_V(IteratePe_) = Source_V(Pe_)
    Defect_V(IterateE_) = Source_V(Energy_)                   + &
         0.5*SourceRho*sum(Old_V(Ux_:Uz_)**2)                 - &
         sum(Old_V(Ux_:Uz_)*Source_V(RhoUx_:RhoUz_))          - &
         (BrBar*SourceBr + sum(SourceBperp_D*BperpOld_D))/Mu0 - &
         0.5*RhoBar*sum(SourceRhoUperp_D**2)

    ! Recover true energy defect from modified one
    Defect_V(IterateE_) = Defect_V(IterateE_) + &
         ReducedPrim_V(Ur_)*Defect_V(RhoUr_)  + &
         sum(ReducedPrim_V(BperpX_:BperpZ_)*Defect_V(BperpX_:BperpZ_))/Mu0

    EnergyDens = RhoBar*0.50*ReducedPrim_V(Ur_)**2      + &
         sum(ReducedPrim_V(BperpX_:BperpZ_)**2)/(2*Mu0) + &
         InvGammaMinus1*ReducedPrim_V(IterateP_)        + &
         Defect_V(IterateE_)*BoostFactor
    ReducedPrim_V(Ur_) = ReducedPrim_V(Ur_) + &
         BoostFactor*Defect_V(RhoUr_)/RhoBar
    ReducedPrim_V(BperpX_:BperpZ_) = ReducedPrim_V(BperpX_:BperpZ_) + &
         Defect_V(BperpX_:BperpZ_)*BoostFactor
    ReducedPrim_V(IteratePe_) = ReducedPrim_V(IteratePe_) + &
         GammaElectronMinus1*Defect_V(IteratePe_)*BoostFactor
    ReducedPrim_V(IterateP_) = GammaMinus1*(EnergyDens  - &
         RhoBar*0.50*ReducedPrim_V(Ur_)**2              - &
         sum(ReducedPrim_V(BperpX_:BperpZ_)**2)/(2*Mu0))
    ! Modify the energy defect
    Defect_V = Defect_V + ConsOld_V
    ! Get modified conserved variables:
    call get_modcons(ReducedPrim_V, ConsOld_V)
    ! Modify the energy defect
    Defect_V = Defect_V - ConsOld_V
    ! Recover modified energy defect from true one
    Defect_V(IterateE_) = Defect_V(IterateE_) - &
         ReducedPrim_V(Ur_)*Defect_V(RhoUr_)  - &
         sum(ReducedPrim_V(BperpX_:BperpZ_)*Defect_V(BperpX_:BperpZ_))/Mu0
    iIter = 0
    ITERATIONS: do
       ! Solve a systen of linear equations
       ! RhoBar*\delta u_R + A_{u_R,Prim}\cdot\delta Prime = Defect_{\rho u_R}
       ! A_{Prim,u_R)\delta u_R + D\cdot \delta Prime = Defect_{Prim},
       ! where D is a diagonal 4x4 matrix
       ReducedPrim_V(IteratePe_:IterateP_) = &
            max(ReducedPrim_V(IteratePe_:IterateP_),0.0)
       ! Get reusable u_R*tOffsetPerR:
       UrNew = tOffsetPerR*ReducedPrim_V(Ur_)

       ! Calculate elements of the invere diagonal D
       InvD_V(Bperp_) = 1/(1 - UrNew - vAlfven2R)
       InvD_V(IteratePe_) = GammaElectronMinus1/(1 - UrNew)
       InvD_V(IterateP_) = GammaMinus1/(1 - UrNew)

       ! Calculate vectors present in the formal solution of the last 4 eqs,
       ! A_{Prim,u_R)\delta u_R \delta Prime = D^{-1}\cdot Defect_{Prim}
       InvDdotDefect_V(BperpX_:BperpZ_) = InvD_V(Bperp_)*&
            Defect_V(BperpX_:BperpZ_)
       InvDdotDefect_V(IteratePe_) = InvD_V(IteratePe_)*Defect_V(IteratePe_)
       InvDdotDefect_V(IterateP_) = InvD_V(IterateP_)*Defect_V(IterateP_)

       InvDdotAprimeUr_V(BperpX_:BperpZ_) = InvD_V(Bperp_)*&
            (-tOffsetPerR*ReducedPrim_V(BperpX_:BperpZ_))
       InvDdotAprimeUr_V(IteratePe_) = InvD_V(IteratePe_)*(&
            -tOffsetPerR*ReducedPrim_V(IteratePe_)*InvGammaElectronMinus1 &
            -tOffsetPerR*PeBar)
       InvDdotAprimeUr_V(IterateP_) = InvD_V(IterateP_)*(-Gamma*&
            tOffsetPerR*ReducedPrim_V(IterateP_)*InvGammaMinus1 &
            -tOffsetPerR*(ReducedPrim_V(IteratePe_) - PeBar))

       ! Vector of derivatives of total pressure over primitives:
       AurPrime_V(BperpX_:BperpZ_) = -tOffsetPerR*&
            ReducedPrim_V(BperpX_:BperpZ_)/Mu0
       AurPrime_V(IteratePe_:IterateP_) = -tOffsetPerR
       Discr = RhoBar - sum(AurPrime_V*InvDdotAprimeUr_V)

       if (Discr <= 0.0) then
          write(*,*)'Old_V=', Old_V
          write(*,*)'Source_V=',Source_V
          write(*,*)'Radial_D=',Radial_D
          write(*,*)'New_V=',New_V
          write(*,'(a,i6)')'At iter=',iIter
          ! Finalize:
          ! Put radial velocity:
          UrNew = ReducedPrim_V(Ur_)
          New_V(Ux_:Uz_) = Radial_D*UrNew
          ! Put perpendicular field:
          New_V(Bx_:Bz_) = ReducedPrim_V(BperpX_:BperpZ_)
          ! Add radial field:
          New_V(Bx_:Bz_) = New_V(Bx_:Bz_) + Radial_D*BrNew
          ! Put pressure:
          New_V(Pe_:P_) = ReducedPrim_V(IteratePe_:IterateP_)
          ! Add perpendicular velocity
          New_V(Ux_:Uz_) = New_V(Ux_:Uz_) + UperpOld_D + SourceRhoUperp_D &
               - tOffsetPerR*BrBar/(Mu0*RhoBar)*&
               (ReducedPrim_V(BperpX_:BperpZ_) - BperpOld_D)
          ! Solve density:
          New_V(Rho_) = RhoBar/(1 - tOffsetPerR*UrNew)
          write(*,*) 'Fast magnetosonic speed=', UrNew + &
               fast_mhd_pe_speed(New_V, Radial_D)
          write(*,*) 'Maximum allowed radial perturbation speed=', 1/&
               tOffsetPerR

          call CON_stop(NameModule//':'//NameSub//&
               ': Negative discriminator, reduce tOffsetPerR')
       end if

       ! Solve \delta u_R
       DeltaUr = (Defect_V(RhoUr_) - sum(AurPrime_V*InvDdotDefect_V))/Discr
       ! Solve all other \delta Prim:
       DeltaPrim_V(BperpX_:) = InvDdotDefect_V - InvDdotAprimeUr_V*DeltaUr
       DeltaPrim_V(Ur_) = DeltaUr

       if(tOffsetPerR*abs(DeltaUr) <= Tolerance) then
          ReducedPrim_V = ReducedPrim_V + DeltaPrim_V
          ReducedPrim_V(IteratePe_:IterateP_) = &
               max(ReducedPrim_V(IteratePe_:IterateP_),0.0)
          EXIT ITERATIONS
       else
          ! Prepare next iteration
          iIter = iIter + 1
          if(iIter > nIter) call CON_stop(NameModule//':'//NameSub//&
               ': No iteration convergence')
          ! Recover true energy defect from modified one
          Defect_V(IterateE_) = Defect_V(IterateE_) + &
               ReducedPrim_V(Ur_)*Defect_V(RhoUr_)  + &
               sum(ReducedPrim_V(BperpX_:BperpZ_)*Defect_V(BperpX_:BperpZ_))&
               /Mu0
          ! Update reduced primitive variables:
          ReducedPrim_V = ReducedPrim_V + DeltaPrim_V
          ReducedPrim_V(IteratePe_:IterateP_) = &
               max(ReducedPrim_V(IteratePe_:IterateP_),0.0)
          if(ReducedPrim_V(IterateP_) < 0.0)call CON_stop(&
               NameSub//': negative pressure')
          ! Update modified conserved variables:
          call get_modcons(ReducedPrim_V, ConsNew_V)
          ! Update defects in the modified conserved variables
          Defect_V = Defect_V + (ConsOld_V - ConsNew_V)
          ! Modify the energy equation and, accordingly, the energy defect,
          ! to facilitate linearization:
          Defect_V(IterateE_) = Defect_V(IterateE_) - &
               ReducedPrim_V(Ur_)*Defect_V(RhoUr_)  - &
               sum(ReducedPrim_V(BperpX_:BperpZ_)*Defect_V(BperpX_:BperpZ_))&
               /Mu0
          ! Store the modified conserved variables:
          ConsOld_V = ConsNew_V
       end if
    end do ITERATIONS
    ! Finalize:
    ! Put radial velocity:
    UrNew = ReducedPrim_V(Ur_)
    New_V(Ux_:Uz_) = Radial_D*UrNew
    ! Put perpendicular field:
    New_V(Bx_:Bz_) = ReducedPrim_V(BperpX_:BperpZ_)
    ! Add radial field:
    New_V(Bx_:Bz_) = New_V(Bx_:Bz_) + Radial_D*BrNew
    ! Put pressure:
    New_V(Pe_:P_) = ReducedPrim_V(IteratePe_:IterateP_)
    ! Add perpendicular velocity
    New_V(Ux_:Uz_) = New_V(Ux_:Uz_) + UperpOld_D + SourceRhoUperp_D &
         - tOffsetPerR*BrBar/(Mu0*RhoBar)*&
         (ReducedPrim_V(BperpX_:BperpZ_) - BperpOld_D)
    ! Solve density:
    New_V(Rho_) = RhoBar/(1 - tOffsetPerR*UrNew)
  contains
    !==========================================================================
    subroutine get_modcons(ReducedPrim_V, ModifiedCons_V)
      ! Solves iteration for the modified conservative variables
      ! in terms of the reduced set of primitive variables:
      ! INPUTS:
      real, intent(in) :: ReducedPrim_V(Ur_:IterateP_)
      ! OUTPUTS:
      real, intent(out) :: ModifiedCons_V(RhoUr_:IterateE_)
      ! Misc:
      real :: Ur
      !------------------------------------------------------------------------

      Ur = ReducedPrim_V(Ur_)
      ModifiedCons_V(RhoUr_) = RhoBar*Ur - tOffsetPerR*(&
           sum(ReducedPrim_V(BperpX_:BperpZ_)**2)/(2*Mu0) + &
           ReducedPrim_V(IteratePe_) + ReducedPrim_V(IterateP_))
      ! Kinetic energy
      ModifiedCons_V(IterateE_) = 0.5*RhoBar*Ur**2
      ! To simplify computations:
      Ur = Ur*tOffsetPerR
      ModifiedCons_V(BperpX_:BperpZ_) = ReducedPrim_V(BperpX_:BperpZ_)*&
           (1 - Ur - vAlfven2R)
      ModifiedCons_V(IteratePe_) = (1 - Ur)*ReducedPrim_V(IteratePe_)* &
           InvGammaElectronMinus1 - Ur*PeBar
      ! Add magnetic and internal energy
      ModifiedCons_V(IterateE_) =  ModifiedCons_V(IterateE_) &
           + (1 - 2*Ur - vAlfven2R)*&
           sum(ReducedPrim_V(BperpX_:BperpZ_)**2)/(2*Mu0) &
           + (1 - Gamma*Ur)*ReducedPrim_V(IterateP_)*InvGammaMinus1 &
           -Ur*(ReducedPrim_V(IteratePe_) - PeBar)
    end subroutine get_modcons
    !==========================================================================
  end subroutine update_mhd_pe_var
  !============================================================================
  real function fast_mhd_speed(Prim_V, Dir_D)
    ! Components of primitive variable vectors
    integer, parameter :: Rho_ = 1, Ux_ = 2, Uz_ = 4, Bx_= 5, Bz_ = 7, P_ = 8
    real, intent(in) :: Prim_V(Rho_:P_), Dir_D(3)
    ! Square of speed of sound, of Alfven wave speed, of it projection on n
    real :: Cs2, Va2, VaN2
    ! Their total; inverse density:
    real :: V2Tot, InvRho
    !--------------------------------------------------------------------------
    InvRho = 1/Prim_V(Rho_)
    ! Squares of speeds:
    Cs2 = Gamma*Prim_V(P_)*InvRho
    Va2 = sum(Prim_V(Bx_:Bz_)**2)*InvRho
    VaN2 = sum(Prim_V(Bx_:Bz_)*Dir_D)**2*InvRho
    V2Tot = Cs2 + Va2
    fast_mhd_speed = sqrt(0.50*(V2Tot + &
         sqrt(max(V2Tot**2 - 4*Cs2*VaN2,0.0))))
  end function fast_mhd_speed
  !============================================================================
  real function fast_mhd_pe_speed(Prim_V, Dir_D)
    ! Components of primitive variable vectors
    integer, parameter :: Rho_ = 1, Ux_ = 2, Uz_ = 4, Bx_= 5, Bz_ = 7, &
         Pe_ = 8, P_ = 9
    real, intent(in) :: Prim_V(Rho_:P_), Dir_D(3)
    ! Square of speed of sound, of Alfven wave speed, of it projection on n
    real :: Cs2, Va2, VaN2
    ! Their total; inverse density:
    real :: V2Tot, InvRho
    !--------------------------------------------------------------------------
    InvRho = 1/Prim_V(Rho_)
    ! Squares of speeds:
    Cs2 = (Gamma*Prim_V(P_) + GammaElectron*Prim_V(Pe_))*InvRho
    Va2 = sum(Prim_V(Bx_:Bz_)**2)*InvRho
    VaN2 = sum(Prim_V(Bx_:Bz_)*Dir_D)**2*InvRho
    V2Tot = Cs2 + Va2
    fast_mhd_pe_speed = sqrt(0.50*(V2Tot + &
         sqrt(max(V2Tot**2 - 4*Cs2*VaN2,0.0))))
  end function fast_mhd_pe_speed
  !============================================================================
end module ModBoostedFrame
!==============================================================================
