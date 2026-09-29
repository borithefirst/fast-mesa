! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

module eos_support

  use const_def, only: dp, ln10, arg_not_provided
  use star_private_def
  use utils_lib, only : is_bad, mesa_error

  implicit none

  private
  public :: get_eos
  public :: solve_eos_given_DE
  public :: solve_eos_given_DEgas
  public :: solve_eos_given_DP
  public :: solve_eos_given_DS
  public :: solve_eos_given_PT
  public :: solve_eos_given_PgasT
  public :: solve_eos_given_PgasT_auto
  public :: get_eos_memo, eos_memo_store, eos_memo_prepare, eos_memo_count, eos_memo_stats_on  ! savethesun eos memo
  public :: eos_fixd_minstep, eos_fixd_minx, eos_dxa_tangent_on, eos_dxa_tangent_partials
  public :: eos_skye_dxa_on, eos_skye_dxa_check, eos_skye_dxa_check_cell

  integer, parameter :: MAX_ITER_FOR_SOLVE = 100

  ! savethesun eos memo: per-cell memo of get_eos keyed on the exact bits of all inputs (opt-in: MESA_EOS_MEMO=1)
  integer, parameter :: memo_slots = 2, memo_ncount = 16
  logical, save :: memo_on = .false., eos_memo_stats_on = .false., memo_checked = .false.
  integer, save :: memo_nz = 0, memo_species = -1, memo_id = -1, memo_handle = -1, memo_last_model = -1
  logical, allocatable, save :: memo_valid(:,:)
  integer(8), allocatable, save :: memo_key(:,:,:)
  real(dp), allocatable, save :: memo_res(:,:,:), memo_dlnd(:,:,:), memo_dlnT(:,:,:)
  integer(8), save :: memo_count(memo_ncount) = 0
  real(dp), save :: eos_fixd_minstep = 0  ! diagnostic, MESA_FIXD_MINSTEP (see hydro_eqns)
  real(dp), save :: eos_fixd_minx = 0  ! diagnostic, MESA_FIXD_MINX: forced FD only if xa_start >= this
  ! lagged tangent composition partials (opt-in: MESA_EOS_DXA_TANGENT=1, see eos_dxa_tangent_partials)
  logical, save :: eos_dxa_tangent_on = .false., tangent_lag = .true.
  ! Skye analytic d_dxa in use (MESA_SKYE_DXA=1); check against central differences (MESA_SKYE_DXA_CHECK=1)
  logical, save :: eos_skye_dxa_on = .false., eos_skye_dxa_check = .false.
  ! check stats per species: max |analytic - FD| for lnE, lnPgas; max |FD|; number of cells
  real(dp), save :: chk_err_E(64) = 0, chk_err_P(64) = 0, chk_max_E(64) = 0, chk_max_P(64) = 0
  integer(8), save :: chk_cells = 0
  real(dp), save :: tangent_minx = 0, tangent_h = 1d-6, tangent_xtol = huge(1d0), tangent_ttol = huge(1d0)
  integer, allocatable, save :: tan_stamp(:)
  logical, allocatable, save :: tan_set(:,:)
  real(dp), allocatable, save :: tan_dlnE(:,:), tan_dlnPgas(:,:), tan_xa(:,:), tan_lnT(:), tan_lnd(:)

contains

  ! savethesun eos memo: called from serial code once per eval_equ; (re)allocates, reads env, prints stats
  subroutine eos_memo_prepare(s)
    use eos_def, only: num_eos_basic_results
    use chem_def, only: chem_isos
    type (star_info), pointer :: s
    character(len=16) :: v
    integer :: st
    if (.not. memo_checked) then
       call get_environment_variable('MESA_EOS_MEMO', v, status=st)
       memo_on = (st == 0 .and. len_trim(v) > 0 .and. trim(v) /= '0')
       call get_environment_variable('MESA_EOS_MEMO_STATS', v, status=st)
       eos_memo_stats_on = (st == 0 .and. len_trim(v) > 0 .and. trim(v) /= '0')
       call get_environment_variable('MESA_FIXD_MINSTEP', v, status=st)
       if (st == 0 .and. len_trim(v) > 0) read(v,*) eos_fixd_minstep
       call get_environment_variable('MESA_FIXD_MINX', v, status=st)
       if (st == 0 .and. len_trim(v) > 0) read(v,*) eos_fixd_minx
       call get_environment_variable('MESA_EOS_DXA_TANGENT', v, status=st)
       eos_dxa_tangent_on = (st == 0 .and. len_trim(v) > 0 .and. trim(v) /= '0')
       call get_environment_variable('MESA_EOS_DXA_MINX', v, status=st)
       if (st == 0 .and. len_trim(v) > 0) read(v,*) tangent_minx
       call get_environment_variable('MESA_EOS_DXA_H', v, status=st)
       if (st == 0 .and. len_trim(v) > 0) read(v,*) tangent_h
       call get_environment_variable('MESA_EOS_DXA_LAG', v, status=st)
       if (st == 0 .and. len_trim(v) > 0) tangent_lag = (trim(v) /= '0')
       call get_environment_variable('MESA_SKYE_DXA', v, status=st)
       eos_skye_dxa_on = (st == 0 .and. len_trim(v) > 0 .and. trim(v) /= '0')
       call get_environment_variable('MESA_SKYE_DXA_CHECK', v, status=st)
       eos_skye_dxa_check = (st == 0 .and. len_trim(v) > 0 .and. trim(v) /= '0')
       call get_environment_variable('MESA_EOS_DXA_XTOL', v, status=st)
       if (st == 0 .and. len_trim(v) > 0) read(v,*) tangent_xtol
       call get_environment_variable('MESA_EOS_DXA_TTOL', v, status=st)
       if (st == 0 .and. len_trim(v) > 0) read(v,*) tangent_ttol
       memo_checked = .true.
    end if
    if ((eos_memo_stats_on .or. eos_skye_dxa_check) .and. s% model_number /= memo_last_model) then
       if (eos_memo_stats_on) write(*,'(a,i8,16i12)') 'EOSMEMO ', s% model_number, memo_count
       if (eos_skye_dxa_check) then
          write(*,'(a,i8,i10)') 'SKYEDXA cells ', s% model_number, chk_cells
          do st = 1, min(s% species, 64)
             write(*,'(a,a8,4es12.3)') 'SKYEDXA ', trim(chem_isos% name(s% chem_id(st))), &
                chk_err_E(st), chk_max_E(st), chk_err_P(st), chk_max_P(st)
          end do
       end if
       memo_last_model = s% model_number
    end if
    if (.not. (memo_on .or. eos_dxa_tangent_on)) return
    if (memo_id /= s% id .or. memo_handle /= s% eos_handle .or. &
        memo_species /= s% species .or. memo_nz < s% nz) then
       if (allocated(memo_valid)) deallocate(memo_valid, memo_key, memo_res, memo_dlnd, memo_dlnT, &
          tan_stamp, tan_set, tan_dlnE, tan_dlnPgas, tan_xa, tan_lnT, tan_lnd)
       memo_nz = s% nz + s% nz/4 + 16
       memo_species = s% species
       memo_id = s% id
       memo_handle = s% eos_handle
       allocate(memo_valid(memo_nz, memo_slots), memo_key(4 + memo_species, memo_nz, memo_slots), &
          memo_res(num_eos_basic_results, memo_nz, memo_slots), &
          memo_dlnd(num_eos_basic_results, memo_nz, memo_slots), &
          memo_dlnT(num_eos_basic_results, memo_nz, memo_slots))
       allocate(tan_stamp(memo_nz), tan_set(memo_species, memo_nz), &
          tan_dlnE(memo_species, memo_nz), tan_dlnPgas(memo_species, memo_nz), &
          tan_xa(memo_species, memo_nz), tan_lnT(memo_nz), tan_lnd(memo_nz))
       memo_valid = .false.
       tan_stamp = -huge(1)
    end if
  end subroutine eos_memo_prepare

  ! Composition partials of lnE and lnPeos for cells where the EOS provides none (Skye/PC/ideal and
  ! their blends), replacing fix_d_eos_dxa_partials' secants (species with |xa - xa_start| >= 1e-4
  ! only, recomputed every iteration). Here: forward-difference tangents d/dX_j of the full blended
  ! EOS result at the current (xa, rho, T), for all species with xa >= MESA_EOS_DXA_MINX, computed at
  ! the first evaluation of each solver call and reused for its later iterations (they only enter
  ! the Jacobian, and rho, T, xa move little within a step). Species below the cutoff keep MESA's values.
  subroutine eos_dxa_tangent_partials(s, k, ierr)
    use eos_def, only: num_eos_basic_results, num_eos_d_dxa_results, i_lnE, i_lnPgas
    type (star_info), pointer :: s
    integer, intent(in) :: k
    integer, intent(out) :: ierr
    real(dp), dimension(num_eos_basic_results) :: res0, d0_dlnd, d0_dlnT, res, d_dlnd, d_dlnT
    real(dp) :: dres_dxa(num_eos_d_dxa_results, s% species), xa1(s% species), h
    real(dp), parameter :: skye_limit = 1d-4  ! eos default mass_fraction_limit_for_Skye
    integer :: j, stamp
    logical :: refresh
    ierr = 0
    if (.not. allocated(tan_stamp)) then
       ierr = -1; return
    end if
    if (k > size(tan_stamp) .or. memo_species /= s% species) then
       ierr = -1; return
    end if
    stamp = s% solver_call_number
    if (.not. tangent_lag) stamp = 1000*s% solver_call_number + s% solver_iter
    refresh = (tan_stamp(k) /= stamp)
    ! lagged mode: also refresh when the cell moved since the tangents were taken
    if (tangent_lag .and. .not. refresh) then
       if (maxval(abs(s% xa(:,k) - tan_xa(:,k))) > tangent_xtol) then
          refresh = .true.
       else if (max(abs(s% lnT(k) - tan_lnT(k)), abs(s% lnd(k) - tan_lnd(k))) > tangent_ttol) then
          refresh = .true.
       end if
    end if
    if (refresh) then
       call eos_memo_count(16)
       tan_xa(:,k) = s% xa(:,k)
       tan_lnT(k) = s% lnT(k)
       tan_lnd(k) = s% lnd(k)
       call get_eos_memo( &
          s, k, 1, s% xa(:,k), &
          s% rho(k), s% lnd(k)/ln10, s% T(k), s% lnT(k)/ln10, &
          res0, d0_dlnd, d0_dlnT, ierr)
       if (ierr /= 0) return
       do j = 1, s% species
          tan_set(j,k) = .false.
          if (s% xa(j,k) < tangent_minx) cycle
          h = tangent_h
          ! stay on one side of Skye's relevant-species cutoff (a step across it is a jump)
          if (s% xa(j,k) < skye_limit .and. s% xa(j,k) + h >= skye_limit) h = -h
          if (s% xa(j,k) + h < 0) cycle
          xa1 = s% xa(:,k)
          xa1(j) = xa1(j) + h
          call eos_memo_count(14)
          call get_eos( &
             s, k, xa1, &
             s% rho(k), s% lnd(k)/ln10, s% T(k), s% lnT(k)/ln10, &
             res, d_dlnd, d_dlnT, dres_dxa, ierr)
          if (ierr /= 0 .or. is_bad(res(i_lnE)) .or. is_bad(res(i_lnPgas))) then
             ierr = 0; cycle  ! keep MESA's value for this species
          end if
          tan_dlnE(j,k) = (res(i_lnE) - res0(i_lnE))/h
          tan_dlnPgas(j,k) = (res(i_lnPgas) - res0(i_lnPgas))/h
          tan_set(j,k) = .true.
       end do
       tan_stamp(k) = stamp
    end if
    do j = 1, s% species
       if (.not. tan_set(j,k)) cycle
       s% dlnE_dxa_for_partials(j,k) = tan_dlnE(j,k)
       s% dlnPeos_dxa_for_partials(j,k) = s% Pgas(k)*tan_dlnPgas(j,k)/s% Peos(k)
    end do
  end subroutine eos_dxa_tangent_partials

  ! diagnostic: compare the composition partials in use (s% dlnE_dxa_for_partials etc., from the normal
  ! EOS call) with central differences of the EOS at the current (xa, rho, T); pure-Skye cells only
  subroutine eos_skye_dxa_check_cell(s, k)
    use eos_def, only: num_eos_basic_results, num_eos_d_dxa_results, i_lnE, i_lnPgas
    type (star_info), pointer :: s
    integer, intent(in) :: k
    real(dp), dimension(num_eos_basic_results) :: resp, resm, d_dlnd, d_dlnT
    real(dp) :: dres_dxa(num_eos_d_dxa_results, s% species), xa1(s% species), fdE, fdP, anE, anP
    real(dp), parameter :: skye_limit = 1d-4, h = 1d-5
    integer :: j, ierr
    do j = 1, min(s% species, 64)
       if (s% xa(j,k) - h < 0) cycle
       if (s% xa(j,k) - h <= skye_limit .and. s% xa(j,k) + h > skye_limit) cycle
       xa1 = s% xa(:,k)
       xa1(j) = s% xa(j,k) + h
       call get_eos(s, k, xa1, s% rho(k), s% lnd(k)/ln10, s% T(k), s% lnT(k)/ln10, &
          resp, d_dlnd, d_dlnT, dres_dxa, ierr)
       if (ierr /= 0) cycle
       xa1(j) = s% xa(j,k) - h
       call get_eos(s, k, xa1, s% rho(k), s% lnd(k)/ln10, s% T(k), s% lnT(k)/ln10, &
          resm, d_dlnd, d_dlnT, dres_dxa, ierr)
       if (ierr /= 0) cycle
       fdE = (resp(i_lnE) - resm(i_lnE))/(2d0*h)
       fdP = (resp(i_lnPgas) - resm(i_lnPgas))/(2d0*h)
       anE = s% dlnE_dxa_for_partials(j,k)
       anP = s% dlnPeos_dxa_for_partials(j,k)*s% Peos(k)/s% Pgas(k)
       if (chk_cells == 0) write(*,'(a,2i6,f8.3,f8.3,es12.4,4es14.6)') 'SKYEDXA_CELL ', k, j, &
          s% lnT(k)/ln10, s% lnd(k)/ln10, s% xa(j,k), anE, fdE, anP, fdP
!$OMP critical (skye_dxa_chk)
       chk_err_E(j) = max(chk_err_E(j), abs(anE - fdE))
       chk_max_E(j) = max(chk_max_E(j), abs(fdE))
       chk_err_P(j) = max(chk_err_P(j), abs(anP - fdP))
       chk_max_P(j) = max(chk_max_P(j), abs(fdP))
!$OMP end critical (skye_dxa_chk)
    end do
!$OMP atomic
    chk_cells = chk_cells + 1
  end subroutine eos_skye_dxa_check_cell

  subroutine eos_memo_count(i)
    integer, intent(in) :: i
    if (.not. eos_memo_stats_on) return
!$OMP atomic
    memo_count(i) = memo_count(i) + 1
  end subroutine eos_memo_count

  logical function memo_usable(s, k, slot, xa)
    type (star_info), pointer :: s
    integer, intent(in) :: k, slot
    real(dp), intent(in) :: xa(:)
    memo_usable = memo_on .and. allocated(memo_valid) .and. k >= 1 .and. k <= memo_nz .and. &
       slot >= 1 .and. slot <= memo_slots .and. memo_id == s% id .and. &
       memo_handle == s% eos_handle .and. size(xa) == memo_species
  end function memo_usable

  subroutine memo_make_key(xa, Rho, logRho, T, logT, key)
    real(dp), intent(in) :: xa(:), Rho, logRho, T, logT
    integer(8), intent(out) :: key(:)
    key(1) = transfer(Rho, 0_8)
    key(2) = transfer(logRho, 0_8)
    key(3) = transfer(T, 0_8)
    key(4) = transfer(logT, 0_8)
    key(5:) = transfer(xa, 0_8, size(xa))
  end subroutine memo_make_key

  subroutine eos_memo_store(s, k, slot, xa, Rho, logRho, T, logT, res, dres_dlnRho, dres_dlnT)
    type (star_info), pointer :: s
    integer, intent(in) :: k, slot
    real(dp), intent(in) :: xa(:), Rho, logRho, T, logT
    real(dp), intent(in) :: res(:), dres_dlnRho(:), dres_dlnT(:)
    if (.not. memo_usable(s, k, slot, xa)) return
    call memo_make_key(xa, Rho, logRho, T, logT, memo_key(:,k,slot))
    memo_res(:,k,slot) = res(1:size(memo_res,1))
    memo_dlnd(:,k,slot) = dres_dlnRho(1:size(memo_dlnd,1))
    memo_dlnT(:,k,slot) = dres_dlnT(1:size(memo_dlnT,1))
    memo_valid(k,slot) = .true.
    call eos_memo_count(3)
  end subroutine eos_memo_store

  ! same as get_eos for res, dres_dlnRho, dres_dlnT (no composition partials), reusing an identical
  ! earlier call for cell k if there is one
  subroutine get_eos_memo( &
       s, k, slot, xa, &
       Rho, logRho, T, logT, &
       res, dres_dlnRho, dres_dlnT, ierr)
    use eos_def, only: num_eos_basic_results, num_eos_d_dxa_results
    type (star_info), pointer :: s
    integer, intent(in) :: k, slot
    real(dp), intent(in) :: xa(:), Rho, logRho, T, logT
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT
    integer, intent(out) :: ierr
    real(dp) :: dres_dxa(num_eos_d_dxa_results,s% species)
    integer(8) :: key(4 + size(xa))
    logical :: usable
    usable = memo_usable(s, k, slot, xa)
    call eos_memo_count(1)
    if (usable) then
       if (memo_valid(k,slot)) then
          call memo_make_key(xa, Rho, logRho, T, logT, key)
          if (all(key == memo_key(:,k,slot))) then
             res = memo_res(:,k,slot)
             dres_dlnRho = memo_dlnd(:,k,slot)
             dres_dlnT = memo_dlnT(:,k,slot)
             ierr = 0
             call eos_memo_count(2)
             return
          end if
       end if
    end if
    call get_eos( &
         s, k, xa, &
         Rho, logRho, T, logT, &
         res, dres_dlnRho, dres_dlnT, &
         dres_dxa, ierr)
    if (ierr == 0 .and. usable) &
       call eos_memo_store(s, k, slot, xa, Rho, logRho, T, logT, res, dres_dlnRho, dres_dlnT)
  end subroutine get_eos_memo

  ! Get eos results data given density & temperature

  subroutine get_eos( &
       s, k, xa, &
       Rho, logRho, T, logT, &
       res, dres_dlnRho, dres_dlnT, &
       dres_dxa, ierr)

    use eos_lib, only: eosDT_get
    use eos_def, only: num_eos_basic_results, num_eos_d_dxa_results, num_helm_results, i_lnE

    type (star_info), pointer :: s
    integer, intent(in) :: k  ! 0 means not being called for a particular cell
    real(dp), intent(in) :: xa(:), Rho, logRho, T, logT
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT
    real(dp), intent(out) :: dres_dxa(num_eos_d_dxa_results,s% species)
    integer, intent(out) :: ierr

    integer :: j

    include 'formats'

    ierr = 0

    if (s% doing_timing) &
       s% timing_num_get_eos_calls = s% timing_num_get_eos_calls + 1

    if(logRho < -25) then
      ! Provide some hard lower limit on what we would even try to evalue the eos at
      ! Going to low causes FPE's when we try to evaluate certain derivatives that need (rho**power)
      s% retry_message = 'eos evaluated at too low a density'
      ierr = -1
      return
    end if

    call eosDT_get( &
       s% eos_handle, s% species, s% chem_id, s% net_iso, xa, &
       Rho, logRho, T, logT, &
       res, dres_dlnRho, dres_dlnT, dres_dxa, ierr)

    if (ierr /= 0) then
       s% retry_message = 'get_eos failed'
       if (s% report_ierr) then
          !$OMP critical (get_eos_critical)
          write(*,*) 'get_eos ierr', ierr
          write(*,2) 'k', k
          do j=1,s% species
             write(*,2) 'xa(j) ' // trim(s% nameofequ(j+s% nvar_hydro)), j, xa(j)
          end do
          write(*,1) 'log10Rho', logRho
          write(*,1) 'log10T', logT
          if (s% stop_for_bad_nums .and. &
               is_bad(logRho+logT)) call mesa_error(__FILE__,__LINE__,'do_eos_for_cell')
          !$OMP end critical (get_eos_critical)
       end if
       return
    end if

  end subroutine get_eos


  ! Solve for temperature & eos results data given density & energy

  subroutine solve_eos_given_DE( &
       s, k, xa, &
       logRho, logE, logT_guess, logT_tol, logE_tol, &
       logT, res, dres_dlnRho, dres_dlnT, dres_dxa, &
       ierr)

    use eos_def
    use eos_lib, only: eosDT_get_T

    type (star_info), pointer :: s
    integer, intent(in) :: k  ! 0 indicates not for a particular cell.
    real(dp), intent(in) :: &
         xa(:), logRho, logE, &
         logT_guess, logT_tol, logE_tol
    real(dp), intent(out) :: logT
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT
    real(dp), intent(out) :: dres_dxa(num_eos_d_dxa_results,s% species)
    integer, intent(out) :: ierr

    integer :: eos_calls

    include 'formats'

    ierr = 0

    call eosDT_get_T( &
       s% eos_handle, &
       s% species, s% chem_id, s% net_iso, xa, &
       logRho, i_lnE, logE*ln10, &
       logT_tol, logE_tol*ln10, MAX_ITER_FOR_SOLVE, logT_guess,  &
       arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
       logT, res, dres_dlnRho, dres_dlnT, dres_dxa, &
       eos_calls, ierr)

    if (s% doing_timing) s% timing_num_solve_eos_calls = s% timing_num_solve_eos_calls + eos_calls

  end subroutine solve_eos_given_DE


  ! Solve for temperature & eos results data given density & gas energy

  subroutine solve_eos_given_DEgas( &
       s, k, xa, &
       logRho, egas, logT_guess, logT_tol, egas_tol, &
       logT, res, dres_dlnRho, dres_dlnT, dres_dxa, &
       ierr)

    use eos_def
    use eos_lib, only: eosDT_get_T

    type (star_info), pointer :: s
    integer, intent(in) :: k  ! 0 indicates not for a particular cell.
    real(dp), intent(in) :: &
         xa(:), logRho, egas, &
         logT_guess, logT_tol, egas_tol
    real(dp), intent(out) :: logT
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT
    real(dp), intent(out) :: dres_dxa(num_eos_d_dxa_results,s% species)
    integer, intent(out) :: ierr

    integer :: eos_calls

    include 'formats'

    ierr = 0

    if (s% doing_timing) s% timing_num_solve_eos_calls = s% timing_num_solve_eos_calls + 1

    call eosDT_get_T( &
       s% eos_handle, &
       s% species, s% chem_id, s% net_iso, xa, &
       logRho, i_egas, egas, logT_tol, egas_tol, MAX_ITER_FOR_SOLVE, logT_guess, &
       arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
       logT, res, dres_dlnRho, dres_dlnT, &
       dres_dxa, eos_calls, ierr)

  end subroutine solve_eos_given_DEgas


  ! Solve for temperature & eos results data given density & pressure

  subroutine solve_eos_given_DP( &
       s, k, xa, &
       logRho, logP, logT_guess, logT_tol, logP_tol, &
       logT, res, dres_dlnRho, dres_dlnT, dres_dxa, &
       ierr)

    use eos_def
    use eos_lib, only: eosDT_get_T

    type (star_info), pointer :: s
    integer, intent(in) :: k  ! 0 indicates not for a particular cell.
    real(dp), intent(in) :: &
         xa(:), logRho, logP, &
         logT_guess, logT_tol, logP_tol
    real(dp), intent(out) :: logT
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT
    real(dp), intent(out) :: dres_dxa(num_eos_d_dxa_results,s% species)
    integer, intent(out) :: ierr

    integer :: eos_calls

    include 'formats'

    ierr = 0

    if (s% doing_timing) s% timing_num_solve_eos_calls = s% timing_num_solve_eos_calls + 1

    call eosDT_get_T( &
       s% eos_handle, &
       s% species, s% chem_id, s% net_iso, xa, &
       logRho, i_logPtot, logP, logT_tol, logP_tol, MAX_ITER_FOR_SOLVE, logT_guess, &
       arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
       logT, res, dres_dlnRho, dres_dlnT, &
       dres_dxa, eos_calls, ierr)

  end subroutine solve_eos_given_DP


  ! Solve for temperature & eos results data for a given density &
  ! entropy

  subroutine solve_eos_given_DS( &
       s, k, xa, &
       logRho, logS, logT_guess, logT_tol, logS_tol, &
       logT, res, dres_dlnRho, dres_dlnT, dres_dxa, &
       ierr)

    use eos_def
    use eos_lib, only: eosDT_get_T

    type (star_info), pointer :: s
    integer, intent(in) :: k  ! 0 indicates not for a particular cell.
    real(dp), intent(in) :: &
         xa(:), logRho, logS, &
         logT_guess, logT_tol, logS_tol
    real(dp), intent(out) :: logT
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT
    real(dp), intent(out) :: dres_dxa(num_eos_d_dxa_results,s% species)
    integer, intent(out) :: ierr

    integer :: eos_calls

    include 'formats'

    ierr = 0

    call eosDT_get_T( &
       s% eos_handle, &
       s% species, s% chem_id, s% net_iso, xa, &
       logRho, i_lnS, logS*ln10, &
       logT_tol, logS_tol*ln10, MAX_ITER_FOR_SOLVE, logT_guess,  &
       arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
       logT, res, dres_dlnRho, dres_dlnT, dres_dxa, &
       eos_calls, ierr)

    if (s% doing_timing) s% timing_num_solve_eos_calls = s% timing_num_solve_eos_calls + eos_calls

  end subroutine solve_eos_given_DS


  ! Solve for density & eos results data given pressure & temperature

  subroutine solve_eos_given_PT( &
       s, k, xa, &
       logT, logP, logRho_guess, logRho_tol, logP_tol, &
       logRho, res, dres_dlnRho, dres_dlnT, dres_dxa, &
       ierr)

    use eos_def
    use eos_lib, only: eosDT_get_Rho

    type (star_info), pointer :: s
    integer, intent(in) :: k  ! 0 indicates not for a particular cell.
    real(dp), intent(in) :: &
         xa(:), logT, logP, &
         logRho_guess, logRho_tol, logP_tol
    real(dp), intent(out) :: logRho
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT
    real(dp), intent(out) :: dres_dxa(num_eos_d_dxa_results,s% species)
    integer, intent(out) :: ierr

    integer :: eos_calls

    include 'formats'

    ierr = 0

    if (s% doing_timing) s% timing_num_solve_eos_calls = s% timing_num_solve_eos_calls + 1

    call eosDT_get_Rho( &
       s% eos_handle, &
       s% species, s% chem_id, s% net_iso, xa, &
       logT, i_logPtot, logP, logRho_tol, logP_tol, MAX_ITER_FOR_SOLVE, logRho_guess, &
       arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
       logRho, res, dres_dlnRho, dres_dlnT, &
       dres_dxa, eos_calls, ierr)

  end subroutine solve_eos_given_PT


  ! Solve for density & eos results data given gas pressure &
  ! temperature

  subroutine solve_eos_given_PgasT( &
       s, k, xa, &
       logT, logPgas, logRho_guess, logRho_tol, logPgas_tol, &
       logRho, res, dres_dlnRho, dres_dlnT, dres_dxa, &
       ierr)

    use eos_def
    use eos_lib, only: eosDT_get_Rho

    type (star_info), pointer :: s
    integer, intent(in) :: k  ! 0 indicates not for a particular cell.
    real(dp), intent(in) :: &
         xa(:), logT, logPgas, &
         logRho_guess, logRho_tol, logPgas_tol
    real(dp), intent(out) :: logRho
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT
    real(dp), intent(out) :: dres_dxa(num_eos_d_dxa_results,s% species)
    integer, intent(out) :: ierr

    integer :: eos_calls

    include 'formats'

    ierr = 0

    call eosDT_get_Rho( &
       s% eos_handle, &
       s% species, s% chem_id, s% net_iso, xa, &
       logT, i_lnPgas, logPgas*ln10, &
       logRho_tol, logPgas_tol*ln10, MAX_ITER_FOR_SOLVE, logRho_guess, &
       arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
       logRho, res, dres_dlnRho, dres_dlnT, &
       dres_dxa, eos_calls, ierr)
    if (ierr /= 0 .and. s% report_ierr) then
       write(*,*) 'Call to eosDT_get_Rho failed in solve_eos_given_PgasT'
       write(*,2) 'logPgas', k, logPgas
       write(*,2) 'logT', k, logT
       write(*,2) 'logRho_guess', k, logRho_guess
    end if

    if (s% doing_timing) s% timing_num_solve_eos_calls = s% timing_num_solve_eos_calls + eos_calls

  end subroutine solve_eos_given_PgasT


  ! Solve for density & eos results data given gas pressure &
  ! temperature, with logRho_guess calculated automatically via an
  ! initial call to eos_gamma_PT_get

  subroutine solve_eos_given_PgasT_auto( &
       s, k, xa, &
       logT, logPgas, logRho_tol, logPgas_tol, &
       logRho, res, dres_dlnRho, dres_dlnT, dres_dxa, &
       ierr)

    use chem_lib, only: basic_composition_info
    use eos_def
    use eos_lib, only: eos_gamma_PT_get

    type (star_info), pointer :: s
    integer, intent(in) :: k  ! 0 indicates not for a particular cell.
    real(dp), intent(in) :: &
         xa(:), logT, logPgas, &
         logRho_tol, logPgas_tol
    real(dp), intent(out) :: logRho
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT
    real(dp), intent(out) :: dres_dxa(num_eos_d_dxa_results,s% species)
    integer, intent(out) :: ierr

    real(dp) :: rho_guess, logRho_guess, gamma

    ! compute composition info
    real(dp) :: Y, Z, X, abar, zbar, z2bar, z53bar, ye, mass_correction, sumx

    call basic_composition_info( &
       s% species, s% chem_id, xa, X, Y, Z, &
       abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)

    gamma = 5d0/3d0
    call eos_gamma_PT_get( &
       s% eos_handle, abar, exp10(logPgas), logPgas, exp10(logT), logT, gamma, &
       rho_guess, logRho_guess, res, dres_dlnRho, dres_dlnT, &
       ierr)
    if (ierr /= 0) then
       ierr = 0
       logRho_guess = arg_not_provided
    end if

    call solve_eos_given_PgasT( &
       s, k, xa, &
       logT, logPgas, logRho_guess, logRho_tol, logPgas_tol, &
       logRho, res, dres_dlnRho, dres_dlnT, dres_dxa, &
       ierr)

  end subroutine solve_eos_given_PgasT_auto

end module eos_support
