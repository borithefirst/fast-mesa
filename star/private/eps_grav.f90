! ***********************************************************************
!
!   Copyright (C) 2010-2021  The MESA Team
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

      module eps_grav

      use star_private_def
      use const_def, only: dp, ln10
      use chem_def, only: chem_isos
      use utils_lib, only: mesa_error, is_bad
      use auto_diff

      implicit none

      ! savethesun epsg lin: check-mode accumulators
      logical, save :: epsg_lin_check = .false.
      integer, save :: epsg_chk_model = -1, epsg_chk_n = 0
      real(dp), save :: epsg_chk_err = 0d0, epsg_chk_mag = 0d0

      private
      public :: eval_eps_grav_and_partials, zero_eps_grav_and_partials
      public :: epsg_lin_report  ! savethesun epsg lin

      contains

      subroutine eval_eps_grav_and_partials(s, k, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: eps_grav
         include 'formats'
         ierr = 0

         if (s% dt <= 0) then
            ierr = -1
            return
         end if

         call zero_eps_grav_and_partials(s, k)

         ! zero composition derivatives
         ! if include_composition_in_eps_grav is true
         ! then these are set in the call to eval_eps_grav_composition
         s% d_eps_grav_dx(:,k) = 0

         call eval1_eps_grav_and_partials(s, k, eps_grav, ierr)
         if (ierr /= 0) return

         s% eps_grav_ad(k) = eps_grav

         if (s% use_other_eps_grav) then
            ! note: call this after 1st doing the standard calculation
            call s% other_eps_grav(s% id, k, s% dt, ierr)
            if (ierr /= 0) return
         end if

         ! apply user-specified scaling factor after hook
         s% eps_grav_ad(k) = s% eps_grav_factor * s% eps_grav_ad(k)

      end subroutine eval_eps_grav_and_partials


      subroutine eval1_eps_grav_and_partials(s, k, eps_grav, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: eps_grav
         integer, intent(out) :: ierr

         type(auto_diff_real_star_order1) :: eps_grav_lnS, eps_grav_std
         real(dp) :: alfa, Gamma
         logical :: using_PC

         include 'formats'
         ierr = 0

         using_PC = (s% eos_frac_PC(k) > 0)

         if (using_PC .and. s% gam_start(k) >= s% Gamma_lnS_eps_grav_full_on) then
            call do_lnS_eps_grav(s, k, eps_grav, ierr)
         else if (using_PC .and. s% gam_start(k) > s% Gamma_lnS_eps_grav_full_off) then
            Gamma = s% gam_start(k)
            alfa = (Gamma - s% Gamma_lnS_eps_grav_full_off) / &
               (s% Gamma_lnS_eps_grav_full_on - s% Gamma_lnS_eps_grav_full_off)
            call do_lnS_eps_grav(s, k, eps_grav_lnS, ierr)
            if (ierr /= 0) return
            call do_std_eps_grav(s, k, eps_grav_std, ierr)
            if (ierr /= 0) return
            ! the derivative of the blending function is missing
            ! but historically we've been able to get away with that
            ! because the two forms should match in the blend region
            eps_grav = alfa * eps_grav_lnS + (1d0 - alfa) * eps_grav_std
         else
            call do_std_eps_grav(s, k, eps_grav, ierr)
         end if


         if (ierr /= 0 .or. is_bad(eps_grav% val)) then
            ierr = -1
            s% retry_message = 'failed in eval_eps_grav_and_partials'
            if (s% report_ierr) then
               write(*,2) &
                  'failed in eval_eps_grav_and_partials', k, eps_grav% val
            end if
            if (s% stop_for_bad_nums) then
               call mesa_error(__FILE__,__LINE__,'eval1_eps_grav_and_partials')
            end if
            return
         end if

      end subroutine eval1_eps_grav_and_partials


      ! this uses the given args to calculate -T*ds/dt
      subroutine do_std_eps_grav(s, k, eps_grav, ierr)
         use auto_diff_support
         use eos_def, only: i_latent_ddlnT, i_latent_ddlnRho
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: eps_grav
         integer, intent(out) :: ierr

         type(auto_diff_real_star_order1) :: T, cp, grada, chiT, chiRho, dlnd_dt, dlnT_dt, &
            latent_ddlnT, latent_ddlnRho, eps_grav_start, eps_grav_composition_term

         logical :: test_partials

         real(dp) :: theta

         include 'formats'
         ierr = 0

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         ! select time-centering
         if (s% use_time_centered_eps_grav .and. .not. s% doing_relax) then
            theta = 0.5_dp
         else
            theta = 1.0_dp
         end if

         ! total time derivatives of (Rho, T) basis
         dlnd_dt = wrap_dxh_lnd(s,k)/s% dt
         dlnT_dt = wrap_dxh_lnT(s,k)/s% dt

         ! thermo quantities
         T = wrap_T_00(s, k)
         Cp = wrap_Cp_00(s, k)
         grada = wrap_grad_ad_00(s, k)
         chiT = wrap_chiT_00(s, k)
         chiRho = wrap_chiRho_00(s, k)

         ! MESA I, Equation (12)
         eps_grav = -T*Cp * ((1d0 - grada*chiT)*dlnT_dt - grada*chiRho*dlnd_dt)

         ! phase transition latent heat
         ! Jermyn et al. (2021), Equation (47)
         latent_ddlnRho = wrap_latent_ddlnRho_00(s,k)
         latent_ddlnT = wrap_latent_ddlnT_00(s,k)
         eps_grav = eps_grav - (dlnd_dt * latent_ddlnRho + dlnT_dt * latent_ddlnT)


         ! for time centered version
         if (s% use_time_centered_eps_grav .and. .not. s% doing_relax) then

            ! start values are constants during Newton iters
            eps_grav_start = -s% T_start(k)*s% cp_start(k) &
                            * ((1d0 - s% grada_start(k)*s% chiT_start(k))*dlnT_dt - s% grada_start(k)*s% chiRho_start(k)*dlnd_dt)

            ! phase transition latent heat
            eps_grav_start = eps_grav_start - (dlnd_dt * s% latent_ddlnRho_start(k) + dlnT_dt * s% latent_ddlnT_start(k))

            eps_grav = theta * eps_grav + (1d0-theta) * eps_grav_start

         end if


         if (s% include_composition_in_eps_grav) then
            call eval_eps_grav_composition(s, k, eps_grav_composition_term, ierr)
            if (ierr /= 0) return
            eps_grav = eps_grav + eps_grav_composition_term
         end if

         if (is_bad(eps_grav% val)) then
            ierr = -1
            s% retry_message = 'do_lnd_eps_grav -- bad value for eps_grav'
            if (s% report_ierr) &
               write(*,2) 'do_lnd_eps_grav -- bad value for eps_grav', k, eps_grav% val
            if (s% stop_for_bad_nums) call mesa_error(__FILE__,__LINE__,'do_lnd_eps_grav')
            return
         end if

         if (test_partials) then
            s% solver_test_partials_val = 0
            s% solver_test_partials_var = 0
            s% solver_test_partials_dval_dx = 0
            write(*,*) 'do_std_eps_grav chiT', s% solver_test_partials_var
         end if

      end subroutine do_std_eps_grav


      subroutine do_lnS_eps_grav(s, k, eps_grav, ierr)
         use auto_diff_support
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: eps_grav
         integer, intent(out) :: ierr

         real(dp) :: entropy_start
         type(auto_diff_real_star_order1) :: entropy, T

         include 'formats'
         ierr = 0

         T = wrap_T_00(s,k)
         entropy = wrap_s_00(s,k)
         entropy_start = exp(s% lnS_start(k))

         eps_grav = -T*(entropy - entropy_start)/s% dt

         ! NOTE: the correct version of the composition term to go with TdS is
         !     -sum_i (\partial e/\partial Y_i)_{s,\rho} dY_i
         ! see for example, MESA IV, equation (59)
         !
         ! When on PC near crystallization (the only place the lnS form is used)
         ! there are typically no composition changes, so we drop this term
         !
         ! The term we normally use comes from expanding in the (rho,T,X) basis
         !     -sum_ (\partial e/\partial X_i)_{T,\rho} dX_i
         ! and in practice can be a reasonable approximation to the above.
         !
         ! If an such approximation is desired, one could use the following code:

         ! if (s% include_composition_in_eps_grav) then
         !    call eval_eps_grav_composition(s, k, eps_grav_composition_term, ierr)
         !    if (ierr /= 0) return
         !    eps_grav = eps_grav + eps_grav_composition_term
         ! end if

         if (is_bad(eps_grav% val)) then
            ierr = -1
            s% retry_message = 'do_lnS_eps_grav -- bad value for eps_grav'
            if (s% report_ierr) &
               write(*,2) 'do_lnS_eps_grav -- bad value for eps_grav', k, eps_grav% val
            if (s% stop_for_bad_nums) call mesa_error(__FILE__,__LINE__,'do_lnS_eps_grav')
            return
         end if

      end subroutine do_lnS_eps_grav


      real(dp) function epsg_lin_threshold()  ! savethesun epsg lin: env MESA_EPSG_LIN=<max |dX|>, 0 = off
         logical, save :: checked = .false.
         real(dp), save :: thr = 0d0
         character(len=64) :: val
         integer :: stat
         if (.not. checked) then
            !$OMP CRITICAL (savethesun_epsg_lin)
            if (.not. checked) then
               call get_environment_variable('MESA_EPSG_LIN', val, status=stat)
               if (stat == 0 .and. len_trim(val) > 0) then
                  read(val, *, iostat=stat) thr
                  if (stat /= 0 .or. thr < 0d0) thr = 0d0
               end if
               call get_environment_variable('MESA_EPSG_LIN_CHECK', val, status=stat)
               epsg_lin_check = (stat == 0 .and. len_trim(val) > 0 .and. trim(val) /= '0')
               if (thr > 0d0) write(*,'(a,es10.2,a,l2)') &
                  'savethesun: linearized eps_grav composition term for 0 < |dX| <', thr, '  check', epsg_lin_check
               checked = .true.
            end if
            !$OMP END CRITICAL (savethesun_epsg_lin)
         end if
         epsg_lin_threshold = thr
      end function epsg_lin_threshold


      subroutine epsg_lin_report(s)  ! savethesun epsg lin: check mode, once per model (called from serial code)
         type (star_info), pointer :: s
         if (.not. epsg_lin_check) return
         if (s% model_number == epsg_chk_model) return
         if (epsg_chk_model >= 0) write(*,'(a,i8,i10,2es12.3)') 'EPSGLIN ', epsg_chk_model, &
            epsg_chk_n, epsg_chk_err, epsg_chk_mag
         epsg_chk_model = s% model_number
         epsg_chk_n = 0; epsg_chk_err = 0d0; epsg_chk_mag = 0d0
      end subroutine epsg_lin_report


      subroutine eval_eps_grav_composition(s, k, eps_grav_composition_term, ierr)
         use auto_diff_support, only: wrap
         use eos_support, only: get_eos, get_eos_memo, eos_memo_count, eos_memo_stats_on  ! savethesun eos memo
         use eos_support, only: eos_skye_dxa_on  ! savethesun epsg lin
         use eos_def, only: num_eos_basic_results, num_eos_d_dxa_results, i_lnE, i_lnPgas

         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: eps_grav_composition_term
         real(dp) :: memo_mx  ! savethesun eos memo
         real(dp) :: lin_thr, lin_mx, lin_S, lin_de, lin_dd, lin_dT  ! savethesun epsg lin
         logical :: lin_ok
         integer, intent(out) :: ierr
         real(dp) :: &
            e, e_start, de, d_de_dlnd, d_de_dlnT, &
            e_with_xa_start, d_e_with_xa_start_dlnd, d_e_with_xa_start_dlnT, &
            e_with_DT_start, Pgas_with_DT_start
         real(dp), dimension(num_eos_basic_results) :: &
            res, dres_dlnd, dres_dlnT
         real(dp) :: dres_dxa(num_eos_d_dxa_results,s% species)
         integer :: j
         logical :: test_partials

         real(dp) :: theta

         include 'formats'
         ierr = 0

         eps_grav_composition_term = 0

         if (s% use_time_centered_eps_grav .and. .not. s% doing_relax) then
            theta = 0.5_dp
         else
            theta = 1.0_dp
         end if

         ! some EOSes have composition partials and some do not
         ! those currently without dx partials are PC & Skye
         ! however, composition partials of lnE & lnP were revised by the call to
         ! fix_d_eos_dxa_partials at the beginning of eval_equ_for_solver

         ! directly fill array with composition derivatives
         ! we only have lnE and lnP derivs, so use
         ! eps_grav = -de/dt + (P/rho) * dlnd/dt
         do j=1, s% species
            s% d_eps_grav_dx(j,k) = -s% energy(k) * s% dlnE_dxa_for_partials(j,k)/s% dt + &
               (s% Peos(k) / s% Rho(k)) * s% dlnPeos_dxa_for_partials(j,k) * s% dxh_lnd(k)/s% dt
         end do

         ! savethesun epsg lin: first-order expansion for tiny composition changes
         lin_ok = .false.
         lin_thr = epsg_lin_threshold()
         if (lin_thr > 0d0) then
            if (s% eos_frac_PC(k) + s% eos_frac_ideal(k) == 0d0 .and. &
                (s% eos_frac_Skye(k) == 0d0 .or. eos_skye_dxa_on)) then
               lin_mx = maxval(abs(s% xa(:,k) - s% xa_start(:,k)))
               lin_ok = (lin_mx > 0d0 .and. lin_mx < lin_thr)
            end if
         end if
         if (lin_ok) then
            lin_S = sum(s% dlnE_dxa_for_partials(:,k)*(s% xa(:,k) - s% xa_start(:,k)))
            if (s% use_time_centered_eps_grav .and. .not. s% doing_relax) then
               lin_de = theta*s% energy(k)*lin_S + (1d0 - theta)*s% energy_start(k)*lin_S
            else
               lin_de = s% energy(k)*lin_S
            end if
            lin_dd = theta*lin_S*s% dE_dRho_for_partials(k)*s% Rho(k)
            lin_dT = theta*lin_S*s% Cv_for_partials(k)*s% T(k)
            if (.not. epsg_lin_check) then
               de = lin_de; d_de_dlnd = lin_dd; d_de_dlnT = lin_dT
               goto 777
            end if
         end if

         e = s% energy(k)
         call get_eos_memo( &  ! savethesun eos memo
            s, k, 1, s% xa_start(:,k), &
            s% rho(k), s% lnd(k)/ln10, s% T(k), s% lnT(k)/ln10, &
            res, dres_dlnd, dres_dlnT, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,2) 'failed in get_eos with xa_start', k
            return
         end if
         if (eos_memo_stats_on) then
            call eos_memo_count(6)
            if (s% eos_frac_PC(k) + s% eos_frac_Skye(k) + s% eos_frac_ideal(k) > 0) call eos_memo_count(15)
            memo_mx = maxval(abs(s% xa(:,k) - s% xa_start(:,k)))
            if (memo_mx == 0) then
               call eos_memo_count(7)
               if (exp(res(i_lnE)) == e) call eos_memo_count(8)
            else if (memo_mx < 1d-12) then
               call eos_memo_count(10)
            else if (memo_mx < 1d-8) then
               call eos_memo_count(11)
            else if (memo_mx < 1d-4) then
               call eos_memo_count(12)
            else
               call eos_memo_count(13)
            end if
         end if

         e_with_xa_start = exp(res(i_lnE))
         d_e_with_xa_start_dlnd = dres_dlnd(i_lnE)*e_with_xa_start
         d_e_with_xa_start_dlnT = dres_dlnT(i_lnE)*e_with_xa_start

         de = (e - e_with_xa_start)
         d_de_dlnd = (s% dE_dRho_for_partials(k)*s% Rho(k) - d_e_with_xa_start_dlnd)
         d_de_dlnT = (s% Cv_for_partials(k)*s% T(k) - d_e_with_xa_start_dlnT)

         if (s% use_time_centered_eps_grav .and. .not. s% doing_relax) then

            e_start = s% energy_start(k)

            call get_eos_memo( &  ! savethesun eos memo
               s, k, 2, s% xa(:,k), &
               s% rho_start(k), s% lnd_start(k)/ln10, s% T_start(k), s% lnT_start(k)/ln10, &
               res, dres_dlnd, dres_dlnT, ierr)
            if (ierr /= 0) then
               if (s% report_ierr) write(*,2) 'failed in get_eos with xa_start', k
               return
            end if

            e_with_DT_start = exp(res(i_lnE))
            if (eos_memo_stats_on) then
               if (memo_mx == 0 .and. e_with_DT_start == e_start) call eos_memo_count(9)
            end if
            de = theta * de + (1d0 - theta) * (e_with_DT_start - e_start)
            d_de_dlnd = theta * d_de_dlnd
            d_de_dlnT = theta * d_de_dlnT

            Pgas_with_DT_start = exp(res(i_lnPgas))

            ! combine with end-of-step composition derivatives

            ! until EOSes always provide composition derivatives, ignore this
            ! the end-of-step term should generally be similar enough to be OK

            ! do j=1, s% species
            !    s% d_eps_grav_dx(j,k) = theta * s% d_eps_grav_dx(j,k) + (1d0 - theta) * &
            !       (-e_with_DT_start*d_dxa(i_lnE,j)/s% dt + &
            !       Pgas_with_DT_start*d_dxa(i_lnPgas,j)/s% rho_start(k)*s% dxh_lnd(k)/s% dt)
            ! end do

         end if

         if (lin_ok .and. epsg_lin_check) then  ! savethesun epsg lin: compare with the stock value (which is kept)
            !$OMP CRITICAL (savethesun_epsg_chk)
            epsg_chk_n = epsg_chk_n + 1
            epsg_chk_err = max(epsg_chk_err, abs(lin_de - de)/s% energy(k))
            epsg_chk_mag = max(epsg_chk_mag, abs(de)/s% energy(k))
            !$OMP END CRITICAL (savethesun_epsg_chk)
         end if
777      continue  ! savethesun epsg lin
         call wrap(eps_grav_composition_term, -de/s% dt, &
            0d0, -d_de_dlnd/s% dt, 0d0, &
            0d0, -d_de_dlnT/s% dt, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0)

         ! add easy access to this quantity in star
         s% eps_grav_composition_term(k) = eps_grav_composition_term% val

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         if (test_partials) then
            s% solver_test_partials_val = de
            s% solver_test_partials_var = s% i_lnT
            s% solver_test_partials_dval_dx = d_de_dlnT
            write(*,*) 'get_dedt', s% solver_test_partials_var
         end if

         if (is_bad(eps_grav_composition_term% val)) then
            if (s% report_ierr) then
               write(*, *) s% retry_message
               write(*,2) 'eps_grav_composition_term', k, eps_grav_composition_term% val
            end if
            if (s% stop_for_bad_nums) then
               write(*,2) 'include_composition_in_eps_grav -- bad value for eps_grav_composition_term', &
                          k, eps_grav_composition_term% val
               call mesa_error(__FILE__,__LINE__,'eval_eps_grav_composition')
            end if
            return
         end if

      end subroutine eval_eps_grav_composition


      subroutine zero_eps_grav_and_partials(s, k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         s% eps_grav_ad(k) = 0
      end subroutine zero_eps_grav_and_partials

      end module eps_grav
