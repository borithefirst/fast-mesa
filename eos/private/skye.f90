! ***********************************************************************
!
!   Copyright (C) 2022  The MESA Team
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

module skye
      use const_def, only: dp, crad, kerg, mp
      use math_lib
      use auto_diff
      use eos_def

      implicit none

      logical, parameter :: dbg = .false.


      private
      public :: Get_Skye_EOS_Results, Get_Skye_alfa, Get_Skye_alfa_simple, get_Skye_for_eosdt

      ! savethesun skye dxa: analytic composition partials, on by default (MESA_SKYE_DXA=0 turns off)
      logical, save :: skye_dxa_checked = .false., skye_dxa_on = .false.

      contains

      logical function skye_dxa_enabled()
         character(len=16) :: v
         integer :: st
         if (.not. skye_dxa_checked) then
!$OMP critical (skye_dxa_init)
            if (.not. skye_dxa_checked) then
               call get_environment_variable('MESA_SKYE_DXA', v, status=st)
               skye_dxa_on = .not. (st == 0 .and. trim(v) == '0')
!$OMP flush
               skye_dxa_checked = .true.
            end if
!$OMP end critical (skye_dxa_init)
         end if
         skye_dxa_enabled = skye_dxa_on
      end function skye_dxa_enabled

      subroutine Get_Skye_alfa( &
            rq, logRho, logT, Z, abar, zbar, &
            alfa, d_alfa_dlogT, d_alfa_dlogRho, &
            ierr)
         use const_def, only: dp
         use eos_blend
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: logRho, logT, Z, abar, zbar
         real(dp), intent(out) :: alfa, d_alfa_dlogT, d_alfa_dlogRho
         integer, intent(out) :: ierr

         logical :: contained
         type(auto_diff_real_2var_order1) :: p(2), blend, dist

         ! Blend parameters
         real(dp) :: skye_blend_width
         integer, parameter :: num_points = 8
         real(dp) :: bounds(8,2)
         type (Helm_Table), pointer :: ht

         ierr = 0
         ht => eos_ht
         skye_blend_width = 0.1d0

         ! Avoid catastrophic loss of precision in HELM tables
         bounds(1,1) = ht% logdlo
         bounds(1,2) = 8.3d0

         ! Rough ionization temperature from Jermyn+2021 Equation 52 (treating denominator as ~1).
         ! We put a lower bound of logT=7.3 to ensure that solar models never use Skye.
         ! This is because the blend even in regions that are 99+% ionized produces noticeable
         ! kinks in the sound speed profile on a scale testable by the observations.
         bounds(2,1) = ht% logdlo
         bounds(2,2) = max(7.3d0,log10(1d5 * pow2(zbar))) + skye_blend_width

         ! Rough ionization density from Jermyn+2021 Equation 53, dividing by 3 so we get closer to Dragons.
         ! Don't let the density get below rho = 200 g/cc so that solar model stays away from blend into Skye.
         bounds(3,1) = max(2.3d0,log10(abar * pow3(zbar))) + skye_blend_width
         bounds(3,2) = max(7.3d0,log10(1d5 * pow2(zbar))) + skye_blend_width

         ! HELM low-T bound
         bounds(4,1) = max(2.3d0,log10(abar * pow3(zbar))) + skye_blend_width
         bounds(4,2) = ht% logtlo

         ! Lower-right of (rho,T) plane
         bounds(5,1) = ht% logdhi
         bounds(5,2) = ht% logtlo

         ! Upper-right of (rho,T) plane
         bounds(6,1) = ht% logdhi
         bounds(6,2) = ht% logthi

         ! Avoid catastrophic loss of precision in HELM tables
         bounds(7,1) = 3d0 * ht% logthi + log10(abar * mp * crad / (3d0 * kerg * (zbar + 1d0))) - 6d0
         bounds(7,2) =  ht% logthi

         ! Avoid catastrophic loss of precision in HELM tables
         bounds(8,1) = 3d0 * 8.3d0 + log10(abar * mp * crad / (3d0 * kerg * (zbar + 1d0))) - 6d0
         bounds(8,2) = 8.3d0

         ! Set up auto_diff point
         p(1) = logRho
         p(1)%d1val1 = 1d0
         p(2) = logT
         p(2)%d1val2 = 1d0

         contained = is_contained(num_points, bounds, p)
         dist = min_distance_to_polygon(num_points, bounds, p)

         if (contained) then  ! Make distance negative for points inside the polygon
            dist = -dist
         end if

         dist = dist / skye_blend_width
         blend = max(dist, 0d0)
         blend = min(blend, 1d0)

         alfa = blend%val
         d_alfa_dlogRho = blend%d1val1
         d_alfa_dlogT = blend%d1val2

      end subroutine Get_Skye_alfa


      subroutine Get_Skye_alfa_simple( &
            rq, logRho, logT, Z, abar, zbar, &
            alfa, d_alfa_dlogT, d_alfa_dlogRho, &
            ierr)
         use eos_blend
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: logRho, logT, Z, abar, zbar
         real(dp), intent(out) :: alfa, d_alfa_dlogT, d_alfa_dlogRho
         integer, intent(out) :: ierr

         type(auto_diff_real_2var_order1) :: logT_auto, logRho_auto
         type(auto_diff_real_2var_order1) :: blend, blend_logT, blend_logRho

         include 'formats'

         ierr = 0

         ! logRho is val1
         logRho_auto% val = logRho
         logRho_auto% d1val1 = 1d0
         logRho_auto% d1val2 = 0d0

         ! logT is val2
         logT_auto% val = logT
         logT_auto% d1val1 = 0d0
         logT_auto% d1val2 = 1d0

         ! logT blend
         if (logT_auto < rq% logT_min_for_any_Skye) then
            blend_logT = 0d0
         else if (logT_auto <= rq% logT_min_for_all_Skye) then
            blend_logT = (logT_auto - rQ% logT_min_for_any_Skye) / (rq% logT_min_for_all_Skye - rq% logT_min_for_any_Skye)
         else if (logT_auto > rq% logT_min_for_all_Skye) then
            blend_logT = 1d0
         end if


         ! logRho blend
         if (logRho_auto < rq% logRho_min_for_any_Skye) then
            blend_logRho = 0d0
         else if (logRho_auto <= rq% logRho_min_for_all_Skye) then
            blend_logRho = (logRho_auto - rQ% logRho_min_for_any_Skye) / (rq% logRho_min_for_all_Skye - rq% logRho_min_for_any_Skye)
         else if (logRho_auto > rq% logRho_min_for_all_Skye) then
            blend_logRho = 1d0
         end if

         ! combine blends
         blend = (1d0 - blend_logRho) * (1d0 - blend_logT)

         alfa = blend% val
         d_alfa_dlogRho = blend% d1val1
         d_alfa_dlogT = blend% d1val2

      end subroutine get_Skye_alfa_simple


      subroutine get_Skye_for_eosdt(handle, dbg, Z, X, abar, zbar, species, chem_id, net_iso, xa, &
                                    rho, logRho, T, logT, remaining_fraction, res, d_dlnd, d_dlnT, d_dxa, skip, ierr)
         integer, intent(in) :: handle
         logical, intent(in) :: dbg
         real(dp), intent(in) :: &
            Z, X, abar, zbar, remaining_fraction
         integer, intent(in) :: species
         integer, pointer :: chem_id(:), net_iso(:)
         real(dp), intent(in) :: xa(:)
         real(dp), intent(in) :: rho, logRho, T, logT
         real(dp), intent(inout), dimension(nv) :: res, d_dlnd, d_dlnT
         real(dp), intent(inout), dimension(nv, species) :: d_dxa
         logical, intent(out) :: skip
         integer, intent(out) :: ierr
         type (EoS_General_Info), pointer :: rq

         rq => eos_handles(handle)

         call Get_Skye_EOS_Results(rq, Z, X, abar, zbar, rho, logRho, T, logT, species, chem_id, xa, &
                                    res, d_dlnd, d_dlnT, d_dxa, ierr)
         skip = .false.

         ! zero all components
         res(i_frac:i_frac+num_eos_frac_results-1) = 0.0d0
         d_dlnd(i_frac:i_frac+num_eos_frac_results-1) = 0.0d0
         d_dlnT(i_frac:i_frac+num_eos_frac_results-1) = 0.0d0

         ! mark this one
         res(i_frac_Skye) = 1.0d0

      end subroutine get_Skye_for_eosdt

      subroutine Get_Skye_EOS_Results( &
               rq, Z, X, abar, zbar, Rho, logRho, T, logT, &
               species, chem_id, xa, res, d_dlnd, d_dlnT, d_dxa, ierr)
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: Z, X, abar, zbar
         real(dp), intent(in) :: Rho, logRho, T, logT
         integer, intent(in) :: species
         integer, pointer :: chem_id(:)
         real(dp), intent(in) :: xa(:)
         integer, intent(out) :: ierr
         real(dp), intent(out), dimension(nv) :: res, d_dlnd, d_dlnT
         real(dp), intent(out), dimension(nv, species) :: d_dxa

         real(dp) :: logT_ion, logT_neutral
         logical :: do_dxa  ! savethesun skye dxa

         include 'formats'

         ierr = 0

         do_dxa = eos_want_skye_dxa
         if (do_dxa) do_dxa = skye_dxa_enabled()
         call skye_eos( &
            T, Rho, X, abar, zbar, &
            rq%Skye_min_gamma_for_solid, rq%Skye_max_gamma_for_liquid, &
            rq%Skye_solid_mixing_rule, rq%mass_fraction_limit_for_Skye, &
            rq%Skye_use_ion_offsets, &
            species, chem_id, xa, &
            res, d_dlnd, d_dlnT, d_dxa, ierr, do_dxa)

         ! composition derivatives not provided  ! savethesun skye dxa: unless MESA_SKYE_DXA=1
         if (.not. do_dxa) d_dxa = 0

         if (ierr /= 0) then
            if (dbg) then
               write(*,*) 'failed in Get_Skye_EOS_Results'
               write(*,1) 'T', T
               write(*,1) 'logT', logT
               write(*,1) 'Rho', Rho
               write(*,1) 'logRho', logRho
               write(*,1) 'abar', abar
               write(*,1) 'zbar', zbar
               write(*,1) 'X', X
               call mesa_error(__FILE__,__LINE__,'Get_Skye_EOS_Results')
            end if
            return
         end if

      end subroutine Get_Skye_EOS_Results


      !>..given a temperature temp [K], density den [g/cm**3], and a composition
      !!..this routine returns most of the other
      !!..thermodynamic quantities. of prime interest is the pressure [erg/cm**3],
      !!..specific thermal energy [erg/gr], the entropy [erg/g/K], along with
      !!..their derivatives with respect to temperature, density, abar, and zbar.
      !!..other quantities such the normalized chemical potential eta (plus its
      !!..derivatives), number density of electrons and positron pair (along
      !!..with their derivatives), adiabatic indices, specific heats, and
      !!..relativistically correct sound speed are also returned.
      !!..
      !!..this routine assumes planckian photons, an ideal gas of ions,
      !!..and an electron-positron gas with an arbitrary degree of relativity
      !!..and degeneracy. interpolation in a table of the helmholtz free energy
      !!..is used to return the electron-positron thermodynamic quantities.
      !!..all other derivatives are analytic.
      !!..
      !!..references: cox & giuli chapter 24 ; timmes & swesty apj 1999

      !!..this routine assumes a call to subroutine read_helm_table has
      !!..been performed prior to calling this routine.
      subroutine skye_eos( &
            temp_in, den_in, Xfrac, abar, zbar,  &
            Skye_min_gamma_for_solid, Skye_max_gamma_for_liquid, &
            Skye_solid_mixing_rule, &
            mass_fraction_limit, use_ion_offsets, &
            species, chem_id, xa, &
            res, d_dlnd, d_dlnT, d_dxa, ierr, do_dxa)

         use eos_def
         use chem_lib, only: composition_info  ! savethesun skye dxa
         use const_def, only: amu
         use utils_lib, only: is_bad
         use chem_def, only: chem_isos
         use ion_offset, only: compute_ion_offset
         use skye_ideal
         use skye_coulomb
         use skye_thermodynamics
         use auto_diff

         integer :: j
         integer, intent(in) :: species
         integer, pointer :: chem_id(:)
         real(dp), intent(in) :: xa(:)
         real(dp), intent(in) :: temp_in, den_in, mass_fraction_limit, Skye_min_gamma_for_solid, Skye_max_gamma_for_liquid
         real(dp), intent(in) :: Xfrac, abar, zbar
         logical, intent(in) :: use_ion_offsets
         character(len=128), intent(in) :: Skye_solid_mixing_rule
         integer, intent(out) :: ierr
         real(dp), intent(out), dimension(nv) :: res, d_dlnd, d_dlnT
         real(dp), intent(out), dimension(nv, species) :: d_dxa

         integer :: relevant_species, lookup(species)
         type(auto_diff_real_2var_order3) :: temp, logtemp, den, logden, din
         real(dp) :: AZION(species), ACMI(species), A(species), select_xa(species), ya(species)
         type (Helm_Table), pointer :: ht
         real(dp) :: ytot1, ye, norm
         type(auto_diff_real_2var_order3) :: etaele, xnefer, phase, latent_ddlnT, latent_ddlnRho
         type(auto_diff_real_2var_order3) :: F_ion_gas, F_rad, F_ideal_ion, F_coul
         type(auto_diff_real_2var_order3) :: F_ele

         ! savethesun skye dxa
         logical, intent(in), optional :: do_dxa
         logical :: want_dxa
         integer :: m
         type(auto_diff_real_2var_order3) :: F_ii0, Kfac, dF_dabar, dF_dzbar, dF_dye, Fybar, G
         type(auto_diff_real_2var_order3) :: dFc_dAY(species), dFii_dya(species), Fy(species)
         real(dp) :: xh_ci, xhe_ci, zz_ci, abar_ci, zbar_ci, z2bar_ci, z53bar_ci, ye_ci, mc_ci, sumx_ci
         real(dp), dimension(species) :: dabar_dx, dzbar_dx, dmc_dx, xpure
         real(dp) :: Srel, norm_all, sumya, F_off, e_tot, pgas_tot

         ht => eos_ht

         want_dxa = .false.
         if (present(do_dxa)) want_dxa = do_dxa
         if (want_dxa) d_dxa = 0

         ierr = 0

         temp = temp_in
         temp%d1val1 = 1d0
         logtemp = log10(temp)

         den = den_in
         den%d1val2 = 1d0
         logden = log10(den)

         ! HELM table lookup uses din rather than den
         ytot1 = 1.0d0 / abar
         ye = ytot1 * zbar
         din = ye*den

         F_rad = 0d0
         F_ion_gas = 0d0
         F_ideal_ion = 0d0
         F_coul = 0d0
         F_ele = 0d0

         ! Radiation free energy, independent of composition
         F_rad = compute_F_rad(temp, den)

         ! Count and pack relevant species for Coulomb corrections. Relevant means mass fraction above limit.
         relevant_species = 0
         norm = 0d0
         do j=1,species
            if (xa(j) > mass_fraction_limit) then
               relevant_species = relevant_species + 1
               AZION(relevant_species) = chem_isos% Z(chem_id(j))
               ACMI(relevant_species) = chem_isos% W(chem_id(j))
               A(relevant_species) = chem_isos% Z_plus_N(chem_id(j))
               select_xa(relevant_species) = xa(j)
               norm = norm + xa(j)
            end if
         end do

         ! Normalize
         do j=1,relevant_species
            select_xa(j) = select_xa(j) / norm
         end do

         ! Compute number fractions
         norm = 0d0
         do j=1,relevant_species
            ya(j) = select_xa(j) / A(j)
            norm = norm + ya(j)
         end do
         do j=1,relevant_species
            ya(j) = ya(j) / norm
         end do

         ! Ideal ion free energy, only depends on abar
         F_ideal_ion = compute_F_ideal_ion(temp, den, abar, relevant_species, ACMI, ya)
         F_ii0 = F_ideal_ion

         if (use_ion_offsets) then
            F_ideal_ion = F_ideal_ion + compute_ion_offset(species, xa, chem_id)  ! Offset so ion ground state energy is zero.
         end if

         ! Ideal electron-positron thermodynamics (s, e, p)
         ! Derivatives are handled by HELM code, so we don't pass *in* any auto_diff types (just get them as return values).
         call compute_ideal_ele(temp%val, den%val, din%val, logtemp%val, logden%val, zbar, ytot1, ye, ht, &
                               F_ele, etaele, xnefer, ierr)

         xnefer = compute_xne(den, ytot1, zbar)

         ! Normalize mass fractions
         do j=1,relevant_species
            select_xa(j) = select_xa(j) / norm
         end do

         ! Compute non-ideal corrections
         if (want_dxa) then  ! savethesun skye dxa
         call nonideal_corrections(relevant_species, ya(1:relevant_species), &
                                     AZION(1:relevant_species), ACMI(1:relevant_species), &
                                     Skye_min_gamma_for_solid, Skye_max_gamma_for_liquid, &
                                     Skye_solid_mixing_rule, den, temp, xnefer, abar, &
                                     F_coul, latent_ddlnT, latent_ddlnRho, phase, &
                                     dF_dAY=dFc_dAY(1:relevant_species))
         else
         call nonideal_corrections(relevant_species, ya(1:relevant_species), &
                                     AZION(1:relevant_species), ACMI(1:relevant_species), &
                                     Skye_min_gamma_for_solid, Skye_max_gamma_for_liquid, &
                                     Skye_solid_mixing_rule, den, temp, xnefer, abar, &
                                     F_coul, latent_ddlnT, latent_ddlnRho, phase)
         end if

         call  pack_for_export(F_ideal_ion, F_coul, F_rad, F_ele, temp, den, xnefer, etaele, abar, zbar, &
                                 phase, latent_ddlnT, latent_ddlnRho, res, d_dlnd, d_dlnT, ierr)
         if(ierr/=0) return

         if (.not. want_dxa) return

         ! savethesun skye dxa: G_j = dF/dxa_j as a function of (T, rho); then
         ! d lnE/dxa_j = (G - T dG/dT)/E and d lnPgas/dxa_j = rho^2 (dG/drho)/Pgas.
         call composition_info( &
            species, chem_id, xa, xh_ci, xhe_ci, zz_ci, &
            abar_ci, zbar_ci, z2bar_ci, z53bar_ci, ye_ci, mc_ci, &
            sumx_ci, dabar_dx, dzbar_dx, dmc_dx)

         Kfac = kerg * temp / (amu * abar)
         call compute_dF_ideal_ion_dya(temp, den, abar, relevant_species, ACMI, ya, dFii_dya)
         sumya = 0d0
         do m=1,relevant_species
            sumya = sumya + ya(m)
         end do

         ! at fixed ye and ya: ideal ions (F_ii0 excludes the offset) and the Coulomb prefactor kT/(abar amu)
         dF_dabar = -(F_ii0 + Kfac*sumya)/abar - F_coul/abar
         ! electrons: F_ele = ye f(T, ye rho); Coulomb depends on ye only through n_e = ye rho N_A
         dF_dye = (F_ele + den*differentiate_2(F_ele))/ye + den*differentiate_2(F_coul)/ye
         dF_dabar = dF_dabar - dF_dye*ye/abar
         dF_dzbar = dF_dye/abar

         Fybar = 0d0
         do m=1,relevant_species
            Fy(m) = dFii_dya(m) + dFc_dAY(m)
            Fybar = Fybar + ya(m)*Fy(m)
         end do

         ! ya_m = (xa_m/A_m)/Srel over relevant species: d ya_i/d xa_m = (delta_im - ya_i)/(A_m Srel)
         Srel = 0d0
         norm_all = 0d0
         do j=1,species
            norm_all = norm_all + xa(j)/chem_isos% Z_plus_N(chem_id(j))
            if (xa(j) > mass_fraction_limit) Srel = Srel + xa(j)/chem_isos% Z_plus_N(chem_id(j))
         end do
         F_off = 0d0
         if (use_ion_offsets) F_off = compute_ion_offset(species, xa, chem_id)

         e_tot = exp(res(i_lnE))
         pgas_tot = exp(res(i_lnPgas))
         m = 0
         do j=1,species
            ! chem_lib: abar = sumx/sum(y), zbar = sum(y Z)*abar, i.e. both carry the sumx factor. Its dabar_dx
            ! includes that factor, its dzbar_dx does not (d of the mean charge); add it back (+ zbar/sumx)
            ! so G is the derivative of the function as evaluated (what finite differences see).
            G = dF_dabar*dabar_dx(j) + dF_dzbar*(dzbar_dx(j) + zbar_ci/sumx_ci)
            if (xa(j) > mass_fraction_limit) then
               m = m + 1
               G = G + (Fy(m) - Fybar)/(A(m)*Srel)
            end if
            if (use_ion_offsets) then
               ! offset = C sum_k I_k y_k / sum_k y_k, y = xa/A: d/dxa_j = (C I_j - offset)/(A_j sum y)
               xpure = 0d0
               xpure(j) = 1d0
               G = G + (compute_ion_offset(species, xpure, chem_id) - F_off) / &
                  (chem_isos% Z_plus_N(chem_id(j))*norm_all)
            end if
            d_dxa(i_lnE,j) = (G%val - temp%val*G%d1val1)/e_tot
            d_dxa(i_lnPgas,j) = den%val*den%val*G%d1val2/pgas_tot
         end do

      end subroutine skye_eos


end module skye
