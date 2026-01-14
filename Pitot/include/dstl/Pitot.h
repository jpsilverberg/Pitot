// Pitot.hpp - High-performance pitot-static + ISA helpers (complete consolidated
// header) Header-only, C++17
//
// Units:
//   Pressure:     lbf/ft^2
//   Temperature:  Kelvin
//   Density:      slug/ft^3
//   Altitude:     feet (geopotential for ISA relations unless stated)
//   Speeds:       knots (kt)
//   Dynamic q:    lbf/ft^2 (consistent with rho in slug/ft^3 and V in ft/s via gc)
//
// Tier-1 naming only:
//   out_from_in(...)
//   out_from_in_and_atm(...)
//   out_from_in_and_hpc_std(...)
//   out_from_in_and_hpc_std_fast(...)
#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

namespace pitot
{
//==========================================================================
// Configuration
//==========================================================================
#ifndef PITOT_DOMAIN_CHECKS
#define PITOT_DOMAIN_CHECKS 1
#endif

// 1 => NaN on invalid domain, 0 => clamp where reasonable
#ifndef PITOT_DOMAIN_POLICY
#define PITOT_DOMAIN_POLICY 1
#endif

// Standard-day fast approximations for ISA ratios:
// 0 => exact ISA (exp/log), 1 => fast series approximations (use with care)
#ifndef PITOT_STD_FAST_MODE
#define PITOT_STD_FAST_MODE 1
#endif

// (gamma-1)/gamma = 2/7 helper:
// 0 => exp(log(x)*2/7)
// 1 => polynomial for x^(2/7)-1 near x=1 on [1..Pt/P@M1], fallback outside
#ifndef PITOT_POW27_MODE
#define PITOT_POW27_MODE 1
#endif

// Include EAS in bundle outputs (minor code-size toggle)
#ifndef PITOT_INCLUDE_EAS
#define PITOT_INCLUDE_EAS 1
#endif

using Real = double;

//==========================================================================
// Branch hints (C++17)
//==========================================================================
#if defined(__clang__) || defined(__GNUC__)
#define PITOT_LIKELY(x) (__builtin_expect(!!(x), 1))
#define PITOT_UNLIKELY(x) (__builtin_expect(!!(x), 0))
#else
#define PITOT_LIKELY(x) (x)
#define PITOT_UNLIKELY(x) (x)
#endif

//==========================================================================
// Fundamental constants (legacy exactness)
//==========================================================================
inline constexpr Real P0 = Real(2116.217);  // lbf/ft^2  sea level standard pressure
inline constexpr Real T0 = Real(288.15);  // K        sea level standard temperature
inline constexpr Real rho0 =
  Real(0.0023769);                          // slug/ft^3 sea level standard density
inline constexpr Real L = Real(0.0019812);  // K/ft     tropospheric lapse rate
inline constexpr Real h_trop =
  Real(36089.0);  // ft       tropopause altitude (geopotential)
inline constexpr Real P_trop = Real(472.6781);  // lbf/ft^2 pressure at tropopause
inline constexpr Real T_trop = Real(216.65);    // K        temperature at tropopause
inline constexpr Real rho_trop =
  Real(0.000706116);  // slug/ft^3 density at tropopause

inline constexpr Real R     = Real(96.034);     // ft*lbf/(lbm*K)
inline constexpr Real g     = Real(32.174049);  // ft/s^2
inline constexpr Real gc    = Real(32.174049);  // lbm/slug
inline constexpr Real gamma = Real(1.4);        // ratio of specific heats

inline constexpr Real kt_to_fps = Real(1.6878098571);
inline constexpr Real fps_to_kt = Real(0.592483801296);

inline constexpr Real lbf_ft_per_slug_to_kt2 =
  Real(0.351037054798);  // lbf-ft/slug -> kt^2 scaling
inline constexpr Real lbf_ft_per_lbm_to_kt2 =
  Real(11.2942832462);  // lbf-ft/lbm  -> kt^2 scaling

inline constexpr Real a0_kt =
  Real(661.477652630007);  // kt speed of sound at SL standard

inline constexpr Real ft_per_nm   = Real(6076.1155);
inline constexpr Real sec_per_min = Real(60);
inline constexpr Real min_per_sec = Real(1) / Real(60);

// Geometric <-> geopotential conversions (approx)
inline constexpr Real earth_radius_ft = Real(20855531.5);

// ISA derived "geopotential scale" constant used by legacy GPS/alt correction
inline constexpr Real C0_GEOPOT_SCALE = Real(145442.15695605625);

//==========================================================================
// Derived constants
//==========================================================================
inline constexpr Real inv_T0        = Real(1) / T0;
inline constexpr Real lapse_over_T0 = L * inv_T0;

// Exponents derived from ISA:
// delta = theta^(g/(gc*L*R))
// sigma = theta^(g/(gc*L*R) - 1)
inline constexpr Real exp_delta_trop = g / (gc * L * R);  // ~5.255863
inline constexpr Real exp_sigma_trop = exp_delta_trop - Real(1);

// Strat exponential constant:
inline constexpr Real k_strat = g / (gc * R * T_trop);

// Helpful thresholds at tropopause (exact ratios from constants)
inline constexpr Real theta_trop = T_trop / T0;      // ~0.75186699
inline constexpr Real delta_trop = P_trop / P0;      // ~0.22336248
inline constexpr Real sigma_trop = rho_trop / rho0;  // ~0.29707711

// Isentropic M=1 thresholds:
inline constexpr Real t_ratio_at_m1 = Real(1) + (gamma - Real(1)) / Real(2);  // 1.2
inline constexpr Real pt_over_p_at_m1 =
  Real(1.892929158737854);  // Pt/P at M=1 for γ=1.4 = (1.2)^(3.5)

inline constexpr Real five = Real(5);  // 2/(gamma-1) for gamma=1.4

//==========================================================================
// Numeric helpers
//==========================================================================
inline Real nan() noexcept
{
   return std::numeric_limits<Real>::quiet_NaN();
}

inline Real safe_sqrt(Real x) noexcept
{
#if PITOT_DOMAIN_CHECKS
   return (x >= Real(0)) ? std::sqrt(x) : (PITOT_DOMAIN_POLICY ? nan() : Real(0));
#else
   return std::sqrt(x);
#endif
}

inline Real safe_log(Real x) noexcept
{
#if PITOT_DOMAIN_CHECKS
   if (PITOT_UNLIKELY(x <= Real(0)))
      return (PITOT_DOMAIN_POLICY ? nan() : -std::numeric_limits<Real>::infinity());
#endif
   return std::log(x);
}

inline Real safe_exp(Real x) noexcept
{
   return std::exp(x);
}

// "Safe pow" for general exponent (rarely needed in hot paths).
// For our library, most exponents are fixed positive rationals and bases are
// positive ratios.
inline Real safe_pow(Real base, Real exp) noexcept
{
#if PITOT_DOMAIN_CHECKS
   if (PITOT_UNLIKELY(base <= Real(0)))
   {
      const Real i    = std::trunc(exp);
      const Real diff = std::abs(exp - i);
      // allow tiny epsilon for "should be integer" exponents
      if (diff > Real(1e-12))
         return (PITOT_DOMAIN_POLICY ? nan() : Real(0));
   }
#endif
   return std::pow(base, exp);
}

//==========================================================================
// Fast exponent helpers specialized to gamma=1.4
//==========================================================================
// x^(3.5) = x^3 * sqrt(x)
inline Real pow_3_5_pos(Real x) noexcept
{
#if PITOT_DOMAIN_CHECKS
   if (PITOT_UNLIKELY(x <= Real(0)))
      return (PITOT_DOMAIN_POLICY ? nan() : Real(0));
#endif
   const Real s  = std::sqrt(x);
   const Real x2 = x * x;
   return (x2 * x) * s;
}

inline Real pow_2_over_7_pos_mode0(Real x) noexcept
{
#if PITOT_DOMAIN_CHECKS
   if (PITOT_UNLIKELY(x <= Real(0)))
      return (PITOT_DOMAIN_POLICY ? nan() : Real(0));
#endif
   return std::exp(std::log(x) * (Real(2) / Real(7)));
}

// Polynomial approximation for x^(2/7)-1 around x=1, good on [1 .. Pt/P@M1]
inline Real term_pow_2_over_7_minus1_subsonic_poly(Real x) noexcept
{
   const Real u = x - Real(1);
   if (PITOT_UNLIKELY(std::abs(u) < Real(1e-14)))
      return (Real(2) / Real(7)) * u;

   constexpr Real c1 = Real(0.285714285714285714285714285714);
   constexpr Real c2 = Real(-0.102040816326530612244897959184);
   constexpr Real c3 = Real(0.0583090379008746355685131195335);
   constexpr Real c4 = Real(-0.0395668471470220741357767596835);
   constexpr Real c5 = Real(0.0293925150235021122151484500506);
   constexpr Real c6 = Real(-0.0230941189470373738833309250397);
   constexpr Real c7 = Real(0.0188523419975815297006783061549);
   constexpr Real c8 = Real(-0.0158225013193987838559264355229);

   return u *
          (c1 +
            u *
              (c2 + u * (c3 + u * (c4 + u * (c5 + u * (c6 + u * (c7 + u * c8)))))));
}

inline Real term_pow_2_over_7_minus1(Real x) noexcept
{
#if PITOT_POW27_MODE == 1
   if (PITOT_LIKELY(x >= Real(1) && x <= pt_over_p_at_m1))
      return term_pow_2_over_7_minus1_subsonic_poly(x);
   return pow_2_over_7_pos_mode0(x) - Real(1);
#else
   return pow_2_over_7_pos_mode0(x) - Real(1);
#endif
}

//==========================================================================
// Original helper layer: dimensionless ratios and conversions
//==========================================================================

// Pressure ratio
inline Real delta_from_P(Real P) noexcept
{
   return P / P0;
}
inline Real P_from_delta(Real d) noexcept
{
   return d * P0;
}

// Temperature ratio
inline Real theta_from_Tk(Real Tk) noexcept
{
   return Tk / T0;
}
inline Real Tk_from_theta(Real th) noexcept
{
   return th * T0;
}

// Density ratio
inline Real sigma_from_rho(Real rho) noexcept
{
   return rho / rho0;
}
inline Real rho_from_sigma(Real s) noexcept
{
   return s * rho0;
}

// delta = sigma * theta (ideal gas)
inline Real delta_from_sigma_theta(Real sigma, Real theta) noexcept
{
   return sigma * theta;
}
inline Real sigma_from_delta_theta(Real delta, Real theta) noexcept
{
   return delta / theta;
}
inline Real theta_from_delta_sigma(Real delta, Real sigma) noexcept
{
   return delta / sigma;
}

// Unit conversions
inline Real fps_from_kt(Real v_kt) noexcept
{
   return v_kt * kt_to_fps;
}
inline Real kt_from_fps(Real v_fps) noexcept
{
   return v_fps * fps_to_kt;
}

//==========================================================================
// ISA: theta/sigma/delta as functions of geopotential altitude (h)
//==========================================================================
inline Real theta_from_h(Real h_ft) noexcept
{
   if (h_ft < h_trop)
      return Real(1) - lapse_over_T0 * h_ft;
   return theta_trop;
}

inline Real sigma_from_h(Real h_ft) noexcept
{
   if (h_ft < h_trop)
   {
      const Real th = theta_from_h(h_ft);
      return std::exp(exp_sigma_trop * std::log(th));
   }
   const Real dh = h_ft - h_trop;
   return sigma_trop * std::exp(-k_strat * dh);
}

inline Real delta_from_h(Real h_ft) noexcept
{
   if (h_ft < h_trop)
   {
      const Real th = theta_from_h(h_ft);
      return std::exp(exp_delta_trop * std::log(th));
   }
   const Real dh = h_ft - h_trop;
   return delta_trop * std::exp(-k_strat * dh);
}

// Inverses: h from theta/delta/sigma (pressure altitude under ISA)
inline Real h_from_theta(Real theta) noexcept
{
   if (theta > theta_trop)
   {
      // theta = 1 - (L/T0) h
      return (Real(1) - theta) / lapse_over_T0;
   }
   // Strat: theta constant (tropopause)
   return h_trop;  // ambiguous; ISA "theta" alone does not encode height above
                   // tropopause
}

inline Real h_from_delta(Real delta) noexcept
{
   if (delta > delta_trop)
   {
      // delta = theta^exp_delta_trop ; theta = delta^(1/exp) ; h = (1-theta)/(L/T0)
      const Real theta = std::exp((Real(1) / exp_delta_trop) * std::log(delta));
      return (Real(1) - theta) / lapse_over_T0;
   }
   // Strat: delta = delta_trop * exp(-k*(h-h_trop))
   const Real ratio = delta / delta_trop;
   return h_trop - (Real(1) / k_strat) * std::log(ratio);
}

inline Real h_from_sigma(Real sigma) noexcept
{
   if (sigma > sigma_trop)
   {
      // sigma = theta^exp_sigma_trop
      const Real theta = std::exp((Real(1) / exp_sigma_trop) * std::log(sigma));
      return (Real(1) - theta) / lapse_over_T0;
   }
   const Real ratio = sigma / sigma_trop;
   return h_trop - (Real(1) / k_strat) * std::log(ratio);
}

// ISA actual values from h
inline Real Tk_from_h(Real h_ft) noexcept
{
   return T0 * theta_from_h(h_ft);
}
inline Real P_from_h(Real h_ft) noexcept
{
   return P0 * delta_from_h(h_ft);
}
inline Real rho_from_h(Real h_ft) noexcept
{
   return rho0 * sigma_from_h(h_ft);
}

//==========================================================================
// Pressure altitude (hpc) helpers under ISA
//==========================================================================
inline Real hpc_from_Ps(Real Ps) noexcept
{
   return h_from_delta(delta_from_P(Ps));
}
inline Real Ps_from_hpc(Real hpc_ft) noexcept
{
   return P_from_h(hpc_ft);
}

//==========================================================================
// Temperature helpers (originals)
//==========================================================================
inline Real Tk_from_degC(Real Tc) noexcept
{
   return Tc + Real(273.15);
}
inline Real degC_from_Tk(Real Tk) noexcept
{
   return Tk - Real(273.15);
}

inline Real theta_from_degC(Real Tc) noexcept
{
   return (Tc + Real(273.15)) / T0;
}

inline Real h_from_Tk_std(Real Tk) noexcept
{
   return h_from_theta(theta_from_Tk(Tk));
}  // ISA only

// Total/impact temperature relations (legacy "k" factor for recovery)
inline Real Tk_from_OATi_degC_k_M(Real OATi_degC, Real k, Real M) noexcept
{
   // TaK = (OATi+273.15) / (1 + (gamma-1)/2 * k * M^2)
   const Real Tt = OATi_degC + Real(273.15);
   return Tt / (Real(1) + (gamma - Real(1)) / Real(2) * k * M * M);
}

inline Real OATiK_from_Tk_M_k(Real Tk, Real M, Real k) noexcept
{
   return Tk * (Real(1) + (gamma - Real(1)) / Real(2) * k * M * M);
}

// Speed of sound from temperature (knots)
inline Real a_from_Tk(Real Tk) noexcept
{
   // a = sqrt(gamma*R*Tk * (lbf*ft/lbm -> kt^2 scaling))
   return safe_sqrt(gamma * R * Tk * lbf_ft_per_lbm_to_kt2);
}

inline Real Tk_from_a(Real a_kt) noexcept
{
   return (a_kt * a_kt) / (gamma * R * lbf_ft_per_lbm_to_kt2);
}

//==========================================================================
// Density relations (original)
//==========================================================================
inline Real h_from_rho_std(Real rho) noexcept
{
   return h_from_sigma(sigma_from_rho(rho));
}
inline Real rho_from_h_std(Real h_ft) noexcept
{
   return rho_from_h(h_ft);
}

//==========================================================================
// Core compressible primitives: qc <-> CAS (Vc) at sea level
//==========================================================================
inline Real qc_from_vc(Real vc_kt) noexcept
{
   const Real term = Real(1) + (vc_kt * vc_kt) * (gamma - Real(1)) /
                                 (Real(2) * gamma) * (rho0 / P0) /
                                 lbf_ft_per_slug_to_kt2;

   // qc = P0 * (term^(gamma/(gamma-1)) - 1) ; gamma/(gamma-1)=3.5
   return (pow_3_5_pos(term) - Real(1)) * P0;
}

inline Real vc_from_qc(Real qc) noexcept
{
   const Real base = Real(1) + qc / P0;               // Pt/P at sea level
   const Real term = term_pow_2_over_7_minus1(base);  // base^(2/7) - 1

   constexpr Real coef =
     (Real(2) * gamma / (gamma - Real(1))) * (P0 / rho0) * lbf_ft_per_slug_to_kt2;

   return safe_sqrt(coef * term);
}

// Indicated airspeed (Vi) uses qci (impact pressure based on Ps)
inline Real Vi_from_qci(Real qci) noexcept
{
   const Real base = Real(1) + qci / P0;
   constexpr Real coef =
     (Real(2) * gamma / (gamma - Real(1))) * (P0 / rho0) * lbf_ft_per_slug_to_kt2;

   const Real term = term_pow_2_over_7_minus1(base);
   return safe_sqrt(coef * term);
}

inline Real qci_from_Vi(Real Vi_kt) noexcept
{
   const Real term = Real(1) + (Vi_kt * Vi_kt) * (gamma - Real(1)) /
                                 (Real(2) * gamma) * (rho0 / P0) /
                                 lbf_ft_per_slug_to_kt2;

   return (pow_3_5_pos(term) - Real(1)) * P0;
}

// Equivalent airspeed Ve from qc and ambient static pressure Pa (using rho0 as
// reference)
inline Real Ve_from_qc_Ps(Real qc, Real Ps) noexcept
{
   const Real base = Real(1) + qc / Ps;
   constexpr Real coef =
     (Real(2) * gamma / (gamma - Real(1))) * lbf_ft_per_slug_to_kt2;

   const Real term = term_pow_2_over_7_minus1(base);
   return safe_sqrt(coef * (Ps / rho0) * term);
}

inline Real qc_from_Ve_Ps(Real Ve_kt, Real Ps) noexcept
{
   const Real term = Real(1) + (Ve_kt * Ve_kt) * (gamma - Real(1)) /
                                 (Real(2) * gamma) * (rho0 / Ps) /
                                 lbf_ft_per_slug_to_kt2;

   return (pow_3_5_pos(term) - Real(1)) * Ps;
}

// Dynamic pressure from equivalent airspeed (SL density)
inline Real q_from_Ve(Real Ve_kt) noexcept
{
   // q = 0.5 * rho0 * V^2 (with V in ft/s)
   const Real Vfps = Ve_kt * kt_to_fps;
   return Real(0.5) * rho0 * (Vfps * Vfps);
}

inline Real Ve_from_q(Real q) noexcept
{
   const Real Vfps = safe_sqrt((Real(2) * q) / rho0);
   return Vfps * fps_to_kt;
}

//==========================================================================
// Mach relations
//==========================================================================
// Mach from qc and static pressure Ps:
inline Real mach_from_qc_Ps(Real qc, Real Ps) noexcept
{
   const Real base = Real(1) + qc / Ps;               // Pt/P
   const Real term = term_pow_2_over_7_minus1(base);  // (Pt/P)^(2/7) - 1
   return safe_sqrt(five * term);                     // M = sqrt( 5*(...))
}

inline Real qc_from_mach_Ps(Real M, Real Ps) noexcept
{
   const Real t = Real(1) + (gamma - Real(1)) / Real(2) * M * M;  // Tt/T
   return Ps * (pow_3_5_pos(t) - Real(1));
}

// Total pressure Pt and static Ps with Mach:
inline Real Pt_from_Ps_M(Real Ps, Real M) noexcept
{
   const Real t = Real(1) + (gamma - Real(1)) / Real(2) * M * M;
   return Ps * pow_3_5_pos(t);
}

inline Real Ps_from_Pt_M(Real Pt, Real M) noexcept
{
   const Real t = Real(1) + (gamma - Real(1)) / Real(2) * M * M;
   return Pt / pow_3_5_pos(t);
}

// Mach from Pt/Ps (isentropic):
inline Real M_from_Pt_over_Ps(Real Pt_over_Ps) noexcept
{
   const Real term = term_pow_2_over_7_minus1(Pt_over_Ps);
   return safe_sqrt(five * term);
}

//==========================================================================
// Performance / kinematic relations (original)
//==========================================================================
inline Real q_from_rho_Vt(Real rho, Real Vt_kt) noexcept
{
   const Real Vfps = Vt_kt * kt_to_fps;
   return Real(0.5) * rho * (Vfps * Vfps);
}

inline Real rho_from_q_Vt(Real q, Real Vt_kt) noexcept
{
   const Real Vfps  = Vt_kt * kt_to_fps;
   const Real denom = Real(0.5) * (Vfps * Vfps);
   return q / denom;
}

inline Real Vt_from_q_rho(Real q, Real rho) noexcept
{
   const Real Vfps = safe_sqrt((Real(2) * q) / rho);
   return Vfps * fps_to_kt;
}

//==========================================================================
// GPS/alt correction helpers (legacy equations preserved)
//==========================================================================
inline Real hcorr_from_gps0_alt0(Real hGPS0, Real hALT0) noexcept
{
   // corr = (1 - k*hALT0)^(5.255863) / (1 - k*hGPS0)^(5.255863)
   constexpr Real k  = Real(6.8755856e-6);
   constexpr Real e1 = Real(5.255863);
   return safe_pow(Real(1) - k * hALT0, e1) / safe_pow(Real(1) - k * hGPS0, e1);
}

inline Real hALT_from_gps1_corr(Real hGPS1, Real corr) noexcept
{
   constexpr Real k  = Real(6.8755856e-6);
   constexpr Real e1 = Real(5.255863);
   constexpr Real e2 = Real(0.1902637112116507);

   const Real tGPS1 = safe_pow(Real(1) - k * hGPS1, e1);
   return C0_GEOPOT_SCALE * (Real(1) - safe_pow(corr * tGPS1, e2));
}

inline Real hALT_from_gps0_alt0_gps1(Real hGPS0, Real hALT0, Real hGPS1) noexcept
{
   constexpr Real k  = Real(6.8755856e-6);
   constexpr Real e1 = Real(5.255863);
   constexpr Real e2 = Real(0.1902637112116507);

   const Real tALT0 = safe_pow(Real(1) - k * hALT0, e1);
   const Real tGPS1 = safe_pow(Real(1) - k * hGPS1, e1);
   const Real tGPS0 = safe_pow(Real(1) - k * hGPS0, e1);

   const Real ratio = (tALT0 * tGPS1) / tGPS0;
   return C0_GEOPOT_SCALE * (Real(1) - safe_pow(ratio, e2));
}

//==========================================================================
// Atmosphere container + builders
//==========================================================================
struct Atm
{
   Real h_ft = Real(0);  // geopotential altitude used for ISA properties

   Real Tk  = T0;
   Real Ps  = P0;
   Real rho = rho0;

   Real theta = Real(1);
   Real delta = Real(1);
   Real sigma = Real(1);

   Real a_kt   = a0_kt;
   Real inv_Ps = Real(1) / P0;

   Real sqrt_sigma     = Real(1);
   Real inv_sqrt_sigma = Real(1);
};

inline Atm standard_day_atm_from_hpc(Real hpc_ft) noexcept
{
   Atm atm;
   atm.h_ft  = hpc_ft;
   atm.theta = theta_from_h(hpc_ft);
   atm.delta = delta_from_h(hpc_ft);
   atm.sigma = sigma_from_h(hpc_ft);

   atm.Tk  = T0 * atm.theta;
   atm.Ps  = P0 * atm.delta;
   atm.rho = rho0 * atm.sigma;

   atm.inv_Ps = (atm.Ps != Real(0)) ? (Real(1) / atm.Ps)
                                    : (PITOT_DOMAIN_POLICY ? nan() : Real(0));
   atm.a_kt   = a_from_Tk(atm.Tk);

   atm.sqrt_sigma     = safe_sqrt(atm.sigma);
   atm.inv_sqrt_sigma = (atm.sqrt_sigma > Real(0))
                          ? (Real(1) / atm.sqrt_sigma)
                          : (PITOT_DOMAIN_POLICY ? nan() : Real(0));
   return atm;
}

// Fast standard-day builder: uses series approximations for theta^exp and exp(-k*dh)
// This is optional speed-over-accuracy; keep PITOT_STD_FAST_MODE=0 for exact.
inline Real log1m_series(Real x) noexcept
{
   const Real x2 = x * x;
   const Real x3 = x2 * x;
   const Real x4 = x3 * x;
   const Real x5 = x4 * x;
   return -(x + Real(0.5) * x2 + (Real(1) / Real(3)) * x3 + Real(0.25) * x4 +
            (Real(1) / Real(5)) * x5);
}

inline Real exp_series_neg(Real y) noexcept
{
   const Real y2 = y * y;
   const Real y3 = y2 * y;
   const Real y4 = y3 * y;
   const Real y5 = y4 * y;
   return Real(1) - y + Real(0.5) * y2 - (Real(1) / Real(6)) * y3 +
          (Real(1) / Real(24)) * y4 - (Real(1) / Real(120)) * y5;
}

inline Real delta_from_h_fast(Real h_ft) noexcept
{
#if PITOT_STD_FAST_MODE == 0
   return delta_from_h(h_ft);
#else
   if (h_ft < h_trop)
   {
      const Real theta    = Real(1) - lapse_over_T0 * h_ft;
      const Real x        = Real(1) - theta;
      const Real ln_theta = log1m_series(x);
      const Real y        = -exp_delta_trop * ln_theta;
      return exp_series_neg(y);
   }
   const Real dh = h_ft - h_trop;
   const Real y  = k_strat * dh;
   return delta_trop * exp_series_neg(y);
#endif
}

inline Real sigma_from_h_fast(Real h_ft) noexcept
{
#if PITOT_STD_FAST_MODE == 0
   return sigma_from_h(h_ft);
#else
   if (h_ft < h_trop)
   {
      const Real theta    = Real(1) - lapse_over_T0 * h_ft;
      const Real x        = Real(1) - theta;
      const Real ln_theta = log1m_series(x);
      const Real y        = -exp_sigma_trop * ln_theta;
      return exp_series_neg(y);
   }
   const Real dh = h_ft - h_trop;
   const Real y  = k_strat * dh;
   return sigma_trop * exp_series_neg(y);
#endif
}

inline Atm standard_day_atm_from_hpc_std_fast(Real hpc_ft) noexcept
{
#if PITOT_STD_FAST_MODE == 0
   return standard_day_atm_from_hpc(hpc_ft);
#else
   Atm atm;
   atm.h_ft  = hpc_ft;
   atm.theta = theta_from_h(hpc_ft);  // exact (already trivial)
   atm.delta = delta_from_h_fast(hpc_ft);
   atm.sigma = sigma_from_h_fast(hpc_ft);

   atm.Tk  = T0 * atm.theta;
   atm.Ps  = P0 * atm.delta;
   atm.rho = rho0 * atm.sigma;

   atm.inv_Ps = (atm.Ps != Real(0)) ? (Real(1) / atm.Ps)
                                    : (PITOT_DOMAIN_POLICY ? nan() : Real(0));
   atm.a_kt   = a_from_Tk(atm.Tk);

   atm.sqrt_sigma     = safe_sqrt(atm.sigma);
   atm.inv_sqrt_sigma = (atm.sqrt_sigma > Real(0))
                          ? (Real(1) / atm.sqrt_sigma)
                          : (PITOT_DOMAIN_POLICY ? nan() : Real(0));
   return atm;
#endif
}

//==========================================================================
// Extended atmosphere models (ISA variants: ΔT, humidity) + Atm solvers/builders
// Append here: right after standard_day_atm_from_hpc_std_fast(...)
//==========================================================================

namespace detail
{
inline Real clamp01(Real x) noexcept
{
   return std::max(Real(0), std::min(Real(1), x));
}

// Convert pressure from lbf/ft^2 to Pa (SI) and back.
// 1 lbf/ft^2 = 47.88025898033584 Pa (exact from lbf and ft definitions)
inline constexpr Real PA_PER_LBF_FT2 = Real(47.88025898033584);
inline constexpr Real LBF_FT2_PER_PA = Real(1) / PA_PER_LBF_FT2;

inline Real Pa_from_lbf_ft2(Real p) noexcept
{
   return p * PA_PER_LBF_FT2;
}
inline Real lbf_ft2_from_Pa(Real p) noexcept
{
   return p * LBF_FT2_PER_PA;
}

// Magnus-Tetens approximation over water, valid roughly -45°C to +60°C.
// For extreme cold/high altitude, consider a different fit.
inline Real saturation_vapor_pressure_Pa_from_Tc(Real Tc) noexcept
{
   // e_s [Pa] = 610.94 * exp(17.625*T / (T + 243.04)) for T in °C
   const Real a = Real(17.625);
   const Real b = Real(243.04);
   return Real(610.94) * std::exp(a * Tc / (Tc + b));
}

inline Real vapor_pressure_Pa_from_Tc_RH(Real Tc, Real RH_0to1) noexcept
{
   const Real rh = clamp01(RH_0to1);
   return rh * saturation_vapor_pressure_Pa_from_Tc(Tc);
}

// Compute moist air density using partial pressures:
// rho = Pd/(Rd*T) + Pv/(Rv*T)
// We do the calculation in SI then convert to slug/ft^3.
// Rd = 287.05 J/(kg*K), Rv = 461.495 J/(kg*K)
// 1 kg/m^3 = 0.00194032033 slug/ft^3
inline constexpr Real RD_SI              = Real(287.05);
inline constexpr Real RV_SI              = Real(461.495);
inline constexpr Real SLUG_FT3_PER_KG_M3 = Real(0.00194032033);

inline Real rho_slugft3_from_P_T_RH(Real P_lbf_ft2, Real Tk, Real RH_0to1) noexcept
{
   const Real P_Pa  = Pa_from_lbf_ft2(P_lbf_ft2);
   const Real Tc    = Tk - Real(273.15);
   const Real Pv_Pa = vapor_pressure_Pa_from_Tc_RH(Tc, RH_0to1);
   const Real Pd_Pa = std::max(Real(0), P_Pa - Pv_Pa);

   // kg/m^3
   const Real rho_si = (Pd_Pa / (RD_SI * Tk)) + (Pv_Pa / (RV_SI * Tk));
   return rho_si * SLUG_FT3_PER_KG_M3;
}

// Dry air density using pitot's imperial gas constant:
// P = rho * (R/gc) * T  => rho = P * gc / (R * T)
inline Real rho_dry_from_P_T(Real P_lbf_ft2, Real Tk) noexcept
{
   // rho [slug/ft^3]
   return (P_lbf_ft2 * gc) / (R * Tk);
}

// Sutherland viscosity (US Std Atmos 1976 typical constants).
// Returns dynamic viscosity mu in slug/(ft*s).
inline Real mu_from_Tk(Real Tk) noexcept
{
   // Sutherland constant (K)
   constexpr Real S = Real(110.4);

   // beta in slug / (ft*s*sqrt(K))
   // (Converted from SI; value chosen for practical accuracy)
   constexpr Real beta_imperial = Real(3.62564354e-7);

   return (beta_imperial * safe_sqrt(Tk) * Tk) / (Tk + S);
}

inline Real reynolds_number_per_ft(Real rho, Real V_kt, Real mu) noexcept
{
   // Re/L = rho*V / mu (V in ft/s)
   return (rho * (V_kt * kt_to_fps)) / mu;
}
}  // namespace detail

//---------------------------------------------------------------------------
// ISA variant builders
//---------------------------------------------------------------------------

// Build Atm under ISA pressure-altitude (hpc) but with measured OAT and optional RH.
// - Ps follows ISA (from hpc)
// - Tk is user-provided (OAT)
// - rho computed either dry or moist (if RH provided)
// Use this when hpc is "pressure altitude under ISA", but you have real
// temperature/humidity.
inline Atm atm_from_hpc_Tk_RH(Real hpc_ft, Real Tk, Real RH_0to1) noexcept
{
   Atm atm = standard_day_atm_from_hpc(hpc_ft);

   // Replace temperature ratio with actual measured temperature
   atm.Tk    = Tk;
   atm.theta = Tk / T0;

   // Pressure remains ISA from hpc (pressure altitude definition)
   // atm.Ps already set, atm.delta already set

   // Density: moist if RH>0 else dry
   if (RH_0to1 > Real(0))
      atm.rho = detail::rho_slugft3_from_P_T_RH(atm.Ps, Tk, RH_0to1);
   else
      atm.rho = detail::rho_dry_from_P_T(atm.Ps, Tk);

   atm.sigma = atm.rho / rho0;

   atm.inv_Ps = (atm.Ps != Real(0)) ? (Real(1) / atm.Ps)
                                    : (PITOT_DOMAIN_POLICY ? nan() : Real(0));
   atm.a_kt   = a_from_Tk(atm.Tk);

   atm.sqrt_sigma     = safe_sqrt(atm.sigma);
   atm.inv_sqrt_sigma = (atm.sqrt_sigma > Real(0))
                          ? (Real(1) / atm.sqrt_sigma)
                          : (PITOT_DOMAIN_POLICY ? nan() : Real(0));
   return atm;
}

// Convenience: OAT in °C + RH
inline Atm atm_from_hpc_OATdegC_RH(Real hpc_ft, Real oat_degC, Real RH_0to1) noexcept
{
   return atm_from_hpc_Tk_RH(hpc_ft, Tk_from_degC(oat_degC), RH_0to1);
}

// ISA + ΔT model (cold/hot day): ΔT is “temperature deviation from ISA” in K at that
// altitude.
// - Ps follows ISA (hpc)
// - Tk = Tk_ISA(hpc) + ΔT
// - rho computed dry (default) or moist if RH provided.
inline Atm atm_from_hpc_dT_RH(Real hpc_ft, Real dT_K, Real RH_0to1) noexcept
{
   const Atm isa = standard_day_atm_from_hpc(hpc_ft);
   return atm_from_hpc_Tk_RH(hpc_ft, isa.Tk + dT_K, RH_0to1);
}

inline Atm atm_from_hpc_dT(Real hpc_ft, Real dT_K) noexcept
{
   return atm_from_hpc_dT_RH(hpc_ft, dT_K, Real(0));
}

//---------------------------------------------------------------------------
// Builders from measured static pressure / temperature (and optionally RH)
//---------------------------------------------------------------------------

// If you have measured Ps directly, you can derive hpc under ISA by inversion.
// Then combine with measured Tk / RH for density and speed of sound.
inline Atm atm_from_Ps_Tk_RH(Real Ps_lbf_ft2, Real Tk, Real RH_0to1) noexcept
{
   const Real hpc = hpc_from_Ps(Ps_lbf_ft2);
   Atm atm        = atm_from_hpc_Tk_RH(hpc, Tk, RH_0to1);

   // But overwrite pressure with the measured one (still consistent with hpc
   // meaning)
   atm.Ps     = Ps_lbf_ft2;
   atm.delta  = atm.Ps / P0;
   atm.inv_Ps = (atm.Ps != Real(0)) ? (Real(1) / atm.Ps)
                                    : (PITOT_DOMAIN_POLICY ? nan() : Real(0));

   // Recompute rho using measured Ps (important)
   if (RH_0to1 > Real(0))
      atm.rho = detail::rho_slugft3_from_P_T_RH(atm.Ps, Tk, RH_0to1);
   else
      atm.rho = detail::rho_dry_from_P_T(atm.Ps, Tk);

   atm.sigma = atm.rho / rho0;

   atm.a_kt           = a_from_Tk(atm.Tk);
   atm.sqrt_sigma     = safe_sqrt(atm.sigma);
   atm.inv_sqrt_sigma = (atm.sqrt_sigma > Real(0))
                          ? (Real(1) / atm.sqrt_sigma)
                          : (PITOT_DOMAIN_POLICY ? nan() : Real(0));
   return atm;
}

inline Atm atm_from_Ps_Tk(Real Ps_lbf_ft2, Real Tk) noexcept
{
   return atm_from_Ps_Tk_RH(Ps_lbf_ft2, Tk, Real(0));
}

inline Atm atm_from_Ps_OATdegC_RH(
  Real Ps_lbf_ft2, Real oat_degC, Real RH_0to1) noexcept
{
   return atm_from_Ps_Tk_RH(Ps_lbf_ft2, Tk_from_degC(oat_degC), RH_0to1);
}

//---------------------------------------------------------------------------
// Pitot-measurement-based “solver” style builders
//---------------------------------------------------------------------------
// Typical measured channels:
//  - Ps (static) (or hpc already computed by avionics)
//  - Pt (total) or qc = Pt - Ps
//  - OAT (Tk)
//  - optionally RH, recovery factor k
//
// These builders create an Atm and also compute Mach/airspeeds if desired.
// Keep it modular: Atm is atmosphere only; use existing
// from_qc_and_atm/from_vt_and_atm bundles.

inline Atm atm_from_hpc_OATdegC_dT_optional_RH(Real hpc_ft,
  Real oat_degC,
  Real dT_override_K,  // use NaN to mean "no override"
  Real RH_0to1) noexcept
{
   const Atm isa = standard_day_atm_from_hpc(hpc_ft);
   Real Tk       = Tk_from_degC(oat_degC);

   if (!std::isnan(dT_override_K))
      Tk = isa.Tk + dT_override_K;

   return atm_from_hpc_Tk_RH(hpc_ft, Tk, RH_0to1);
}

//---------------------------------------------------------------------------
// Optional extras commonly needed in aero computations
//---------------------------------------------------------------------------

inline Real mu_from_Tk(Real Tk) noexcept
{
   return detail::mu_from_Tk(Tk);
}

inline Real reynolds_number_per_ft(Real rho, Real V_kt, Real mu) noexcept
{
   return detail::reynolds_number_per_ft(rho, V_kt, mu);
}

// Crossover (CAS/Mach transition): delta at which constant CAS equals constant Mach.
// Uses your core qc and isentropic ratio.
inline Real delta_crossover_from_vc_mach(Real vc_kt, Real mach) noexcept
{
   // qc for CAS target (sea-level definition)
   const Real qc_cas = qc_from_vc(vc_kt);

   // qc/P for Mach target: (1 + 0.2 M^2)^(3.5) - 1
   const Real t_mach     = Real(1) + (gamma - Real(1)) / Real(2) * mach * mach;
   const Real ratio_mach = pow_3_5_pos(t_mach) - Real(1);

   // P_crossover = qc_cas / ratio_mach ; delta = P/P0
   return (qc_cas / ratio_mach) / P0;
}

inline Real hpc_crossover_from_vc_mach(Real vc_kt, Real mach) noexcept
{
   return h_from_delta(delta_crossover_from_vc_mach(vc_kt, mach));
}

//==========================================================================
// Tier-1 combined helpers: *_and_atm
//==========================================================================
inline Real mach_from_qc_and_atm(Real qc, const Atm& atm) noexcept
{
   return mach_from_qc_Ps(qc, atm.Ps);
}

inline Real qc_from_mach_and_atm(Real M, const Atm& atm) noexcept
{
   return qc_from_mach_Ps(M, atm.Ps);
}

inline Real vt_from_mach_and_atm(Real M, const Atm& atm) noexcept
{
   return M * atm.a_kt;
}

inline Real mach_from_vt_and_atm(Real vt_kt, const Atm& atm) noexcept
{
   return vt_kt / atm.a_kt;
}

inline Real vt_from_vc_and_atm(Real vc_kt, const Atm& atm) noexcept
{
   const Real qc = qc_from_vc(vc_kt);
   const Real M  = mach_from_qc_and_atm(qc, atm);
   return vt_from_mach_and_atm(M, atm);
}

inline Real vc_from_vt_and_atm(Real vt_kt, const Atm& atm) noexcept
{
   const Real M  = mach_from_vt_and_atm(vt_kt, atm);
   const Real qc = qc_from_mach_and_atm(M, atm);
   return vc_from_qc(qc);
}

inline Real vc_from_mach_and_atm(Real M, const Atm& atm) noexcept
{
   const Real qc = qc_from_mach_and_atm(M, atm);
   return vc_from_qc(qc);
}

inline Real mach_from_vc_and_atm(Real vc_kt, const Atm& atm) noexcept
{
   const Real qc = qc_from_vc(vc_kt);
   return mach_from_qc_and_atm(qc, atm);
}

//==========================================================================
// Bundled helpers (combined; no redundant recomputation)
//==========================================================================
struct FromVc
{
   Real qc    = Real(0);
   Real mach  = Real(0);
   Real vt_kt = Real(0);
#if PITOT_INCLUDE_EAS
   Real eas_kt = Real(0);
#endif
};

struct FromVt
{
   Real mach  = Real(0);
   Real qc    = Real(0);
   Real vc_kt = Real(0);
#if PITOT_INCLUDE_EAS
   Real eas_kt = Real(0);
#endif
};

struct FromQc
{
   Real mach  = Real(0);
   Real vt_kt = Real(0);
   Real vc_kt = Real(0);
#if PITOT_INCLUDE_EAS
   Real eas_kt = Real(0);
#endif
};

inline FromVc from_vc_and_atm(Real vc_kt, const Atm& atm) noexcept
{
   FromVc out;
   out.qc    = qc_from_vc(vc_kt);
   out.mach  = mach_from_qc_and_atm(out.qc, atm);
   out.vt_kt = out.mach * atm.a_kt;
#if PITOT_INCLUDE_EAS
   out.eas_kt = out.vt_kt * atm.sqrt_sigma;
#endif
   return out;
}

inline FromVt from_vt_and_atm(Real vt_kt, const Atm& atm) noexcept
{
   FromVt out;
   out.mach  = vt_kt / atm.a_kt;
   out.qc    = qc_from_mach_and_atm(out.mach, atm);
   out.vc_kt = vc_from_qc(out.qc);
#if PITOT_INCLUDE_EAS
   out.eas_kt = vt_kt * atm.sqrt_sigma;
#endif
   return out;
}

inline FromQc from_qc_and_atm(Real qc, const Atm& atm) noexcept
{
   FromQc out;
   out.mach  = mach_from_qc_and_atm(qc, atm);
   out.vt_kt = out.mach * atm.a_kt;
   out.vc_kt = vc_from_qc(qc);
#if PITOT_INCLUDE_EAS
   out.eas_kt = out.vt_kt * atm.sqrt_sigma;
#endif
   return out;
}

//==========================================================================
// Standard-day wrappers: *_and_hpc_std and *_and_hpc_std_fast
//==========================================================================
inline FromVc from_vc_and_hpc_std(Real vc_kt, Real hpc_ft) noexcept
{
   const Atm atm = standard_day_atm_from_hpc(hpc_ft);
   return from_vc_and_atm(vc_kt, atm);
}

inline FromVt from_vt_and_hpc_std(Real vt_kt, Real hpc_ft) noexcept
{
   const Atm atm = standard_day_atm_from_hpc(hpc_ft);
   return from_vt_and_atm(vt_kt, atm);
}

inline FromQc from_qc_and_hpc_std(Real qc, Real hpc_ft) noexcept
{
   const Atm atm = standard_day_atm_from_hpc(hpc_ft);
   return from_qc_and_atm(qc, atm);
}

inline FromVc from_vc_and_hpc_std_fast(Real vc_kt, Real hpc_ft) noexcept
{
   const Atm atm = standard_day_atm_from_hpc_std_fast(hpc_ft);
   return from_vc_and_atm(vc_kt, atm);
}

inline FromVt from_vt_and_hpc_std_fast(Real vt_kt, Real hpc_ft) noexcept
{
   const Atm atm = standard_day_atm_from_hpc_std_fast(hpc_ft);
   return from_vt_and_atm(vt_kt, atm);
}

inline FromQc from_qc_and_hpc_std_fast(Real qc, Real hpc_ft) noexcept
{
   const Atm atm = standard_day_atm_from_hpc_std_fast(hpc_ft);
   return from_qc_and_atm(qc, atm);
}

// Single-output wrappers (implemented through bundles)
inline Real mach_from_vc_and_hpc_std(Real vc_kt, Real hpc_ft) noexcept
{
   return from_vc_and_hpc_std(vc_kt, hpc_ft).mach;
}

inline Real vt_from_vc_and_hpc_std(Real vc_kt, Real hpc_ft) noexcept
{
   return from_vc_and_hpc_std(vc_kt, hpc_ft).vt_kt;
}

inline Real vc_from_vt_and_hpc_std(Real vt_kt, Real hpc_ft) noexcept
{
   return from_vt_and_hpc_std(vt_kt, hpc_ft).vc_kt;
}

inline Real mach_from_vt_and_hpc_std(Real vt_kt, Real hpc_ft) noexcept
{
   return from_vt_and_hpc_std(vt_kt, hpc_ft).mach;
}

inline Real mach_from_vc_and_hpc_std_fast(Real vc_kt, Real hpc_ft) noexcept
{
   return from_vc_and_hpc_std_fast(vc_kt, hpc_ft).mach;
}

inline Real vt_from_vc_and_hpc_std_fast(Real vc_kt, Real hpc_ft) noexcept
{
   return from_vc_and_hpc_std_fast(vc_kt, hpc_ft).vt_kt;
}

inline Real vc_from_vt_and_hpc_std_fast(Real vt_kt, Real hpc_ft) noexcept
{
   return from_vt_and_hpc_std_fast(vt_kt, hpc_ft).vc_kt;
}

inline Real mach_from_vt_and_hpc_std_fast(Real vt_kt, Real hpc_ft) noexcept
{
   return from_vt_and_hpc_std_fast(vt_kt, hpc_ft).mach;
}

//==========================================================================
// Additional convenience relations you asked about (Tier-1 naming)
//==========================================================================
inline Real vc_from_mach_and_hpc_std(Real M, Real hpc_ft) noexcept
{
   const Atm atm = standard_day_atm_from_hpc(hpc_ft);
   return vc_from_mach_and_atm(M, atm);
}

inline Real vc_from_mach_and_hpc_std_fast(Real M, Real hpc_ft) noexcept
{
   const Atm atm = standard_day_atm_from_hpc_std_fast(hpc_ft);
   return vc_from_mach_and_atm(M, atm);
}

inline Real vt_from_mach_and_hpc_std(Real M, Real hpc_ft) noexcept
{
   const Atm atm = standard_day_atm_from_hpc(hpc_ft);
   return vt_from_mach_and_atm(M, atm);
}

inline Real vt_from_mach_and_hpc_std_fast(Real M, Real hpc_ft) noexcept
{
   const Atm atm = standard_day_atm_from_hpc_std_fast(hpc_ft);
   return vt_from_mach_and_atm(M, atm);
}

// EAS convenience
inline Real eas_from_vt_and_atm(Real vt_kt, const Atm& atm) noexcept
{
   return vt_kt * atm.sqrt_sigma;
}

//==========================================================================
// Geometric <-> Geopotential altitude
//==========================================================================
inline Real h_geopot_from_h_geom(Real z_ft) noexcept
{
   return (earth_radius_ft * z_ft) / (earth_radius_ft + z_ft);
}

inline Real h_geom_from_h_geopot(Real h_ft) noexcept
{
   return (earth_radius_ft * h_ft) / (earth_radius_ft - h_ft);
}

//==========================================================================
// Cleanup macros
//==========================================================================
#undef PITOT_LIKELY
#undef PITOT_UNLIKELY
}  // namespace pitot
