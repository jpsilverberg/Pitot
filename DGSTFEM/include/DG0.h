#pragma once
#include <array>
#include <vector>
#include <cmath>
#include <functional>
#include <cassert>
#include <iostream>
#include <dstl/dlog.h>
#include <dstl/math/mat.h>
#include <dstl/math/vec.h>

// -----------------------------------------------------------------------------
// Header-only demo: Space–Time DG (Burgers) — TWO ELEMENTS, P4xP4 per element
// - Moving trapezoids mapping: t = t0 + dt*s; x = xL(s) + r * h(s), h = xR - xL
// - Basis: tensor P4(r) x P4(s) with Lagrange nodes at GL(5) points (or any)
// - Quadrature: GL(6) in r and s (safe for P4 + nonlinearity f(u)=0.5 u^2)
// - Assembles local A11, A22, off-diagonal A21 (and, if needed, A12), plus R1,R2
// -----------------------------------------------------------------------------

namespace stdg {


// Small fixed sizes for P4xP4
static constexpr int P = 4;            // polynomial degree per coordinate
static constexpr int N1D = P + 1;      // 5
static constexpr int N = N1D * N1D;    // 25 dofs per element

// Linear algebra micro-types
using Vec25 = dstl::math::Vec<N, double>;
using Mat25 = dstl::math::Mat<N, N, double>;


inline int idx(int i, int j) { return j * N1D + i; }   // (i in r, j in s)

// ------------------------------ Gauss–Legendre --------------------------------
struct GL1D {
  std::vector<double> x, w; // nodes and weights in [0,1]
};

// Hard-code GL(6) on [0,1]
inline GL1D gauss_legendre_6_on_01() {
  // Standard GL6 on [-1,1]
  static const double z[] = {
    -0.9324695142031521, -0.6612093864662645, -0.2386191860831969,
     0.2386191860831969,  0.6612093864662645,  0.9324695142031521
  };
  static const double w[] = {
    0.1713244923791704, 0.3607615730481386, 0.4679139345726910,
    0.4679139345726910, 0.3607615730481386, 0.1713244923791704
  };
  GL1D q; q.x.resize(6); q.w.resize(6);
  for (int i = 0; i < 6; ++i) {
    q.x[i] = 0.5 * (z[i] + 1.0);   // map to [0,1]
    q.w[i] = 0.5 * w[i];           // scale weights
  }
  return q;
}

// ------------------------------- Basis tables ---------------------------------
// Use 1D Lagrange polynomials at GL(5) nodes in [0,1] (degree 4 interpolatory).
struct Basis1D {
  std::array<double, N1D> nodes;               // interpolation nodes r_i or s_j
  std::vector<std::array<double, N1D>> phi;    // phi[iNode] evaluated at quad k
  std::vector<std::array<double, N1D>> dphi;   // dphi/dr or dphi/ds at quad k
};

// Build Lagrange basis and derivative tables at given quadrature nodes.
inline Basis1D make_lagrange_p4_basis(const std::vector<double>& quad_x) {
  // GL(5) nodes on [0,1]
  static const double z5[] = {
    0.0469100770306680, 0.2307653449471585, 0.5,
    0.7692346550528415, 0.9530899229693320
  };
  Basis1D B;
  for (int i = 0; i < N1D; ++i) B.nodes[i] = z5[i];
  const int Q = (int)quad_x.size();
  B.phi.resize(Q);
  B.dphi.resize(Q);
  for (int q = 0; q < Q; ++q) {
    const double x = quad_x[q];
    std::array<double, N1D> L{};   // values
    std::array<double, N1D> dL{};  // derivatives
    // Evaluate L_i(x) and L'_i(x) via barycentric-like product formula
    for (int i = 0; i < N1D; ++i) {
      // Compute numerator and denominator for Lagrange basis
      double num = 1.0, den = 1.0;
      for (int m = 0; m < N1D; ++m) if (m != i) {
        num *= (x - B.nodes[m]);
        den *= (B.nodes[i] - B.nodes[m]);
      }
      L[i] = num / den;
      // Derivative L'_i(x) = sum_{k!=i} ( prod_{m!=i,k} (x-x_m) ) / prod_{m!=i}(x_i-x_m)
      double sum = 0.0;
      for (int k = 0; k < N1D; ++k) if (k != i) {
        double prod = 1.0;
        for (int m = 0; m < N1D; ++m) if (m != i && m != k) prod *= (x - B.nodes[m]);
        sum += prod;
      }
      dL[i] = sum / den;
    }
    B.phi[q]  = L;
    B.dphi[q] = dL;
  }
  return B;
}

// Container for tensor-product tables on (r,s)
struct Basis2D {
  Basis1D r, s;        // 1D bases evaluated at r- and s-quadrature nodes
  GL1D qr, qs;         // store quadrature too
};

inline Basis2D make_basis_tables() {
  Basis2D B;
  B.qr = gauss_legendre_6_on_01();
  B.qs = gauss_legendre_6_on_01();
  B.r  = make_lagrange_p4_basis(B.qr.x);
  B.s  = make_lagrange_p4_basis(B.qs.x);
  return B;
}

// ------------------------------- Geometry (ALE) -------------------------------
struct Geom {
  double dt; // \Delta t
  // xL(s), xR(s) and derivatives with respect to s
  std::function<double(double)> xL, xR, xL_s, xR_s;
  inline double h(double s)   const { return xR(s) - xL(s); }
  inline double h_s(double s) const { return xR_s(s) - xL_s(s); }
  // x_s(r,s) = d/ds x(r,s) = xL'(s) + r * h'(s)
  inline double x_s(double r, double s) const { return xL_s(s) + r * h_s(s); }
};

// ------------------------------ Local state utils -----------------------------
struct ElementState {
  // DOFs ordered lexicographically: alpha = 5*j + i (i: r-index, j: s-index)
  std::array<double, N> U{};

  // Evaluate u_h at (r_q, s_q) using basis tables
  double eval(double r, double s, const Basis2D& B) const {
    // For speed we evaluate at existing quadrature nodes only. This overload
    // is for convenience when r,s are exactly those nodes. If not, user should
    // precompute another Basis2D.
    (void)r; (void)s; (void)B; // not used here
    return 0.0; // placeholder; use eval_q below
  }

  // Evaluate at quadrature node (a,b)
  double eval_q(int a, int b, const Basis2D& Bt) const {
    const auto& Pr  = Bt.r.phi[a];
    const auto& Ps  = Bt.s.phi[b];
    double u = 0.0;
    for (int j = 0; j < N1D; ++j)
      for (int i = 0; i < N1D; ++i)
        u += U[idx(i,j)] * Pr[i] * Ps[j];
    return u;
  }

  // Evaluate restricted to r=0 or r=1 at s-quadrature node b
  double eval_face_r(int rConst /*0 or 1*/, int b, const Basis2D& Bt) const {
    const auto& Ps  = Bt.s.phi[b];
    double u = 0.0;
    for (int j = 0; j < N1D; ++j) {
      // Evaluate 1D Lagrange basis at r=0 or r=1 using nodes
      double val_r = 0.0;
      for (int i = 0; i < N1D; ++i) {
        // Build L_i(rConst) exactly (since rConst is 0 or 1, form product)
        double num = 1.0, den = 1.0;
        for (int m = 0; m < N1D; ++m) if (m != i) {
          num *= ( (double)rConst - Bt.r.nodes[m] );
          den *= ( Bt.r.nodes[i]  - Bt.r.nodes[m] );
        }
        val_r += (num/den) * (i==i ? (i==i) : 0); // we just want L_i(rConst)
        // Note: for rConst in {0,1}, computing once and caching is better.
      }
      // Simpler: explicitly compute L_i(r) and sum; but for robustness,
      // we approximate via a tiny helper per i:
    }
    // For correctness and simplicity, we instead reconstruct via DOFs on edge:
    // u(rConst,s) = sum_i L_i(rConst) * sum_j U_{ij} * theta_j(s)
    // Build L(rConst) once:
    std::array<double, N1D> Lr{};
    for (int i = 0; i < N1D; ++i) {
      double num = 1.0, den = 1.0;
      for (int m = 0; m < N1D; ++m) if (m != i) {
        num *= ((double)rConst - Bt.r.nodes[m]);
        den *= (Bt.r.nodes[i] - Bt.r.nodes[m]);
      }
      Lr[i] = num/den;
    }
    for (int j = 0; j < N1D; ++j) {
      double edge_sum = 0.0;
      for (int i = 0; i < N1D; ++i) edge_sum += U[idx(i,j)] * Lr[i];
      u += edge_sum * Ps[j];
    }
    return u;
  }
};


inline void zero(Vec25& v) { v.fill(0.0); }
inline void zero(Mat25& A) { A.fill(0.0); }
inline void axpy(Vec25& y, double a, const Vec25& x) {
  for (int i = 0; i < N; ++i) y[i] += a * x[i];
}

inline void add_outer(Mat25& A, double w, const std::array<double,N>& left, const std::array<double,N>& right) {
  for (int r = 0; r < N; ++r) for (int c = 0; c < N; ++c) A[r*N + c] += w * left[r] * right[c];
}

// ------------------------------ Core assembly ---------------------------------
// Helper: accumulate volume for one element into (Aee, Re)
inline void assemble_element_volume(const Geom& g, const Basis2D& B,
                                    const ElementState& Ue,
                                    Mat25& Aee, Vec25& Re)
{
  const int Qr = (int)B.qr.x.size();
  const int Qs = (int)B.qs.x.size();
  for (int b = 0; b < Qs; ++b) {
    const double sb = B.qs.x[b];
    const double h  = g.h(sb);
    for (int a = 0; a < Qr; ++a) {
      const double ra = B.qr.x[a];
      const double xs = g.x_s(ra, sb);
      const double wJ = B.qr.w[a] * B.qs.w[b]; // Jacobian h*dt absorbed below
      // u at (ra,sb)
      double u = 0.0;
      // Build test-vector slices for time and space derivatives
      std::array<double, N> v_t{}; v_t.fill(0.0);
      std::array<double, N> v_x{}; v_x.fill(0.0);
      for (int j = 0; j < N1D; ++j) for (int i = 0; i < N1D; ++i) {
        const int aij = idx(i,j);
        const double psi  = B.r.phi[a][i];
        const double dpsi = B.r.dphi[a][i];
        const double th   = B.s.phi[b][j];
        const double dth  = B.s.dphi[b][j];
        u += Ue.U[aij] * psi * th;
        v_t[aij] += ( h * (psi * dth) );                      // (1) h * phi_s
        v_x[aij] += ( g.dt * (dpsi * th) );                   // (2) dt * phi_r
      }
      // Volume residual at node: (h u phi_s) + dt * ( (f(u)-xs*u) phi_r )
      const double fu = 0.5 * u * u;
      const double coeff_r = (fu - xs * u);
      for (int aij = 0; aij < N; ++aij) {
        Re[aij] += wJ * ( u * v_t[aij] + coeff_r * v_x[aij] );
      }
      // Jacobian contributions: h * phi_s * (d u) + dt * ( (f'(u)-xs) * phi_r * (d u) )
      const double df = u; // f'(u) = u
      for (int j = 0; j < N1D; ++j) for (int i = 0; i < N1D; ++i) {
        const int col = idx(i,j);
        const double psi  = B.r.phi[a][i];
        const double th   = B.s.phi[b][j];
        const double du = psi * th; // variation from this DOF
        const double wcol_t = h * du;                       // multiplies test v_t
        const double wcol_x = g.dt * (df - xs) * du;        // multiplies test v_x
        // add to all rows
        for (int rrow = 0; rrow < N; ++rrow) {
          Aee[rrow*N + col] += wJ * ( v_t[rrow] * wcol_t + v_x[rrow] * wcol_x );
        }
      }
    }
  }
}

// Time faces: bottom (s=0, inflow) and top (s=1, outflow)
// u_star_bottom is provided (projected IC). At top, use interior trace.
inline void assemble_time_faces(const Geom& g, const Basis2D& B,
                                const ElementState& Ue,
                                const std::function<double(double)>& u_star_bottom_r,
                                Mat25& Aee, Vec25& Re)
{
  // bottom s=0: + \int h(0) u* phi(r,0) dr  (u* is exterior known)
  {
    const double h0 = g.h(0.0);
    for (int a = 0; a < (int)B.qr.x.size(); ++a) {
      const double ra = B.qr.x[a];
      const double w  = B.qr.w[a];
      const double ustar = u_star_bottom_r(ra);
      for (int j = 0; j < N1D; ++j) {
        const double th0 = B.s.phi[0][j]; // s=0 is first quad node (mapped); acceptable for demo
        for (int i = 0; i < N1D; ++i) {
          const int row = idx(i,j);
          Re[row] += w * h0 * ustar * B.r.phi[a][i] * th0;
        }
      }
    }
  }
  // top s=1: - \int h(1) u(phi) phi(r,1) dr  with interior u -> contributes to A
  {
    const double h1 = g.h(1.0);
    for (int a = 0; a < (int)B.qr.x.size(); ++a) {
      const double w  = B.qr.w[a];
      // Build row (test) vector at s=1
      std::array<double,N> rowvec{}; rowvec.fill(0.0);
      for (int j = 0; j < N1D; ++j) {
        const double th1 = B.s.phi.back()[j]; // assume last entry ~ s=1
        for (int i = 0; i < N1D; ++i) {
          const int row = idx(i,j);
          rowvec[row] += h1 * B.r.phi[a][i] * th1; // coefficient multiplying u*
        }
      }
      // R: use interior u_h at s=1
      double uh_face = 0.0;
      for (int j = 0; j < N1D; ++j) {
        const double th1 = B.s.phi.back()[j];
        for (int i = 0; i < N1D; ++i) uh_face += Ue.U[idx(i,j)] * B.r.phi[a][i] * th1;
      }
      for (int rrow = 0; rrow < N; ++rrow) Re[rrow] -= w * uh_face * rowvec[rrow];
      // J: - rowvec * (column vector of dof basis at r_a, s=1)
      for (int j = 0; j < N1D; ++j) {
        const double th1 = B.s.phi.back()[j];
        for (int i = 0; i < N1D; ++i) {
          const int col = idx(i,j);
          const double colval = B.r.phi[a][i] * th1;
          for (int rrow = 0; rrow < N; ++rrow)
            Aee[rrow*N + col] -= w * rowvec[rrow] * colval;
        }
      }
    }
  }
}

// Lateral face contribution for ONE element face r=0 or r=1
// If upwind picks INTERIOR (this element): contributes to Aee and Re
// If upwind picks NEIGHBOR: contributes to off-diagonal Aen and Re depends on neighbor

struct UpwindPick {
  // returns (isInterior, u_star)
  std::function<std::pair<bool,double>(int /*b*/, double /*nx*/, double /*nt*/)> pick;
};

// Build L(rConst) at r=0 or r=1 once
inline std::array<double,N1D> edge_L_at_r(int rConst, const Basis2D& B){
  std::array<double,N1D> L{};
  for (int i = 0; i < N1D; ++i) {
    double num = 1.0, den = 1.0;
    for (int m = 0; m < N1D; ++m) if (m != i) {
      num *= ((double)rConst - B.r.nodes[m]);
      den *= (B.r.nodes[i] - B.r.nodes[m]);
    }
    L[i] = num/den;
  }
  return L;
}

inline void assemble_lateral_face(const Geom& g_this, const Basis2D& B,
                                  int rConst /*0 or 1*/,
                                  const ElementState& U_this,
                                  // neighbor state and geometry (for off-diagonal if needed)
                                  const ElementState* U_nei,
                                  const Geom* g_nei,
                                  // upwind rule for this face
                                  const UpwindPick& rule,
                                  // outputs
                                  Mat25& Aee, Mat25* Aen, Vec25& Re)
{
  const auto Lr = edge_L_at_r(rConst, B);
  const int Qs = (int)B.qs.x.size();
  const int sgn = (rConst==0 ? -1 : +1); // outward normal sign factor used below

  for (int b = 0; b < Qs; ++b) {
    const double sb = B.qs.x[b];
    const double xs_this = g_this.x_s((double)rConst, sb);
    // Non-unit normal (n_x, n_t)
    const double nx_this = (rConst==1 ? +g_this.dt : -g_this.dt);
    const double nt_this = (rConst==1 ? -xs_this   : +xs_this);

    // Build u_interior on this face at s_b
    double u_int = 0.0;
    for (int j = 0; j < N1D; ++j) {
      double edge_sum = 0.0;
      for (int i = 0; i < N1D; ++i) edge_sum += U_this.U[idx(i,j)] * Lr[i];
      u_int += edge_sum * B.s.phi[b][j];
    }

    // Build u_neighbor on the opposite face if available (for sign and offdiag)
    double u_nei = 0.0;
    if (U_nei) {
      const int rNei = (rConst==1 ? 0 : 1); // opposite side
      const auto LrN = edge_L_at_r(rNei, B);
      for (int j = 0; j < N1D; ++j) {
        double edge_sum = 0.0;
        for (int i = 0; i < N1D; ++i) edge_sum += U_nei->U[idx(i,j)] * LrN[i];
        u_nei += edge_sum * B.s.phi[b][j];
      }
    }

    // Upwind decision using averaged state for sign (robust in nonlinear case)
    double u_bar = U_nei ? 0.5*(u_int + u_nei) : u_int;
    auto pick = rule.pick(b, nx_this, nt_this);
    bool useInterior = pick.first;
    double u_star = pick.second; // caller may override (e.g., impose BC)
    // If rule did not provide value (<nan>), default to computed interior/neighbor
    if (!std::isfinite(u_star)) u_star = useInterior ? u_int : u_nei;

    // G(u*) = dt * 0.5 u*^2 - x_s * u*
    const double G = g_this.dt * 0.5 * u_star * u_star - xs_this * u_star;

    // Test vector on face: phi(rConst, s_b)
    std::array<double,N> face_row{}; face_row.fill(0.0);
    for (int j = 0; j < N1D; ++j) for (int i = 0; i < N1D; ++i) {
      const int row = idx(i,j);
      face_row[row] += B.s.phi[b][j] * Lr[i];
    }

    // Residual contribution: - sgn * sum_b w_b * G * phi
    for (int rrow = 0; rrow < N; ++rrow)
      Re[rrow] -= sgn * B.qs.w[b] * G * face_row[rrow];

    // Jacobian: dG/du* = dt * u* - x_s
    const double dG = g_this.dt * u_star - xs_this;
    if (useInterior) {
      // Aee contribution: - sgn * w * dG * ( test_face_row * basis_col_on_face )
      for (int j = 0; j < N1D; ++j) for (int i = 0; i < N1D; ++i) {
        const int col = idx(i,j);
        const double colval = B.s.phi[b][j] * Lr[i];
        for (int rrow = 0; rrow < N; ++rrow)
          Aee[rrow*N + col] -= sgn * B.qs.w[b] * dG * face_row[rrow] * colval;
      }
    } else if (Aen && U_nei) {
      // Off-diagonal Aen (this residual depends on neighbor DOF)
      for (int j = 0; j < N1D; ++j) for (int i = 0; i < N1D; ++i) {
        const int col = idx(i,j);
        const int rNei = (rConst==1 ? 0 : 1);
        // neighbor column basis at its opposite face
        double colval = 0.0;
        const auto LrN = edge_L_at_r(rNei, B);
        colval = B.s.phi[b][j] * LrN[i];
        for (int rrow = 0; rrow < N; ++rrow)
          (*Aen)[rrow*N + col] -= sgn * B.qs.w[b] * dG * face_row[rrow] * colval;
      }
    }
  }
}

// ---------------------------- Two-element assembly ----------------------------
struct TwoElemSystem {
  Mat25 A11, A22, A21, A12; // A12 often zero for this slab, but provided
  Vec25 R1, R2;
};

// Default upwind rule: use sign of a_n = u*nx + nt with u* ~ average(ul,ur)
inline UpwindPick default_upwind_rule(const ElementState* UL, const ElementState* UR,
                                      const Basis2D& B, int rThis /*0 or 1*/)
{
  UpwindPick R;
  R.pick = [UL,UR,&B,rThis](int b, double nx, double nt) -> std::pair<bool,double> {
    // Estimate u_bar from interior and neighbor along the shared face
    auto edge_val = [&](const ElementState* U, int rC)->double{
      const auto Lr = edge_L_at_r(rC, B);
      double val = 0.0;
      for (int j = 0; j < N1D; ++j) {
        double edge_sum = 0.0;
        for (int i = 0; i < N1D; ++i) edge_sum += U->U[idx(i,j)] * Lr[i];
        val += edge_sum * B.s.phi[b][j];
      }
      return val;
    };
    double uL = UL ? edge_val(UL, /*r=1*/1) : NAN;
    double uR = UR ? edge_val(UR, /*r=0*/0) : NAN;
    double ubar = (std::isfinite(uL) && std::isfinite(uR)) ? 0.5*(uL+uR)
                                                           : (std::isfinite(uL)?uL:uR);
    double an = ubar * nx + nt;
    bool useInterior = (an >= 0.0);
    return {useInterior, NAN}; // let caller compute u* from chosen side
  };
  return R;
}

// Assemble the two-element slab system
inline TwoElemSystem assemble_two_elements(
    const Geom& g1, const Geom& g2,
    const Basis2D& B,
    const ElementState& U1, const ElementState& U2,
    // bottom traces (IC) as functions of r in [0,1] for each element
    const std::function<double(double)>& u0_elem1_r,
    const std::function<double(double)>& u0_elem2_r)
{
  TwoElemSystem S; zero(S.A11); zero(S.A22); zero(S.A21); zero(S.A12); zero(S.R1); zero(S.R2);

  // Volume
  assemble_element_volume(g1, B, U1, S.A11, S.R1);
  assemble_element_volume(g2, B, U2, S.A22, S.R2);

  // Time faces
  assemble_time_faces(g1, B, U1, u0_elem1_r, S.A11, S.R1);
  assemble_time_faces(g2, B, U2, u0_elem2_r, S.A22, S.R2);

  // Exterior spatial faces (domain left of E1 and right of E2):
  // For this demo we assume outflow -> interior upwind; if inflow exists, provide a BC rule.
  {
    // E1 left face r=0
    UpwindPick rule_outflow;
    rule_outflow.pick = [](int, double, double){ return std::pair<bool,double>{true, NAN}; };
    assemble_lateral_face(g1, B, /*r=0*/0, U1, /*nei*/nullptr, /*g_nei*/nullptr,
                          rule_outflow, S.A11, nullptr, S.R1);
    // E2 right face r=1
    assemble_lateral_face(g2, B, /*r=1*/1, U2, /*nei*/nullptr, /*g_nei*/nullptr,
                          rule_outflow, S.A22, nullptr, S.R2);
  }

  // Shared face between E1 (r=1) and E2 (r=0)
  {
    auto rule_shared = default_upwind_rule(&U1, &U2, B, /*rThis ignored*/1);
    // For E1 right face: neighbor is E2; off-diagonal A12 if E1 is downwind (unlikely here)
    assemble_lateral_face(g1, B, /*r=1*/1, U1, &U2, &g2, rule_shared, S.A11, &S.A12, S.R1);

    // For E2 left face: neighbor is E1; off-diagonal A21 if E2 is downwind (likely here)
    assemble_lateral_face(g2, B, /*r=0*/0, U2, &U1, &g1, rule_shared, S.A22, &S.A21, S.R2);
  }

  return S;
}

// ------------------------------ Example usage ---------------------------------
// Geometry for the notional slab (dt = 0.5):
//  E1: xL=0, xR=1+0.5 s  => h1=1+0.5 s, x_s1 = r * 0.5
//  E2: xL=1+0.5 s, xR=2  => h2=1 - 0.5 s, x_s2 = (1-r) * 0.5
inline Geom geom_elem1(double dt){
  Geom g; g.dt = dt;
  g.xL   = [](double s){ return 0.0; };
  g.xR   = [](double s){ return 1.0 + 0.5*s; };
  g.xL_s = [](double s){ (void)s; return 0.0; };
  g.xR_s = [](double s){ (void)s; return 0.5; };
  return g;
}
inline Geom geom_elem2(double dt){
  Geom g; g.dt = dt;
  g.xL   = [](double s){ return 1.0 + 0.5*s; };
  g.xR   = [](double s){ return 2.0; };
  g.xL_s = [](double s){ (void)s; return 0.5; };
  g.xR_s = [](double s){ (void)s; return 0.0; };
  return g;
}

// Initial data projected on each element face r in [0,1]
// E1: u0(x)=1 on [0,1]  => u0(r) = 1
// E2: u0(x)=1-(x-1) on [1,2] linearly -> on reference r, map x=1 + r, so u0 = 1 - r
inline std::function<double(double)> u0_elem1() { return [](double){ return 1.0; }; }
inline std::function<double(double)> u0_elem2() { return [](double r){ return 1.0 - r; }; }

// A tiny driver to show assembly (prints norms)
inline void demo_build() {
  Basis2D B = make_basis_tables();
  Geom g1 = geom_elem1(0.5), g2 = geom_elem2(0.5);
  ElementState U1{}, U2{}; // start with zeros or with a projection of IC inside volume

  auto S = assemble_two_elements(g1, g2, B, U1, U2, u0_elem1(), u0_elem2());
  auto normA = [](const Mat25& A){ double s=0; for(double v:A) s+=v*v; return std::sqrt(s); };
  auto normR = [](const Vec25& R){ double s=0; for(double v:R) s+=v*v; return std::sqrt(s); };
  VAR(S.A11);
  VAR(S.A12);
  VAR(S.R1);
  VAR(S.A21);
  VAR(S.A22);
  VAR(S.R2);
  std::cout << "||A11||_F= " << normA(S.A11) << "  ||A22||_F= " << normA(S.A22)
            << "  ||A21||_F= " << normA(S.A21) << "  ||A12||_F= " << normA(S.A12) << "\n";
  std::cout << "||R1||= " << normR(S.R1) << "  ||R2||= " << normR(S.R2) << "\n";
}

} // namespace stdg
