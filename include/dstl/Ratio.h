#pragma once
#include <cstdint>
#include <limits>
#include <numeric>   // std::gcd (C++17)
#include <type_traits>
#include <sstream>
#include <string>
#include <cmath>
#include <stdexcept>

namespace dstl {

class Ratio {
public:
  using int_t = std::int64_t;

  // ----- ctors -----
  Ratio() : num_(0), den_(1) {}
  ~Ratio() = default;

  Ratio(const Ratio&) = default;
  Ratio(Ratio&&) noexcept = default;
  Ratio& operator=(const Ratio&) = default;
  Ratio& operator=(Ratio&&) noexcept = default;

  // “backend Set(...)” – required by Number<>
  void Set(const double& x)            { from_double(x); }
  void Set(const float& x)             { from_double(static_cast<double>(x)); }
  void Set(const std::int32_t& x)      { num_ = x; den_ = 1; }
  void Set(const std::int64_t& x)      { num_ = x; den_ = 1; }

  // Optional: allow direct exact set
  void Set(const Ratio& r)             { num_ = r.num_; den_ = r.den_; }

  // ----- value access -----
  double Double() const {
    // convert safely
    return static_cast<double>(num_) / static_cast<double>(den_);
  }

  // “backend String(precision, scientific)” – required by Number<>
  std::string String(unsigned precision = 6, bool scientific = false) const {
    if (den_ == 1) {
      std::ostringstream os;
      os.setf(std::ios::fmtflags(0), std::ios::floatfield);
      os << num_;
      return os.str();
    }
    if (!scientific) {
      std::ostringstream os;
      os << num_ << "/" << den_;
      return os.str();
    }
    std::ostringstream os;
    os.setf(scientific ? std::ios::scientific : std::ios::fmtflags(0), std::ios::floatfield);
    os.precision(precision);
    os << Double();
    return os.str();
  }

  // ----- unary -----
  void Neg() { num_ = -num_; }

  // recip()/inv() – required by Number<>
  Ratio recip() const {
    if (num_ == 0) throw std::domain_error("Ratio::recip: divide by zero");
    Ratio r(den_, num_);
    r.normalize_();
    return r;
  }
  void inv() {
    if (num_ == 0) throw std::domain_error("Ratio::inv: divide by zero");
    std::swap(num_, den_);
    normalize_();
  }

  // ----- arithmetic returning new backend (Plus/Minus/Times/Divide) -----
  Ratio Plus(const Ratio& rhs) const {
    return add_(rhs);
  }
  Ratio Minus(const Ratio& rhs) const {
    return add_(Ratio(-rhs.num_, rhs.den_));
  }
  Ratio Times(const Ratio& rhs) const {
    return mul_(rhs);
  }
  Ratio Divide(const Ratio& rhs) const {
    if (rhs.num_ == 0) throw std::domain_error("Ratio::Divide: divide by zero");
    return mul_(Ratio(rhs.den_, rhs.num_));
  }

  // ----- in-place arithmetic (Add/Sub/Mul/Div) -----
  void Add(const Ratio& rhs) { *this = add_(rhs); }
  void Sub(const Ratio& rhs) { *this = add_(Ratio(-rhs.num_, rhs.den_)); }
  void Mul(const Ratio& rhs) { *this = mul_(rhs); }
  void Div(const Ratio& rhs) {
    if (rhs.num_ == 0) throw std::domain_error("Ratio::Div: divide by zero");
    *this = mul_(Ratio(rhs.den_, rhs.num_));
  }

  // ----- convenience ctors for backend internals -----
  Ratio(int_t n, int_t d) : num_(n), den_(d) { normalize_(); }

private:
  int_t num_;
  int_t den_;

  // normalize sign and reduce
  void normalize_() {
    if (den_ == 0) throw std::domain_error("Ratio: zero denominator");
    if (den_ < 0) { den_ = -den_; num_ = -num_; }
    if (num_ == 0) { den_ = 1; return; }
    auto g = std::gcd(num_ < 0 ? -num_ : num_, den_);
    if (g > 1) { num_ /= g; den_ /= g; }
  }

  // safe add with cross-cancellation:
  // a/b + c/d = (a*(d/g2) + c*(b/g1)) / lcm(b,d) where g1=gcd(a,c), g2=gcd(b,d)
  Ratio add_(const Ratio& r) const {
    if (num_ == 0) return r;
    if (r.num_ == 0) return *this;

    auto g1 = std::gcd(num_ < 0 ? -num_ : num_, r.num_ < 0 ? -r.num_ : r.num_);
    auto g2 = std::gcd(den_, r.den_);

    // reduce before multiply to limit overflow
    // a' = a/g1, c' = c/g1, b' = b/g2, d' = d/g2
    __int128 a1 = (__int128)num_ / g1;
    __int128 c1 = (__int128)r.num_ / g1;
    __int128 b1 = (__int128)den_ / g2;
    __int128 d1 = (__int128)r.den_ / g2;

    __int128 n = a1 * d1 + c1 * b1;
    __int128 d = (__int128)den_ / g2 * r.den_; // lcm-like denominator

    // normalize back to int64_t with extra gcd to reduce growth
    int_t N, D;
    to_i64_reduce_(n, d, N, D);
    return Ratio(N, D);
  }

  // safe mul with cross-cancellation:
  // (a/b)*(c/d) = ( (a/g1)*(c/g2) ) / ( (b/g2)*(d/g1) ), where g1=gcd(a,d), g2=gcd(c,b)
  Ratio mul_(const Ratio& r) const {
    if (num_ == 0 || r.num_ == 0) return Ratio(0,1);

    auto g1 = std::gcd(num_ < 0 ? -num_ : num_, r.den_);
    auto g2 = std::gcd(r.num_ < 0 ? -r.num_ : r.num_, den_);

    __int128 a = (__int128)num_ / g1;
    __int128 d = (__int128)r.den_ / g1;
    __int128 c = (__int128)r.num_ / g2;
    __int128 b = (__int128)den_ / g2;

    __int128 n = a * c;
    __int128 m = b * d;

    int_t N, D;
    to_i64_reduce_(n, m, N, D);
    return Ratio(N, D);
  }

  // convert 128-bit frac to int64_t with reduction
  static void to_i64_reduce_(__int128 n128, __int128 d128, int_t& N, int_t& D) {
    if (d128 < 0) { d128 = -d128; n128 = -n128; }
    if (n128 == 0) { N = 0; D = 1; return; }

    auto g = gcd128_(n128 < 0 ? -n128 : n128, d128);
    n128 /= g; d128 /= g;

    // clamp if still out of range (very rare unless pathological)
    const __int128 min64 = (__int128)std::numeric_limits<int_t>::min();
    const __int128 max64 = (__int128)std::numeric_limits<int_t>::max();
    if (n128 < min64 || n128 > max64 || d128 < 1 || d128 > max64) {
      // scale down proportionally (lossy but prevents UB)
      long double val = (long double)n128 / (long double)d128;
      // best-effort bounded denominator
      approx_from_long_double_(val, 1000000000LL, N, D);
      return;
    }
    N = static_cast<int_t>(n128);
    D = static_cast<int_t>(d128);
  }

  static __int128 gcd128_(__int128 a, __int128 b) {
    while (b != 0) { __int128 t = a % b; a = b; b = t; }
    return a < 0 ? -a : a;
  }

  // Continued fraction approximation for Set(double) and clamping fallbacks
  void from_double(double x, std::int64_t max_den = 1000000) {
    if (!std::isfinite(x)) { // NaN/inf → store as 0/1 to avoid UB
      num_ = 0; den_ = 1; return;
    }
    int sign = (x < 0) ? -1 : 1;
    long double v = std::fabs((long double)x);

    // continued fraction expansion
    std::int64_t n0 = 0, d0 = 1;
    std::int64_t n1 = 1, d1 = 0;

    while (true) {
      long double a = std::floor(v);
      long double tmp = n1; n1 = (std::int64_t)(a) * n1 + n0; n0 = (std::int64_t)(tmp);
      tmp = d1; d1 = (std::int64_t)(a) * d1 + d0; d0 = (std::int64_t)(tmp);

      if (d1 > max_den || std::fabs(v - a) < 1e-18L) break;
      v = 1.0L / (v - a);
    }

    num_ = sign * n1; den_ = d1 == 0 ? 1 : d1;
    normalize_();
  }

  static void approx_from_long_double_(long double x, std::int64_t max_den, int_t& N, int_t& D) {
    Ratio tmp; tmp.from_double((double)x, max_den); N = tmp.num_; D = tmp.den_;
  }
};

} // namespace dstl
