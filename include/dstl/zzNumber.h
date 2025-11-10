#pragma once
#include <type_traits>
#include <utility>
#include <string>
#include <ostream>
#include <cmath>
#include <limits>
#include <sstream>
#include <iomanip>

// ----------------------------------------
// Simple arithmetic backend for primitive types
// ----------------------------------------
template<typename T>
class ArithmeticBackend {
public:
  ArithmeticBackend() : value_(T{}) {}
  ArithmeticBackend(T val) : value_(val) {}
  
  void Set(T val) { value_ = val; }
  double Double() const { return static_cast<double>(value_); }
  
  std::string String(unsigned precision = 6, bool scientific = false) const {
    std::ostringstream oss;
    if (scientific) {
      oss << std::scientific;
    }
    oss << std::setprecision(precision) << value_;
    return oss.str();
  }
  
  void Add(const ArithmeticBackend& other) { value_ += other.value_; }
  void Sub(const ArithmeticBackend& other) { value_ -= other.value_; }
  void Mul(const ArithmeticBackend& other) { value_ *= other.value_; }
  void Div(const ArithmeticBackend& other) { value_ /= other.value_; }
  void Neg() { value_ = -value_; }
  
  ArithmeticBackend Plus(const ArithmeticBackend& other) const { return ArithmeticBackend(value_ + other.value_); }
  ArithmeticBackend Minus(const ArithmeticBackend& other) const { return ArithmeticBackend(value_ - other.value_); }
  ArithmeticBackend Times(const ArithmeticBackend& other) const { return ArithmeticBackend(value_ * other.value_); }
  ArithmeticBackend Divide(const ArithmeticBackend& other) const { return ArithmeticBackend(value_ / other.value_); }
  
  T value() const { return value_; }
  
private:
  T value_;
};

// ----------------------------------------
// Utilities: void_t for C++17 (and C++14)
// ----------------------------------------
template<class...> using void_t = void;

// ----------------------------------------
// Backend detection (SFINAE)
// ----------------------------------------
namespace nb_detail {
  template<class, class = void> struct has_Set_double   : std::false_type {};
  template<class T> struct has_Set_double<T, void_t<decltype(std::declval<T&>().Set(std::declval<double>()))>> : std::true_type {};

  template<class, class = void> struct has_sin : std::false_type {};
  template<class T> struct has_sin<T, void_t<decltype(std::declval<const T&>().sin())>> : std::true_type {};

  template<class, class = void> struct has_cos : std::false_type {};
  template<class T> struct has_cos<T, void_t<decltype(std::declval<const T&>().cos())>> : std::true_type {};

  template<class, class = void> struct has_sqrt : std::false_type {};
  template<class T> struct has_sqrt<T, void_t<decltype(std::declval<const T&>().sqrt())>> : std::true_type {};

  template<class, class = void> struct has_abs : std::false_type {};
  template<class T> struct has_abs<T, void_t<decltype(std::declval<const T&>().abs())>> : std::true_type {};

  template<class, class = void> struct has_exp : std::false_type {};
  template<class T> struct has_exp<T, void_t<decltype(std::declval<const T&>().exp())>> : std::true_type {};

  template<class, class = void> struct has_log : std::false_type {};
  template<class T> struct has_log<T, void_t<decltype(std::declval<const T&>().log())>> : std::true_type {};

  template<class, class = void> struct has_fma : std::false_type {};
  template<class T> struct has_fma<T, void_t<decltype(std::declval<const T&>().fma(std::declval<const T&>(), std::declval<const T&>()))>> : std::true_type {};

  template<class T> struct is_arith_or_enum : std::integral_constant<bool,
      std::is_arithmetic<T>::value || std::is_enum<T>::value> {};
}

// ----------------------------------------
// Promotion trait (customize as you like)
// Specialize PromoteBackend<A,B> to control mixed ops.
// Default: if same type -> A, else A (or pick a common type).
// ----------------------------------------
template<class A, class B> struct PromoteBackend { using type = A; };
template<class A> struct PromoteBackend<A,A> { using type = A; };
template<class A, class B>
using PromoteBackendT = typename PromoteBackend<A,B>::type;

// ----------------------------------------
// Number wrapper (no concepts; C++17)
// Backend must provide:
//  Set(double/float/int32_t/int64_t), Plus/Minus/Times/Divide, Add/Sub/Mul/Div,
//  Neg(), recip()/inv(), Double(), String(unsigned,bool)
// ----------------------------------------
// Number class
// ----------------------------------------
namespace dstl {

// Type alias for common arithmetic types
template<typename T>
using NumberBackend = std::conditional_t<
  std::is_arithmetic_v<T>,
  ArithmeticBackend<T>,
  T
>;

template<class Backend>
class Number {
public:
  using backend_type = NumberBackend<Backend>;

  // ctors
  Number() = default;
  Number(const Number&) = default;
  Number(Number&&) noexcept = default;

  explicit Number(const backend_type& b) : m_(b) {}
  explicit Number(backend_type&& b) noexcept : m_(std::move(b)) {}

  template<class T, typename std::enable_if<nb_detail::is_arith_or_enum<T>::value, int>::type = 0>
  explicit Number(T v) {
    assign_from_arith(v);
  }

  Number& operator=(const Number&) = default;
  Number& operator=(Number&&) noexcept = default;

  template<class T, typename std::enable_if<nb_detail::is_arith_or_enum<T>::value, int>::type = 0>
  Number& operator=(T v) { assign_from_arith(v); return *this; }

  backend_type&       backend()       noexcept { return m_; }
  const backend_type& backend() const noexcept { return m_; }

  double to_double() const { return m_.Double(); }
  explicit operator double() const { return m_.Double(); }

  std::string string(unsigned precision = 6, bool scientific = false) const {
    return m_.String(precision, scientific);
  }

  // Additional utility methods
  template<class T, typename std::enable_if<nb_detail::is_arith_or_enum<T>::value, int>::type = 0>
  void Add(T value) {
    *this += Number(value);
  }

  // unary
  Number  operator+() const { return *this; }
  Number  operator-() const { Number t(*this); t.m_.Neg(); return t; }

  // compound
  Number& operator+=(const Number& r) { m_.Add(r.m_); return *this; }
  Number& operator-=(const Number& r) { m_.Sub(r.m_); return *this; }
  Number& operator*=(const Number& r) { m_.Mul(r.m_); return *this; }
  Number& operator/=(const Number& r) { m_.Div(r.m_); return *this; }

  template<class T, typename std::enable_if<nb_detail::is_arith_or_enum<T>::value, int>::type = 0>
  Number& operator+=(T v) { return (*this += Number(v)); }
  template<class T, typename std::enable_if<nb_detail::is_arith_or_enum<T>::value, int>::type = 0>
  Number& operator-=(T v) { return (*this -= Number(v)); }
  template<class T, typename std::enable_if<nb_detail::is_arith_or_enum<T>::value, int>::type = 0>
  Number& operator*=(T v) { return (*this *= Number(v)); }
  template<class T, typename std::enable_if<nb_detail::is_arith_or_enum<T>::value, int>::type = 0>
  Number& operator/=(T v) { return (*this /= Number(v)); }

  // binary
  friend Number operator+(const Number& a, const Number& b) { return Number(a.m_.Plus(b.m_)); }
  friend Number operator-(const Number& a, const Number& b) { return Number(a.m_.Minus(b.m_)); }
  friend Number operator*(const Number& a, const Number& b) { return Number(a.m_.Times(b.m_)); }
  friend Number operator/(const Number& a, const Number& b) { return Number(a.m_.Divide(b.m_)); }

  // mixed Number<A> op Number<B> → Number<PromoteBackendT<A,B>>
  template<class B2>
  friend Number<PromoteBackendT<Backend,B2>>
  operator+(const Number& a, const Number<B2>& b) {
    using OutB = PromoteBackendT<Backend,B2>;
    OutB x = OutB(a.m_.Plus(b.backend()));
    return Number<OutB>(std::move(x));
  }
  template<class B2>
  friend Number<PromoteBackendT<Backend,B2>>
  operator-(const Number& a, const Number<B2>& b) {
    using OutB = PromoteBackendT<Backend,B2>;
    OutB x = OutB(a.m_.Minus(b.backend()));
    return Number<OutB>(std::move(x));
  }
  template<class B2>
  friend Number<PromoteBackendT<Backend,B2>>
  operator*(const Number& a, const Number<B2>& b) {
    using OutB = PromoteBackendT<Backend,B2>;
    OutB x = OutB(a.m_.Times(b.backend()));
    return Number<OutB>(std::move(x));
  }
  template<class B2>
  friend Number<PromoteBackendT<Backend,B2>>
  operator/(const Number& a, const Number<B2>& b) {
    using OutB = PromoteBackendT<Backend,B2>;
    OutB x = OutB(a.m_.Divide(b.backend()));
    return Number<OutB>(std::move(x));
  }

  // arithmetic with builtins
  template<class T, typename std::enable_if<nb_detail::is_arith_or_enum<T>::value, int>::type = 0>
  friend Number operator+(const Number& a, T b) { return a + Number(b); }
  template<class T, typename std::enable_if<nb_detail::is_arith_or_enum<T>::value, int>::type = 0>
  friend Number operator+(T a, const Number& b) { return Number(a) + b; }

  template<class T, typename std::enable_if<nb_detail::is_arith_or_enum<T>::value, int>::type = 0>
  friend Number operator-(const Number& a, T b) { return a - Number(b); }
  template<class T, typename std::enable_if<nb_detail::is_arith_or_enum<T>::value, int>::type = 0>
  friend Number operator-(T a, const Number& b) { return Number(a) - b; }

  template<class T, typename std::enable_if<nb_detail::is_arith_or_enum<T>::value, int>::type = 0>
  friend Number operator*(const Number& a, T b) { return a * Number(b); }
  template<class T, typename std::enable_if<nb_detail::is_arith_or_enum<T>::value, int>::type = 0>
  friend Number operator*(T a, const Number& b) { return Number(a) * b; }

  template<class T, typename std::enable_if<nb_detail::is_arith_or_enum<T>::value, int>::type = 0>
  friend Number operator/(const Number& a, T b) { return a / Number(b); }
  template<class T, typename std::enable_if<nb_detail::is_arith_or_enum<T>::value, int>::type = 0>
  friend Number operator/(T a, const Number& b) { return Number(a) / b; }

  // comparisons (by double unless you add a backend Compare)
  friend bool operator==(const Number& a, const Number& b) { return a.to_double() == b.to_double(); }
  friend bool operator!=(const Number& a, const Number& b) { return !(a == b); }
  friend bool operator<(const Number& a, const Number& b)  { return a.to_double() <  b.to_double(); }
  friend bool operator>(const Number& a, const Number& b)  { return b < a; }
  friend bool operator<=(const Number& a, const Number& b) { return !(b < a); }
  friend bool operator>=(const Number& a, const Number& b) { return !(a < b); }

private:
  template<class T>
  void assign_from_arith(T v) {
    // C++17 path with if constexpr (for C++14 replace with overload set)
#if __cplusplus >= 201703L
    if constexpr (std::is_same<T,double>::value) m_.Set(v);
    else if constexpr (std::is_same<T,float>::value) m_.Set(v);
    else if constexpr (std::is_same<T,std::int32_t>::value) m_.Set(v);
    else if constexpr (std::is_same<T,std::int64_t>::value) m_.Set(v);
    else m_.Set(static_cast<double>(v));
#else
    m_.Set(static_cast<double>(v));
#endif
  }

  backend_type m_{};
};

// ostream
template<class B>
inline std::ostream& operator<<(std::ostream& os, const Number<B>& x) {
  return os << x.string();
}

// ----------------------------------------
// ADL math overloads: backend-first, std:: fallback
// ----------------------------------------
template<class B>
inline Number<B> sin(const Number<B>& x) {
  if (nb_detail::has_sin<B>::value) return Number<B>( x.backend().sin() );
  return Number<B>( Number<B>( std::sin(x.to_double()) ).backend() );
}
template<class B>
inline Number<B> cos(const Number<B>& x) {
  if (nb_detail::has_cos<B>::value) return Number<B>( x.backend().cos() );
  return Number<B>( Number<B>( std::cos(x.to_double()) ).backend() );
}
template<class B>
inline Number<B> sqrt(const Number<B>& x) {
  if (nb_detail::has_sqrt<B>::value) return Number<B>( x.backend().sqrt() );
  return Number<B>( Number<B>( std::sqrt(x.to_double()) ).backend() );
}
template<class B>
inline Number<B> abs(const Number<B>& x) {
  if (nb_detail::has_abs<B>::value) return Number<B>( x.backend().abs() );
  return Number<B>( Number<B>( std::fabs(x.to_double()) ).backend() );
}
template<class B>
inline Number<B> exp(const Number<B>& x) {
  if (nb_detail::has_exp<B>::value) return Number<B>( x.backend().exp() );
  return Number<B>( Number<B>( std::exp(x.to_double()) ).backend() );
}
template<class B>
inline Number<B> log(const Number<B>& x) {
  if (nb_detail::has_log<B>::value) return Number<B>( x.backend().log() );
  return Number<B>( Number<B>( std::log(x.to_double()) ).backend() );
}
template<class B>
inline Number<B> fma(const Number<B>& a, const Number<B>& b, const Number<B>& c) {
  if (nb_detail::has_fma<B>::value) return Number<B>( a.backend().fma(b.backend(), c.backend()) );
  return Number<B>( Number<B>( std::fma(a.to_double(), b.to_double(), c.to_double()) ).backend() );
}

} // namespace dstl

// ----------------------------------------
// std::numeric_limits passthrough
// ----------------------------------------
namespace std {
  template<class B>
  struct numeric_limits<dstl::Number<B>> : numeric_limits<double> {};
}
