#include <dstl/dlog.h>
#include <dstl/dmutex.h>
#include <dstl/BitFlag.h>
#include <dstl/Numbers.h>
#include <dstl/Quadrature.h>

enum class Flags : unsigned
{
  None = 0,
  One = 1u << 0
};

int main()
{
  dstl::BitFlag<Flags> flag;
  flag.SetFlag(Flags::One);

  dstl::MultiThreadMutex mutex;
  {
    dstl::LockGuard<dstl::MultiThreadMutex> guard(mutex);
    (void)guard;
  }

  auto level = dstl::LogLevel::INFO;

  dstl::numbers::Float64 scalar{1.0};
  auto params = dstl::quadrature::AdaptParams{};
  (void)scalar;

  const bool ok = flag.HasFlag(Flags::One) && level == dstl::LogLevel::INFO && params.max_depth == 20;
  return ok ? 0 : 1;
}
