#ifdef MAC
#include <functional>

namespace std
{
  void __throw_bad_function_call()
  {
    throw bad_function_call();
  }
}
#endif /* MAC */