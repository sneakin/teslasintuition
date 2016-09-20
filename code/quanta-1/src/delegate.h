#define DELEGATE(type, method, member) \
  type method() const { return member.method(); }
