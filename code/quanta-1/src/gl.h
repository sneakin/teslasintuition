#ifndef _GL_H_
#define _GL_H_

#include <GL/glew.h>

#define GL_GLEXT_PROTOTYPES
#if MAC
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <iostream>
#include <string>

#define ACCESSOR(kind, attr)		       \
  private:				       \
  kind _##attr;				       \
public:					       \
  const kind &attr() const { return _##attr; } \
  void attr(const kind &v) { _##attr = v; }

class GLError
{
  ACCESSOR(std::string, file);
  ACCESSOR(std::string, error_msg);
  ACCESSOR(GLenum, error);
  ACCESSOR(int, line);

 public:
  GLError(const char *file, int line, GLenum err)
    : _file(file), _line(line), _error(err), _error_msg((const char *)gluErrorString(err))
  {
  }
};

#define ASSERT_GL_ERROR { \
    GLenum last_error = GL_NO_ERROR, error = GL_NO_ERROR;           \
    do {                                                            \
      error = glGetError();                                      \
      if(error != GL_NO_ERROR) {                                        \
        std::cerr << "GL error:" << __FILE__ << ":" << __LINE__ << ": " << gluErrorString(error) << std::endl; \
        throw GLError(__FILE__, __LINE__, error);			\
      }                                                                 \
      last_error = error;                                               \
    } while(error != GL_NO_ERROR);                                      \
  }

#endif /* _GL_H_ */
