#ifndef _GL_H_
#define _GL_H_

#define GL_GLEXT_PROTOTYPES
#if MAC
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#define ASSERT_GL_ERROR { \
    GLenum error = glGetError(); \
    if(error != GL_NO_ERROR) { \
      std::cerr << "GL error:" << __LINE__ << ": " << gluErrorString(error) << std::endl; \
    } \
  }

#endif /* _GL_H_ */