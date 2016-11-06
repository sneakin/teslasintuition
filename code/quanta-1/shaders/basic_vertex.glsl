#version 330
layout(location = 0) in vec4 inPosition;
layout(location = 1) in vec4 inColor;
layout(location = 2) in vec2 inTexture;

uniform mat4 ModelView, Projection;

void main()
{
  gl_Position = Projection * ModelView * inPosition;
}