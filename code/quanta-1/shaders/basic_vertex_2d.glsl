#version 330
layout(location = 0) in vec4 inPosition;
layout(location = 1) in vec4 inColor;
layout(location = 2) in vec4 inTexture;

uniform mat4 ModelView, Projection, TS;

out vec4 Color;
out vec2 Texture;

void main()
{
  Color = inColor;
  Texture = vec2(TS * inTexture);
  gl_Position = Projection * ModelView * inPosition;
}