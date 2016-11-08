#version 330
layout(location = 0) in vec4 inPosition;
layout(location = 1) in vec4 inColor;
layout(location = 2) in vec4 inNormal;
layout(location = 3) in vec4 inTexture;

uniform mat4 mModelView, mProjection, mTexture, mColor;
uniform vec4 uColor;

out vec4 Color;
out vec2 Texture;
out vec3 Normal;
out vec3 Vert;

void main()
{
  Color = mColor * (uColor.a * uColor + (1.0 - uColor.a) * inColor);
  Normal = vec3(inNormal);
  Vert = vec3(inPosition);
  Texture = vec2(mTexture * inTexture);
  gl_Position = mProjection * mModelView * inPosition;
}