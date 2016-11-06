#version 330
layout(location = 0) in vec4 inPosition;
layout(location = 1) in vec4 inColor;
layout(location = 2) in vec4 inTexture;

uniform mat4 mModelView, mProjection, mTexture, mColor;

out vec4 Color;
out vec2 Texture;

void main()
{
  Color = mColor * inColor;
  Texture = vec2(mTexture * inTexture);
  gl_Position = mProjection * mModelView * inPosition;
}