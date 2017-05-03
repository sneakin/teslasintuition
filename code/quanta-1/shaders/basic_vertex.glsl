#version 330
layout(location = 0) in vec4 inPosition;
layout(location = 1) in vec4 inColor;
layout(location = 2) in vec4 inNormal;
layout(location = 3) in vec4 inTexture;

uniform mat4 mModelView, mProjection, mTexture, mColor;
uniform vec4 uColor = vec4(1.0, 1.0, 1.0, 0.0);
uniform vec3 uVelocity = vec3(0.0, 0.0, 0.0);

out vec4 Color;
out vec2 Texture;
out vec3 Normal;
out vec3 Vert;
out vec4 Vert_screen;
out vec3 Velocity;
out vec3 Vert_Lastframe;
out vec4 Vert_screen_lastframe;

void main()
{
  Color = mColor * mix(inColor, uColor, uColor.a);
  Normal = mat3(mModelView) * vec3(inNormal);
  Vert = vec3(mModelView * inPosition);
  Vert_Lastframe = vec3(mModelView * inPosition - vec4(uVelocity, 0.0));
  Vert_screen = mProjection * mModelView * inPosition;
  Vert_screen_lastframe = mProjection * (mModelView * inPosition - vec4(uVelocity, 0.0));
  Velocity = uVelocity;
  Texture = vec2(mTexture * inTexture);
  gl_Position = Vert_screen;
}
