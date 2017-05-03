#version 330
uniform sampler2D tex;
//uniform sampler2D depth;
//uniform sampler2D velocity;

uniform mat4 mModelView, mProjection, mColor;
uniform int uFragmentId;

smooth in vec4 Color;
smooth in vec2 Texture;
smooth in vec3 Normal;
smooth in vec3 Vert;
smooth in vec4 Vert_screen;
in vec3 Velocity;
smooth in vec3 Vert_lastframe;
smooth in vec4 Vert_screen_lastframe;

out vec4 outColor;
out vec4 outVelocity;
out vec4 outDepth;

int MAX_INT = 0xFFFFFFF;

void main()
{
  outColor = mColor * (vec4(texture(tex, Texture).rgb, 1.0) * Color);

  outVelocity = (Vert_screen - Vert_screen_lastframe);
  outVelocity = outVelocity * 0.5 + 0.5;
  float far = 200.0, near = -200.0;
  vec4 d = vec4(1.0 / Vert_screen.z, float(uFragmentId) / float(10), 0, 0);
  outDepth = d; //0.5 + vec3(Vert) / (far - near);
}
