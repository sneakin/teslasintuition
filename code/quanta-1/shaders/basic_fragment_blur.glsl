#version 330
uniform sampler2D tex;
uniform sampler2D depth;
uniform sampler2D velocity;
uniform sampler2D fragment_data;
uniform sampler2D last_frame;

uniform mat4 mProjection, mColor;
uniform int uFragmentId;

uniform float uFPS = 60.0, uTargetFPS = 60.0;

smooth in vec4 Color;
smooth in vec2 Texture;
smooth in vec4 Normal;
smooth in vec3 Vert;
smooth in vec4 Vert_screen;
in vec3 Velocity;

out vec4 outColor;
out vec3 outVelocity;
out vec4 outDepth;


vec4 powv(vec4 v, float e)
{
  return vec4(pow(v.x, e),
              pow(v.y, e),
              pow(v.z, e),
              pow(v.w, e));
}

const int MAX_SAMPLES = 16;
const int MAX_VELOCITY = 32;

void main()
{
  vec2 texel = 1.0 / vec2(textureSize(tex, 0));
  vec4 vv = texture(velocity, Texture);;
  //vv = powv(vv, 1.0 / 3.0);
  //vv = (vv - 0.5) * 2.0;
  float d = texture(depth, Texture).r;
  float scale = uTargetFPS / uFPS;
  vec2 v = clamp(vec2(vv.x, vv.y) * scale, float(-MAX_VELOCITY) * texel.x, float(MAX_VELOCITY) * texel.x);
  int num_samples = int(clamp(length(v) / texel.x, 0, MAX_SAMPLES));
  int samples = 1;
  outColor = texture(tex, Texture);
  for(int i = 0; i < num_samples; i++) {
    float t = float(i) / float(num_samples - 1);
    vec2 offset = Texture + v * (t - 0.5);
    vec2 coffset = clamp(offset, 0.0, 1.0);
    if(offset == coffset) {
      float nd = texture(depth, coffset).r;
      if(d <= nd) {
        //outColor += mix(texture(last_frame, coffset), texture(tex, coffset), t);
        outColor += texture(tex, coffset);
        samples++;
      }
    }
  }
  outColor /= samples;
  
  outVelocity = Velocity;
  outDepth = vec4(1.0 / Vert_screen.z, float(uFragmentId) / float(10), 0, 0);
}
