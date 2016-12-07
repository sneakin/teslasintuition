#version 330
uniform sampler2D inColors;
uniform sampler2D inBloom;
uniform sampler2D inBloomSamples;

uniform int inBloomNumSamples = 0;
uniform float inBloomBrightness = 1.0;

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

vec4 downsample(sampler2D tex, vec2 pos, int samples, int start = 0)
{
  vec4 s;
  for(int i = start; i < samples; i++) {
    float scale = pow(0.5, i + 1);
    s += texture(tex, vec2(1.0 - scale * 2.0, 1.0 - scale) + pos * scale);
  }
  s /= float(samples);

  return s;
}

void main()
{
  vec4 bloom;// = texture(inBloom, Texture); //blur(inBloom, Texture, size);
  bloom = downsample(inBloomSamples, Texture, inBloomNumSamples, 0);
  bloom *= inBloomBrightness;
  outColor = texture(inColors, Texture) + bloom;
}
