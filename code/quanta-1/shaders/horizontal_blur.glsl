#version 330
uniform sampler2D inColors;
uniform int inBlurSize = 0;
uniform float inBlurFilter[32];// = float[](0.125, 0.25, 0.25, 0.125);

smooth in vec4 Color;
smooth in vec2 Texture;

out vec4 outColor;

vec4 fourPointSample(sampler2D tex, vec2 x, vec2 offset)
{
  return (texture(tex, x + offset * vec2(1.0, 0.0))
          + texture(tex, x - offset * vec2(1.0, 0.0))
          + texture(tex, x + offset * vec2(0.0, 1.0))
          + texture(tex, x - offset * vec2(0.0, 1.0)))
    * 0.25;
}

vec4 hblur(sampler2D tex, vec2 pos)
{
  vec4 pixel;
  vec2 texel = 1.0 / textureSize(tex, 0);

  vec2 vv = inBlurSize * texel;
  vec2 v = vv * vec2(1.0, 0.0);
  
  for(int n = 0; n < inBlurSize; n++) {
    float t = float(n) / float(inBlurSize - 1) - 0.5;
    vec2 x = pos + v * t;
    pixel += inBlurFilter[n] * fourPointSample(tex, x, vv * 0.25);
  }

  return pixel;
}

void main()
{
  outColor = hblur(inColors, Texture);
}
