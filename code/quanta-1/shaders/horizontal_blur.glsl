#version 330
uniform sampler2D inColors;
uniform int inBlurSize = 4;
uniform float inBlurFilter[16];// = float[](0.125, 0.25, 0.25, 0.125);

smooth in vec4 Color;
smooth in vec2 Texture;

out vec4 outColor;

vec4 hblur(sampler2D tex, vec2 pos)
{
  vec4 pixel;
  vec2 texel = 1.0 / textureSize(tex, 0);
  int size = inBlurSize;
  float stepsize = 1.0 / float(size);
  for(int n = 0; n < size; n++) {
    pixel += inBlurFilter[n] * texture(tex, pos + vec2(float(n) - 0.5 * float(size), 0) * texel);
    //pixel += stepsize * texture(tex, pos + vec2(n - 0.5 * float(size), 0) * texel);
  }

  return pixel;
}

void main()
{
  outColor = hblur(inColors, Texture);
}
