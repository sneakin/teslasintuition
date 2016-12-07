#version 330
uniform sampler2D inColors;
uniform int inBlurSize = 4;
uniform float inBlurFilter[16];

smooth in vec4 Color;
smooth in vec2 Texture;

out vec4 outColor;

vec4 vblur(sampler2D tex, vec2 pos)
{
  vec4 pixel;
  vec2 texel = 1.0 / textureSize(tex, 0);
  int size = inBlurSize;
  float stepsize = 1.0 / float(size);
  for(int n = 0; n < inBlurSize; n++) {
    pixel += inBlurFilter[n] * texture(tex, pos + vec2(0.0, float(n) - 0.5 * float(size)) * texel);
    //pixel += stepsize * texture(tex, pos + vec2(0.0, n - 0.5 * float(size)) * texel);
  }
  
  return pixel;
}

void main()
{
  outColor = vblur(inColors, Texture);
}
