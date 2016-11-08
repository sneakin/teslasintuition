#version 330
uniform sampler2D tex;

smooth in vec4 Color;
smooth in vec2 Texture;
smooth in vec4 Normal;

out vec4 outColor;

void main()
{
  outColor = texture(tex, Texture) * Color;
}
