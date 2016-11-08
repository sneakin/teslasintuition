#version 330
uniform sampler2D tex;
uniform mat4 mModelView, mColor;
uniform vec3 uForward, uCamera;
uniform vec3 uLight0Color, uAmbientLight;
uniform vec4 uLight0Position, uLight0Specular;
uniform float uLight0Intensity;

smooth in vec4 Color;
smooth in vec2 Texture;
smooth in vec3 Normal;
smooth in vec3 Vert;
smooth in float Light0Distance;
smooth in vec3 Light0Position;

out vec4 outColor;

void main()
{
  vec3 normal = normalize(Normal);
  if(!gl_FrontFacing) {
    normal = -normal;
  }
  vec3 light_dir = normalize(Light0Position);
  float d = Light0Distance * Light0Distance;
  float brightness = max(dot(light_dir, normal), 0.0);
  vec3 diffuse = brightness * uLight0Intensity * vec3(uLight0Color) / d;

  vec3 surf_to_cam = -normalize(Vert);
  vec3 r = reflect(-light_dir, normal);
  float spec_angle = max(dot(r, surf_to_cam), 0.0);
  vec3 specular = uLight0Intensity * pow(spec_angle, uLight0Specular.a) * vec3(uLight0Specular) / d;

  vec4 tex_color = texture(tex, Texture);
  vec3 surf_color = mix(vec3(tex_color), vec3(Color), Color.a);
  vec3 c = (uAmbientLight + diffuse + specular) * vec3(surf_color);
  outColor = vec4(c.r, c.g, c.b, tex_color.a);
}
