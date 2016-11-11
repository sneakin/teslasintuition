#version 330
uniform sampler2D tex;
uniform mat4 mModelView, mColor;
uniform vec3 uForward, uCamera;

struct Light {
  vec3 color;
  vec4 position;
  vec4 specular;
  float intensity;
};

uniform Light uLight = {
  vec3(1.0, 1.0, 1.0),
  vec4(0.0, 0.0, 1.0, 0.0),
  vec4(1.0, 1.0, 1.0, 0.5),
  1.0
};
uniform vec3 uAmbientLight = vec3(0.1, 0.1, 0.1);

struct Material {
  float shine;
  vec4 specular;
  vec4 emission;
  vec4 diffuse;
  vec4 ambient;
};
uniform Material uMaterial = {
  48.0,
  vec4(1.0, 1.0, 1.0, 1.0),
  vec4(1.0, 1.0, 1.0, 1.0),
  vec4(1.0, 1.0, 1.0, 1.0),
  vec4(1.0, 1.0, 1.0, 1.0)
};

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
  float d = Light0Distance;
  float falloff = uLight.intensity / d;
  float brightness = max(dot(light_dir, normal), 0.0);
  vec3 diffuse = vec3(uMaterial.diffuse) * brightness * vec3(uLight.color) * falloff;

  vec3 surf_to_cam = -normalize(Vert);
  vec3 r = reflect(-light_dir, normal);
  float spec_angle = max(dot(r, surf_to_cam), 0.0);
  vec3 base_spec = uMaterial.specular.a * vec3(uMaterial.specular) * vec3(uLight.color) * falloff;
  vec3 specular = base_spec * pow(spec_angle, 0.321 * uMaterial.shine);
  vec3 super_specular = base_spec * pow(spec_angle, 0.321 * uMaterial.shine * 1.5);

  vec3 ambient = vec3(uMaterial.ambient) * uAmbientLight;
  
  vec4 tex_color = texture(tex, Texture);
  vec3 surf_color = mix(vec3(tex_color), vec3(Color), Color.a);
  vec3 emission = vec3(uMaterial.emission) * surf_color;
  vec3 c = clamp(ambient + diffuse + specular, 0.0, 1.0) * vec3(surf_color) + super_specular;
  c = max(c, emission);
  outColor = vec4(c.r, c.g, c.b, tex_color.a);
}
