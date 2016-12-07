#version 330
uniform sampler2D tex;

uniform mat4 mModelView, mColor;
uniform vec3 uForward, uCamera;
uniform int uFragmentId;

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
smooth in vec4 Vert_screen;
smooth in vec3 Vert_lastframe;
smooth in vec4 Vert_screen_lastframe;
in vec3 Velocity;
smooth in float Light0Distance;
smooth in vec3 Light0Position;

out vec4 outColor;
out vec4 outVelocity;
out vec4 outDepth;
out vec4 outBloom;

vec4 powv(vec4 v, float e)
{
  return vec4(pow(v.x, e),
              pow(v.y, e),
              pow(v.z, e),
              pow(v.w, e));
}

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

  vec3 surf_to_cam = -normalize(vec3(Vert.x, Vert.y, Vert.z));
  vec3 r = reflect(-light_dir, normal);
  float spec_angle = max(dot(r, surf_to_cam), 0.0);
  vec3 base_spec = uMaterial.specular.a * vec3(uMaterial.specular) * vec3(uLight.color) * falloff;
  vec3 specular = max(vec3(0.0, 0.0, 0.0), base_spec * pow(spec_angle, 0.321 * uMaterial.shine));
  vec3 super_specular = base_spec * pow(spec_angle, 0.321 * uMaterial.shine * 1.5);

  vec3 ambient = vec3(uMaterial.ambient) * uAmbientLight;
  
  vec4 tex_color = texture(tex, Texture);
  vec3 surf_color = mix(vec3(tex_color), vec3(Color), Color.a);
  vec3 emission = vec3(uMaterial.emission) * surf_color;
  vec3 c = clamp(ambient + diffuse + specular, 0.0, 1.0) * vec3(surf_color) + super_specular;
  c = max(c, emission);

  outColor = mColor * vec4(c.r, c.g, c.b, 1.0);
  outColor.a = tex_color.a;
  //outVelocity = powv((Vert_screen - Vert_screen_lastframe) * 0.5 + 0.5, 3.0);
  outVelocity = (Vert_screen - Vert_screen_lastframe) * 0.5 + 0.5;
  //outVelocity = (Vert_screen - Vert_screen_lastframe);
  outDepth = vec4(1.0 / Vert_screen.z, float(uFragmentId) / float(10), 0, 0);
  //outBloom = vec4(clamp(specular + emission, 0.0, 1.0), 1.0);
  outBloom = vec4(specular + emission, 1.0);
}
