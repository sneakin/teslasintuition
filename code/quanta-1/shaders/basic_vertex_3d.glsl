#version 330
layout(location = 0) in vec4 inPosition;
layout(location = 1) in vec4 inColor;
layout(location = 2) in vec4 inNormal;
layout(location = 3) in vec4 inTexture;

uniform mat4 mModelView, mProjection, mTexture, mColor;
uniform vec4 uColor = vec4(1.0, 1.0, 1.0, 0.0);

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

out vec4 Color;
out vec2 Texture;
out vec3 Normal;
out vec3 Vert;
out float Light0Distance;
out vec3 Light0Position;

void main()
{
  Color = mColor * uColor; //mix(inColor, uColor, uColor.a);
  //Normal = mat3(mModelView) * normalize(vec3(inNormal));
  Normal = normalize(mat3(mModelView) * vec3(inNormal));
  vec4 pos = vec4(inPosition.x, inPosition.y, inPosition.z, 1.0);
  Vert = vec3(mModelView * pos);
  Texture = vec2(mTexture * inTexture);
  Light0Position = vec3(uLight.position) - Vert;
  Light0Distance = length(Light0Position);
  gl_Position = mProjection * mModelView * pos;
}
