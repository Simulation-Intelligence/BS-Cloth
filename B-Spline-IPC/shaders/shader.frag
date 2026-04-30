#version 450

layout(location = 0) in vec3 fragColor;
layout(location = 1) in vec3 pos;
layout(location = 2) in vec3 normal;

layout(location = 0) out vec4 outColor;

const vec3 lightPos = vec3(1.0, 1.0, 1.0);

void main()
{
    vec3 n = normalize(normal);
    vec3 l = normalize(pos - lightPos);
    float lNorm = length(lightPos - pos);
    float intensity = abs(dot(n, l)) * (1.0 / (lNorm * lNorm));
    outColor = vec4(fragColor * intensity, 1.0);

    // outColor = vec4(fragColor, 1.0);
}