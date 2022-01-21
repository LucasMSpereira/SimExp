#version 330 core

layout (location = 0) in vec2 aPos;
layout (location = 1) in float dens;

out float densAlpha;

uniform mat4 model;
//uniform mat4 view;
//uniform mat4 projection;

void main()
{
    
    densAlpha = dens;
    gl_Position = model*vec4(aPos, 0.0f, 1.0f);

}