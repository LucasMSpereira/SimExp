#version 330 core

layout (location = 0) in vec2 aPos;
layout (location = 1) in vec2 disp;

out vec2 dispColor;

uniform int deformScale;
uniform int stepsCurrent;
uniform int stepsTot;
uniform mat4 model;
//uniform mat4 view;
//uniform mat4 projection;

void main()
{
    
    dispColor = disp;
    gl_Position = model*(vec4(deformScale*float(stepsCurrent)/float(stepsTot)*disp, 0.0f, 1.0f) + vec4(aPos, 0.0f, 1.0f));

}