#version 330 core

uniform float maxDisp;
uniform int stepsCurrent;
uniform int stepsTot;

in vec2 dispColor;

out vec4 FragColor;

void main()
{
    float totalDisp;
    float percent = float(stepsCurrent)/float(stepsTot);
    totalDisp = sqrt(pow(dispColor.x, 2) + pow(dispColor.y, 2));
    FragColor = vec4(1.0f, 1.0f - percent*(totalDisp/maxDisp), 1.0f - percent*(totalDisp/maxDisp), 1.0f);
}