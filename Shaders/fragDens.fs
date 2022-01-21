#version 330 core

in float densAlpha;

out vec4 FragColor;

void main()
{
    FragColor = vec4(1.0f, 1.0f, 1.0f, densAlpha);
}