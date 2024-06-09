#version 300 es

layout(location=0) in vec4 position;
layout(location=1) in vec4 color;

uniform float seconds;
out vec4 vColor;

void main() {
	float seed = float(gl_VertexID);
    vColor = color;
    float jitterAmount = 0.06; 
	
	float randomX = (fract(sin( seed + seconds / 4.5) * 20.0));
    float randomY = (fract(cos(	seed + seconds / 4.5) * 20.0));

    
    randomX *= jitterAmount;
    randomY *= jitterAmount;

    vec4 jitteredPosition = position + vec4(randomX, randomY, 0.0, 0.0);

    gl_Position = jitteredPosition;
}

