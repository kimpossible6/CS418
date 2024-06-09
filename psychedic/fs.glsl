#version 300 es
precision highp float;

in vec4 position_copy;

uniform float seconds;

out vec4 fragColor;

void main() {
	vec2 uv = position_copy.xy;
	float color = cos( uv.x * uv.x * sin( seconds / 4.0 ) * 30.0) + sin( uv.x * uv.x * cos( seconds / 4.0 ) * 20.0 ) + cos( uv.y * uv.y * sin( seconds / 4.0 ) * 30.0 ) + sin( uv.y * uv.y * cos( seconds / 4.0 ) * 20.0 );
	
	fragColor = vec4(cos( color + seconds / 2.5 ) * 0.6
	,   sin( color + seconds / 2.5 ) * 0.6, sin( color 
	+ seconds / 1.5 ) * 0.8, 1.0);
}