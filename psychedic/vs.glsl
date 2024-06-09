#version 300 es

layout(location=0) in vec4 position;

out vec4 position_copy;

void main() {

  gl_Position = gl_VertexID == 0 || gl_VertexID == 4 ? vec4(1,1,0,1) : (gl_VertexID == 1 ? vec4(-1,1,0,1) : (gl_VertexID == 2 || gl_VertexID == 5 ? vec4(-1,-1,0,1) : vec4(1,-1,0,1)));
  position_copy = gl_Position;
}
