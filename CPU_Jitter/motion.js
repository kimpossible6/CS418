function compileShader(vs_source, fs_source) {
    const vs = gl.createShader(gl.VERTEX_SHADER)
    gl.shaderSource(vs, vs_source)
    gl.compileShader(vs)
    if (!gl.getShaderParameter(vs, gl.COMPILE_STATUS)) {
        console.error(gl.getShaderInfoLog(vs))
        throw Error("Vertex shader compilation failed")
    }

    const fs = gl.createShader(gl.FRAGMENT_SHADER)
    gl.shaderSource(fs, fs_source)
    gl.compileShader(fs)
    if (!gl.getShaderParameter(fs, gl.COMPILE_STATUS)) {
        console.error(gl.getShaderInfoLog(fs))
        throw Error("Fragment shader compilation failed")
    }

    const program = gl.createProgram()
    gl.attachShader(program, vs)
    gl.attachShader(program, fs)
    gl.linkProgram(program)
    if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
        console.error(gl.getProgramInfoLog(program))
        throw Error("Linking failed")
    }

    // loop through all uniforms in the shader source code
    // get their locations and store them in the GLSL program object for later use
    const uniforms = {}
    for(let i=0; i<gl.getProgramParameter(program, gl.ACTIVE_UNIFORMS); i+=1) {
        let info = gl.getActiveUniform(program, i)
        uniforms[info.name] = gl.getUniformLocation(program, info.name)
    }
    program.uniforms = uniforms

    return program
}

//cpu jittering randomize the position
function createRandomnessToData(randomData) {
    for(let i=0; i<randomData.attributes[0].length; i+=1) {
        randomData.attributes[0][i] = [randomData.attributes[0][i][0] + ((Math.random() * 0.004) - 0.002), randomData.attributes[0][i][1] + ((Math.random() * 0.004) - 0.002)];
    }
    return randomData;
}

let globalBuffer;
let colorBuffer;
let f32;
let randomData;
let data;
function setupGeomery(geom) {
    var triangleArray = gl.createVertexArray()
    gl.bindVertexArray(triangleArray)

    // for(let i=0; i<geom.attributes.length; i+=1) {
        globalBuffer = gl.createBuffer()
        colorBuffer = gl.createBuffer()
        gl.bindBuffer(gl.ARRAY_BUFFER, globalBuffer)
        f32 = new Float32Array(geom.attributes[0].flat())
        gl.bufferData(gl.ARRAY_BUFFER, f32, gl.DYNAMIC_DRAW)
        gl.vertexAttribPointer(0, geom.attributes[0][0].length, gl.FLOAT, false, 0, 0)
        gl.enableVertexAttribArray(0)

        gl.bindBuffer(gl.ARRAY_BUFFER, colorBuffer)
        let colorf32 = new Float32Array(geom.attributes[1].flat())
        gl.bufferData(gl.ARRAY_BUFFER, colorf32, gl.DYNAMIC_DRAW)
        gl.vertexAttribPointer(1, geom.attributes[1][0].length, gl.FLOAT, false, 0, 0)
        gl.enableVertexAttribArray(1)

    // }

    var indices = new Uint16Array(geom.triangles.flat())
    var indexBuffer = gl.createBuffer()
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, indexBuffer)
    gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, indices, gl.STATIC_DRAW)

    return {
        mode: gl.TRIANGLES,
        count: indices.length,
        type: gl.UNSIGNED_SHORT,
        vao: triangleArray
    }
}

function draw(milliseconds) {
    gl.clear(gl.COLOR_BUFFER_BIT) 
    gl.useProgram(program)
    // values that do not vary between vertexes or fragments are called "uniforms"
    gl.uniform1f(program.uniforms.seconds, milliseconds/1000);
    gl.bindVertexArray(geom.vao)
    gl.drawElements(geom.mode, geom.count, geom.type, 0)
}

function tick(milliseconds) {
    draw(milliseconds)
    randomData = createRandomnessToData(data)
    gl.bindBuffer(gl.ARRAY_BUFFER, globalBuffer)
    f32 = new Float32Array(randomData.attributes[0].flat())
    gl.bufferData(gl.ARRAY_BUFFER, f32, gl.DYNAMIC_DRAW)
    requestAnimationFrame(tick) // asks browser to call tick before next frame
}

window.addEventListener('load', async (event) => {
    window.gl = document.querySelector('canvas').getContext('webgl2')
    let vs = await fetch('vertex.glsl').then(res => res.text())
    let fs = await fetch('fragment.glsl').then(res => res.text())
    window.program = compileShader(vs,fs)
    data = await fetch('geometry.json').then(r=>r.json())
    window.geom = setupGeomery(data)
    requestAnimationFrame(tick) // asks browser to call tick before first frame
})
