"use strict";

define(function (){

    var exports = {};

    /**
        @param test - property value to test, may be undefined
        @param control - must be either string, array, or object
    */
    exports.validate_property = function (test, control)  {

        if (typeof test === 'undefined') {

            // nothing will happen if control is an array, and the first element is undefined.
            if (!Array.isArray(control) || typeof control[0] !== 'undefined') {
                // all errors will be caught and redefined so just throw nothing
                throw new Error(" <- Non-optional property is undefined!");
            }
        }else{
            /*
                The test property is defined, so now we just have to see if it's defined properly.

                The most likely ideal scenario is the type of property is simply equal
                to the control. But it could fail for different reasons.
            */
            var testType = typeof test;

            if (testType !== control) {
                if (typeof control === 'object') {
                    if (Array.isArray(control)) {

                        /*
                            Ok, so the control property is an array, which means that
                            the test property has to pass only one of the tests.
                            Each element can either be undefined, a string, or an object.
                        */

                        var i = 0;
                        var failed = true;

                        while(failed && i < control.length) {
                            try {
                                if (typeof control[i] !== 'undefined') {
                                    exports.validate_property(test, control[i]);
                                    failed = false;
                                }
                            }finally{
                                i++;
                            }
                        }

                        if (failed) {
                            throw new Error(" <- Property does not match any given possible types!");
                        }

                    }else if (testType === 'object'){
                        /*
                            The control is a non-array object, which means the test property
                            must also be an object. The is a call for recursion.
                        */

                            exports.validate_object(test, control);
                    }else{
                        throw new Error(" <- Property does not match any given possible types!");
                    }
                }else{
                    throw new Error(" <- Property does not match any given possible types!");
                }
            }
        }
    }

    /**
        Validates a spec input according to some form object (which can be recursive)

        The form works like this. For each property the form object has, it checks the
        spec object to make sure it has it. Also, that the type of the property in
        the spec object is what is listed in the from property.

        E.G.

        control = {
            property1 : 'string'
        };

        This says spec must have a property called 'property1' and that
        typeof spec.property1 === 'string' must be true.

        A property can be optional if the form property is an array where the
        first element is left empty

        control = {
            property1 : [,'string']
        };

        This means spec.property1 could be undefined, or be a 'string'.

        control = {
            property1 : ['number','string']
        };

        As expected, more than one type can be specified to broaden the options.

        A control property can also itself be an object, in order to specify sub-properties

        control = {
            property1 : {
                nestedProp : 'number'
            }
        };

        This requires typeof spec.property1.nestedProp === 'number' be true.

        @param test - object to validate
        @param control - object with properties test object must have
     */
    exports.validate_object = function(f, test, control) {

        for(var prop in control) {
            try {
                exports.validate_property(test[prop], control[prop]);
            }catch(error) {
                throw new Error("." + prop + error.message)
            }
        }
    };


    /**
     * This is a wrapper for the webgl interface.
     */
    exports.webGL = function(canvas) {

        // the out object will hold the created wrapper interface.
        var out = {};

        // if we can't get a canvas then abandon since we can't do anything
        if (!canvas) {
            throw new Error("canvas context undefined.");
        }

        // accessor for getting canvas object we're using for webgl
        Object.defineProperty(out, 'canvas', {
            configurable: true,
            get: function() {
                return canvas;
            }
        });

        // get the webgl interface object from the canvas context.
        var gl = canvas.getContext("webgl") || canvas.getContext("experimental-webgl");

        // if that fails then abandon
        if (!gl) {
            throw new Error("webgl context undefined.");
        }

        // accessor for gl object, although using it to alter anything could cause problems
        // but I want to allow access for now at least for testing purposes
        Object.defineProperty(out, 'gl', {
            configurable: true,
            get: function() {
                return gl;
            }
        });


        /**
            @typedef {Object} LinkProgramSpec
            @property {string} vertexShaderSource - vertex shader source code
            @property {string} fragmentShaderSource - fragment shader source code
        */
        /**
            add/link a program source code to webgl

            @param {LinkProgramSpec} spec
        */
        out.linkProgram = function(spec) {
            // validate inputs
            exports.validate_object(spec, {
                vertexShaderSource : 'string',
                fragmentShaderSource : 'string'
            });

            var error;
            var compiled

            // create the vertex shader
            var vshader = gl.createShader(gl.VERTEX_SHADER);

            // add vertex shader source code
            gl.shaderSource(vshader, spec.vertexShaderSource);

            // try to compile shader soruce code
            gl.compileShader(vshader);

            compiled = gl.getShaderParameter(vshader, gl.COMPILE_STATUS);

            if (!compiled) {

              error = gl.getShaderInfoLog(vshader);
              gl.deleteShader(vshader);

              throw new Error("Vertex shader compilation failed: " + error);
            }

            // create fragment shader
            var fshader = gl.createShader(gl.FRAGMENT_SHADER);;

            // add fragment shader source code
            gl.shaderSource(fshader, spec.fragmentShaderSource);

            // try to compile fragment shader
            gl.compileShader(fshader);

            compiled = gl.getShaderParameter(fshader, gl.COMPILE_STATUS);

            if (!compiled) {

              error = gl.getShaderInfoLog(fshader);
              gl.deleteShader(fshader);

              throw new Error("Fragment shader compilation failed: " + error);
            }

            // now that the shaders are compiled, we can create a program and
            // attach the shaders to it. It is assumed there is only 1 vertex
            // shader and 1 fragment shader per program.
            var program = gl.createProgram();

            gl.attachShader(program, vshader);
            gl.attachShader(program, fshader);

            // try to link program.
            gl.linkProgram(program);


            var linked = gl.getProgramParameter(program, gl.LINK_STATUS);

            if (!linked) {
                var error = gl.getProgramInfoLog(program);

                gl.deleteProgram(program);

                throw new Error("Program linking failed: " + error);
            }

            // create an interface to use the program. The webgl interface is
            // exposed for testing purposes but using it to alter anything may
            // lead to undefined behaviour.
            var prog = {
                program : program,
                bindings : []
            };

            /**
                Sets the value of a uniform float variable.

                The value can be any number, or an array of numbers up to length 4.

                Arrays assume the variable is of type vec2, vec3, or vec4 in the source.

                If longer arrays are needed, one must use textures.

                @param {string} uniform_name - The variable name as in the source code.
                @param value - The value to assign to the uniform variable in the next draw

            */
            prog.setUniformFloat = function(uniform_name, value) {
                gl.useProgram(program);

                // get location of variable
                var uniformLocation = gl.getUniformLocation(program, uniform_name);

                if (uniformLocation === null) {
                    throw new Error("Could not find uniform: " + uniform_name);
                }

                // assign value to variable
                if (typeof value === 'number') {

                    gl.uniform1f(uniformLocation, value);

                }else if (Array.isArray(value)){
                    if (value.length === 2) {

                        gl.uniform2f(uniformLocation, value[0], value[1]);

                    }else if (value.length === 3) {

                        gl.uniform3f(uniformLocation, value[0], value[1], value[2]);

                    }else if (value.length === 4) {

                        gl.uniform4f(uniformLocation, value[0], value[1], value[2], value[3]);
                    }else{
                        throw new Error("Cannot add uniform value: " + value);
                    }
                }else{
                    throw new Error("Cannot add uniform value type: " + value);
                }

                // revert back to what was specified as program to use for drawing
                gl.useProgram(null);
            };

            prog.drawTriangles = function(first, count, target, resolution_name) {
                gl.useProgram(program);


                for(var i = 0; i < prog.bindings.length; i++) {
                    prog.bindings[i]();
                }

                if (target) {
                    // make this the framebuffer we are rendering to.
                    gl.bindFramebuffer(gl.FRAMEBUFFER, target.frameBuffer);

                    if (resolution_name) {
                        // Tell the shader the resolution of the framebuffer.
                        var resolutionLocation = gl.getUniformLocation(program, resolution_name);
                        gl.uniform2f(resolutionLocation, target.width, target.height);
                    }

                    // Tell webgl the viewport setting needed for framebuffer.
                    gl.viewport(0, 0, target.width, target.height);
                }else{
                    gl.bindFramebuffer(gl.FRAMEBUFFER, null);

                    if (resolution_name) {
                        // Tell the shader the resolution of the framebuffer.
                        var resolutionLocation = gl.getUniformLocation(program, resolution_name);
                        gl.uniform2f(resolutionLocation, canvas.width, canvas.height);
                    }

                    // Tell webgl the viewport setting needed for framebuffer.
                    gl.viewport(0, 0, canvas.width, canvas.height);
                }

                gl.drawArrays(gl.TRIANGLES, first, count);

                // revert back to what was specified as program to use for drawing
                gl.useProgram(null);
            };


            prog.drawPoints = function(first, count, target, resolution_name) {
                gl.useProgram(program);


                for(var i = 0; i < prog.bindings.length; i++) {
                    prog.bindings[i]();
                }

                if (target) {
                    // make this the framebuffer we are rendering to.
                    gl.bindFramebuffer(gl.FRAMEBUFFER, target.frameBuffer);

                    if (resolution_name) {
                        // Tell the shader the resolution of the framebuffer.
                        var resolutionLocation = gl.getUniformLocation(program, resolution_name);
                        gl.uniform2f(resolutionLocation, target.width, target.height);
                    }

                    // Tell webgl the viewport setting needed for framebuffer.
                    gl.viewport(0, 0, target.width, target.height);
                }else{
                    gl.bindFramebuffer(gl.FRAMEBUFFER, null);

                    if (resolution_name) {
                        // Tell the shader the resolution of the framebuffer.
                        var resolutionLocation = gl.getUniformLocation(program, resolution_name);
                        gl.uniform2f(resolutionLocation, canvas.width, canvas.height);
                    }

                    // Tell webgl the viewport setting needed for framebuffer.
                    gl.viewport(0, 0, canvas.width, canvas.height);
                }

                gl.drawArrays(gl.POINTS, first, count);

                // revert back to what was specified as program to use for drawing
                gl.useProgram(null);
            };

            return prog;
        };

        /**
            Adds vertex data in the form of 2D array. First index is vertex, second
            index are components of vertex data per vertex.

            The number of components per vertex can be 1, 2, 3, or 4 and can bind
            to float, vec2, vec3, and vec4 varying variables.

            @param array - 2D array of vertex data
        */
        out.addVertexData = function(array) {
            var numComponents;
            var monolithicArray;

            if (Array.isArray(array[0]) === false) {
                numComponents = 1;
                monolithicArray = new Float32Array(array.length);

                for(var i = 0; i < array.length; i++) {
                    monolithicArray[i] = array[i];
                }
            }else{
                var numComponents = array[0].length;

                monolithicArray = new Float32Array(numComponents*array.length);
                var mIndex = 0;

                for(var i = 0; i < array.length; i++) {
                    for(var j = 0; j < numComponents; j++){
                        monolithicArray[mIndex] = array[i][j];
                        mIndex++;
                    }
                }
            }


            var attributeBuffer = gl.createBuffer();
            gl.bindBuffer(gl.ARRAY_BUFFER, attributeBuffer);
            gl.bufferData(gl.ARRAY_BUFFER, monolithicArray, gl.STATIC_DRAW);
            gl.bindBuffer(gl.ARRAY_BUFFER, null);

            var buff ={};

            buff.bind = function(program, attribute_name) {
                gl.useProgram(program.program);

                // look up where the texture coordinates need to go.
                var attributeLocation = gl.getAttribLocation(program.program, attribute_name);

                if (attributeLocation === -1) {
                    throw new Error("Could not get attribute: " + attribute_name);
                }

                gl.useProgram(null);

                program.bindings.push(function(){
                    gl.bindBuffer(gl.ARRAY_BUFFER, attributeBuffer);
                    gl.enableVertexAttribArray(attributeLocation);
                    gl.vertexAttribPointer(attributeLocation, numComponents, gl.FLOAT, false, 0, 0);
                    gl.bindBuffer(gl.ARRAY_BUFFER, null);
                });

                return buff;
            };

            return buff;
        };

        // keep a reference to all textures that have been created in webgl
        var textures = [];

        /**
            Enables the use of floating point textures FLOAT, instead of UNSIGNED_BYTE.
            But only if it's supported by webgl
        */
        out.enableFloatTexture = function(){
            var OES_texture_float = gl.getExtension('OES_texture_float');

            if (!OES_texture_float) {
                throw new Error("No support for OES_texture_float");
            }
        };

        out.additiveBlending = function() {
            gl.enable(gl.BLEND);
            gl.blendFunc(gl.ONE, gl.ONE);

            gl.disable(gl.DEPTH_TEST);
        };

        out.disableBlending = function() {
            gl.disable(gl.BLEND);
            gl.enable(gl.DEPTH_TEST);
        };

        /**
            Stores an array as a texture unit.

            If you wish to store it as a FLOAT texture unit, you must successfully
            call enableFloatTexture(), and pass true into useFloat paramter. The array
            must also be of type Float32Array, with rbga in sequence.

            [r,g,b,a,r,g,b,a,....]

            The length of the array must be 4*width*height for rgba values.

            @param {number} width - width of texture in pixels
            @param {number} height - height of texture in pixels
            @param array - data to assign to texture unit
        */
        out.addTextureArray = function(width, height, array, useFloat) {

            // use the next available texture unit to store the array
            gl.activeTexture(gl.TEXTURE0 + textures.length);

            // create a new texture in the given texture unit
            var texture = gl.createTexture();
            gl.bindTexture(gl.TEXTURE_2D, texture);
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);


            if (useFloat) {
                gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, width, height, 0, gl.RGBA, gl.FLOAT, array);
            }else{
                gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, width, height, 0, gl.RGBA, gl.UNSIGNED_BYTE, array);
            }

            // object reference to texture created
            var tex = {
                index : textures.length,
                texture : texture,
                width : width,
                height : height,
                array : array
            };

            /**
                Assigns the texture to a variable in a program.

                @param {object} program - The program to search for the variable
                @param {string} tex_name - the variable name of type sampler2D
                @param {string} tex_size_name - optional name of variable to store the pixel size of the texture as type uniform vec2
            */
            tex.bind = function(program, tex_name, tex_size_name) {

                gl.useProgram(program.program);

                var u_textureLocation = gl.getUniformLocation(program.program, tex_name);
                gl.uniform1i(u_textureLocation, tex.index);

                // if a variable name is specified for texture size, bind that too
                if (tex_size_name) {
                    var textureSizeLocation = gl.getUniformLocation(program.program, tex_size_name);
                    gl.uniform2f(textureSizeLocation, tex.width, tex.height);
                }

                gl.useProgram(null);

                return tex;
            };

            /**
                re-Pushes the current local array data to the texture unit.
            */
            tex.update = function() {
                gl.activeTexture(gl.TEXTURE0 + tex.index);
                gl.bindTexture(gl.TEXTURE_2D, tex.texture);

                if (useFloat) {
                    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, width, height, 0, gl.RGBA, gl.FLOAT, array);
                }else{
                    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, width, height, 0, gl.RGBA, gl.UNSIGNED_BYTE, array);
                }
            };

            textures.push(tex);


            return tex;

        };

        out.addTextureImage = function(image) {

            // use the next available texture unit to store the image
            gl.activeTexture(gl.TEXTURE0 + textures.length);

            var texture = gl.createTexture();
            gl.bindTexture(gl.TEXTURE_2D, texture);
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);

            // add the image
            gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, image);

            gl.bindTexture(gl.TEXTURE_2D, null);

            var tex = {
                index : textures.length,
                texture : texture,
                image : image,
                webgl : out
            };

            tex.bind = function(program, tex_name, tex_size_name) {

                gl.useProgram(program.program);

                var u_textureLocation = gl.getUniformLocation(program.program, tex_name);

                gl.uniform1i(u_textureLocation, tex.index);


                if (tex_size_name) {
                    var textureSizeLocation = gl.getUniformLocation(program.program, tex_size_name);

                    if (textureSizeLocation === null) {
                        throw new Error("Could not get uniform: " + tex_size_name);
                    }

                    gl.uniform2f(textureSizeLocation, tex.width, tex.height);
                }

                gl.useProgram(null);

                return tex;
            };

            tex.update = function() {
                gl.activeTexture(gl.TEXTURE0 + tex.index);
                gl.bindTexture(gl.TEXTURE_2D, tex.texture);

                gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, image);
            };

            textures.push(tex);

            return tex;

        };

        var framebuffers = [];

        /**
            A frame buffer can be used to take the output of one program and
            use it as data for another program. A texture unit is automatically
            used to bind the output of the program it is bound to.

        */
        out.addFrameBuffer = function(width, height, useFloat) {

            var texture = out.addTextureArray(width, height, null, useFloat);

            //out.gl.bindTexture(gl.TEXTURE_2D, texture);

            var fbo = gl.createFramebuffer();


            gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);

            // Attach a texture to it.
            gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, texture.texture, 0);


            if (gl.checkFramebufferStatus(gl.FRAMEBUFFER) !== gl.FRAMEBUFFER_COMPLETE) {
                    throw new Error("Frame buffer creation failed");
            }

            gl.bindFramebuffer(gl.FRAMEBUFFER, null);

            var fb = {
                index : framebuffers.length,
                frameBuffer: fbo,
                texture : texture,
                width : width,
                height : height,
                webgl : out
            };

            /**
                Use this framebuffer as the output of the program.
            */
            fb.bind = function(program, resolution_name) {

                gl.useProgram(program.program);

                // make this the framebuffer we are rendering to.
                gl.bindFramebuffer(gl.FRAMEBUFFER, fb.frameBuffer);

                if (resolution_name) {
                    // Tell the shader the resolution of the framebuffer.
                    var resolutionLocation = gl.getUniformLocation(program.program, resolution_name);
                    gl.uniform2f(resolutionLocation, fb.width, fb.height);
                }

                // Tell webgl the viewport setting needed for framebuffer.
                gl.viewport(0, 0, fb.width, fb.height);

                gl.useProgram(null);
            };

            fb.readPixels = function(array) {
                gl.bindFramebuffer(gl.FRAMEBUFFER, fb.frameBuffer);

                if (useFloat) {
                    gl.readPixels(0, 0, width, height, gl.RGBA, gl.FLOAT, array);
                }else{
                    gl.readPixels(0, 0, width, height, gl.RGBA, gl.UNSIGNED_BYTE, array);
                }

                gl.bindFramebuffer(gl.FRAMEBUFFER, null);
            }

            fb.clear = function (r, g, b, a) {
                gl.bindFramebuffer(gl.FRAMEBUFFER, fb.frameBuffer);
                gl.clearColor(r, g, b, a);
                gl.clear(gl.COLOR_BUFFER_BIT);
                gl.bindFramebuffer(gl.FRAMEBUFFER, null);
            };

            framebuffers.push(fb);

            return fb;
        };

        out.clear = function(r, g, b, a) {
            gl.bindFramebuffer(gl.FRAMEBUFFER, null);
            gl.clearColor(r, g, b, a);
            gl.clear(gl.COLOR_BUFFER_BIT);
        };

        /**
            Use the canvas as the output buffer of the program
        */
        out.bindCanvas = function (program, resolution_name) {
            gl.useProgram(program.program);

            gl.bindFramebuffer(gl.FRAMEBUFFER, null);

            if (resolution_name) {
                // Tell the shader the resolution of the framebuffer.
                var resolutionLocation = gl.getUniformLocation(program.program, resolution_name);
                gl.uniform2f(resolutionLocation, canvas.width, canvas.height);
            }

            // Tell webgl the viewport setting needed for framebuffer.
            gl.viewport(0, 0, canvas.width, canvas.height);

            gl.useProgram(null);
        };

        out.getCanvasPixels = function(array) {

            gl.bindFramebuffer(gl.FRAMEBUFFER, null);
            gl.readPixels(0, 0, canvas.width, canvas.height, gl.RGBA, gl.UNSIGNED_BYTE, array);

        };

        return out;
    };

    exports.webgl_max = function(input, validated){

        if (!validated) {
            exports.validate_object(input, {
                in : 'object',
                out : 'object'
            });
        }

        if (input.in.webgl !== input.out.webgl) {
            throw new Error("Buffers must be in same webgl context.");
        }

        var webgl = input.in.webgl;

        var output = {};

        var source_buffer = input.in;
        var target_buffer = input.out;

        var cells_i = source_buffer.width / target_buffer.width;
        var cells_j = source_buffer.height / target_buffer.height;

        if (cells_i % 1 !== 0 || cells_j % 1 !== 0) {
            throw new Error("source dimensions do not match integer target dimensions.");
        }

        var render_vert = function() {
            var src_arr = [
                "attribute vec2 a_position;",
                "attribute vec2 a_texCoord;",
                //"uniform vec2 u_resolution;",
                "varying vec2 v_texCoord;",

                "void main() {",
                "    gl_Position = vec4(a_position, 0, 1);",

                "    v_texCoord = a_texCoord;",
                "}",
            ];


            return src_arr.join('\n');
        } // render_vert()

        var N = function(num) {
            return num.toFixed(20);
        };

        var max_frag = function() {

            var max_arr = [];
            var startj = 1;
            for(var i = 0; i < cells_i; i++) {
                for(var j = startj; j < cells_j; j++) {
                    max_arr.push("    gl_FragColor = max(gl_FragColor, texture2D(u_source, v_texCoord + vec2(" + N(i/source_buffer.width) + "," + N(j/source_buffer.height) + ")));");
                }
                startj = 0;
            }


            var src_arr = [
                "precision mediump float;",

                "uniform sampler2D u_source;",

                // the texCoords passed in from the vertex shader.
                "varying vec2 v_texCoord;",

                "void main() {",
                "    gl_FragColor = texture2D(u_source, v_texCoord);",
                max_arr.join("\n"),
                "}"
            ];
            var src = src_arr.join('\n');

            //console.log(src)


            return src;
        }; // max_frag()

        var vertex_positions = webgl.addVertexData([
            [-1, 1],
            [1, 1],
            [-1, -1],
            [-1, -1],
            [1, 1],
            [1, -1]
        ]);

        // texture coordinets for vertices
        var texture_coordinates = webgl.addVertexData([
            [0.0,  1.0],
            [1.0,  1.0],
            [0.0,  0.0],
            [0.0,  0.0],
            [1.0,  1.0],
            [1.0,  0.0]
        ]);

        var program = source_buffer.webgl.linkProgram({
            vertexShaderSource : render_vert(),
            fragmentShaderSource : max_frag()
        });

        source_buffer.texture.bind(program, "u_source");
        vertex_positions.bind(program, "a_position");
        texture_coordinates.bind(program, "a_texCoord");

        output.compute = function () {

            program.drawTriangles(0, 6, target_buffer);

        };


        return output;
    };


    exports.webgl_avg = function(input, validated){

        if (!validated) {
            exports.validate_object(input, {
                in : 'object',
                out : 'object'
            });
        }

        if (input.in.webgl !== input.out.webgl) {
            throw new Error("Buffers must be in same webgl context.");
        }

        var webgl = input.in.webgl;

        var output = {};

        var source_buffer = input.in;
        var target_buffer = input.out;

        var cells_i = source_buffer.width / target_buffer.width;
        var cells_j = source_buffer.height / target_buffer.height;

        if (cells_i % 1 !== 0 || cells_j % 1 !== 0) {
            throw new Error("source dimensions do not match integer target dimensions.");
        }

        var render_vert = function() {
            var src_arr = [
                "attribute vec2 a_position;",
                "attribute vec2 a_texCoord;",
                //"uniform vec2 u_resolution;",
                "varying vec2 v_texCoord;",

                "void main() {",
                "    gl_Position = vec4(a_position, 0, 1);",

                "    v_texCoord = a_texCoord;",
                "}",
            ];


            return src_arr.join('\n');
        } // render_vert()

        var N = function(num) {
            return num.toFixed(20);
        };

        var avg_frag = function() {

            var avg_arr = [];

            var startj = 1;
            for(var i = 0; i < cells_i; i++) {
                for(var j = startj; j < cells_j; j++) {
                    avg_arr.push("        texture2D(u_source, v_texCoord + vec2(" + N(i/source_buffer.width) + "," + N(j/source_buffer.height) + "))");
                }

                startj = 0;
            }


            var src_arr = [
                "precision mediump float;",

                "uniform sampler2D u_source;",

                // the texCoords passed in from the vertex shader.
                "varying vec2 v_texCoord;",

                "void main() {",
                "    gl_FragColor = " + N(1.0/(cells_i*cells_j)) + " * (texture2D(u_source, v_texCoord) + ",
                avg_arr.join(" + \n"),
                "    );",
                "}"
            ];
            var src = src_arr.join('\n');

            //console.log(src)


            return src;
        }; // avg_frag()

        var vertex_positions = webgl.addVertexData([
            [-1, 1],
            [1, 1],
            [-1, -1],
            [-1, -1],
            [1, 1],
            [1, -1]
        ]);

        // texture coordinets for vertices
        var texture_coordinates = webgl.addVertexData([
            [0.0,  1.0],
            [1.0,  1.0],
            [0.0,  0.0],
            [0.0,  0.0],
            [1.0,  1.0],
            [1.0,  0.0]
        ]);

        var program = source_buffer.webgl.linkProgram({
            vertexShaderSource : render_vert(),
            fragmentShaderSource : avg_frag()
        });

        source_buffer.texture.bind(program, "u_source");
        vertex_positions.bind(program, "a_position");
        texture_coordinates.bind(program, "a_texCoord");

        output.compute = function () {

            program.drawTriangles(0, 6, target_buffer);

        };


        return output;
    }; // webgl_avg

    /**
     *
     * @param {type} input
     {
     max:,
     min:,
     }
     * @returns {newRange.output}
     */
    exports.makeRange = function (input, validated)
    {
        if (!validated) {
            exports.validate_object(input, {
                max: 'number',
                min: 'number'
            });
        }

        var slope;
        var intercept;

        var output = {};

        var recalc = function(){
            slope = 1 / (input.max - input.min);
            intercept = -slope * input.min;
        };

        recalc();

        Object.defineProperty(output, 'max', {

            get: function () {
                return input.max;
            },

            set: function (value) {
                input.max = value;
                recalc();
            }
        });

        Object.defineProperty(output, 'min', {
            get: function () {
                return input.min;
            },

            set: function (value) {
                input.min = value;
                recalc();
            }
        });



        output.norm = function (x)
        {
            return Math.max(0, Math.min(1, slope * x + intercept));
        };

        return output;
    };

    /**
     * @param min
     * @param max
     * @param {type} n
     * @param {type} params structure
     * {
     *     // piecewise function for each color x \in [0,1], y \in [0,1]
     *     r : [[x1, x2, y(x1), y(x2)],...],
     *     g : [[x1, x2, y(x1), y(x2)],...],
     *     b : [[x1, x2, y(x1), y(x2)],...]
     * }
     * @returns {ColorMap}
     */
    exports.makeColorMap = function (input, validated){
        if(!validated){
            exports.validate_object(input, {
                max: 'number',
                min: 'number',
                n: 'number',
                params : 'object'
            });
        }

        var output = exports.makeRange(input, true);

        var mapr = new Uint8ClampedArray(input.n);
        var mapg = new Uint8ClampedArray(input.n);
        var mapb = new Uint8ClampedArray(input.n);

        var x = 0;
        var s = 0;
        var y = 0;

        // build the color map based on the color function parameters
        for (var i = 0; i < input.n; i++)
        {
            // inputs are always x \in [0, 1.0]
            x = i / (input.n - 1);

            // red component

            for (var j = 0; j < input.params.r.length; j++)
            {
                if (x >= input.params.r[j][0] && x <= input.params.r[j][1])
                {
                    s = (x - input.params.r[j][0]) / (input.params.r[j][1] - input.params.r[j][0]);
                    y = (1 - s) * input.params.r[j][2] + s * input.params.r[j][3];

                    mapr[i] = Math.floor(255 * y);
                }
            }


            // green component

            for (var j = 0; j < input.params.g.length; j++)
            {
                if (x >= input.params.g[j][0] && x <= input.params.g[j][1])
                {
                    s = (x - input.params.g[j][0]) / (input.params.g[j][1] - input.params.g[j][0]);
                    y = (1 - s) * input.params.g[j][2] + s * input.params.g[j][3];

                    mapg[i] = Math.floor(255 * y);
                }
            }

            // blue component
            for (var j = 0; j < input.params.b.length; j++)
            {
                if (x >= input.params.b[j][0] && x <= input.params.b[j][1])
                {
                    s = (x - input.params.b[j][0]) / (input.params.b[j][1] - input.params.b[j][0]);
                    y = (1 - s) * input.params.b[j][2] + s * input.params.b[j][3];

                    mapb[i] = Math.floor(255 * y);
                }
            }
        }

        Object.defineProperty(output, 'mapr', {
            get: function () {
                return mapr;
            }
        });

        Object.defineProperty(output, 'mapg', {
            get: function () {
                return mapg;
            }
        });

        Object.defineProperty(output, 'mapb', {
            get: function () {
                return mapb;
            }
        });

        output.n = input.n;

        output.r = function (x)
        {

            var index = Math.floor((input.n - 1) * output.norm(x));

            return mapr[index];
        };

        output.g = function (x)
        {
            var index = Math.floor((input.n - 1) * output.norm(x));

            return mapg[index];
        };

        output.b = function (x)
        {
            var index = Math.floor((input.n - 1) * output.norm(x));

            return mapb[index];
        };

        return output;
    };


    exports.colormap_presets = {
        jet: {
            r: [[0.4, 0.6, 0, 1], [0.6, 0.9, 1, 1], [0.9, 1, 1, 0.5]],
            g: [[0.1, 0.4, 0, 1], [0.4, 0.6, 1, 1], [0.6, 0.9, 1, 0]],
            b: [[0, 0.1, 0.5, 1], [0.1, 0.4, 1, 1], [0.4, 0.6, 1, 0]]
        },
        hot: {
            r: [[0, 0.35, 0, 1], [0.35, 1, 1, 1]],
            g: [[0.35, 0.65, 0, 1], [0.65, 1, 1, 1]],
            b: [[0.65, 1, 0, 1]]
        },
        rainbow: {
            r: [[0, 0.2, 1, 1], [0.2, 0.4, 1, 0], [0.8, 1, 0, 1]],
            g: [[0, 0.2, 0, 1], [0.2, 0.6, 1, 1], [0.6, 0.8, 1, 0]],
            b: [[0.4, 0.6, 0, 1], [0.6, 1, 1, 1]]
        },
        gray: {
            r: [[0, 1, 0, 1]],
            g: [[0, 1, 0, 1]],
            b: [[0, 1, 0, 1]]
        },
        bone: {
            r: [[0, 1, 0, 1]],
            g: [[0, 1, 0, 1]],
            b: [[0, 0.5, 0, 0.65], [0.5, 1, 0.65, 1]]
        },
        violet: {
            r: [[0, 0.5, 0, 1], [0.5, 1, 1, 1]],
            g: [[0.5, 1, 0, 1]],
            b: [[0, 0.5, 0, 1], [0.5, 1, 1, 1]]
        },
        yellow: {
            r: [[0, 0.5, 0, 1], [0.5, 1, 1, 1]],
            g: [[0, 0.5, 0, 1], [0.5, 1, 1, 1]],
            b: [[0.5, 1, 0, 1]]
        },
        cyan: {
            r: [[0.5, 1, 0, 1]],
            g: [[0, 0.5, 0, 1], [0.5, 1, 1, 1]],
            b: [[0, 0.5, 0, 1], [0.5, 1, 1, 1]]
        },
        red_violet: {
            r: [[0, 0.33, 0, 1], [0.33, 1, 1, 1]],
            g: [[0.66, 1, 0, 1]],
            b: [[0.33, 0.66, 0, 1], [0.66, 1, 1, 1]]
        },
        green_cyan: {
            r: [[0.66, 1, 0, 1]],
            g: [[0, 0.33, 0, 1], [0.33, 1, 1, 1]],
            b: [[0.33, 0.66, 0, 1], [0.66, 1, 1, 1]]
        },
        green: {
            r: [[0.5, 1, 0, 1]],
            g: [[0, 0.5, 0, 1], [0.5, 1, 1, 1]],
            b: [[0.5, 1, 0, 1]]
        },
        red: {
            r: [[0, 0.5, 0, 1], [0.5, 1, 1, 1]],
            g: [[0.5, 1, 0, 1]],
            b: [[0.5, 1, 0, 1]]
        },
        blue: {
            r: [[0.5, 1, 0, 1]],
            g: [[0.5, 1, 0, 1]],
            b: [[0, 0.5, 0, 1], [0.5, 1, 1, 1]]
        },
        blue_cyan: {
            r: [[0.66, 1, 0, 1]],
            g: [[0.33, 0.66, 0, 1], [0.66, 1, 1, 1]],
            b: [[0, 0.33, 0, 1], [0.33, 1, 1, 1]]
        },
        mud: {
            r: [[0, 1, 0, 1]],
            g: [[0.33, 1, 0, 1]],
            b: [[0.66, 1, 0, 1]]
        },
        grass: {
            r: [[0.33, 1, 0, 1]],
            g: [[0, 1, 0, 1]],
            b: [[0.66, 1, 0, 1]]
        },
        purplehaze: {
            r: [[0.33, 1, 0, 1]],
            g: [[0.66, 1, 0, 1]],
            b: [[0, 1, 0, 1]]
        },
        atmosphere: {
            r: [[0.66, 1, 0, 1]],
            g: [[0.33, 1, 0, 1]],
            b: [[0, 1, 0, 1]]
        },
        pond: {
            r: [[0.66, 1, 0, 1]],
            g: [[0, 1, 0, 1]],
            b: [[0.33, 1, 0, 1]]
        },
        berry: {
            r: [[0, 1, 0, 1]],
            g: [[0.66, 1, 0, 1]],
            b: [[0.33, 1, 0, 1]]
        },
        doppler: {
            r: [[0, 0.5, 1, 1], [0.5, 0.75, 1, 0]],
            g: [[0, 0.5, 0, 1], [0.5, 1, 1, 0]],
            b: [[0.25, 0.5, 0, 1], [0.5, 1, 1, 1]]
        },
        autumn: {
            r: [[0, 1, 1, 1]],
            g: [[0, 1, 0, 1]],
            b: []
        },
        spring: {
            r: [[0, 1, 1, 1]],
            g: [[0, 1, 0, 1]],
            b: [[0, 1, 1, 0]]
        },
        winter: {
            r: [[0, 1, 0, 0.3]],
            g: [[0, 1, 0, 1]],
            b: [[0, 1, 1, 0.3]]
        },
        anime: {
            r: [[0, 0.33, 0.7, 0.9], [0.33, 0.66, 0.9, 0], [0.95, 1, 0, 1]],
            g: [[0, 0.33, 0.1, 0.7], [0.33, 0.66, 0.7, 0], [0.66, 1, 0, 1]],
            b: [[0, 0.33, 0.15, 0.3], [0.33, 0.66, 0.3, 0.7], [0.66, 1, 0.7, 1]]
        }
    };

    exports.makePlot2DArea = function (input) {
        exports.validate_object(input, {
            x: 'number',
            y: 'number',
            width: 'number',
            height: 'number',
            data_array : [,'object'],
            data_canvas : [, 'object'],
            sx : [,'number'],
            sy : [,'number'],
            swidth : [,'number'],
            sheight : [,'number'],
            colormap : [,'object']
        });

        var output = {};


        if (input.data_array) {

            var canvas = document.createElement("canvas");
            var ctx;
            var imgData;



            if (canvas.getContext !== null) {
                canvas.width = input.data_array.num_i;
                canvas.height = input.data_array.num_j;

                ctx = canvas.getContext('2d');

                imgData = ctx.createImageData(canvas.width, canvas.height);

            } else {
                throw new Error("could not get context of canvas");
            }
        }else if (input.data_canvas) {
            // no initialization for canvas data
        }else{
            throw new Error("No data to bind to!");
        }

        /**
         * Updates the plot using current contents of the input.data field
         *
         * @returns {undefined}
         */
        if (input.data_array) {
            output.redraw = function (viewport, t)
            {

                var k = 0;
                var arr = input.data_array.arr;


                for (var j = 0; j < input.data_array.num_j; j++)
                {

                    for (var i = 0; i < input.data_array.num_i; i++)
                    {

                        imgData.data[k * 4] = input.colormap.r(arr[k]);
                        imgData.data[k * 4 + 1] = input.colormap.g(arr[k]);
                        imgData.data[k * 4 + 2] = input.colormap.b(arr[k]);
                        imgData.data[k * 4 + 3] = 255;


                        k++;
                    }
                }

                ctx.putImageData(imgData, 0, 0);

                viewport.ctx.drawImage(canvas, input.x, input.y, input.width, input.height);
            }
        }else if (input.data_canvas) {

            if (typeof input.sx !== 'undefined' && typeof input.sy !== 'undefined' && typeof input.swidth !== 'undefined' && typeof input.sheight !== 'undefined'){
                output.redraw = function (viewport, t)
                {
                    viewport.ctx.drawImage(
                        input.data_canvas,
                        input.sx,
                        input.sy,
                        input.swidth,
                        input.sheight,
                        input.x,
                        input.y,
                        input.width,
                        input.height);
                };
            }else{
                output.redraw = function (viewport, t)
                {
                    viewport.ctx.drawImage(input.data_canvas, input.x, input.y, input.width, input.height);
                };
            }
        }



        return output;
    };

    /**
     *
     * @param {type} input
     {
     x:,
     y:,
     width:,
     height:,
     colormap : ColorMap object
     }
     * @returns {unresolved}
     */
    exports.makeColorBar = function (input) {

        exports.validate_object(input, {
            x: 'number',
            y: 'number',
            width: 'number',
            height: 'number',
            colormap : 'object'
        });

        var output = {};

        var canvas = document.createElement("canvas");
        var ctx;
        var imgData;



        if (canvas.getContext !== null) {
            canvas.width = input.width;
            canvas.height = input.height;

            ctx = canvas.getContext('2d');

            imgData = ctx.createImageData(canvas.width, canvas.height);

        } else {
            throw "could not get context of canvas";
        }

        output.redraw = function (viewport, t){
            var k = 0;
            var x = 0;

            for (var j = 0; j < input.height; j++)
            {

                for (var i = 0; i < input.width; i++)
                {

                    x = input.colormap.min + (input.colormap.max - input.colormap.min) * (input.height - 1 - j) / (input.height - 1);

                    imgData.data[k * 4] = input.colormap.r(x);
                    imgData.data[k * 4 + 1] = input.colormap.g(x);
                    imgData.data[k * 4 + 2] = input.colormap.b(x);
                    imgData.data[k * 4 + 3] = 255;


                    k++;
                }
            }

            ctx.putImageData(imgData, 0, 0);

            viewport.ctx.drawImage(canvas, input.x, input.y, input.width, input.height);
        };

        return output;
    };


    exports.makeCanvasFigure = function (input) {

        exports.validate_object(input, {
            id: 'string'
        });

        var output = {};

        output.input = input;



        output.canvas = document.getElementById(input.id);
        output.ctx = output.canvas.getContext('2d');
        output.element = $("#" + input.id);

        output.beginFrame = null;
        output.endFrame = null;

        output.layers = [];

        /**
         *
         * @param {type} layer
         {
         redraw: function(viewport, t)
         }
         * @returns {undefined}
         */
        output.addLayer = function (layer) {
            if (typeof layer.redraw === 'undefined')
            {
                throw "layer.redraw(figure, t) must be defined";
            }

            output.layers.push(layer);
        };


        /**
         * Redraws the canvas drawing each layer in the order that they were added
         * @param
         * @param {type} t time since the beginning of the animation.
         * @returns {undefined}
         */
        output.redraw = function (t) {


            output.ctx.setTransform(1, 0, 0, 1, 0, 0);
            output.ctx.clearRect(0, 0, output.canvas.width, output.canvas.height);

            if (output.beginFrame) {
                output.beginFrame(output, t);
            }

            var l = output.layers.length;

            for (var i = 0; i < l; i++)
            {
                output.ctx.setTransform(1, 0, 0, 1, 0, 0);
                output.layers[i].redraw(output, t);
            }

            if (output.endFrame) {
                output.endFrame(output, t);
            }

        };

        var selectables = [];
        var selected = [];

        /**
            @typedef {Object} Selectable
            @property {function} in - test if (x,y) is inside selectable area
            @property {function} selected - callback when selected
            @property {function} unselected - callback when un-selected
        */
        /**
         *
         * @param {Selectable} selectable
         * @returns {undefined}
         */
        output.addSelectable = function (selectable) {

            if (typeof selectable.selected === 'undefined' || typeof selectable.unselected === 'undefined')
            {
                throw new Error("selected() and unselected() must be defined.");
            }

            selectables.push(selectable);
            selected.push(false);
        };

        output.getSelections = function () {
            var selections = [];

            var l = selectables.length;

            for (var i = 0; i < l; i++)
            {
                if (selected[i]) {
                    selections.push(selectables[i]);
                }
            }

            return selections;
        };


        var mouse = {
            over: false,
            down: false,
            downx: 0,
            downy: 0,
            upx: 0,
            upy: 0,
            curx: 0,
            cury: 0
        };


        output.element.mousedown(function (e) {
            mouse.down = true;
            mouse.downx = output.canvas.width*(e.pageX - output.element.offset().left)/output.element.width();
            mouse.downy = output.canvas.height*(e.pageY - output.element.offset().top)/output.element.height();


            var l = selectables.length;
            var selectable;

            for (var i = 0; i < l; i++)
            {
                selectable = selectables[i];

                if (selected[i])
                {

                    if (selectable.in(mouse.downx, mouse.downy))
                    {
                        selected[i] = false;
                        selectable.unselected();
                    } else if (!e.ctrlKey) {
                        selected[i] = false;
                        selectable.unselected();
                    }

                } else {
                    if (selectable.in(mouse.downx, mouse.downy)) {

                        selected[i] = true;
                        selectable.selected(output);

                    }
                }

            }
        });

        output.element.mouseup(function (e) {
            if (mouse.down === true)
            {
                mouse.down = false;
                mouse.upx = e.pageX - output.element.offset().left;
                mouse.upy = e.pageY - output.element.offset().top;
            }

            //console.log("@up");
        });

        output.element.mouseout(function (e) {
            mouse.over = false;
            mouse.down = false;

            //console.log("@out");
        });

        output.element.mouseover(function (e) {
            mouse.over = true;

            //console.log("@over");
        });

        output.element.mousemove(function (e) {
            mouse.curx = e.pageX - output.element.offset().left;
            mouse.cury = e.pageY - output.element.offset().top;

        });



        return output;
    };

    /**
        @typedef {Object} SquareClickAreaInput
        @property {object} viewport - the viewport to put selectable object
        @property {number} x1 - left edge
        @property {number} x2 - right edge
        @property {number} y1 - top edge
        @property {number} y2 - bottom edge
    */
    /**
     *
     * @param {SquareClickAreaInput} input
     * @returns {object}
     */
    exports.makeSquareClickArea = function (input) {

        exports.validate_object(input, {
            x1: 'number',
            x2: 'number',
            y1: 'number',
            y2: 'number',
        });

        var output = {};


        var selected = false;

        output.in = function (x, y) {
            if (x > input.x1 && x < input.x2) {

                if (y > input.y1 && y < input.y2)
                {
                    return true;
                }
            }

            return false;
        };



        output.selected = function () {
            selected = true;
        };

        output.unselected = function () {
            selected = false;
        };

        output.redraw = function (viewport) {
            var ctx = viewport.ctx;

            ctx.beginPath();
            ctx.moveTo(input.x1, input.y1);
            ctx.lineTo(input.x2, input.y1);
            ctx.lineTo(input.x2, input.y2);
            ctx.lineTo(input.x1, input.y2);
            ctx.closePath();
            ctx.stroke();

            if (selected) {

                ctx.fillStyle = "rgba(0, 0, 200, 0.5)";
                ctx.fillRect(input.x1, input.y1, input.x2 - input.x1, input.y2 - input.y1);
            } else {

            }
        };

        return output;
    };

    /**
        @typedef {Object} ImageClickAreaInput
        @property {object} viewport - the viewport to put selectable object
        @property {number} x - left
        @property {number} y - top
        @property {number} width -
        @property {number} height -
    */

    /**
     *
     * @param {type} input
     {
        x:
        y:
        width:
        height:
        image:
        callback:
     }
     * @returns {}
     */
    exports.makeImageClickArea = function (input) {
        exports.validate_object(input, {
            x: 'number',
            y: 'number',
            width: 'number',
            height: 'number',
            image: 'object',
            callback: 'function'
        });

        var output = {};

        var selected = false;

        output.in = function (x, y) {
            if (x > input.x && x < input.x + input.width) {

                if (y > input.y && y < input.y + input.width) {
                    return true;
                }
            }

            return false;
        };



        output.selected = function (figure) {
            selected = true;

            input.callback(figure);
        };

        output.unselected = function () {
            selected = false;
        };

        output.redraw = function (figure) {
            var ctx = figure.ctx;

            ctx.drawImage(input.image, input.x, input.y, input.width, input.height);

            ctx.beginPath();
            ctx.moveTo(input.x, input.y);
            ctx.lineTo(input.x + input.width, input.y);
            ctx.lineTo(input.x + input.width, input.y + input.height);
            ctx.lineTo(input.x, input.y + input.height);
            ctx.closePath();

            ctx.lineWidth = 1;
            ctx.strokeStyle = 'black';
            ctx.stroke();

        };

        return output;
    };

    /**
        An animation is a collection of figures which are redrawn together in time
    */
    exports.makeAnimation = function () {

        var output = {};

        output.beginFrame = null;
        output.endFrame = null;

        var figures = [];

        output.addFigure = function (figure) {

            if (typeof figure.redraw === 'undefined')
            {
                throw new Error("layer.redraw(t) must be defined");
            }

            figures.push(figure);
        };

        /**
         * Redraws the canvas drawing each layer in the order that they were added
         * @param
         * @param {type} t time since the beginning of the animation.
         * @returns {undefined}
         */
        output.redraw = function (t) {

            if (output.beginFrame) {
                output.beginFrame(t);
            }


            for (var i=0, l = figures.length; i < l; i++) {
                figures[i].redraw(t);
            }

            if (output.endFrame) {
                output.endFrame(t);
            }

        };

        var animate = false;
        output.fpsCallback = null;


        /**
         * Begins an animation loop of drawing the canvas
         * @param {type} duration number of milliseconds to run animation
         * @returns {undefined}
         */
        output.startAnimation = function (duration) {

            animate = true;

            if (typeof duration === 'undefined') {
                duration = 0;
            }

            window.requestAnimationFrame = window.requestAnimationFrame || window.mozRequestAnimationFrame ||
                    window.webkitRequestAnimationFrame || window.oRequestAnimationFrame;

            var start = null;
            var fpsstart, fpstime;
            var framecount = 0;
            var t, dt, last_t=0;

            function step(timestamp) {

                if (animate === true)
                {


                    if (start === null) {
                        start = timestamp;
                        fpsstart = timestamp;
                        last_t = timestamp;
                    }

                    t = timestamp - start;
                    dt = last_t - timestamp;
                    last_t = timestamp;

                    if (duration === 0 || t < duration)
                    {

                        output.redraw(t);



                        fpstime = timestamp - fpsstart;
                        framecount++;

                        if (fpstime >= 1000) {
                            if (output.fpsCallback)
                            {
                                output.fpsCallback(Math.round(1000 * framecount / fpstime), fpstime/1000.0);
                            }

                            fpsstart = timestamp;
                            framecount = 0;
                        }

                        requestAnimationFrame(step);
                    } else {
                        if (output.fpsCallback)
                        {
                            output.fpsCallback(0);
                        }
                    }


                } else {
                    if (output.fpsCallback)
                    {
                        output.fpsCallback(0);
                    }
                }
            }

            requestAnimationFrame(step);
        };

        output.stopAnimation = function () {
            animate = false;

        };



        Object.defineProperty(output, 'animate', {
            configurable: true,
            set: function (value) {
                if (value === true && animate === false) {

                    output.startAnimation();

                } else if (value === false) {

                    output.stopAnimation();
                }
            },
            get: function () {
                return animate;
            }
        });

        return output;
    };

    return exports;

});
