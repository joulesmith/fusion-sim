/*
    Copyright (C) 2015  Joulesmith Energy Technologies, LLC

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

define(['utilities'], function (util){
    "use strict";

    var exports = {};

    var N = function(num) {
        return num.toFixed(20);
    };

    var speed_of_light = 2.998e8; // m/s


    exports.makeCylindricalParticlePusher = function(spec) {
        util.validate_object(spec, {
            radius : 'number', // meters
            height : 'number',
            nr : 'number',
            nz : 'number',
            dt : 'number',
            nparticles : 'number',
            particle_mass : 'number',
            particle_charge : 'number'
        });

        // physical quantities
        var h = spec.particle_charge * spec.dt / (2 * spec.particle_mass);
        var factor_r = 1/spec.radius;
        var factor_z = 1/spec.height;

        var out = {};

        var i, j, k;

        // create canvas element for webgl to work on
        var canvas = document.createElement("CANVAS");
        canvas.id = "webgl_canvas_CylindricalParticlePusher";
        canvas.width = spec.nr;
        canvas.height = spec.nz;
        canvas.style.display = "none";

        document.body.appendChild(canvas);
        out.canvas = canvas;

        var webgl = util.webGL(canvas);

        webgl.enableFloatTexture();

        // so that particle rendering adds densities
        webgl.additiveBlending();

        var vertex_positions = webgl.addVertexData([
            [-1, 1],
            [1, 1],
            [-1, -1],
            [-1, -1],
            [1, 1],
            [1, -1]
        ]);

        var render_vertex_positions = webgl.addVertexData([
            [-1, -1],
            [1, -1],
            [-1, 1],
            [-1, 1],
            [1, -1],
            [1, 1]
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

        var render_texture_coordinates = webgl.addVertexData([
            [0.0,  0.0],
            [1.0,  0.0],
            [0.0,  1.0],
            [0.0,  1.0],
            [1.0,  0.0],
            [1.0,  1.0]
        ]);

        //
        // The particles
        //

        var nparticles_h = spec.nparticles;
        var nparticles_w = spec.nparticles;
        var nparticles = nparticles_h * nparticles_w;

        var position_arr = new Float32Array(4 * nparticles);

        var position_tex = webgl.addTextureArray(
            nparticles_w,
            nparticles_h,
            position_arr,
            true
        );

        var velocity_arr = new Float32Array(4 * nparticles);

        var velocity_tex = webgl.addTextureArray(
            nparticles_w,
            nparticles_h,
            velocity_arr,
            true
        );

        //
        // The electromagnetic fields
        //

        var E = webgl.addFrameBuffer(spec.nr, spec.nz, true);

        var E_arr = new Float32Array(4 * spec.nr * spec.nz);

        var E_tex = webgl.addTextureArray(
            spec.nr,
            spec.nz,
            E_arr,
            true
        );

        var B = webgl.addFrameBuffer(spec.nr, spec.nz, true);

        var B_arr = new Float32Array(4 * spec.nr * spec.nz);

        var B_tex = webgl.addTextureArray(
            spec.nr,
            spec.nz,
            B_arr,
            true
        );


        //
        // General shader to compute a texture frame buffer.
        //
        var precompute_vert = function() {
            var src_arr = [
                "attribute vec2 a_position;",
                "attribute vec2 a_texCoord;",
                "varying vec2 v_texCoord;",

                "void main() {",
                "    gl_Position = vec4(a_position, 0, 1);",

                "    v_texCoord = a_texCoord;",
                "}",
            ];

            return src_arr.join('\n');
        } // precompute_vert()


        //
        // Program to calculate magnetic field from a current loop.
        //

        var B_loop_half = webgl.addFrameBuffer(spec.nr, spec.nz, true);
        var B_loop_tenth = webgl.addFrameBuffer(spec.nr, spec.nz, true);

        // ---------------------------------------------------------------------
        var programCurrentLoopShape = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform float u_R;",
                    "uniform float u_a;",

                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",

                        "float Bx = 0.0;",
                        "float Bz = 0.0;",
                        "float r = 0.0;",
                        "float factor = 0.0;",
                        "float constant = u_R * 0.001 * 1.25663706e-6 / (4.0 * 3.14159265359);",
                        "for(float i = 0.0; i < 1000.0; i++){",
                            "r = sqrt(u_R * u_R + v_texCoord.x * v_texCoord.x + v_texCoord.y * v_texCoord.y - 2.0 * v_texCoord.x * u_R * cos(3.14159265359 * i / 999.0));",
                            "factor = (r > 0.0) ? constant * (1.0 - exp(- r * r / (u_a * u_a))) / (r * r * r) : 0.0;",
                            "Bx += v_texCoord.y * factor * cos(3.14159265359 * i / 999.0);",
                            "Bz += factor * (u_R - v_texCoord.x * cos(3.14159265359 * i / 999.0));",
                        "}",
                        "gl_FragColor = vec4(Bx, 0.0, Bz, 1.0);",

                    "}"
                ];

                return src_arr.join('\n');
            })()
        });

        // triangle vertices
        vertex_positions.bind(programCurrentLoopShape, "a_position");
        // texture coordinets for vertices
        texture_coordinates.bind(programCurrentLoopShape, "a_texCoord");

        programCurrentLoopShape.setUniformFloat("u_a", 50.0 / spec.nr);

        programCurrentLoopShape.setUniformFloat("u_R", 0.5);

        programCurrentLoopShape.drawTriangles(0, 6, B_loop_half);

        programCurrentLoopShape.setUniformFloat("u_R", 0.1);

        programCurrentLoopShape.drawTriangles(0, 6, B_loop_tenth);


        // ---------------------------------------------------------------------
        var programCurrentLoop = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform float u_R;",
                    "uniform float u_Z;",
                    "uniform float u_I;",

                    "uniform sampler2D u_shape_half;",
                    "uniform sampler2D u_shape_tenth;",

                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "float a = v_texCoord.x / u_R;",
                        "float b = (v_texCoord.y - u_Z) / u_R;",
                        "vec4 field;",
                        "if(a > 2.0 || b > 2.0){",
                            "field = u_I * vec4(sign(b), 1.0, 1.0, 1.0) * texture2D(u_shape_tenth, vec2(a / 10.0, abs(b) / 10.0));",
                        "}else{",
                            "field = u_I * vec4(sign(b), 1.0, 1.0, 1.0) * texture2D(u_shape_half, vec2(a / 2.0, abs(b) / 2.0));",
                        "}",
                        "gl_FragColor += field;",

                    "}"
                ];

                return src_arr.join('\n');
            })()
        });

        B_loop_half.texture.bind(programCurrentLoop, "u_shape_half");
        B_loop_tenth.texture.bind(programCurrentLoop, "u_shape_tenth");


        // ---------------------------------------------------------------------
        var programCurrentZ = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform float u_I;",

                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "gl_FragColor += vec4(0.0, u_I * 1.25663706e-6 / (2.0 * 3.14159265359 * v_texCoord.x), 0.0, 1.0);",

                    "}"
                ];

                return src_arr.join('\n');
            })()
        });

        // ---------------------------------------------------------------------
        var programBZ = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform float u_Bz;",

                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "gl_FragColor += vec4(0.0, 0.0, u_Bz, 1.0);",

                    "}"
                ];

                return src_arr.join('\n');
            })()
        });

        // ---------------------------------------------------------------------
        var programBTheta = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform float u_Btheta;",

                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "gl_FragColor += vec4(0.0, u_Btheta, 0.0, 1.0);",

                    "}"
                ];

                return src_arr.join('\n');
            })()
        });

        //
        // programs for pre-calculations of the fields
        //

        var R1 = webgl.addFrameBuffer(spec.nr, spec.nz, true);
        var R2 = webgl.addFrameBuffer(spec.nr, spec.nz, true);
        var R3 = webgl.addFrameBuffer(spec.nr, spec.nz, true);
        var A = webgl.addFrameBuffer(spec.nr, spec.nz, true);


        // ---------------------------------------------------------------------
        var programPre1 = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform float u_h;",

                    "uniform sampler2D u_B;",
                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "vec3 B = texture2D(u_B, v_texCoord).xyz;",
                        "float Bmag = length(B);",
                        "float hB2 = u_h * u_h * Bmag * Bmag;",
                        "float factor = 2.0 / (1.0 + hB2);",

                        "float R11 = (1.0 - hB2 * factor) + factor * u_h * u_h * B.x * B.x;",
                        "float R12 = factor * u_h * (B.z + u_h * B.x * B.y);",
                        "float R13 = factor * u_h * (-B.y + u_h * B.x * B.z);",
                        "vec3 R1 = vec3(R11, R12, R13 * " + N(factor_r/factor_z) + ");",
                        //"vec3 R1 = vec3(R11 , R12, R13);",

                        // R1
                        "gl_FragColor = vec4(R1, 1.0);",
                    "}"
                ];

                return src_arr.join('\n');
            })()
        });

        // triangle vertices
        vertex_positions.bind(programPre1, "a_position");
        // texture coordinets for vertices
        texture_coordinates.bind(programPre1, "a_texCoord");

        B.texture.bind(programPre1, "u_B");
        programPre1.setUniformFloat("u_h", h);

        // ---------------------------------------------------------------------
        var programPre2 = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform float u_h;",

                    "uniform sampler2D u_B;",
                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "vec3 B = texture2D(u_B, v_texCoord).xyz;",
                        "float Bmag = length(B);",
                        "float hB2 = u_h * u_h * Bmag * Bmag;",
                        "float factor = 2.0 / (1.0 + hB2);",

                        "float R21 = factor * u_h * (-B.z + u_h * B.y * B.x);",
                        "float R22 = (1.0 - hB2 * factor) + factor * u_h * u_h * B.y * B.y;",
                        "float R23 = factor * u_h * (B.x + u_h * B.y * B.z);",
                        "vec3 R2 = vec3(R21, R22, R23 * " + N(factor_r/factor_z) + ");",
                        //"vec3 R2 = vec3(R21 , R22, R23);",

                        // R2
                        "gl_FragColor = vec4(R2, 1.0);",
                    "}"
                ];

                return src_arr.join('\n');
            })()
        });

        // triangle vertices
        vertex_positions.bind(programPre2, "a_position");
        // texture coordinets for vertices
        texture_coordinates.bind(programPre2, "a_texCoord");

        B.texture.bind(programPre2, "u_B");
        programPre2.setUniformFloat("u_h", h);


        // ---------------------------------------------------------------------
        var programPre3 = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform float u_h;",

                    "uniform sampler2D u_B;",
                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "vec3 B = texture2D(u_B, v_texCoord).xyz;",
                        "float Bmag = length(B);",
                        "float hB2 = u_h * u_h * Bmag * Bmag;",
                        "float factor = 2.0 / (1.0 + hB2);",

                        "float R31 = factor * u_h * (B.y + u_h * B.z * B.x);",
                        "float R32 = factor * u_h * (-B.x + u_h * B.z * B.y);",
                        "float R33 = (1.0 - hB2 * factor) + factor * u_h * u_h * B.z * B.z;",
                        "vec3 R3 = vec3(R31 * " + N(factor_z/factor_r) + ", R32 * " + N(factor_z/factor_r) + ", R33);",
                        //"vec3 R3 = vec3(R31 , R32, R33);",

                        // R2
                        "gl_FragColor = vec4(R3, 1.0);",
                    "}"
                ];

                return src_arr.join('\n');
            })()
        });

        // triangle vertices
        vertex_positions.bind(programPre3, "a_position");
        // texture coordinets for vertices
        texture_coordinates.bind(programPre3, "a_texCoord");

        B.texture.bind(programPre3, "u_B");
        programPre3.setUniformFloat("u_h", h);


        // ---------------------------------------------------------------------
        var programPreA = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform float u_h;",

                    "uniform sampler2D u_E;",
                    "uniform sampler2D u_B;",

                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "vec3 B = texture2D(u_B, v_texCoord).xyz;",
                        "vec3 E = texture2D(u_E, v_texCoord).xyz;",
                        "float Bmag = length(B);",
                        "float hB2 = u_h * u_h * Bmag * Bmag;",
                        "float factor = 2.0/(1.0 + hB2);",
                        "vec3 A = (u_h * (2.0 - hB2 * factor) * E + u_h * u_h * factor * (cross(E,B) + u_h * dot(E,B)))/2.998e8;",
                        // divide by speed of light because of normalized velocity
                        "gl_FragColor = vec4(A.x * " + N(factor_r) + ", A.y * " + N(factor_r) + ", A.z * " + N(factor_z) + ", 1.0);",
                    "}"
                ];

                return src_arr.join('\n');
            })()
        });

        // triangle vertices
        vertex_positions.bind(programPreA, "a_position");
        // texture coordinets for vertices
        texture_coordinates.bind(programPreA, "a_texCoord");

        E.texture.bind(programPreA, "u_E");
        B.texture.bind(programPreA, "u_B");
        programPreA.setUniformFloat("u_h", h);


        //
        // step programs in leap-frog
        //

        var position_A = webgl.addFrameBuffer(nparticles_w, nparticles_h, true);
        var velocity_A = webgl.addFrameBuffer(nparticles_w, nparticles_h, true);
        var position_B = webgl.addFrameBuffer(nparticles_w, nparticles_h, true);
        var velocity_B = webgl.addFrameBuffer(nparticles_w, nparticles_h, true);

        var step_vert = function() {
            var src_arr = [
                "attribute vec2 a_position;",
                "attribute vec2 a_texCoord;",
                "varying vec2 v_texCoord;",

                "void main() {",
                "    gl_Position = vec4(a_position, 0, 1);",

                "    v_texCoord = a_texCoord;",
                "}",
            ];

            return src_arr.join('\n');
        } // step_vert()

        var step_position_frag = function() {
            var src_arr = [
                "precision highp float;",

                "uniform float u_step_factor;",

                // current phase
                "uniform sampler2D u_position;",
                "uniform sampler2D u_velocity;",

                // the texCoords passed in from the vertex shader.
                "varying vec2 v_texCoord;",

                "void main() {",

                    "vec4 next_position = texture2D(u_position, v_texCoord) + u_step_factor * texture2D(u_velocity, v_texCoord);",
                    "if (next_position.x * next_position.x + next_position.y + next_position.y > 1.0 || next_position.z > 1.0 || next_position.z < 0.0) {",
                        //"gl_FragColor = vec4(0.0, 0.0, 0.5, 1.0);",
                        "gl_FragColor = next_position;",
                    "}else{",
                        "gl_FragColor = next_position;",
                    "}",
                "}"
            ];


            return src_arr.join('\n');
        }; // step_position_frag()

        var step_velocity_frag = function() {
            var src_arr = [
                "precision highp float;",

                // current phase
                "uniform sampler2D u_position;",
                "uniform sampler2D u_velocity;",

                // fields
                "uniform sampler2D u_R_1;",
                "uniform sampler2D u_R_2;",
                "uniform sampler2D u_R_3;",
                "uniform sampler2D u_A;",

                // the texCoords passed in from the vertex shader.
                "varying vec2 v_texCoord;",

                "void main() {",
                    "vec3 position = texture2D(u_position, v_texCoord).xyz;",
                    "vec3 velocity = texture2D(u_velocity, v_texCoord).xyz;",

                    "float r = sqrt(position.x * position.x + position.y * position.y);",
                    "vec2 direction = vec2(position.x/r, position.y/r);",

                    "float vr = velocity.x * direction.x + velocity.y * direction.y;",
                    "float va = velocity.y * direction.x - velocity.x * direction.y;",
                    "vec2 cylindrical_position = vec2(r, position.z);",
                    "vec3 cylindrical_velocity = vec3(vr, va, velocity.z);",

                    "vec3 R1 = texture2D(u_R_1, cylindrical_position).xyz;",
                    "vec3 R2 = texture2D(u_R_2, cylindrical_position).xyz;",
                    "vec3 R3 = texture2D(u_R_3, cylindrical_position).xyz;",
                    "vec3 A = texture2D(u_A, cylindrical_position).xyz;",

                    "cylindrical_velocity = vec3(dot(R1, cylindrical_velocity), dot(R2, cylindrical_velocity), dot(R3, cylindrical_velocity)) + A;",
                    "vec3 next_velocity = vec3(cylindrical_velocity.x * direction.x - cylindrical_velocity.y * direction.y, cylindrical_velocity.x * direction.y + cylindrical_velocity.y * direction.x, cylindrical_velocity.z);",

                    "gl_FragColor = vec4(next_velocity, 1.0);",
                "}"
            ];


            return src_arr.join('\n');
        }; // step_velocity_frag()



        // ---------------------------------------------------------------------
        // computes value of velocity_B
        var programStepVelocityB = webgl.linkProgram({
            vertexShaderSource : step_vert(),
            fragmentShaderSource : step_velocity_frag()
        });

        // triangle vertices
        vertex_positions.bind(programStepVelocityB, "a_position");
        // texture coordinets for vertices
        texture_coordinates.bind(programStepVelocityB, "a_texCoord");

        position_A.texture.bind(programStepVelocityB, "u_position");
        velocity_A.texture.bind(programStepVelocityB, "u_velocity");
        R1.texture.bind(programStepVelocityB, "u_R_1");
        R2.texture.bind(programStepVelocityB, "u_R_2");
        R3.texture.bind(programStepVelocityB, "u_R_3");
        A.texture.bind(programStepVelocityB, "u_A");

        // ---------------------------------------------------------------------
        // computes value of position_B
        var programStepPositionB = webgl.linkProgram({
            vertexShaderSource : step_vert(),
            fragmentShaderSource : step_position_frag()
        });

        // triangle vertices
        vertex_positions.bind(programStepPositionB, "a_position");
        // texture coordinets for vertices
        texture_coordinates.bind(programStepPositionB, "a_texCoord");

        programStepPositionB.setUniformFloat("u_step_factor", spec.dt * speed_of_light);

        position_A.texture.bind(programStepPositionB, "u_position");
        velocity_B.texture.bind(programStepPositionB, "u_velocity");

        // ---------------------------------------------------------------------
        // computes value of velocity_A
        var programStepVelocityA = webgl.linkProgram({
            vertexShaderSource : step_vert(),
            fragmentShaderSource : step_velocity_frag()
        });

        // triangle vertices
        vertex_positions.bind(programStepVelocityA, "a_position");
        // texture coordinets for vertices
        texture_coordinates.bind(programStepVelocityA, "a_texCoord");

        position_B.texture.bind(programStepVelocityA, "u_position");
        velocity_B.texture.bind(programStepVelocityA, "u_velocity");
        R1.texture.bind(programStepVelocityA, "u_R_1");
        R2.texture.bind(programStepVelocityA, "u_R_2");
        R3.texture.bind(programStepVelocityA, "u_R_3");
        A.texture.bind(programStepVelocityA, "u_A");

        // ---------------------------------------------------------------------
        // computes value of position_A
        var programStepPositionA = webgl.linkProgram({
            vertexShaderSource : step_vert(),
            fragmentShaderSource : step_position_frag()
        });

        // triangle vertices
        vertex_positions.bind(programStepPositionA, "a_position");
        // texture coordinets for vertices
        texture_coordinates.bind(programStepPositionA, "a_texCoord");

        programStepPositionA.setUniformFloat("u_step_factor", spec.dt * speed_of_light);

        position_B.texture.bind(programStepPositionA, "u_position");
        velocity_A.texture.bind(programStepPositionA, "u_velocity");

        // ---------------------------------------------------------------------
        // program for rendering particles into a density function
        var programRenderDensity = webgl.linkProgram({
            vertexShaderSource : (function() {
                var src_arr = [
                    "attribute vec2 a_particleTexCoord;",

                    "uniform sampler2D u_position;",
                    "uniform sampler2D u_velocity;",

                    "uniform float u_pointsize;",

                    "varying vec4 v_color;",

                    "void main() {",
                        "vec3 position = texture2D(u_position, a_particleTexCoord).xyz;",
                        "vec3 velocity = texture2D(u_velocity, a_particleTexCoord).xyz;",
                        "float r = sqrt(position.x * position.x + position.y * position.y);",
                        "gl_Position = vec4(2.0*r-1.0, 2.0*position.z-1.0, 0, 1.0);",
                        //"gl_Position = vec4(position.x, position.y, 0, 1.0);",
                        "gl_PointSize = u_pointsize;",

                        "vec2 direction = vec2(position.x/r, position.y/r);",

                        "float v = length(velocity);",
                        "float va = velocity.y * direction.x - velocity.x * direction.y;",
                        "v_color = vec4((va < 0.0) ? (1.0 + va / v) : 1.0, (1.0 + va / v) * (1.0 - va / v), (va > 0.0) ? (1.0 - va / v) : 1.0, 1.0);",
                    "}"
                ];

                return src_arr.join('\n');
            })(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform float u_weight;",
                    "uniform sampler2D u_shape;",

                    "varying vec4 v_color;",

                    "void main() {",
                        "gl_FragColor += u_weight * v_color * texture2D(u_shape, gl_PointCoord);",
                    "}"
                ];


                return src_arr.join('\n');
            })()
        });

        // texture coordinates of particles
        var particles_texcoord_arr = [];

        for(j = 0; j < nparticles_h; j++) {
            for(i = 0; i < nparticles_w; i++) {
                particles_texcoord_arr[i + nparticles_w * j] = [(i+0.5)/(nparticles_w), (j+0.5)/(nparticles_h)];
            }
        }


        var particles_texcoord = webgl.addVertexData(particles_texcoord_arr);
        particles_texcoord.bind(programRenderDensity, "a_particleTexCoord");


        position_A.texture.bind(programRenderDensity, "u_position");
        velocity_A.texture.bind(programRenderDensity, "u_velocity");

        programRenderDensity.setUniformFloat("u_weight", 0.01);

        var nshape = 11;
        var shape_arr = new Float32Array(4 * nshape * nshape);

        var mid = (nshape-1)/2;

        for(j = 0; j < nshape; j++) {
            for(i = 0; i < nshape; i++) {
                var d = Math.sqrt(Math.pow(i - mid, 2) + Math.pow(j - mid, 2));
                shape_arr[4*(i + nshape*j)] = Math.pow(Math.max(0.0, Math.cos(0.5*Math.PI*d/mid)), 2);
                shape_arr[4*(i + nshape*j)+1] = shape_arr[4*(i + nshape*j)];
                shape_arr[4*(i + nshape*j)+2] = shape_arr[4*(i + nshape*j)];
                shape_arr[4*(i + nshape*j)+3] = 1.0;
            }
        }

        var shape_tex = webgl.addTextureArray(
            nshape,
            nshape,
            shape_arr,
            true
        );

        shape_tex.bind(programRenderDensity, "u_shape");
        programRenderDensity.setUniformFloat("u_pointsize", nshape);

        // ---------------------------------------------------------------------
        // program for setting value of particles
        //
        var programSet = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision mediump float;",

                    "uniform sampler2D u_value;",
                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "gl_FragColor = texture2D(u_value, v_texCoord);",
                    "}"
                ];

                return src_arr.join('\n');
            })()
        });

        // triangle vertices
        vertex_positions.bind(programSet, "a_position");
        // texture coordinets for vertices
        texture_coordinates.bind(programSet, "a_texCoord");

        //
        // Output interface for simulation programs
        //

        out.set = function(value) {
            if (value.E) {
                for(i = 0; i < spec.nr; i++) {
                    for(j = 0; j < spec.nz; j++){
                        E_arr[4*(i + j*spec.nr)] = value.E[i][j][0];
                        E_arr[4*(i + j*spec.nr)+1] = value.E[i][j][1];
                        E_arr[4*(i + j*spec.nr)+2] = value.E[i][j][2];
                        E_arr[4*(i + j*spec.nr)+3] = 1.0;
                    }
                }

                E_tex.update();
                E_tex.bind(programSet, "u_value");

                programSet.drawTriangles(0, 6, E);
            }

            if (value.B) {
                for(i = 0; i < spec.nr; i++) {
                    for(j = 0; j < spec.nz; j++){
                        B_arr[4*(i + j*spec.nr)] = value.B[i][j][0];
                        B_arr[4*(i + j*spec.nr)+1] = value.B[i][j][1];
                        B_arr[4*(i + j*spec.nr)+2] = value.B[i][j][2];
                        B_arr[4*(i + j*spec.nr)+3] = 1.0;
                    }
                }

                B_tex.update();
                B_tex.bind(programSet, "u_value");
                programSet.drawTriangles(0, 6, B);
            }

            if (value.position) {

                for(i = 0; i < nparticles; i++) {
                    position_arr[4*i] = value.position[i][0] * factor_r;
                    position_arr[4*i+1] = value.position[i][1] * factor_r;
                    position_arr[4*i+2] = value.position[i][2] * factor_z;
                    position_arr[4*i+3] = 1.0;
                }

                // send values to card
                position_tex.update();
                position_tex.bind(programSet, "u_value");
                programSet.drawTriangles(0, 6, position_A);
                programSet.drawTriangles(0, 6, position_B);
            }

            if (value.velocity) {

                for(i = 0; i < nparticles; i++) {
                    velocity_arr[4*i] = value.velocity[i][0] * factor_r;
                    velocity_arr[4*i+1] = value.velocity[i][1] * factor_r;
                    velocity_arr[4*i+2] = value.velocity[i][2] * factor_z;
                    velocity_arr[4*i+3] = 1.0;
                }

                // send values to card
                velocity_tex.update();
                velocity_tex.bind(programSet, "u_value");
                programSet.drawTriangles(0, 6, velocity_A);
                programSet.drawTriangles(0, 6, velocity_B);
            }
        };

        out.addCurrentLoop = function(r, z, I) {
            programCurrentLoop.setUniformFloat("u_R", r * factor_r);
            programCurrentLoop.setUniformFloat("u_Z", z * factor_z);
            programCurrentLoop.setUniformFloat("u_I", I);
            programCurrentLoop.drawTriangles(0, 6, B);
        };

        out.addCurrentZ = function(I) {
            programCurrentZ.setUniformFloat("u_I", I);
            programCurrentZ.drawTriangles(0, 6, B);
        };

        out.addBZ = function(Bz) {
            programBZ.setUniformFloat("u_Bz", Bz);
            programBZ.drawTriangles(0, 6, B);
        };

        out.addBTheta = function(Btheta) {
            programBTheta.setUniformFloat("u_Btheta", Btheta);
            programBTheta.drawTriangles(0, 6, B);
        };

        out.precalc = function () {
            R1.clear(0.0, 0.0, 0.0, 1.0);
            R2.clear(0.0, 0.0, 0.0, 1.0);
            R3.clear(0.0, 0.0, 0.0, 1.0);
            A.clear(0.0, 0.0, 0.0, 1.0);

            programPre1.drawTriangles(0, 6, R1);
            programPre2.drawTriangles(0, 6, R2);
            programPre3.drawTriangles(0, 6, R3);
            programPreA.drawTriangles(0, 6, A);
        };

        out.step = function () {
            // have to clear one at a time otherwise data is wiped out
            velocity_B.clear(0.0, 0.0, 0.0, 1.0);
            programStepVelocityB.drawTriangles(0, 6, velocity_B);

            position_B.clear(0.0, 0.0, 0.0, 1.0);
            programStepPositionB.drawTriangles(0, 6, position_B);

            velocity_A.clear(0.0, 0.0, 0.0, 1.0);
            programStepVelocityA.drawTriangles(0, 6, velocity_A);

            position_A.clear(0.0, 0.0, 0.0, 1.0);
            programStepPositionA.drawTriangles(0, 6, position_A);
        };

        out.density = function() {
            webgl.clear();
            programRenderDensity.drawPoints(0, nparticles);

            //B.texture.bind(programSet, "u_value");
            //programSet.drawTriangles(0, 6);

        };

        return out;
    };

    return exports;
});
