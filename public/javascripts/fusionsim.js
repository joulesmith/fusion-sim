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



define([
    'angular',
    'utilities',
    'empic_webgl'
], function (angular, utilities, empic_webgl){
    "use strict";

    var app = angular.module('fusionsim', []);

    app.controller('simulation', [
        '$scope',
        '$interval',
        function($scope, $interval){

            var nparticles = 160000;

            var simulation = empic_webgl.makeCylindricalParticlePusher({
                radius : 1,
                height : 2,
                nr : 400,
                nz : 800,
                dt : 2E-8,
                nparticles : 400,
                particle_mass : 1.67e-27, //kg : proton
                particle_charge : 1.602e-19 // C
            });

            var init_position = [];
            var init_velocity = [];

            for(var i = 0; i < nparticles; i++) {
                init_position[i] = [0.2 * Math.random() + 0.4, 0, 0.2*Math.random() + 0.9];
                init_velocity[i] = [0.001 * (Math.random()-0.5), 0.001 * (Math.random()-0.5), 0.001 * (Math.random()-0.5)];
            }



            simulation.set({
                position : init_position,
                velocity : init_velocity,
            });

            //simulation.addCurrentLoop(1.0, 2.0, -10000000);
            //simulation.addCurrentLoop(1.0, 0.0, 10000000);

            simulation.addCurrentLoop(0.5, 1.0, 10000000);
            simulation.addCurrentZ(5000000);
            simulation.addBZ(0.01);

            //simulation.addBTheta(0.01);

            simulation.precalc();

            var plot = document.getElementById("plot");
            var plot_ctx = plot.getContext("2d");

            simulation.density();
            plot_ctx.drawImage(simulation.canvas, 0, 0);

            $scope.run = false;

            console.log("initialized");

            window.requestAnimationFrame = window.requestAnimationFrame || window.mozRequestAnimationFrame ||
                    window.webkitRequestAnimationFrame || window.oRequestAnimationFrame;


            $scope.start = function() {
                console.log("started");

                $scope.run = true;

                var step = function(){

                    simulation.step();

                    simulation.density();

                    plot_ctx.drawImage(simulation.canvas, 0, 0);

                    if ($scope.run) {
                        requestAnimationFrame(step);
                    }
                };

                requestAnimationFrame(step);
            };

            $scope.stop = function() {
                console.log("stopped");
                $scope.run = false;
            };
	}]);

    app.controller('MainCtrl', [
        '$scope',
        '$interval',
        function($scope, $interval){

	}]);

});
