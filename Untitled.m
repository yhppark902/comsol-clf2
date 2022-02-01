function out = model
%
% Untitled.m
%
% Model exported on Feb 1 2022, 10:08 by COMSOL 6.0.0.318.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('C:\usr');

model.param.set('M0', '0[N*m]');
model.param.set('M1', '0[N*m]');
model.param.set('dtFem', '0.02[s]');
model.param.set('phi', '0[rad]');
model.param.set('phi_t', '0[rad/s]');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);
model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('r1').set('pos', [-0.025 0]);
model.component('comp1').geom('geom1').feature('r1').set('size', [0.05 1]);
model.component('comp1').geom('geom1').create('pt1', 'Point');
model.component('comp1').geom('geom1').create('r2', 'Rectangle');
model.component('comp1').geom('geom1').feature('r2').set('pos', [-0.025 1]);
model.component('comp1').geom('geom1').feature('r2').set('size', [0.05 1]);
model.component('comp1').geom('geom1').create('pt2', 'Point');
model.component('comp1').geom('geom1').feature('pt2').set('p', [0 1]);
model.component('comp1').geom('geom1').run;

model.component('comp1').physics.create('mbd', 'MultibodyDynamics', 'geom1');
model.component('comp1').physics('mbd').create('gr1', 'Gravity', 2);
model.component('comp1').physics('mbd').feature('gr1').selection.set([1 2]);
model.component('comp1').physics('mbd').create('rd1', 'RigidDomain', 2);
model.component('comp1').physics('mbd').feature('rd1').selection.set([1]);
model.component('comp1').physics('mbd').feature('rd1').set('CenterOfRotationType', 'CentroidOfSelectedEntities');
model.component('comp1').physics('mbd').feature('rd1').feature('crp1').selection.set([4]);
model.component('comp1').physics('mbd').create('rd2', 'RigidDomain', 2);
model.component('comp1').physics('mbd').feature('rd2').selection.set([2]);
model.component('comp1').physics('mbd').feature('rd2').set('CenterOfRotationType', 'CentroidOfSelectedEntities');
model.component('comp1').physics('mbd').feature('rd2').feature('crp1').selection.set([5]);
model.component('comp1').physics('mbd').create('hgj1', 'HingeJoint', -1);
model.component('comp1').physics('mbd').feature('hgj1').set('Source', 'fixed');
model.component('comp1').physics('mbd').feature('hgj1').set('Destination', 'rd1');
model.component('comp1').physics('mbd').feature('hgj1').create('afm1', 'AppliedForceAndMoment', -1);
model.component('comp1').physics('mbd').create('hgj2', 'HingeJoint', -1);
model.component('comp1').physics('mbd').feature('hgj2').set('Source', 'rd1');
model.component('comp1').physics('mbd').feature('hgj2').set('Destination', 'rd2');
model.component('comp1').physics('mbd').feature('hgj2').create('afm1', 'AppliedForceAndMoment', -1);

model.component('comp1').mesh.create('mesh1');
model.component('comp1').mesh('mesh1').create('map1', 'Map');

model.component('comp1').physics('mbd').feature('rd1').set('EntityLevel', 'Point');
model.component('comp1').physics('mbd').feature('rd2').set('EntityLevel', 'Point');

out = model;
