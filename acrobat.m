function out = model
%
% acrobat.m
%
% Model exported on Feb 1 2022, 09:48 by COMSOL 6.0.0.318.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('C:\usr\repo - Copy');

model.label('acrobat.mph');

model.param.set('M0', '10[N*m]');
model.param.set('M1', '10[N*m]');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);

model.result.table.create('tbl1', 'Table');
model.result.table.create('tbl2', 'Table');

model.component('comp1').mesh.create('mesh1');

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

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material('mat1').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
model.component('comp1').material('mat1').propertyGroup.create('Murnaghan', 'Murnaghan');
model.component('comp1').material('mat1').propertyGroup.create('Lame', ['Lam' native2unicode(hex2dec({'00' 'e9'}), 'unicode') ' parameters']);

model.component('comp1').physics.create('mbd', 'MultibodyDynamics', 'geom1');
model.component('comp1').physics('mbd').create('gr1', 'Gravity', 2);
model.component('comp1').physics('mbd').feature('gr1').selection.set([1]);
model.component('comp1').physics('mbd').create('rd1', 'RigidDomain', 2);
model.component('comp1').physics('mbd').feature('rd1').selection.set([1]);
model.component('comp1').physics('mbd').feature('rd1').feature('crp1').selection.set([4]);
model.component('comp1').physics('mbd').create('rd2', 'RigidDomain', 2);
model.component('comp1').physics('mbd').feature('rd2').selection.set([2]);
model.component('comp1').physics('mbd').feature('rd2').feature('crp1').selection.set([5]);
model.component('comp1').physics('mbd').create('hgj1', 'HingeJoint', -1);
model.component('comp1').physics('mbd').feature('hgj1').set('Source', 'fixed');
model.component('comp1').physics('mbd').feature('hgj1').set('Destination', 'rd1');
model.component('comp1').physics('mbd').feature('hgj1').create('afm1', 'AppliedForceAndMoment', -1);
model.component('comp1').physics('mbd').create('hgj2', 'HingeJoint', -1);
model.component('comp1').physics('mbd').feature('hgj2').set('Source', 'rd1');
model.component('comp1').physics('mbd').feature('hgj2').set('Destination', 'rd2');
model.component('comp1').physics('mbd').feature('hgj2').create('afm1', 'AppliedForceAndMoment', -1);

model.component('comp1').probe.create('var1', 'GlobalVariable');
model.component('comp1').probe.create('var2', 'GlobalVariable');
model.component('comp1').probe.create('var3', 'GlobalVariable');
model.component('comp1').probe.create('var4', 'GlobalVariable');

model.result.table('tbl1').label('Probe Table 1');
model.result.table('tbl2').comments('System Matrix 1');

model.component('comp1').view('view1').axis.set('xmin', -1.8480901718139648);
model.component('comp1').view('view1').axis.set('xmax', 2.1926355361938477);
model.component('comp1').view('view1').axis.set('ymin', -2.212296485900879);
model.component('comp1').view('view1').axis.set('ymax', 1.8736404180526733);

model.component('comp1').material('mat1').label('Structural steel');
model.component('comp1').material('mat1').set('family', 'custom');
model.component('comp1').material('mat1').set('customspecular', [0.7843137254901961 0.7843137254901961 0.7843137254901961]);
model.component('comp1').material('mat1').set('customdiffuse', [0.6666666666666666 0.6666666666666666 0.6666666666666666]);
model.component('comp1').material('mat1').set('customambient', [0.6666666666666666 0.6666666666666666 0.6666666666666666]);
model.component('comp1').material('mat1').set('fresnel', 0.9);
model.component('comp1').material('mat1').set('roughness', 0.3);
model.component('comp1').material('mat1').set('metallic', 0);
model.component('comp1').material('mat1').set('pearl', 0);
model.component('comp1').material('mat1').set('diffusewrap', 0);
model.component('comp1').material('mat1').set('clearcoat', 0);
model.component('comp1').material('mat1').set('reflectance', 0);
model.component('comp1').material('mat1').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', '475[J/(kg*K)]');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', {'44.5[W/(m*K)]' '0' '0' '0' '44.5[W/(m*K)]' '0' '0' '0' '44.5[W/(m*K)]'});
model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'4.032e6[S/m]' '0' '0' '0' '4.032e6[S/m]' '0' '0' '0' '4.032e6[S/m]'});
model.component('comp1').material('mat1').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', {'12.3e-6[1/K]' '0' '0' '0' '12.3e-6[1/K]' '0' '0' '0' '12.3e-6[1/K]'});
model.component('comp1').material('mat1').propertyGroup('def').set('density', '7850[kg/m^3]');
model.component('comp1').material('mat1').propertyGroup('Enu').set('E', '');
model.component('comp1').material('mat1').propertyGroup('Enu').set('nu', '');
model.component('comp1').material('mat1').propertyGroup('Enu').set('E', '');
model.component('comp1').material('mat1').propertyGroup('Enu').set('nu', '');
model.component('comp1').material('mat1').propertyGroup('Enu').set('E', '');
model.component('comp1').material('mat1').propertyGroup('Enu').set('nu', '');
model.component('comp1').material('mat1').propertyGroup('Enu').set('E', '200e9[Pa]');
model.component('comp1').material('mat1').propertyGroup('Enu').set('nu', '0.30');
model.component('comp1').material('mat1').propertyGroup('Murnaghan').set('l', '');
model.component('comp1').material('mat1').propertyGroup('Murnaghan').set('m', '');
model.component('comp1').material('mat1').propertyGroup('Murnaghan').set('n', '');
model.component('comp1').material('mat1').propertyGroup('Murnaghan').set('l', '');
model.component('comp1').material('mat1').propertyGroup('Murnaghan').set('m', '');
model.component('comp1').material('mat1').propertyGroup('Murnaghan').set('n', '');
model.component('comp1').material('mat1').propertyGroup('Murnaghan').set('l', '');
model.component('comp1').material('mat1').propertyGroup('Murnaghan').set('m', '');
model.component('comp1').material('mat1').propertyGroup('Murnaghan').set('n', '');
model.component('comp1').material('mat1').propertyGroup('Murnaghan').set('l', '-3.0e11[Pa]');
model.component('comp1').material('mat1').propertyGroup('Murnaghan').set('m', '-6.2e11[Pa]');
model.component('comp1').material('mat1').propertyGroup('Murnaghan').set('n', '-7.2e11[Pa]');
model.component('comp1').material('mat1').propertyGroup('Lame').set('lambLame', '');
model.component('comp1').material('mat1').propertyGroup('Lame').set('muLame', '');
model.component('comp1').material('mat1').propertyGroup('Lame').set('lambLame', '');
model.component('comp1').material('mat1').propertyGroup('Lame').set('muLame', '');
model.component('comp1').material('mat1').propertyGroup('Lame').set('lambLame', '');
model.component('comp1').material('mat1').propertyGroup('Lame').set('muLame', '');
model.component('comp1').material('mat1').propertyGroup('Lame').set('lambLame', '1.15e11[Pa]');
model.component('comp1').material('mat1').propertyGroup('Lame').set('muLame', '7.69e10[Pa]');

model.component('comp1').physics('mbd').prop('PhysicsSymbols').set('PhysicsSymbols', true);
model.component('comp1').physics('mbd').feature('rd1').set('CenterOfRotationType', 'CentroidOfSelectedEntities');
model.component('comp1').physics('mbd').feature('rd1').set('EntityLevel', 'Point');
model.component('comp1').physics('mbd').feature('rd2').set('CenterOfRotationType', 'CentroidOfSelectedEntities');
model.component('comp1').physics('mbd').feature('rd2').set('EntityLevel', 'Point');
model.component('comp1').physics('mbd').feature('hgj1').set('CenterOfJointType', 'UserDefined');
model.component('comp1').physics('mbd').feature('hgj1').feature('afm1').set('AppliedOnSelectedAttachment', 'Joint');
model.component('comp1').physics('mbd').feature('hgj1').feature('afm1').set('Ms', 'M0');
model.component('comp1').physics('mbd').feature('hgj2').set('CenterOfJointType', 'UserDefined');
model.component('comp1').physics('mbd').feature('hgj2').set('xc', [0; 1; 0]);
model.component('comp1').physics('mbd').feature('hgj2').feature('afm1').set('Mz', 'M1');

model.component('comp1').probe('var1').set('expr', 'mbd.rd1.thz');
model.component('comp1').probe('var1').set('descr', 'Rigid body rotation, z component');
model.component('comp1').probe('var1').set('table', 'tbl1');
model.component('comp1').probe('var1').set('window', 'window1');
model.component('comp1').probe('var2').set('expr', 'mbd.rd1.th_tz');
model.component('comp1').probe('var2').set('unit', 'rad/s');
model.component('comp1').probe('var2').set('descr', 'Rigid body angular velocity, z component');
model.component('comp1').probe('var2').set('table', 'tbl1');
model.component('comp1').probe('var2').set('window', 'window1');
model.component('comp1').probe('var3').set('expr', 'mbd.rd2.thz');
model.component('comp1').probe('var3').set('descr', 'Rigid body rotation, z component');
model.component('comp1').probe('var3').set('table', 'tbl1');
model.component('comp1').probe('var3').set('window', 'window1');
model.component('comp1').probe('var4').set('expr', 'mbd.rd2.th_tz');
model.component('comp1').probe('var4').set('unit', 'rad/s');
model.component('comp1').probe('var4').set('descr', 'Rigid body angular velocity, z component');
model.component('comp1').probe('var4').set('table', 'tbl1');
model.component('comp1').probe('var4').set('window', 'window1');

model.study.create('std1');
model.study('std1').create('time', 'Transient');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('t1', 'Time');
model.sol('sol1').create('sp1', 'StateSpace');
model.sol('sol1').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').feature.remove('fcDef');

model.result.dataset.create('dset2', 'Solution');
model.result.dataset('dset2').set('probetag', 'var4');
model.result.numerical.create('gev1', 'EvalGlobal');
model.result.numerical.create('sys1', 'SystemMatrix');
model.result.numerical.create('gev2', 'EvalGlobal');
model.result.numerical.create('gev3', 'EvalGlobal');
model.result.numerical.create('gev4', 'EvalGlobal');
model.result.numerical('gev1').set('data', 'dset2');
model.result.numerical('gev1').set('probetag', 'var1');
model.result.numerical('gev2').set('data', 'dset2');
model.result.numerical('gev2').set('probetag', 'var2');
model.result.numerical('gev3').set('data', 'dset2');
model.result.numerical('gev3').set('probetag', 'var3');
model.result.numerical('gev4').set('data', 'dset2');
model.result.numerical('gev4').set('probetag', 'var4');
model.result.create('pg1', 'PlotGroup2D');
model.result.create('pg2', 'PlotGroup2D');
model.result.create('pg3', 'PlotGroup2D');
model.result.create('pg4', 'PlotGroup1D');
model.result('pg1').create('surf1', 'Surface');
model.result('pg1').feature('surf1').create('def1', 'Deform');
model.result('pg2').create('surf1', 'Surface');
model.result('pg2').create('arwl1', 'ArrowLine');
model.result('pg2').feature('surf1').set('expr', 'mod(dom,10)');
model.result('pg2').feature('surf1').create('def1', 'Deform');
model.result('pg2').feature('arwl1').create('def1', 'Deform');
model.result('pg3').create('arws1', 'ArrowSurface');
model.result('pg3').create('surf1', 'Surface');
model.result('pg3').feature('arws1').create('col', 'Color');
model.result('pg3').feature('arws1').create('def', 'Deform');
model.result('pg3').feature('arws1').feature('col').set('expr', 'mbd.gr1.F_V_Mag');
model.result('pg3').feature('surf1').set('expr', '1');
model.result('pg3').feature('surf1').create('def', 'Deform');
model.result('pg3').feature('surf1').create('sel1', 'Selection');
model.result('pg3').feature('surf1').feature('sel1').selection.set([1]);
model.result('pg4').set('probetag', 'window1_default');
model.result('pg4').create('tblp1', 'Table');
model.result('pg4').feature('tblp1').set('probetag', 'var1,var2,var3,var4');
model.result.export.create('anim1', 'Animation');

model.component('comp1').probe('var1').genResult([]);
model.component('comp1').probe('var2').genResult([]);
model.component('comp1').probe('var3').genResult([]);
model.component('comp1').probe('var4').genResult([]);

model.nodeGroup.create('dset1mbdlgrp', 'Results');
model.nodeGroup('dset1mbdlgrp').set('type', 'plotgroup');
model.nodeGroup('dset1mbdlgrp').placeAfter('plotgroup', 'pg2');

model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').label('Compile Equations: Time Dependent');
model.sol('sol1').feature('v1').label('Dependent Variables 1.1');
model.sol('sol1').feature('v1').feature('comp1_u').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_u').set('scaleval', 1);
model.sol('sol1').feature('v1').feature('comp1_mbd_rd1_u').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_mbd_rd1_u').set('scaleval', 1);
model.sol('sol1').feature('v1').feature('comp1_mbd_rd1_phi').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_mbd_rd1_phi').set('scaleval', 1);
model.sol('sol1').feature('v1').feature('comp1_mbd_hgj1_th').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_mbd_hgj1_th').set('scaleval', 1);
model.sol('sol1').feature('v1').feature('comp1_mbd_rd2_u').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_mbd_rd2_phi').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_mbd_hgj2_th').set('scalemethod', 'manual');
model.sol('sol1').feature('t1').label('Time-Dependent Solver 1.1');
model.sol('sol1').feature('t1').set('rtol', 0.001);
model.sol('sol1').feature('t1').set('tstepsbdf', 'strict');
model.sol('sol1').feature('t1').set('maxorder', 2);
model.sol('sol1').feature('t1').feature('dDef').label('Direct 1');
model.sol('sol1').feature('t1').feature('aDef').label('Advanced 1');
model.sol('sol1').feature('t1').feature('aDef').set('cachepattern', true);
model.sol('sol1').feature('t1').feature('fc1').label('Fully Coupled 1.1');
model.sol('sol1').feature('t1').feature('fc1').set('maxiter', 50);
model.sol('sol1').feature('t1').feature('fc1').set('jtech', 'onevery');
model.sol('sol1').feature('sp1').label('State Space 1.1');
model.sol('sol1').feature('sp1').set('input', {'M0' 'M1'});
model.sol('sol1').feature('sp1').set('output', {'comp1.var1' 'comp1.var2' 'comp1.var3' 'comp1.var4'});
model.sol('sol1').feature('sp1').set('MA', true);
model.sol('sol1').feature('sp1').set('MB', true);
model.sol('sol1').feature('sp1').set('D', true);
model.sol('sol1').feature('sp1').set('C', true);
model.sol('sol1').feature('sp1').set('static', false);
model.sol('sol1').feature('sp1').set('Mc', true);
model.sol('sol1').feature('sp1').set('Null', true);
model.sol('sol1').feature('sp1').set('ud', true);
model.sol('sol1').feature('sp1').set('x0', true);
model.sol('sol1').feature('sp1').feature('aDef').label('Advanced 1');
model.sol('sol1').feature('sp1').feature('aDef').set('complexfun', true);
model.sol('sol1').feature('sp1').feature('aDef').set('checkmatherr', true);
model.sol('sol1').runAll;

model.result.dataset('dset2').label('Probe Solution 2');
model.result.numerical('gev1').set('const', {'mbd.refpntx' '0' 'Reference point for moment computation, x component'; 'mbd.refpnty' '0' 'Reference point for moment computation, y component'; 'mbd.refpntz' '0' 'Reference point for moment computation, z component'});
model.result.numerical('sys1').set('table', 'tbl2');
model.result.numerical('sys1').set('matrixstatespace', 'ssc');
model.result.numerical('sys1').set('format', 'filled');
model.result.numerical('gev2').set('descr', {''});
model.result.numerical('gev2').set('const', {'mbd.refpntx' '0' 'Reference point for moment computation, x component'; 'mbd.refpnty' '0' 'Reference point for moment computation, y component'; 'mbd.refpntz' '0' 'Reference point for moment computation, z component'});
model.result.numerical('gev3').set('const', {'mbd.refpntx' '0' 'Reference point for moment computation, x component'; 'mbd.refpnty' '0' 'Reference point for moment computation, y component'; 'mbd.refpntz' '0' 'Reference point for moment computation, z component'});
model.result.numerical('gev4').set('const', {'mbd.refpntx' '0' 'Reference point for moment computation, x component'; 'mbd.refpnty' '0' 'Reference point for moment computation, y component'; 'mbd.refpntz' '0' 'Reference point for moment computation, z component'});
model.result.numerical('gev1').setResult;
model.result.numerical('sys1').setResult;
model.result('pg1').label('Displacement (mbd)');
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').feature('surf1').label('Surface');
model.result('pg1').feature('surf1').set('const', {'mbd.refpntx' '0' 'Reference point for moment computation, x component'; 'mbd.refpnty' '0' 'Reference point for moment computation, y component'; 'mbd.refpntz' '0' 'Reference point for moment computation, z component'});
model.result('pg1').feature('surf1').set('colortable', 'SpectrumLight');
model.result('pg1').feature('surf1').set('resolution', 'normal');
model.result('pg1').feature('surf1').feature('def1').label('Deformation');
model.result('pg1').feature('surf1').feature('def1').set('scaleactive', true);
model.result('pg2').label('Velocity (mbd)');
model.result('pg2').set('frametype', 'spatial');
model.result('pg2').feature('surf1').label('Surface');
model.result('pg2').feature('surf1').set('const', {'mbd.refpntx' '0' 'Reference point for moment computation, x component'; 'mbd.refpnty' '0' 'Reference point for moment computation, y component'; 'mbd.refpntz' '0' 'Reference point for moment computation, z component'});
model.result('pg2').feature('surf1').set('colortable', 'Cyclic');
model.result('pg2').feature('surf1').set('colorlegend', false);
model.result('pg2').feature('surf1').set('resolution', 'normal');
model.result('pg2').feature('surf1').feature('def1').label('Deformation');
model.result('pg2').feature('arwl1').label('Arrow Line');
model.result('pg2').feature('arwl1').set('expr', {'mbd.u_tX' 'mbd.u_tY'});
model.result('pg2').feature('arwl1').set('descr', 'Velocity');
model.result('pg2').feature('arwl1').set('const', {'mbd.refpntx' '0' 'Reference point for moment computation, x component'; 'mbd.refpnty' '0' 'Reference point for moment computation, y component'; 'mbd.refpntz' '0' 'Reference point for moment computation, z component'});
model.result('pg2').feature('arwl1').set('placement', 'elements');
model.result('pg2').feature('arwl1').feature('def1').label('Deformation');
model.result('pg3').label('Volume Loads (mbd)');
model.result('pg3').set('titletype', 'label');
model.result('pg3').set('frametype', 'spatial');
model.result('pg3').set('showlegendsunit', true);
model.result('pg3').feature('arws1').label('Gravity 1');
model.result('pg3').feature('arws1').set('expr', {'mbd.gr1.F_Vx' 'mbd.gr1.F_Vy'});
model.result('pg3').feature('arws1').set('descr', 'Load (spatial frame)');
model.result('pg3').feature('arws1').set('const', {'mbd.refpntx' '0' 'Reference point for moment computation, x component'; 'mbd.refpnty' '0' 'Reference point for moment computation, y component'; 'mbd.refpntz' '0' 'Reference point for moment computation, z component'});
model.result('pg3').feature('arws1').set('placement', 'gausspoints');
model.result('pg3').feature('arws1').feature('col').set('coloring', 'gradient');
model.result('pg3').feature('arws1').feature('col').set('topcolor', 'red');
model.result('pg3').feature('arws1').feature('col').set('bottomcolor', 'custom');
model.result('pg3').feature('arws1').feature('col').set('custombottomcolor', [0.5882353186607361 0.5137255191802979 0.5176470875740051]);
model.result('pg3').feature('arws1').feature('def').set('scale', 0);
model.result('pg3').feature('arws1').feature('def').set('scaleactive', true);
model.result('pg3').feature('surf1').label('Gray Surfaces');
model.result('pg3').feature('surf1').set('const', {'mbd.refpntx' '0' 'Reference point for moment computation, x component'; 'mbd.refpnty' '0' 'Reference point for moment computation, y component'; 'mbd.refpntz' '0' 'Reference point for moment computation, z component'});
model.result('pg3').feature('surf1').set('coloring', 'uniform');
model.result('pg3').feature('surf1').set('color', 'gray');
model.result('pg3').feature('surf1').set('resolution', 'normal');
model.result('pg3').feature('surf1').feature('def').set('scale', 0);
model.result('pg3').feature('surf1').feature('def').set('scaleactive', true);
model.result('pg4').label('Probe Plot Group 4');
model.result('pg4').set('xlabel', 'Time (s)');
model.result('pg4').set('windowtitle', 'Probe Plot 1');
model.result('pg4').set('xlabelactive', false);
model.result('pg4').feature('tblp1').label('Probe Table Graph 1');
model.result('pg4').feature('tblp1').set('legend', true);
model.result.export('anim1').set('target', 'player');
model.result.export('anim1').set('maxframes', 11);
model.result.export('anim1').set('showframe', 11);
model.result.export('anim1').set('shownparameter', '1');
model.result.export('anim1').set('fontsize', '9');
model.result.export('anim1').set('colortheme', 'globaltheme');
model.result.export('anim1').set('customcolor', [1 1 1]);
model.result.export('anim1').set('background', 'color');
model.result.export('anim1').set('gltfincludelines', 'on');
model.result.export('anim1').set('title1d', 'on');
model.result.export('anim1').set('legend1d', 'on');
model.result.export('anim1').set('logo1d', 'on');
model.result.export('anim1').set('options1d', 'on');
model.result.export('anim1').set('title2d', 'on');
model.result.export('anim1').set('legend2d', 'on');
model.result.export('anim1').set('logo2d', 'on');
model.result.export('anim1').set('options2d', 'off');
model.result.export('anim1').set('title3d', 'on');
model.result.export('anim1').set('legend3d', 'on');
model.result.export('anim1').set('logo3d', 'on');
model.result.export('anim1').set('options3d', 'off');
model.result.export('anim1').set('axisorientation', 'on');
model.result.export('anim1').set('grid', 'on');
model.result.export('anim1').set('axes1d', 'on');
model.result.export('anim1').set('axes2d', 'on');
model.result.export('anim1').set('showgrid', 'on');

model.nodeGroup('dset1mbdlgrp').label('Applied Loads (mbd)');
model.nodeGroup('dset1mbdlgrp').add('plotgroup', 'pg3');

out = model;
