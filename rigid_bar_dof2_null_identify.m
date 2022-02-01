clear all
close all
clc
% Model Import
import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model')
model.modelPath(['C:\usr']);

%% FEM parameters

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
model.component('comp1').physics('mbd').feature('gr1').selection.set([1,2]);
model.component('comp1').physics('mbd').create('rd1', 'RigidDomain', 2);
model.component('comp1').physics('mbd').feature('rd1').selection.set([1]);
model.component('comp1').physics('mbd').feature('rd1').set('CenterOfRotationType', 'CentroidOfSelectedEntities');
model.component('comp1').physics('mbd').feature('rd1').set('EntityLevel', 'Point');
model.component('comp1').physics('mbd').feature('rd1').feature('crp1').selection.set([4]);
model.component('comp1').physics('mbd').create('rd2', 'RigidDomain', 2);
model.component('comp1').physics('mbd').feature('rd2').selection.set([2]);
model.component('comp1').physics('mbd').feature('rd2').set('CenterOfRotationType', 'CentroidOfSelectedEntities');
model.component('comp1').physics('mbd').feature('rd2').set('EntityLevel', 'Point');
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
model.component('comp1').probe.create('var1', 'GlobalVariable');
model.component('comp1').probe.create('var2', 'GlobalVariable');
model.component('comp1').probe.create('var3', 'GlobalVariable');
model.component('comp1').probe.create('var4', 'GlobalVariable');

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material('mat1').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
model.component('comp1').material('mat1').label('Structural steel');
model.component('comp1').material('mat1').propertyGroup('def').set('density', '7850e-3[kg/m^3]');
model.component('comp1').material('mat1').propertyGroup('Enu').set('E', '200e9[Pa]');
model.component('comp1').material('mat1').propertyGroup('Enu').set('nu', '0.30');

model.component('comp1').physics('mbd').prop('PhysicsSymbols').set('PhysicsSymbols', true);
model.component('comp1').physics('mbd').prop('d').set('d', 0.05);
model.component('comp1').physics('mbd').feature('rd1').set('ConsistentInitialization', 'ForceInitialValues');
model.component('comp1').physics('mbd').feature('hgj1').set('CenterOfJointType', 'UserDefined');
model.component('comp1').physics('mbd').feature('hgj1').feature('afm1').set('AppliedOnSelectedAttachment', 'Joint');
model.component('comp1').physics('mbd').feature('hgj1').feature('afm1').set('Ms', 'M0');
model.component('comp1').physics('mbd').feature('hgj2').set('CenterOfJointType', 'UserDefined');
model.component('comp1').physics('mbd').feature('hgj2').set('xc', [0; 1; 0]);
model.component('comp1').physics('mbd').feature('hgj2').feature('afm1').set('Mz', 'M1');

% model.component('comp1').physics('mbd').feature('rd1').set('InitialValueType', 'locallyDefined');
% model.component('comp1').physics('mbd').feature('rd1').feature('init1').set('phi', 'phi');
% model.component('comp1').physics('mbd').feature('rd1').feature('init1').set('phit', 'phi_t');

model.component('comp1').physics('mbd').feature('hgj1').create('fric1', 'Friction', -1);
model.component('comp1').physics('mbd').feature('hgj1').feature('fric1').set('mu', 0.2);
model.component('comp1').physics('mbd').feature('hgj1').feature('fric1').set('v0', 0);
model.component('comp1').physics('mbd').feature('hgj2').create('fric1', 'Friction', -1);
model.component('comp1').physics('mbd').feature('hgj2').feature('fric1').set('mu', 0.2);
model.component('comp1').physics('mbd').feature('hgj2').feature('fric1').set('v0', 0);

model.component('comp1').mesh('mesh1').run;
model.component('comp1').probe('var1').set('expr', 'mbd.rd1.thz');
model.component('comp1').probe('var2').set('expr', 'mbd.rd1.th_tz');
model.component('comp1').probe('var2').set('unit', 'rad/s');
model.component('comp1').probe('var3').set('expr', 'mbd.rd2.thz');
model.component('comp1').probe('var4').set('expr', 'mbd.rd2.th_tz');
model.component('comp1').probe('var4').set('unit', 'rad/s');

model.study.create('std1');
model.study('std1').create('time', 'Transient');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('t1', 'Time');
model.sol('sol1').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').feature.remove('fcDef');


model.study('std1').feature('time').set('tlist', 'range(0,dtFem/3,dtFem)');


model.sol('sol1').attach('std1');

%% run
model.sol('sol1').feature('st1').label('Compile Equations: Time Dependent');
model.sol('sol1').feature('v1').label('Dependent Variables 1.1');
model.sol('sol1').feature('v1').set('control', 'user');
model.sol('sol1').feature('v1').set('initsol', 'zero');
model.sol('sol1').feature('v1').set('initmethod', 'init');
model.sol('sol1').feature('t1').set('control', 'time');
model.sol('sol1').feature('v1').set('scalemethod', 'manual');
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
model.sol('sol1').feature('t1').set('control', 'user');
model.sol('sol1').feature('t1').set('rtol', 0.001);
model.sol('sol1').feature('t1').set('tstepsbdf', 'strict');
model.sol('sol1').feature('t1').set('maxorder', 2);
model.sol('sol1').feature('t1').feature('dDef').label('Direct 1');
model.sol('sol1').feature('t1').feature('aDef').label('Advanced 1');
model.sol('sol1').feature('t1').feature('aDef').set('cachepattern', true);
model.sol('sol1').feature('t1').feature('fc1').label('Fully Coupled 1.1');
model.sol('sol1').feature('t1').feature('fc1').set('maxiter', 50);
model.sol('sol1').feature('t1').feature('fc1').set('jtech', 'onevery');


%% Assumed feedback gian, CLF
% kp =6e1;
% kd =5e2;
% clf_rate=3;
% slack = 1e5;
% umax=100;
% umin=-umax;
% syms x1 x2;
% x = [x1; x2];
% 
%% Init condition
u=[1;1];
phi=0;
phi_t=0;

%% System Identification
t=[0];
array_t=[];
th=[phi];
array_th=th;
th_t=[0];
array_th_t=th_t;
array_t_ode23=[];
array_rd1_th_ode23=[];
array_rd1_thz_ode23=[];
array_rd2_th_ode23=[];
array_rd2_thz_ode23=[];
array_rd1_th=[];
array_rd1_th_t=[];
array_rd2_th=[];
array_rd2_th_t=[];

x_ode_init=[0,0,0,0];
model.sol('sol1').runAll;
model.sol('sol1').feature('v1').set('initmethod', 'sol');
model.sol('sol1').feature('v1').set('initsol', 'sol1');
model.sol('sol1').feature('v1').set('solnum', 'last');

for k=1:1
%% Comsol Paramset
model.param.set('M0', strcat(num2str(u(1)),'[N*m]'));
model.param.set('M1', strcat(num2str(u(2)),'[N*m]'));
model.sol('sol1').runAll;
M = mphstate(model,'sol1','out',{'Mc' 'MA' 'MB' 'A' 'B' 'C' 'D' 'x0', 'Null', 'ud'},'input', {'M0','M1'}, 'output', {'comp1.var1','comp1.var2','comp1.var3','comp1.var4'}, 'sparse', 'off', 'initmethod','sol','solnum','first');

%% ODE23s Param set(DAE is stiff. So, manipulate with jacobian and Mc)
func = @(tt,xx)M.MA*xx+M.MB*u;
opt=odeset('mass',M.Mc,'jacobian',M.MA);
[t_ode,x_ode]=ode23s(func,[0, 0.02],x_ode_init,opt);
x_ode_init= x_ode(end,:);
M.C(2,:)=M.Null(9,:);
M.C(4,:)=M.Null(10,:);
y=M.C*x_ode';

array_t_ode23=[array_t_ode23;0.02*(k-1)+t_ode];
array_rd1_th_ode23=[array_rd1_th_ode23;y(1,:)'];
array_rd1_thz_ode23=[array_rd1_thz_ode23;y(2,:)'];
array_rd2_th_ode23=[array_rd2_th_ode23;y(3,:)'];
array_rd2_thz_ode23=[array_rd2_thz_ode23;y(4,:)'];

t=mphglobal(model,'root.t')+t(end);
array_t=[array_t;t];
rd1_th=mphglobal(model,'mbd.hgj1.th');
array_rd1_th=[array_rd1_th;rd1_th];
rd1_th_t=mphglobal(model,'mbd.hgj1.th_t');
array_rd1_th_t=[array_rd1_th_t;rd1_th_t];

rd2_th=mphglobal(model,'mbd.hgj2.th');
array_rd2_th=[array_rd2_th;rd2_th];
rd2_th_t=mphglobal(model,'mbd.hgj2.th_t');
array_rd2_th_t=[array_rd2_th_t;rd2_th_t];

end

figure(1)
plot(array_t,array_rd1_th)
hold on
plot(array_t_ode23,array_rd1_th_ode23,'--')
legend('comsol','ode23s');
xlabel('time(s)');
ylabel('theta(rad)');

figure(2)
plot(array_t,array_rd1_th_t)
hold on
plot(array_t_ode23,array_rd1_thz_ode23,'--')
legend('comsol','ode23s');
xlabel('time(s)');
ylabel('dtheta(rad/s)');

figure(3)
plot(array_t,array_rd2_th)
hold on
plot(array_t_ode23,array_rd2_th_ode23,'--')
legend('comsol','ode23s');
xlabel('time(s)');
ylabel('theta(rad)');

figure(4)
plot(array_t,array_rd2_th_t)
hold on
plot(array_t_ode23,array_rd2_thz_ode23,'--')
legend('comsol','ode23s');
xlabel('time(s)');
ylabel('theta(rad/s)');