clear all
clc
% Model Import
import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model')
model.modelPath(['D:\usr']);

%% FEM parameters

model.param.set('M0', '0[N*m]');
model.param.set('dtFem', '0.02[s]');
model.param.set('phi', '0[rad]');
model.param.set('phi_t', '0[rad/s]');

model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 2);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('r1').set('pos', [-0.025 0]);
model.component('comp1').geom('geom1').feature('r1').set('size', [0.05 1]);
model.component('comp1').geom('geom1').create('pt1', 'Point');
model.component('comp1').geom('geom1').run('pt1');
model.component('comp1').geom('geom1').create('rot1', 'Rotate');
model.component('comp1').geom('geom1').feature('rot1').selection('input').set({'r1'});
model.component('comp1').geom('geom1').feature('rot1').set('rot', 0);
model.component('comp1').geom('geom1').run;

model.component('comp1').physics.create('mbd', 'MultibodyDynamics', 'geom1');
model.component('comp1').physics('mbd').create('rd1', 'RigidDomain', 2);
model.component('comp1').physics('mbd').feature('rd1').selection.set([1]);
model.component('comp1').physics('mbd').feature('rd1').feature('crp1').selection.set([3]);
model.component('comp1').physics('mbd').create('hgj1', 'HingeJoint', -1);
model.component('comp1').physics('mbd').feature('hgj1').set('Source', 'fixed');
model.component('comp1').physics('mbd').feature('hgj1').set('Destination', 'rd1');
model.component('comp1').physics('mbd').feature('hgj1').create('afm1', 'AppliedForceAndMoment', -1);
model.component('comp1').physics('mbd').create('gr1', 'Gravity', 2);
model.component('comp1').physics('mbd').feature('gr1').selection.set([1]);

model.component('comp1').mesh('mesh1').create('map1', 'Map');
model.component('comp1').probe.create('dom1', 'Domain');
model.component('comp1').probe.create('dom2', 'Domain');

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material('mat1').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
model.component('comp1').material('mat1').label('Structural steel');
model.component('comp1').material('mat1').propertyGroup('def').set('density', '7850e-3[kg/m^3]');
model.component('comp1').material('mat1').propertyGroup('Enu').set('E', '200e9[Pa]');
model.component('comp1').material('mat1').propertyGroup('Enu').set('nu', '0.30');

model.component('comp1').physics('mbd').prop('PhysicsSymbols').set('PhysicsSymbols', true);
model.component('comp1').physics('mbd').prop('d').set('d', 0.05);
% model.component('comp1').physics('mbd').feature('rd1').set('CenterOfRotationType', 'CentroidOfSelectedEntities');
% model.component('comp1').physics('mbd').feature('rd1').set('EntityLevel', 'Point');
model.component('comp1').physics('mbd').feature('rd1').set('ConsistentInitialization', 'ForceInitialValues');
model.component('comp1').physics('mbd').feature('hgj1').set('CenterOfJointType', 'UserDefined');
model.component('comp1').physics('mbd').feature('hgj1').feature('afm1').set('AppliedOnSelectedAttachment', 'Joint');
model.component('comp1').physics('mbd').feature('hgj1').feature('afm1').set('PointOfApplicationType', 'CenterOfJoint');
model.component('comp1').physics('mbd').feature('hgj1').feature('afm1').set('Mz', 'input1');
model.component('comp1').physics('mbd').feature('hgj1').feature('afm1').set('Ms', 'M0');
model.component('comp1').physics('mbd').feature('rd1').set('InitialValueType', 'locallyDefined');
model.component('comp1').physics('mbd').feature('rd1').feature('init1').set('phi', 'phi');
model.component('comp1').physics('mbd').feature('rd1').feature('init1').set('phit', 'phi_t');

model.component('comp1').physics('mbd').feature('hgj1').create('fric1', 'Friction', -1);
model.component('comp1').physics('mbd').feature('hgj1').feature('fric1').set('mu', 0.2);
model.component('comp1').physics('mbd').feature('hgj1').feature('fric1').set('v0', 0);

model.component('comp1').mesh('mesh1').run;
model.component('comp1').probe('dom1').set('expr', 'mbd.rd1.phi');
model.component('comp1').probe('dom2').set('intsurface', true);
model.component('comp1').probe('dom2').set('intvolume', true);
model.component('comp1').probe('dom2').set('expr', 'mbd.rd1.th_tz');

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

model.result.dataset.create('dset2', 'Solution');
model.result.dataset.create('avh1', 'Average');
model.result.dataset('dset2').set('probetag', 'dom1');
model.result.dataset('avh1').set('probetag', 'dom1');
model.result.dataset('avh1').set('data', 'dset2');
model.result.dataset('avh1').selection.geom('geom1', 2);
model.result.dataset('avh1').selection.set([1]);
model.result.numerical.create('pev1', 'EvalPoint');
model.result.numerical('pev1').set('probetag', 'dom1');
model.result.create('pg1', 'PlotGroup1D');
model.result('pg1').create('tblp1', 'Table');
model.result('pg1').feature('tblp1').set('probetag', 'dom1');
model.component('comp1').probe('dom1').genResult([]);

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
model.sol('sol1').feature('v1').feature('comp1_mbd_rd1_u').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_mbd_rd1_phi').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_mbd_hgj1_th').set('scalemethod', 'manual');
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
kp =6e1;
kd =5e2;
clf_rate=3;
slack = 1e5;
umax=100;
umin=-umax;
syms x1 x2;
x = [x1; x2];

%% Init condition
x_init=[0;0];
u=0;
phi=0.1;
phi_t=0;

%% System Identification
t=[0];
array_t=t;
th=[phi];
array_th=th;
th_t=[0];
array_th_t=th_t;
array_t_ode23=[];
array_th_ode23=[];

for k=1:5e3
%% Comsol Paramset
model.param.set('M0', strcat(num2str(u),'[N*m]'));
model.param.set('phi', phi,'[rad]');
model.param.set('phi_t', phi_t,'[rad/s]');
model.component('comp1').physics('mbd').feature('rd1').feature('init1').set('phi', 'phi');
model.component('comp1').physics('mbd').feature('rd1').feature('init1').set('phit', 'phi_t');
model.sol('sol1').runAll;
M = mphstate(model,'sol1','out',{'Mc' 'MA' 'MB' 'A' 'B' 'C' 'D' 'x0', 'Null', 'ud'},'input', {'M0'}, 'output', {'comp1.dom1','comp1.dom2'}, 'sparse', 'off', 'initmethod','sol','solnum','first');

%% ODE23s Param set(DAE is stiff. So, manipulate with jacobian and Mc)
x_init=[phi/M.Null(1,1); phi_t/M.Null(5,2)];
x1=x_init(1)
x2=x_init(2);
func = @(tt,xx)M.MA*xx+M.MB*u;
opt=odeset('mass',M.Mc,'jacobian',M.MA);
[t_ode,x_ode]=ode23s(func,[0, 0.02],x_init,opt);
y=M.C*x_ode';
y(2,:)=[];

array_t_ode23=[array_t_ode23;0.02*(k-1)+t_ode];
array_th_ode23=[array_th_ode23;y'];

t=mphglobal(model,'root.t')+t(end);
array_t=[array_t;t];
th=mphglobal(model,'mbd.hgj1.th');
array_th=[array_th;th];
phi=array_th(end)
th_t=mphglobal(model,'mbd.hgj1.th_t');
array_th_t=[array_th_t;th_t];
phi_t=array_th_t(end);

%% Constraints : A[u; slack] <= b
% Determine control input
A=[M.A(1,1),M.A(1,2); M.A(2,1)-kp/M.Null(1,1),M.A(2,2)-kd/M.Null(5,2)];
Q=clf_rate*eye(size(A,1));
P=lyap(A',Q);
clf= x'*P*x;
dclf=simplify(jacobian(clf,x));
xdim=size(A,1);
udim=size(M.B,2);

AA=[double(subs(dclf*M.B)),-1];
b=-double(subs(dclf*A*x))-clf_rate*double(subs(clf));
%  And max input constraint
AA=[AA; eye(udim),zeros(udim,1)];
b=[b;umax];
% And min inpu contstraint
AA =[AA;-eye(udim),zeros(udim,1)];
b=[b;-umin];
%% Cost
H=[eye(udim),zeros(udim,1);
    zeros(1,udim),slack];
u_ref=zeros(1,udim);
ff=[u_ref; 0];
u=quadprog(H,ff,AA,b);
u=u(1)
end

figure(1)
plot(array_t,array_th)


hold on

plot(array_t_ode23,array_th_ode23,'--')
legend('comsol','ode23s');
xlabel('time(s)');
ylabel('theat(rad)');