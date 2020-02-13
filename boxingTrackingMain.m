function [Xopt,fval,exitflag,output] = boxingTrackingMain(varargin)

% Main function for Boxing Project

% Import the OpenSim modeling classes
import org.opensim.modeling.*


%% Define options
Options.N = 30; % number of mesh intervals

% Path of main OCP folder 
OCPpath = pwd;


%% Check function inputs  - IMU data

% % Create paths for data
% pathData = [erase(OCPpath,'\OCP'),'\ACC_DATA\'];
% addpath(genpath(pathData));
% 
% load('Ax_raw.mat'); % punches_x2
% load('Ay_raw.mat'); % punches_y2
% load('Az_raw.mat'); % punches_z2
% 
% if length(varargin) == 0
% 
%     accX = punches_x2{1,1}{12,1};
%     accY = punches_y2{1,1}{12,1};
%     accZ = punches_z2{1,1}{12,1};
%     
% end


%% Spline IMU data 
% 
% % Accelerometer fs
% accFs = 250;
% 
% % Trial duration
% time = ([0:length(accX)-1]*(accFs^-1))';
% 
% % DC time
% timeDC = linspace(0,time(end),Options.N+1);
% 
% % Create Piecewise Polynomials (PP)
% pp_accX = spline(time,accX);
% pp_accY = spline(time,accY);
% pp_accZ = spline(time,accZ);
% 
% % Evaluate PP at collocation points
% accX_s = ppval(pp_accX,timeDC);
% accY_s = ppval(pp_accY,timeDC);
% accZ_s = ppval(pp_accZ,timeDC);
% 

%% Osim Model
pathModel = [OCPpath,'\Model\'];
addpath(genpath(pathModel));

% Read in the osim model
osimModel = Model([pathModel,'Boxing_model_3DOF.osim']);

% Initialize the model
osimState = osimModel.initSystem();

% Get the number of states, coordinates, actuators and controls of the
% model; in this case number of controls == number of actuators
nStates    = osimModel.getNumStateVariables();
stateNames = osimModel.getStateVariableNames();
nControls  = osimModel.getNumControls();
nCoords    = osimModel.getNumCoordinates();
nActuators = osimModel.getActuators().getSize(); 

% Get the names of the states from the model
states_all = cell(nStates,1);
for i = 1:nStates
    states_all(i,1) = cell(osimModel.getStateVariableNames().getitem(i-1));
end

% Get the names of controls/actuators from the model
controls_all = cell(nControls,1);
for i = 1:nControls
    controls_all(i,1) = cell(osimModel.getActuators().get(i-1).getName());
end


%% DC Time
duration = 0.2;
h        = duration/(Options.N-1);
timeDC   = h*(0:Options.N-1);


%% Auxiliary data to be passed to the optimiser
auxdata.model_obj  = osimModel;
auxdata.model_con  = osimModel;
auxdata.time       = timeDC;
auxdata.N          = Options.N;
auxdata.h          = h;
auxdata.nStates    = nStates;
auxdata.nControls  = nControls;
auxdata.nCoords    = nCoords;
auxdata.nActuators = nActuators;


%% Toy tracking data
rShoAng1 = deg2rad(25);  rShoAng2 = deg2rad(90);
rElbAng1 = deg2rad(150); rElbAng2 = deg2rad(0);
rSupAng1 = deg2rad(0);   rSupAng2 = deg2rad(130);

coeSho = polyfit([min(timeDC) max(timeDC)],[rShoAng1 rShoAng2],1);
coeElb = polyfit([min(timeDC) max(timeDC)],[rElbAng1 rElbAng2],1);
coeSup = polyfit([min(timeDC) max(timeDC)],[rSupAng1 rSupAng2],1);

rShoAng = polyval(coeSho,timeDC);
rElbAng = polyval(coeElb,timeDC);
rSupAng = polyval(coeSup,timeDC);

rShoAngV = polyval(coeSho(1),timeDC); 
rElbAngV = polyval(coeElb(1),timeDC);
rSupAngV = polyval(coeSup(1),timeDC);

x0_temp = [rShoAng' rElbAng' rSupAng' rShoAngV' rElbAngV' rSupAngV'];

% Begin populating vector of design variables (DVs)
for i = 1:nStates
    X0(Options.N*(i-1)+1:Options.N*i,1) = x0_temp(:,i);
end

% Initial guess for torque actuator controls
u0_temp = zeros(Options.N,nActuators);
u0_temp(:,[1,3]) = 0.05;
u0_temp(:,2) = -0.05;

% Populate vector of DVs with control variables
for i = 1:nControls
    X0(nStates*Options.N+Options.N*(i-1)+1:nStates*Options.N+Options.N*i,1) = u0_temp(:,i);
end

%% Lower & Upper Bounds

Pos_LB(1:nCoords*Options.N)       = -0.5; 
Vel_LB(1:Options.N)               = -0.5;
Vel_LB(Options.N+1:Options.N*2)   = -25;
Vel_LB(Options.N*2+1:Options.N*3) = -2;
Con_LB(1:nActuators*Options.N)    = -1.0;

Pos_UB(1:nCoords*Options.N)       = 3.0; 
Vel_UB(1:Options.N)               = 20;
Vel_UB(Options.N+1:Options.N*2)   = 1;
Vel_UB(Options.N*2+1:Options.N*3) = 25;
Con_UB(1:nActuators*Options.N)    = 1.0;

lb = [Pos_LB Vel_LB Con_LB]';
ub = [Pos_UB Vel_UB Con_UB]';

%% Define anonymous functions for the objective and constraint functions

objfun = @(x)boxingObjFun(x,auxdata);
confun = @(x)boxingConFun(x,auxdata);


%% Define fmincon settings and solve problem

options = optimset('algorithm','interior-point','TolFun',1e-5,'TolX',1e-6,...
    'TolCon',1e-3,'FinDiffType','forward','MaxFunEvals',1e5,'Hessian','bfgs',...
    'display','iter');

[Xopt,fval,exitflag,output] = fmincon(objfun,X0,[],[],[],[],lb,ub,confun,options);

%% Print the optimal states and controls to STO files

X_state_opt = zeros(Options.N,nStates);
for i = 1:nStates
    X_state_opt(:,i) = Xopt(Options.N*(i-1)+1:Options.N*i,1);
end

X_controls_opt = zeros(Options.N,nControls);
for i = 1:nControls
    X_controls_opt(:,i) = Xopt(nStates*Options.N+Options.N*(i-1)+1:nStates*Options.N+Options.N*i,1);
end

% Create data structure for the states file
StatesData = struct();
StatesData.name = [char(osimModel.getName()), '_Optimal_States_', char(date)];
StatesData.nRows = size(timeDC, 2);
StatesData.nColumns = nStates+1; %All the states + time
StatesData.inDegrees = false;
StatesData.labels = cell(1,StatesData.nColumns); 
StatesData.labels{1}= 'time';
for j = 2:1:StatesData.nColumns
   StatesData.labels{j} = char(states_all(j-1));
end
StatesData.data = [timeDC', X_state_opt];
writeOpenSimStatesFile(StatesData);

% Create data structure for the controls file
ControlData = struct();
ControlData.name = [char(osimModel.getName()), '_Optimal_Controls_', char(date)];
ControlData.nRows = size(dc_time, 2);
ControlData.nColumns = Ncontrols+1; %All the controls + time
ControlData.inDegrees = false;
ControlData.labels = cell(1,ControlData.nColumns); 
ControlData.labels{1}= 'time';
for j = 2:1:ControlData.nColumns
   ControlData.labels{j} = char(controls_all(j-1));
end
ControlData.data = [timeDC', X_controls_opt];
writeOpenSimControlFile(ControlData)


numVar = osimState.getNY();
% for i = 0:numVar-1
%     osimState.updY().set(i,0.5);
% end

states = [0.2 0 -0.2 5 4 2];

for i = 0:nStates-1
    osimModel.setStateVariable(osimState,stateNames.get(i),states(i+1));
end

osimModel.computeStateVariableDerivatives(osimState);

modelControls = osimModel.updControls(osimState);
actuatorControls = Vector(1,0.0);

cont = [1; 0.0; 0];

for i = 0:nActuators-1
    actuatorControls.set(0,cont(i+1));
    osimModel.updActuators().get(i).addInControls(actuatorControls,modelControls);
end

% Set the control (excitation) values
osimModel.setControls(osimState,modelControls);

osimModel.computeStateVariableDerivatives(osimState);

% Get Simbody Engine
simbodyEngine = osimModel.getSimbodyEngine();
% Hand body
hand = osimModel.getBodySet().get('hand_r');
upperarm = osimModel.getBodySet().get('r_humerus');
ulna   = osimModel.getBodySet().get('ulna_r');

% Transform
transform = simbodyEngine.getTransform(osimState,hand); % Returns transform wrt Joint
tformR    = simbodyEngine.getTransform(osimState,upperarm);

% Body wrt Ground
rot = tformR.R();
trans = tformR.T();

rotR = transform.R();

rotMat = zeros(3,3);
rotMatHand = zeros(3,3);
for i = 0:2
    for j = 0:2
        rotMat(i+1,j+1) = rot.get(i,j); % Local wrt Global
        rotMatHand(i+1,j+1) = rotR.get(i,j);
    end
end

UA_mc = Vec3(0);
upperarm.getMassCenter(UA_mc);

R_mc = Vec3(0);
ulna.getMassCenter(R_mc);

H_mc = Vec3(0);
hand.getMassCenter(H_mc);

origin = Vec3(0);

UA_pos = Vec3(0);
simbodyEngine.getPosition(osimState,upperarm,UA_mc,UA_pos);

UA_vel = Vec3(0);
simbodyEngine.getVelocity(osimState,upperarm,UA_mc,UA_vel);

UA_o_vel = Vec3(0);
simbodyEngine.getVelocity(osimState,upperarm,origin,UA_o_vel);

UA_o_acc = Vec3(0);
simbodyEngine.getAcceleration(osimState,upperarm,origin,UA_o_acc);

UA_acc = Vec3(0);
simbodyEngine.getAcceleration(osimState,upperarm,UA_mc,UA_acc);

hand_acc = Vec3(0);
simbodyEngine.getAcceleration(osimState,hand,origin,hand_acc);

handCOM_acc = Vec3(0);
simbodyEngine.getAcceleration(osimState,hand,H_mc,handCOM_acc);

% upperarm_acc = Vec3(0);
% simbodyEngine.getAcceleration(osimState,upperarm,mc,upperarm_acc);

% UA_w_accG = Vec3(0);
% simbodyEngine.getAngularAcceleration(osimState,upperarm,UA_w_accG);

UA_w_accL = Vec3(0);
simbodyEngine.getAngularAccelerationBodyLocal(osimState,upperarm,UA_w_accL);

UA_w_accG = Vec3(0);
simbodyEngine.getAngularAcceleration(osimState,upperarm,UA_w_accG);

UA_w_velL = Vec3(0);
simbodyEngine.getAngularVelocityBodyLocal(osimState,upperarm,UA_w_velL);

UA_w_velG = Vec3(0);
simbodyEngine.getAngularVelocity(osimState,upperarm,UA_w_velG);


R_w_velL = Vec3(0);
simbodyEngine.getAngularVelocity(osimState,ulna,R_w_velL);

H_w_velG = Vec3(0);
simbodyEngine.getAngularVelocity(osimState,hand,H_w_velG);

H_w_velL = Vec3(0);
simbodyEngine.getAngularVelocityBodyLocal(osimState,hand,H_w_velL);


R_w_accL = Vec3(0);
simbodyEngine.getAngularAcceleration(osimState,ulna,R_w_accL);

H_w_accL = Vec3(0);
simbodyEngine.getAngularAcceleration(osimState,hand,H_w_accL);

derivs = osimState.getYDot();


acc = eval(UA_w_accL.toString().substring(1))';
pos = eval(UA_mc.toString().substring(1))';
vel = eval(UA_w_velL.toString().substring(1))';
acc_hand = eval(hand_acc.toString().substring(1))';
pos_UA = eval(UA_pos.toString().substring(1))';
pos_mc = eval(UA_mc.toString().substring(1))';

COM_accHand = eval(handCOM_acc.toString().substring(1))';

pos_H_mc = eval(H_mc.toString().substring(1))';
vel_H_w = eval(H_w_velL.toString().substring(1))';
acc_H_w = eval(H_w_accL.toString().substring(1))';

s = acc_hand+cross(acc_H_w,rotMatHand*pos_H_mc)+cross(vel_H_w,cross(vel_H_w,rotMatHand*pos_H_mc));
ss = rotMatHand'*COM_accHand;
sss = rotMatHand'*s;
ssss = rotMatHand'*acc_hand+cross(acc_H_w,rotMatHand*pos_H_mc)+cross(vel_H_w,cross(vel_H_w,rotMatHand*pos_H_mc));
sssss = rotMatHand'*acc_hand+cross(rotMatHand'*acc_H_w,pos_H_mc)+cross(rotMatHand'*vel_H_w,cross(rotMatHand'*vel_H_w,pos_H_mc)); % Linear Acceleration in local
% s = rotMat'*acc_hand+cross(acc,pos)+cross(vel,cross(vel,pos));

t=1;
