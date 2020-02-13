function [f] = boxingObjFun(x,auxdata)

% Import the OpenSim modeling classes
import org.opensim.modeling.*

% Extract the necessary auxiliary data
N           = auxdata.N;
h           = auxdata.h;
nStates     = auxdata.nStates;
nActuators  = auxdata.nActuators;
nCoords     = auxdata.nCoords;
nControls   = auxdata.nControls;
osimModel   = auxdata.model_obj;
time        = auxdata.time;

osimState   = osimModel.updWorkingState();

states = zeros(N,nStates);
for i = 1:nStates
    states(:,i) = x(N*(i-1)+1:N*i,1); % column: state; row: nodes
end

controls = zeros(N,nControls);
for i = 1:nActuators
    controls(:,i) = x(nStates*N+N*(i-1)+1:nStates*N+N*i,1); % column: control
end

% Get the state derivatives and simulated sensor acceleration (SSA) in local
% frame 
x_dot = zeros(N,nStates);
SSA   = zeros(3,N);
for i = 1:N
    [x_dot(i,:),SSA(:,i)] = computeXdotObj(states(i,:)',controls(i,:)',osimModel,osimState);
end

% Calculate the integral of squared anterior-posterior hand accelerations
integral = 0;

integral = integral + (h*trapz(SSA(1,:)).^2);

f = -integral./(3/(time(end)^2));

end