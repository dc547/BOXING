function [c,ceq] = boxingConFun(x,auxdata)

% Import the OpenSim modeling classes
import org.opensim.modeling.*

% Extract the necessary auxiliary data
N           = auxdata.N;
h           = auxdata.h;
nStates     = auxdata.nStates;
nActuators  = auxdata.nActuators;
nCoords     = auxdata.nCoords;
nControls   = auxdata.nControls;
osimModel   = auxdata.model_con;
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

% Compute equality constraint violation using Backward Euler method
states_dot = zeros(N-1,nStates);
for i = 1:nStates
    for j = 1:N-1
        states_dot(j,i) = (states(j+1,i)-states(j,i))/h;
    end
end

% Get the state derivatives from OpenSim
x_dot = zeros(N-1,nStates);
for i = 1:N-1
    x_dot(i,:) = computeXdotCon(states(i+1,:)',controls(i+1,:)',osimModel,osimState)';
end

% Evaluate constraint violations
ceq_temp = zeros(N-1,nStates);
for i = 1:nStates
    ceq_temp(:,i) = states_dot(1:end,i)-x_dot(1:end,i);
end

% Re-arrange the constraint violations in one long vector
for i = 1:nStates
    ceq((N-1)*(i-1)+1:(N-1)*i,:) = ceq_temp(:,i);
end

% Additional task constraints

% R Sho q_ti
ceq(end+1,1) = states(1,1) - deg2rad(20);

% R Sho q_tf
c(1,1) = deg2rad(60) - states(N,1);
c(2,1) = states(N,1) - deg2rad(110);

% R Elb q_ti
ceq(end+1,1) = states(1,2) - deg2rad(150);

% R Elb q_tf
ceq(end+1,1) = states(N,2) - 0;

end