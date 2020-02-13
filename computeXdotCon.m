function x_dot = computeXdotCon(states,controls,osimModel,osimState)

% Import the OpenSim modeling classes
import org.opensim.modeling.*

% Update model state with current values
numVar = osimState.getNY();
stateNames = osimModel.getStateVariableNames();
for i = 0:numVar-1
    osimModel.setStateVariable(osimState,stateNames.get(i),states(i+1));
end

osimModel.computeStateVariableDerivatives(osimState);

modelControls = osimModel.updControls(osimState);
actuatorControls = Vector(1,0.0);

nActuators = osimModel.getActuators().getSize();

for i = 0:nActuators-1
    actuatorControls.set(0,controls(i+1,1));
    osimModel.updActuators().get(i).addInControls(actuatorControls,modelControls);
end

osimModel.setControls(osimState,modelControls);

osimModel.computeStateVariableDerivatives(osimState);

x_dot = zeros(numVar,1);

for i = 0:numVar-1
    x_dot(i+1,1) = osimState.getYDot().get(i);
end


end