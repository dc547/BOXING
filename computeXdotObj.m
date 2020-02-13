function [x_dot,SSA] = computeXdotObj(states,controls,osimModel,osimState)

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

x_dot = x_dot';

% Get Simbody Engine
simbodyEngine = osimModel.getSimbodyEngine();

hand = osimModel.getBodySet().get('hand_r');

% Transform
transform = simbodyEngine.getTransform(osimState,hand); % Returns transform wrt Joint
rotHand = transform.R();

rotMatHand = zeros(3,3);
for i = 0:2
    for j = 0:2
        rotMatHand(i+1,j+1) = rotHand.get(i,j); % Local wrt Global
    end
end

hand_mc = Vec3(0);
hand.getMassCenter(hand_mc);

origin = Vec3(0);

% Get acceleration of Hand body origin in Global frame
hand_acc = Vec3(0);
simbodyEngine.getAcceleration(osimState,hand,origin,hand_acc);

hand_w_velL = Vec3(0);
simbodyEngine.getAngularVelocity(osimState,hand,hand_w_velL);

hand_a_accL = Vec3(0);
simbodyEngine.getAngularAcceleration(osimState,hand,hand_a_accL);

pos_hand_mc = eval(hand_mc.toString().substring(1))'; % MC local
vel_hand_w  = eval(hand_w_velL.toString().substring(1))'; 
acc_hand_a  = eval(hand_a_accL.toString().substring(1))';

acc_hand    = eval(hand_acc.toString().substring(1))';

gravity = [0;-9.81;0];

% SSA global
ssa_G = acc_hand-gravity+cross(acc_hand_a,rotMatHand*pos_hand_mc)+cross(vel_hand_w,cross(vel_hand_w,rotMatHand*pos_hand_mc));

% SSA expressed in body's local frame
SSA = rotMatHand'*ssa_G;

end