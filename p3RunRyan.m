% Run profile
% Badgerloop Pod 3
% Justin Williams
% Started 8/28/2017
% V0.1.2
clc, clear

%% Input Parameters
% Tube details
tube_length = 1250; %m
tube_pressure = 70; %torr

% General pod details
pod_mass = 0;               %kg--defined by systems
pod_mass_structure = 30;
pod_mass_stability = 15;

% Braking details
braking_pressure_1 = 120; %psi
braking_pressure_2 = 120; %psi
braking_mass = 20;

% Cold Gas Thruster Details
propulsion_pressure_1 = 4000; %psi
propulsion_pressure_2 = 4000; %psi
propulsion_mass = 15; %10 kg estimate for the wheel
propulsion_mass_thrusters = 60; %
% Electrical Details
battery_internal_resistance = 10; %mOhm per cell
battery_cell_series = 75; 
battery_cell_parallel = 1;
battery_cell_voltage = 3.6;
battery_cell_amphours = 10;
battery_cell_Crating = 18;

% Motor Details
motor_current = 200; %Amps
motor_volt = 240; %Volts
motor_rpm_per_vdc = 25;

motor_rpm = 6000;%battery_cell_series * battery_cell_voltage * motor_rpm_per_vdc;  
motor_torque = 100; %Nm
motor_power = motor_rpm * motor_torque / 9550; %kW

motor_wheel_diameter = 6; %Inches
motor_wheel_torque_ratio = .6;
motor_mass = 7;
motor_mass_controller = 8;
motor_mass_battery = 20; 
motor_heat_capacity = 0.896; %J/g-C
motor_efficiency = 0.90; % percent of electrical power to mechanical power



%% Calculations
pod_mass = motor_mass + motor_mass_controller + motor_mass_battery + propulsion_mass + braking_mass+pod_mass_structure+pod_mass_stability + propulsion_mass_thrusters;

% Pod Kinematics
mWheelDia = motor_wheel_diameter*2.54/100;
mPower = motor_torque*6.28/60*motor_rpm;
vMax = motor_rpm / 60*3.14*mWheelDia / motor_wheel_torque_ratio;
fMax = 2*motor_torque * motor_wheel_torque_ratio / mWheelDia;
aMax = fMax / pod_mass;
xMax = vMax^2 / (2*aMax);
aBraking = vMax^2 / 2 / (tube_length - xMax);
fBraking = aBraking*pod_mass;
tMax = vMax/aMax;

% Motor Temperature 
mPowerLoss = (1-motor_efficiency)*mPower;
mHeatCapacity = motor_heat_capacity*motor_mass;
mTempIncreasePS = mPowerLoss/mHeatCapacity/1000;
mFinalT = tMax*mTempIncreasePS + 25;

% CGTs
force_thrusters = 1775; %N
thrust_time = 3; % s
aMax_thrusters = force_thrusters / pod_mass; 
vMax_thrusters = vMax + aMax_thrusters * thrust_time;

xMax_thrusters = vMax_thrusters^2 / (2*aMax_thrusters);
aBraking_thrusters = vMax_thrusters^2 / 2 / (tube_length - xMax_thrusters);
fBraking_thrusters = aBraking_thrusters*pod_mass;
tMax = vMax_thrusters/aMax_thrusters;



%% Results
fprintf('Pod mass %3.0f kg\n', pod_mass);
fprintf('Motor acceleration: %2.1f m/s^2 \t Braking acceleration: %2.1f m/s^2 \t Top speed: %3.1f m/s \n', aMax, aBraking, vMax);
fprintf('Motor Temp: %3.1fC\n', mFinalT);
fprintf('Motor acceleration with thrusters: %2.1f m/s^2 \t Braking acceleration with thrusters: %2.1f m/s^2 \t Top speed with thrusters: %3.1f m/s \n', aMax_thrusters, aBraking_thrusters, vMax_thrusters);
