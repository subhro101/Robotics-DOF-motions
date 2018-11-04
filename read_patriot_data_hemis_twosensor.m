function [A0_all_trackable, pos_all_trackable, nframe, ntrackable, trackableID] ...
    = read_patriot_data_hemis_twosensor(filename)
% derived from read_patriot_data_hemis()
% read patriot data file with up to two sensors.
% the file may have one OR two sensor data. This code should be able to 
% process both.
%
% 07/11/2016

M = csvread(filename); 
%%!!!!!!!!!!!!!!!!!!get the first column for the times M(:,1)
times = M(:,1);
% determine number of trackables (actually, number of sensors)
ncol = size(M, 2);
if ncol == 15
    ntrackable = 2;
elseif ncol == 8
    ntrackable = 1;
else
    error('number of columns of patriot data file is neither 8 nor 15.');
end

% get trackable ID
trackableID = zeros(1, ntrackable);
for i = 0:ntrackable-1
    trackableID(i+1) = M(2, 2 + i*7);
end

nframe = size(M, 1) - 1;

A0_all_trackable  = cell(1, ntrackable); % includes all 6dof data
pos_all_trackable = cell(1, ntrackable); % all position of trackable

angular_velocity = [];

for i = 0:ntrackable-1

    %% TA's notes: get position and orientation from the patriot.csv file. No need to modify.
    pos = M(2:end, (3 + i*7):(5 + i*7))' / 100; % convert from centimeter to meter
    ori = M(2:end, (6 + i*7):(8 + i*7))' * pi / 180; % convert from degree to radians
    
    A0_all_trackable{i+1} = cell(1, nframe);
    pos_all_trackable{i+1} = zeros(3, nframe);
  n=0;  
    for f = 1:nframe
        %% TA's Notes: the function 'euler2rotm' takes the orientation data from the patriot.csv file
        %%      and converts it into the rotation matrix based on Eq. 2.65 in the textbook. 
        
        n=n+1;    
        [R,D] = euler2rotm(ori(:, f)); % rotation matrix in patriot frame
       if n==1
            Rn = R;
        else
            R = Rn * R;
        end 
        %% TA's Notes: the position data (x,y,z) is taken from the file and put into a vector form here:
        d = pos(:, f);

        %% TA's Notes: the rotation R and the position d are taken together to form a homogeneous transformation matrix 
        %%  (based on Eq. 2.91 in the textbook)
        A2 = [R            d;
              zeros(1, 3)   1];

        %% TA's Notes: you'll have to do some coordinate system transformation here to get it with respect of the tool.

        %% Question 6 a) done here.6c You should also calculate the linear and angular velocities here since you
        %%  will have access to the (x,y,z) and orientation.
       
        skew = D * transpose(R);
        ang =[  skew(1, end) skew(2) skew(end, 2)       ];
        angular_velocity=[angular_velocity; ang];
        %% -- perform your functions here before inputting it into patriot_frame2_to_0_transfm(...) function.

        %% TA's Notes: the Patriot's frame is different to that of the viewer. This part 
        %%      changes its coordinate system to that of the interface.
        %% -- remember that A0 and A2 are homogeneous transformation matrices. You will need to compute inverse of it 
        %%      in a different way than usual.
        A0 = patriot_frame2_to_0_transfm(A2);

        %% TA's notes: use the angular velocity equation to obtain the skew matrix and then take the values for x, y, z:
        %%      You can probably save everything into a matrix and use it to generate plots.
        %% ang = ......; 
        %% angular_velocity = [angular_velocity; ang];
        
        A0_all_trackable{i+1}{f} = A0;
        pos_all_trackable{i+1}(:, f) = A0(1:3, 4);

    end
    
end


%% TA's Notes: the position and orientation stays within the scope of this file. You should perform your operations here
%%      while you have easier access to it; otherwise, you will need to access the data stored within cell structures 'A0_all_trackable'.
linear_velocity = [];
linear_velocity = [];
%% velocity = pos_x - pos_x-1 / tendlinear_velocity = [linear_velocity; temp;];end

for i = 2 : nframe
    
    %% calculate the velocity from one point to another using the positions (they are stored in pos)
    %%  -- linear_velocity = [linear_velocity; temp;];
end
linear_velocity = [];
pos1 = transpose(pos_all_trackable);
% for i = 2 : nframe
%     
%     linear_velocity_x= [linear_velocity,pos_all_trackable{1,i}-pos_all_trackable{1,i-1}/times(i)];
%     linear_velocity_y= [linear_velocity,pos_all_trackable{2,i}-pos_all_trackable{2,i-1}/times(i)];
%     linear_velocity_z= [linear_velocity,pos_all_trackable{3,i}-pos_all_trackable{3,i-1}/times(i)];
%     %for j = 1 : 3
%         % velocity = pos_x - pos_x-1 / t
%        linear_velocity = [linear_velocity_x;linear_velocity_y;linear_velocity_z];
% end
%  
for i = 2 : nframe
linear_velocity = [linear_velocity; ((pos(:,i) - pos(:,i-1))/(M(i,1) -M(i-1,1)))'];
end




%% 2. creating linear velocity plots: %%
% % figure;
% % title('Angle Representation: Theta');
% % plot(1:size(thetas,1),thetas(:,1));

figure;
subplot(3,1,1)
plot(1:size(angular_velocity,1),angular_velocity(:,1));
title('Angular Velocity for X axis');
subplot(3,1,2)
plot(1:size(angular_velocity,1),angular_velocity(:,2));
title('Angular Velocity for Y axis');
subplot(3,1,3)
plot(1:size(angular_velocity,1),angular_velocity(:,3));
title('Angular Velocity for Z axis');

figure;
subplot(3,1,1)
plot(1:size(linear_velocity,1),linear_velocity(:,1));
title('Linear Velocity for X axis');
subplot(3,1,2)
plot(1:size(linear_velocity,1),linear_velocity(:,2));
title('Linear Velocity for Y axis');
subplot(3,1,3)
plot(1:size(linear_velocity,1),linear_velocity(:,3));
title('Linear Velocity for Z axis');

figure;
scatter3(linear_velocity(:,1),linear_velocity(:,2),linear_velocity(:,3),'.');
xlabel('x')
ylabel('y')
zlabel('z')

plot(x,y);
%% TA's Notes: Try using 3D plots such as scatter3, comet3, etc. for those which require 3D visualization.
%%      Otherwise, use regular subplot/plot for anything that can be visualized in 1D/2D.
