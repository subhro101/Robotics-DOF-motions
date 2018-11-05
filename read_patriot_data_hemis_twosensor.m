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

   
    pos = M(2:end, (3 + i*7):(5 + i*7))' / 100; % convert from centimeter to meter
    ori = M(2:end, (6 + i*7):(8 + i*7))' * pi / 180; % convert from degree to radians
    
    A0_all_trackable{i+1} = cell(1, nframe);
    pos_all_trackable{i+1} = zeros(3, nframe);
  n=0;  
    for f = 1:nframe
        
        n=n+1;    
        [R,D] = euler2rotm(ori(:, f)); % rotation matrix in patriot frame
       if n==1
            Rn = R;
        else
            R = Rn * R;
        end 
        
        d = pos(:, f);

       
        A2 = [R            d;
              zeros(1, 3)   1];

       
        skew = D * transpose(R);
        ang =[  skew(1, end) skew(2) skew(end, 2)       ];
        angular_velocity=[angular_velocity; ang];
        %% -- perform your functions here before inputting it into patriot_frame2_to_0_transfm(...) function.

       
        A0 = patriot_frame2_to_0_transfm(A2);

        
        
        A0_all_trackable{i+1}{f} = A0;
        pos_all_trackable{i+1}(:, f) = A0(1:3, 4);

    end
    
end



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
