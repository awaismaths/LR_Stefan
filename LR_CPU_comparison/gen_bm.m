%% Brownian motions for the finest mesh to be used for all other coarse meshes
% Reminder: we are using dt = h^2;
% Constructing Brownian motion (Financial toolbox is required)
clear all

% Parameters
T = 1; % Total time
nbm = 25; % Number of Brownian motions
dt_initial = [0.1; 0.01; 0.001; 0.0001]; % Initial time steps
meshes = {'mesh1_1.mat'; 'mesh1_2.mat'; 'mesh1_3.mat'; 'mesh1_4.mat'; 'mesh1_5.mat'; 'mesh1_6.mat'}; % Mesh files
nbmeshes = size(meshes, 1); % Number of meshes

% Loop over each initial time step
for dti = 1:size(dt_initial, 1)
    % Loop over each mesh
    for imesh = 1:nbmeshes
        % Load fine mesh data
        fine_mesh = matfile(strcat('../HHO-Lapl-OM-master/matlab_meshes/', meshes{imesh}));
        h(imesh) = max(abs(fine_mesh.diam)); % Calculate maximum diameter of the finest mesh
        Ndt(imesh)=2*round(0.5*T/h(imesh)^2);
        %Ndt(imesh) = T / dt_initial(dti); % Calculate number of time steps
        dt = T / Ndt(imesh); % Calculate time step size
        WF = zeros(Ndt(imesh), nbm); % Initialize array to store Brownian motion

        % Generate Brownian motions
        for bmm = 1:nbm
            [Path, Time, dW] = simulate(bm(0, 1), Ndt(imesh)); % Simulate Brownian motion
            WF(:, bmm) = sqrt(dt) * dW; % Scale the Brownian increments
        end
        
        % Save the Brownian motions to a file
        save(strcat('BM/BM_mesh', num2str(meshes{imesh}(1:7)), 'Ndt', num2str(Ndt(imesh))), 'WF');
        clear WF; % Clear WF to free up memory
    end
end



% %% Brownian motions for the finest mesh to be used for all other coarse meshes
% %Reminder: we are using dt=h^2;
% %Constructing Brownian motion (Finacial toolbox is  required)
% clear all
% T=1;
% nbm=100;
% dt_initial = [0.1; 0.01; 0.001; 0.0001];
% meshes={'mesh1_1.mat';'mesh1_2.mat';'mesh1_3.mat';'mesh1_4.mat';'mesh1_5.mat';'mesh1_6.mat'};
% nbmeshes=size(meshes,1);
% for dti=1:size(dt_initial,1)
%     for imesh=1:nbmeshes
%         fine_mesh=matfile(strcat('../../HHO-Lapl-OM-master/matlab_meshes/',meshes{imesh}));
%         h(imesh)=max(abs(fine_mesh.diam));
%         Ndt(imesh)=T/dt_initial(dti);%2*round(0.5*T/h(imesh)^2);
%         dt=T/Ndt(imesh);
%         WF=zeros(Ndt(imesh),nbm);
%         for bmm=1:nbm
%             [Path,Time,dW]=simulate(bm(0,1),Ndt(imesh));
%             WF(:,bmm)=sqrt(dt)*dW;
%         end
%         save(strcat('BM/BM_mesh',num2str(meshes{imesh}(1:7)),'Ndt',num2str(Ndt(imesh))),'WF');
%         clear WF;
%     end
% end