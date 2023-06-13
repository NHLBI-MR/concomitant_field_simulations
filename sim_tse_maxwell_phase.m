function maxwell_results = sim_tse_maxwell_phase(tse_gradient_s)
% function maxwell_results = sim_tse_maxwell_phase(tse_gradient_s)
%
% Simulations for https://pubmed.ncbi.nlm.nih.gov/37306784/ : figure 3
%
% This function assumes:
%   Spiral Gradients and Concomitant Correction Gradients have been
%   calculated and imported in the following format:
%
%   tse_gradient_s.gradients_tse_train [ samples  x  tse train  x  axes ] 
%   - 2D/3D imaging gradients between each refocussing pulse
%
%   tse_gradient_s.gradients_90_180 [ samples  x  1  x  axes ] 
%   - 2D/3D imaging gradients between excitation and refocussing pulse
%
%   These gradients have been sampled at 10 us, the gradient raster time of
%   Siemens systems
%
% multi-echo TSE simulations - June 2023
% R Ramasawmy, NHLBI

% return input to output
maxwell_results = tse_gradient_s;

%% Quick-peek: TSE Gradient assessment

% Extract number of echoes in train
tse_etl = size(tse_gradient_s.gradients_tse_train,2);

% Visualise x-gradients over TSE TRAIN
figure('Name','Gradient Assessment'), 
color_ax = [.3 .3 .3];

subplot(2,1,1); 
colorful_plots( tse_gradient_s.gradients_tse_train(:,:,1));
title('Gradients over echo train (mT/m)')

% self-squared concomitant field assessment 

% HARD CODED ECHO TIMES
echo_center = [ 3796 7647 11498 15349 19200 23051];


self_squared_grads = cumsum((tse_gradient_s.gradients_tse_train).^2, 1);
ss_running_phase_t = cumsum(squeeze(tse_gradient_s.gradients_90_180).^2);

for i = 1:tse_etl
    ss_running_phase_t = cat(1, ss_running_phase_t, -1*ss_running_phase_t(end,:) +  squeeze(self_squared_grads(:,i,:))   );
end

subplot(2,1,2), hold on,
plot(ss_running_phase_t); legend({'x','y'}); 
plot(echo_center, ss_running_phase_t(echo_center), 'ko')
yline(0,'--', 'Color',color_ax)
title('Gradient self-squared integral')

maxwell_results.self_squared_grads = ss_running_phase_t;

%% PANEL 1 - trajectory visualisation


% Integrate gradients to trajectory (no scaling)
trajectory = cumsum(tse_gradient_s.gradients_tse_train,1);
% Normalize for not properly scaled trajectories..
trajectory = trajectory./max(trajectory(:));

view_1 = .2;            % Squeeze to 20% to view center of k-space

figure, 
subplot(1,8,1); hold on; xline(0,'--', 'Color',color_ax);yline(0,'--', 'Color',color_ax)

% Choose center of trajectories to avoid pre/rewinders (HARD CODED for example)
colorful_plots(trajectory(1000:2844,:,1),trajectory(1000:2844,:,2));
axis image; xlim(view_1*[-1 1]);ylim(view_1*[-1 1]);xticklabels(''); yticklabels('');
title('Trajectory');

%% PANEL 2 - TSE simulation

    %% Physical Sim // Input parameters

    B_0 = 0.55;                             % T
    % R = [1 0 0; 0 1 0; 0 0 1];            % axial
    R = [0 0 1; 0 1 0; 1 0 0];              % sagittal
    % R = [1 0 0; 0 0 1; 0 1 0];            % coronal

    fov   = [24 24 0];                      % cm
    fov   = fov/100;                        % convert to m

    steps = [160 160 1];                    % matrix
    res   = fov./steps;                     %

    offset = [0 0 -0.5];                    % cm [R P S]
    offset = offset/100;                    % convert to m

    %% Map to physical coordinates

    phys_coords     = R*(fov)';
    phys_steps      = R*steps';
    phys_res        = R*res';
    phys_offset     = R*offset';

    phys_coords_x    = linspace(-1*(phys_coords(1)*0.5-phys_res(1)*0.5)+ phys_offset(1),(phys_coords(1)*0.5-phys_res(1)*0.5)+ phys_offset(1), phys_steps(1) );
    phys_coords_y    = linspace(-1*(phys_coords(2)*0.5-phys_res(2)*0.5)+ phys_offset(2),(phys_coords(2)*0.5-phys_res(2)*0.5)+ phys_offset(2), phys_steps(2) );
    phys_coords_z    = linspace(-1*(phys_coords(3)*0.5-phys_res(3)*0.5)+ phys_offset(3),(phys_coords(3)*0.5-phys_res(3)*0.5)+ phys_offset(3), phys_steps(3) );

    % whos phys_coords_x phys_coords_y phys_coords_z
    %     [min(phys_coords_x) max(phys_coords_x);
    %      min(phys_coords_y) max(phys_coords_y);
    %      min(phys_coords_z) max(phys_coords_z)]

    [Y, X, Z] = meshgrid(phys_coords_y, phys_coords_x, phys_coords_z);
    % whos X Y Z
    %   [min(X(:)) max(X(:));
    %    min(Y(:)) max(Y(:));
    %    min(Z(:)) max(Z(:))]

    xyz = zeros(length(phys_coords_x), length(phys_coords_y), length(phys_coords_z), 3, 'single');
    xyz(:,:,:,1) = X; xyz(:,:,:,2) = Y; xyz(:,:,:,3) = Z;

    %% Concomitant Phase - 90-180

    pre_comp.gradients_nominal           = tse_gradient_s.gradients_90_180;
    pre_comp.gradients_nominal(:,:,3)    = 0;                               % PAD 2D gradients to 3D
    pre_comp.ints                        = 1; 

    disp('Simulating 90-180 gradient concomitant phase')
    [pre_comp.phase_plot_ints, pre_comp.phase_plot] = calc_Bc_ppm_phase(pre_comp, xyz, B_0, R);


    %% Concomitant Phase - Echo Train
    tse_train.gradients_nominal           = tse_gradient_s.gradients_tse_train; 
    tse_train.gradients_nominal(:,:,3)    = 0;                                  % PAD 2D gradients to 3D    
    tse_train.ints                        = tse_etl;

    disp('Simulating Echo Train gradient concomitant phase')
    [tse_train.phase_plot_ints, tse_train.phase_plot] = calc_Bc_ppm_phase(tse_train, xyz, B_0, R);

%% echo spacing with history

running_phase_t = pre_comp.phase_plot;

for i = 1:tse_etl
    running_phase_t = cat(3, running_phase_t, -1*running_phase_t(:,:,end) +  squeeze(tse_train.phase_plot(:,:,:,i))   );
end

running_phase_t = single(running_phase_t);

phase_echo_time = running_phase_t(:,:,echo_center);

for i = 1:tse_etl
    subplot(1,8,i+1); 
    imshow(phase_echo_time(:,:,i),[-360 360]); colormap('parula')
    title(['Echo ' num2str(i)])
end

%% PANEL 3 - TSE RSS

phase_rss = sqrt(sum((phase_echo_time - phase_echo_time(:,:,1)).^2,3));

sample_roi = circ_roi(40, [80 80],[160 160]);

subplot(1,8,8); hold on;
imshow(phase_rss,[-360 360]); colormap('parula')
contour(sample_roi, 'r-')
title('RSS phase')

mean_phase_rss = mean(phase_rss(sample_roi));
disp(['Mean phase RSS: ' num2str(mean_phase_rss) ' std phase RSS: ' num2str(std(phase_rss(sample_roi))) ' deg' ]);

%% EXPORT

maxwell_results.phase_echo_time     = phase_echo_time;
maxwell_results.phase_rss           = phase_rss;
maxwell_results.mean_phase_rss      = mean_phase_rss;


end


function [phase_plot_ints, phase_plot_store] = calc_Bc_ppm_phase(s_s, xyz, B_0, R)
%%
gamma_c = 42.57e6;          % Hz/T

time_step = 1e-5;           % s

Bc_pc = @(xyz, gradients_physical, B_0, t) (1/(2*B_0))*...
    (gradients_physical(t,1).^2.*xyz(:,:,:,3).^2 + gradients_physical(t,2).^2.*xyz(:,:,:,3).^2 + ...
    gradients_physical(t,3).^2.*(xyz(:,:,:,1).^2 + xyz(:,:,:,2).^2).*0.25 ...
    - gradients_physical(t,1).*gradients_physical(t,3).*xyz(:,:,:,1).*xyz(:,:,:,3) ...
    - gradients_physical(t,2).*gradients_physical(t,3).*xyz(:,:,:,2).*xyz(:,:,:,3)     );

%% Rotating spirals

Bc_plot     = zeros([size(xyz,[1:3]),length(s_s.gradients_nominal), s_s.ints], 'single');
phase_plot  = zeros([size(xyz,[1:3]),length(s_s.gradients_nominal), s_s.ints], 'single');
%         size(phase_plot_s_ints)
    
for solid_int= 1:s_s.ints
    
    disp(['echo ' num2str(solid_int) ' / ' num2str(s_s.ints)])
    
    % Extract grads
    gradients_nominal       = [s_s.gradients_nominal(:,solid_int,1) s_s.gradients_nominal(:,solid_int,2) s_s.gradients_nominal(:,solid_int,3)];
    gradients_physical      = (R*gradients_nominal')';
    
    % Run through gradient   
    
    %% Calc & Sum
    
    for i = 1:length(gradients_physical)
        if i == 1
            Bc_plot(:,:,:,i, solid_int)      = (Bc_pc(xyz, gradients_physical, B_0, i));
            phase_plot(:,:,:,i, solid_int)   = (Bc_pc(xyz, gradients_physical, B_0, i))*time_step*2*pi*gamma_c;
        else
            Bc_plot(:,:,:,i, solid_int)      = (Bc_pc(xyz, gradients_physical, B_0, i));
            phase_plot(:,:,:,i, solid_int)   = (Bc_pc(xyz, gradients_physical, B_0, i))*time_step*2*pi*gamma_c + phase_plot(:,:,:, i-1, solid_int);
        end
    end
    
end

% rotate for scanner orientation
phase_plot_store= rot90(squeeze(phase_plot)*(180/pi)); % convert to degrees
phase_plot_ints = rot90(squeeze(phase_plot(:,:,:,end,:))*(180/pi)); % convert to degrees


end

function [c_mask] = circ_roi(r, centers, msize)
if nargin == 1
    cx=64;  cy=64;
    ix=128; iy=128;
else
    cx=centers(1);  cy=centers(2);
    ix=msize(1);    iy=msize(2);
end
[x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
c_mask=((x.^2+y.^2)<=r^2);

end