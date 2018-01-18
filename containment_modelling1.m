%% Spatial and time parameters - 1,1-dichloroethylene
x = linspace(1,500,1000); % [ft]
 y = linspace(-50,50,1000); % [ft]
 z = (0,20,100); % [ft]
z= (0,20,100); % [ft]
 z = linspace(0,20,100); % [ft]
 T = [30,365,1825] ; % [days]

 [X,Y,Z] = meshgrid(x,y,z);

%% 1,1-dichloroethylene
 for t=1:length(T)
Rf = 1.96; %[ ]
M = 12500/Rf; % [g]
ne = 0.5; % []
up = 103.68/Rf; % [ft/day]
%1.054*(10^-8); % [ft2/s]
Dx = up*(0.83*(log10(X).^2.414)); % [ft2/day]
Dy = Dx/10; % [ft2/day]
Dz = Dx/10; % [ft2/day]
k = 0.004/Rf; % [1/day]
Term1 = M/((8*ne)*((pi*t)^(3/2))*((Dx.*Dy.*Dz).^(1/2)));
Term3 = ((X-(up*t)).^2)./(Dx*t);
Term4 = (Y.^2)./(Dy*t);
Term5 = (Z.^2)./(Dz*t);
Term2 = exp((-k*t) - 0.25*(Term3 + Term4 + Term5));
Conc(:,:,:,t) = abs(Term1.*Term2);

 Conc_t30 = Conc(:,:,:,1);
Conc_tyear = Conc(:,:,:,2);
 Conc_t5year = Conc(:,:,:,3);

 %% plot surf and contour plots
 [~,idx_z] = min(abs(z-0)); % finds z-values closest to 2 ft
 figure;
 subplot(311); hold on;

surf(X(:,:,idx_z),Y(:,:,idx_z),Conc_t30(:,:,idx_z),'EdgeColor','none');
 xlim([min(x) max(x)]);
 ylin([min(y) max(y)]);

 ylim([min(y) max(y)]);
 view([0 90]);

 title([num2str(30) ' days after spill']);
 subplot(312); hold on;
surf(X(:,:,idx_z),Y(:,:,idx_z),Conc_tyear(:,:,idx_z),'EdgeColor','none');
 ylabel('y [ft]');
 xlim([min(x) max(x)]);
 ylim([min(y) max(y)]);
 view([0 90]);
 title([num2str(365) ' days after spill']);
 subplot(313); hold on;
surf(X(:,:,idx_z),Y(:,:,idx_z),Conc_t5year(:,:,idx_z),'EdgeColor','none');
 xlabel( 'x [ft]');
 xlim([min(x) max(x)]);
 ylin([min(y) max(y)]);
 ylim([min(y) max(y)]);
 view([0 90]);
 title([num2str(1825) ' days after spill']);
 figure;
 subplot(311);
surf(X(:,:,idx_z),Y(:,:,idx_z),Conc_t30(:,:,idx_z),'EdgeColor','none')'
 figure;
 subplot(311); hold on;
surf(X(:,:,idx_z),Y(:,:,idx_z),Conc_t30(:,:,idx_z),'EdgeColor','none');
 xlim([min(x) max(x)]);
 ylim([min(y) max(y)]);
 view([0 90]);
 title([num2str(30) ' days after spill']);
 title([num2str(30) ' days after spill']);
 subplot(312); hold on;
surf(X(:,:,idx_z),Y(:,;,idx_z),Conc_tyear(:,:,idx_z),'EdgeColor','none');
surf(X(:,:,idx_z),Y(:,:,idx_z),Conc_tyear(:,:,idx_z),'EdgeColor','none');
 ylabel('y [ft]');
 xlim([min(x) max(x)]);
 ylim([min(y) max(y)]);
 view([0 90]);
 title([num2str(365) ' days after spill']);

 %% Spatial and time parameters - TCE
 x = linspace(1,200,600); % [ft]
 y = linspace(-25,25,500); % [ft]
 z = linspace(0,20,100); % [ft]
 [X,Y,Z] = meshgrid(x,y,z);

 %% TCE
 for t=1:length(T)
Rf = 5.48; % []
ne = 0.5; % []
M = 17500/Rf; % [g]
up = 103.68/Rf; % [ft/day]
%1.054*(10^-8); %[ ft2/s]
Dx = up*(0.83*(log10(X).^2.414)); % [ft2/day]
Dy = Dx/10; % [ft2/day]
Dz = Dx/10; % [ft2/day]
k = 0.173/Rf; % [1/day]
Term1 = M/((8*ne)*((pi*t)^(3/2))*((Dx.*Dy.*Dz).^(1/2)));
Term3 = ((X-(up*t)).^2)./(Dx*t);
Term4 = (Y.^2)./(Dy*t);
Term5 = (Z.^2)./(Dz*t);
Term2 = exp((-k*t) - 0.25*(Term3 + Term4 + Term5));
Conc(:,:,:,t) = abs(Term1.*Term2);

 T = [30,365,1825]; % [days]
 for t=1:length(T)

 for t=1:length(T)
Rf = 5.48; % []
ne = 0.5; % []

M = 17500/Rf; % [g]
up = 103.68/Rf; % [ft/day]
%1.054*(10^-8); %[ ft2/s]
Dx = up*(0.83*(log10(X).^2.414)); % [ft2/day]
Dy = Dx/10; % [ft2/day]
Dz = Dx/10; % [ft2/day]
k = 0.173/Rf; % [1/day]
Term1 = M/((8*ne)*((pi*t)^(3/2))*((Dx.*Dy.*Dz).^(1/2)));
Term3 = ((X-(up*t)).^2)./(Dx*t);
Term4 = (Y.^2)./(Dy*t);
Term5 = (Z.^2)./(Dz*t);
Term2 = exp((-k*t) - 0.25*(Term3 + Term4 + Term5));
Conc(:,:,:,t) = abs(Term1.*Term2);

 Conc_t30 = Conc(:,:,:,1);
 Conc_tyear = Conc(:,:,:,2);
 Conc_t5year = Conc(:,:,:,3);
 figure;
 subplot(311); hold on;

surf(X(:,:,idx_z),Y(:,:,idx_z),Conc_t30(:,:,idx_z),'EdgeColor','none')
 %plot surf and contour plots
[~,idx_z] = min(abs(z-0)); % finds z-value closest to 2 ft
 figure;

 subplot(311); hold on;

surf(X(:,:,idx_z),Y(:,:,idx_z),Conc_t30(:,:,idx_z),'EdgeColor','none');
 xlim([min(x) max(x)]);
 ylim([min(y) max(y)]);
 view([0 90]);
 title([num2str(30) ' days after spill']);
 subplot(312); hold on;

surf(X(:,:,idx_z),Y(:,:,idx_z),Conc_tyear(:,:,idx_z),'EdgeColor','none');
 ylabel('y [ft]');
 xlim([min(x) max(x)]);
 ylim([min(y) max(y)]);
 view([0 90]);
 title([num2str(365) 'days after spill']);
 subplot(313); hold on;

surf(X(:,:,idx_z),Y(:,:,idx_z),Conc_t5year(:,:,idx_z),'EdgeColor','none');
 xlabel('x [ft]');
 xlim([min(x) max(x)]);
 ylim([min(y) max(y)]);

 view([0 90]);
 title([num2str(1825) ' days after spill']);
 clear all;close all;
 %% Spatial and time parameters - PCE
 x = linspace(1,300,500); % [ft]
 y = linspace(-50,50,1000); % [ft]
 z = linspace(0,20,100); % [ft]
 T = [30,365,1825]; % [days]
 [X,Y,Z] = meshgrid(x,y,z);
 %% PCE
 for t=1:length(T)
Rf = 5.18; % []
M = 12500/Rf; % [g]
ne = 0.5; % []
up = 103.68/Rf; %[ft/day]
%1.054*(10^-8); % [ft2/s]
Dx = up*(0.83*(log10(X).^2.414)); % [ft2/day]
Dy = Dx/10; % [ft2/day]
Dz = Dx/10; % [ft2/day]
k = 0.009/Rf; % [1/day]
Term1 = M/((8*ne)*((pi*t)^(3/2))*((Dx.*Dy.*Dz).^(1/2)));
Term3 = ((X-(up*t)).^2)./(Dx*t);
Term4 = (Y.^2)./(Dy*t);

Term5 = (Z.^2)./(Dz*t);
Term2 = exp((-k*t) - 0.25*(Term3 + Term4 + Term5));
Conc(:,:,:,t) = abs(Term1.*Term2);

 Conc_t30 = Conc(:,:,:,1);
 Conc_tyear = Conc(:,:,:,2);
 Conc_t5year = Conc(:,:,:,3);
 % plot surf and contour plots
 [~,idx_z] = min(abs(z-0)); % finds z-value closest to 2 ft
 figure;
 subplot(311); hold on;

surf(X(:,:,idx_z),Y(:,:,idx_z),Conc_t30(:,:,idx_z),'EdgeColor','none');
 xlim([min(x) max(x)]);
 ylim([min(y) max(y)]);
 view([0 90]);
 title([num2str(30) ' days after spill']);
 subplot(312); hold on;

surf(x(:,:,idx_z),Y(:,:,idx_z),Conc_tyear(:,:,idx_z),'EdgeColor','none');
 ylabel('y [ft]');
 xlim([min(x) max(x)]);

 ylim([min(x) max(x)]);
 view([0 90]);
 title([num2str(365) ' days after spill']);
 subplot(313); hold on;
surf(X(:,:,idx_z),Y(:,:,idx_z),Conc_t5year(:,:,idx_z),'EdgeColor','none');
 xlabel('x [ft]');
 xlim([min(x) max(x)]);
 ylim([min(y) max(y)]);
 view([0 90]);
 title([num2str(1825) ' days after spill']);
 subplot(312); hold on;
surf(x(:,:,idx_z),Y(:,:,idx_z),Conc_tyear(:,:,idx_z),'EdgeColor','none');
 ylabel('y [ft]');
 xlim([min(x) max(x)]);
 ylim([min(y) max(y)]);

 %% Spatial and time parameters - Vinyl Chloride
 x = linspace(1,300,500); % [ft]
 y = linspace(-50,50,500); % [ft]
 z = linspace(0,20,100); % [ft]
 T = [30,365,1825]; % [days]
 [X,Y,Z] = meshgrid(x,y,z);
 %% Vinyl Chloride
 for t=1:length(T)
Rf = 1.50; % []
M = 7500/Rf; % [g]
ne = 0.5; % []
up = 103.68/Rf; % [ft/day]
%1.054*(10^-8); %[ft2/s]
Dx = up*(0.83*(log10(X).^2.414)); % [ft2/day]
Dy = Dx/10; % [ft2/day]
Dz = Dx/10; % [ft2/day]
k = 0.229/Rf; % [1/day]
Term1 = M/((8*ne)*((pi*t)^(3/2))*((Dx.*Dy.Dz).^(1/2)));
Term3 = ((X-(up*t)).^2)/(Dx*t);
Term4 = (Y.^2)./(Dy*t);
Term5 = (Z.^2)./(Dz*t);
Term2 = exp((-k*t) - 0.25*(Term3 + Term4 + Term5));
Conc(:,:,:,t) = abs(Term1.*Term2);
 for t=1:length(T)
Term1 = M/((8*ne)*((pi*t)^(3/2))*((Dx.*Dy.*Dz).^(1/2)));
Term3 = ((X-(up*t)).^2)/(Dx*t);
Term4 = (Y.^2)./(Dy*t);
Term5 = (Z.^2)./(Dz*t);
Term2 = exp((-k*t) - 0.25*(Term3 + Term4 + Term5));
Conc(:,:,:,t) = abs(Term1.*Term2);

 for t=1:length(T)
Term1 = M/((8*ne)*((pi*t)^(3/2))*((Dx.*Dy.*Dz).^(1/2)));
Term3 = ((X-(up*t)).^2)./(Dx*t);
Term4 = (Y.^2)./(Dy*t);
Term5 = (Z.^2)./(Dz*t);
Term2 = exp((-k*t) - 0.25*(Term3 + Term4 + Term5));
Conc(:,:,:,t) = abs(Term1.*Term2);
end
 Conc_t30 = Conc(:,:,:,1);
 Conc_tyear = Conc(:,:,:,2);
 Conc_t5year = Conc(:,:,:,3);
 % plot surf plots
 [~,idx_z] = min(abs(z-0)); % finds z-value closest to 2 ft
 figure;
 subplot(311); hold on;


surf(X(:,:,idx_z),Y(:,:,idx_z),Conc_t30(:,:,idx_z),'EdgeColor','none');
 xlim([min(x) max(x)]);
 ylim([min(y) max(y)]);
 view([0 90])l
 view([0 90]);
 title([num2str(30) ' days after spill']);
 subplot(312); hold on;
surf(X(:,:,idx_z),Y(:,:,idx_z),Conc_tyear(:,:,idx_z),'EdgeColor','none');
 ylabel('y [ft]');
 xlim([min(x) max(x)]);
 ylim([min(y) max(y)]);
 view([0 90]);
 title([num2str(365) ' days after spill']);
 subplot(313); hold on;
surf(X(:,:,idx_z),Y(:,:,idx_z),Conc_t5year(:,:,idx_z),'EdgeColor','none');
 xlabel('x [ft]');
 xlim([min(x) max(x)]);
 ylim([min(y) max(y)]);
 view([0 90]);
 title([num2str(1825) ' days after spill']);