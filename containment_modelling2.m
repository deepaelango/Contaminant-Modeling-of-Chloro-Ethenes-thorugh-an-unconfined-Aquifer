%% Spatial and time parameters - 1,1-dichloroethylene
 x = linspace(1,500,1000); % [ft]
 y = linspace(-50,50,1000); % [ft]
 z = (0,20,100); % [ft]
z = (0,20,100); % [ft]
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
%% Plot depth z chart
[~,idx_max_t30] = max(Conc_t30(idx_y0,:,1));
[~,idx_max_tyear] = max(Conc_tyear(idx_y0,:,1));
[~,idx_max_t5year] = max(Conc_t5year(idx_y0,:,1));
figure;
subplot(311);
plot(squeeze(Conc_t30(idx_y0,idx_max_t30),:)),z);
set(gca, ‘Ydir’,’reverse’)
title([num2str(30) ‘ days after spill’]);
subplot(312);
plot(squeeze(Conc_tyear(idx_y0,idx_max_tyear),:)),z);
set(gca, ‘Ydir’,’reverse’)
title([num2str(365) ‘ days after spill’]);
subplot(312);
plot(squeeze(Conc_t5year(idx_y0,idx_max_t5year),:)),z);

set(gca, ‘Ydir’,’reverse’)
title([num2str(1825) ‘ days after spill’]);
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

%% Plot depth z chart
[~,idx_max_t30] = max(Conc_t30(idx_y0,:,1));
[~,idx_max_tyear] = max(Conc_tyear(idx_y0,:,1));
[~,idx_max_t5year] = max(Conc_t5year(idx_y0,:,1));
figure;
subplot(311);
plot(squeeze(Conc_t30(idx_y0,idx_max_t30),:)),z);
set(gca, ‘Ydir’,’reverse’)
title([num2str(30) ‘ days after spill’]);
subplot(312);
plot(squeeze(Conc_tyear(idx_y0,idx_max_tyear),:)),z);
set(gca, ‘Ydir’,’reverse’)
title([num2str(365) ‘ days after spill’]);
subplot(312);
plot(squeeze(Conc_t5year(idx_y0,idx_max_t5year),:)),z);
set(gca, ‘Ydir’,’reverse’)
title([num2str(1825) ‘ days after spill’]);
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

%% Plot depth z chart
[~,idx_max_t30] = max(Conc_t30(idx_y0,:,1));
[~,idx_max_tyear] = max(Conc_tyear(idx_y0,:,1));
[~,idx_max_t5year] = max(Conc_t5year(idx_y0,:,1));
figure;
subplot(311);
plot(squeeze(Conc_t30(idx_y0,idx_max_t30),:)),z);
set(gca, ‘Ydir’,’reverse’)
title([num2str(30) ‘ days after spill’]);
subplot(312);
plot(squeeze(Conc_tyear(idx_y0,idx_max_tyear),:)),z);
set(gca, ‘Ydir’,’reverse’)
title([num2str(365) ‘ days after spill’]);
subplot(312);
plot(squeeze(Conc_t5year(idx_y0,idx_max_t5year),:)),z);
set(gca, ‘Ydir’,’reverse’)
title([num2str(1825) ‘ days after spill’]);
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
 for t=1:length(T)
Term1 = M/((8*ne)*((pi*t)^(3/2))*((Dx.*Dy.*Dz).^(1/2)));
Term3 = ((X-(up*t)).^2)./(Dx*t);
Term4 = (Y.^2)./(Dy*t);
Term5 = (Z.^2)./(Dz*t);
Term2 = exp((-k*t) - 0.25*(Term3 + Term4 + Term5));
Conc(:,:,:,t) = abs(Term1.*Term2);

 Conc_t30 = Conc(:,:,:,1);
 Conc_tyear = Conc(:,:,:,2);
 Conc_t5year = Conc(:,:,:,3);
%% Plot depth z chart
[~,idx_max_t30] = max(Conc_t30(idx_y0,:,1));

[~,idx_max_tyear] = max(Conc_tyear(idx_y0,:,1));
[~,idx_max_t5year] = max(Conc_t5year(idx_y0,:,1));
figure;
subplot(311);
plot(squeeze(Conc_t30(idx_y0,idx_max_t30),:)),z);
set(gca, ‘Ydir’,’reverse’)
title([num2str(30) ‘ days after spill’]);
subplot(312);
plot(squeeze(Conc_tyear(idx_y0,idx_max_tyear),:)),z);
set(gca, ‘Ydir’,’reverse’)
title([num2str(365) ‘ days after spill’]);
subplot(312);
plot(squeeze(Conc_t5year(idx_y0,idx_max_t5year),:)),z);
set(gca, ‘Ydir’,’reverse’)
title([num2str(1825) ‘ days after spill’]);
