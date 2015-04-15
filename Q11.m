% explosion simulation

close all;
clear all;

profiles = [ -10 5 -5 0 -50 ;
             -6  0 -5 15 -50;
             -6 5 -5 5 -50]*10^3;

s = 1:5; % spline range
queryint = 0.1;
sq = 1:queryint:5; % spline query points

identifiers = {'H_2' 'O_2' 'O_3' 'H_2O' 'OH^.' 'HOOH' 'OOH^.' 'H^.', 'O'};
arrow = '->';
rxns = {strsplit('O_2 -> O O');% 1
        strsplit('O_2 O -> O_3');
        strsplit('O H_2 -> OH^. H^.');
        strsplit('H^. O_2 -> OH^. O');
        strsplit('OH^. H_2 -> H_2O H^.');% 5
        strsplit('OH^. O_2 -> OOH^. O');
        strsplit('H^. OH^. -> H_2O');
        strsplit('OH^. OH^. -> HOOH');
        strsplit('H^. OOH^. -> HOOH');
        strsplit('H_2 O_2 -> OH^. OH^.')};% 10



forward =[1e-4, 1.0e-3,	1.0e-3,	1.0e-2,	1.0e-2,	1.0e-2,	1.0e-1,...
            1.0e-1,	1.0e-1,	1.0e-4];
backward =[1e-2, 1.0e-5, 1.0e-2, 1.0e-4	, 1.0e-3, 1.0e-2, 1.0e-6, ...
            1.0e-6, 1.0e-6, 1.0e-1];

model = modelbuilder(identifiers, arrow, rxns, forward, backward);

initial_y = [ 20 10 0 0 0 0 0 0 0 ];

initialOx = 10;

% preequilib
mask = eye(length(identifiers));
mask = mask(:, ismember(identifiers, 'O_2'));
preequilib_y = mask*initialOx;
tspan = [ 0 50000 ];

[t, Y] = ode45(model, tspan, preequilib_y);

%plot sim results
figure
plot(t, Y);
legend(identifiers);
title('Preequilibrium');

% actual sim time
initial_y = initial_y + Y(end, :);

tspan = [ 0 1000 ];

[t, Y] = ode45(model, tspan, initial_y);

%plot sim results
figure
plot(t, Y);
legend(identifiers);
title('Concentrations');

% moar plots
% add all radicals together
radidentifiers = {'OH^.' 'OOH^.' 'O' 'H^.'};
radY = Y(:, ismember(identifiers, radidentifiers));
radTotY = sum(radY, 2);
figure
plot(t, radTotY);
hold on
plot(t, radY);
legend(['All Radicals' radidentifiers]);
title('Radical Concentrations');

% do the log
figure
radTotYScaled = (-sqrt(-log(radTotY)));
orig_state = warning('off','all'); 
plot(t, radTotYScaled);
warning(orig_state);
title('Looking for a BOOM');
