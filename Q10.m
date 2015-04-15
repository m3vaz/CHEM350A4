% energy profile!

close all;
clear all;
R = 8.314;
T = 250;

profiles = [ -10 5 -5 0 -50 ;
             -6  0 -5 15 -50;
             -6 5 -5 5 -50]*10^3;

s = 1:5; % spline range
queryint = 0.1;
sq = 1:queryint:5; % spline query points

identifiers = {'a' 'i' 'p'};
arrow = '->';
rxns = {strsplit('a -> i');
        strsplit('i -> p')};


for index = 1:size(profiles, 1)
    profile = profiles(index, :);
    
    forward = [ exp(-(profile(2)-profile(1))/(R*T)) ...
                exp(-(profile(4)-profile(3))/(R*T)) ];
    backward = [ exp(-(profile(2)-profile(3))/(R*T)) ...
                exp(-(profile(4)-profile(5))/(R*T)) ];
    
    model = modelbuilder(identifiers, arrow, rxns, forward, backward);

    initial_y = [ 10 0 0 ];

    tspan = [ 0 1000 ];

    [t, Y] = ode45(model, tspan, initial_y);
    
    %plot sim results
    figure
    subplot(1, 2, 1)
    plot(t, Y)
    legend(identifiers)
    title(strcat('Case  ', num2str(index)));
    
    % plot energy diagram
    subplot(1, 2, 2)
    sq1 = interp1(s,profile,sq,'pchip'); % the smooth surface
    sq2 = interp1(s,profile,sq,'nearest'); % the levels
    
    plot(sq, sq1);
    title(strcat('Case ',num2str(index)));
    
    
    hold on
    plot(sq(1:5), sq2(1:5));
    plot(sq(6:15), sq2(6:15));
    plot(sq(16:25), sq2(16:25));
    plot(sq(26:35), sq2(26:35));
    plot(sq(36:41), sq2(36:41));
    hold off
    
    figure
    semilogy(t, Y(:, 1));
    title(strcat('Case ',num2str(index)));

end