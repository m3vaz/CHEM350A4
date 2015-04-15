% Rate constant consistency

close all;
clear all;
%format long;

% top line is forward rate, bottom line is backward rate
% column n is rates for reaction n
forward = [1.0 1.0 1.0];
backward = [ 0.5, 0.2, 2.0];

identifiers = {'a' 'b' 'c'};
arrow = '->';
rxns = {strsplit('a -> b');
        strsplit('a -> c');
        strsplit('b -> c')};
rxnidentifiers = cellfun(@strjoin, rxns, 'UniformOutput', false);
    
model = modelbuilder(identifiers, arrow, rxns, forward, backward);
    
initial_y = [ 1 0 0 ];

tspan = [ 0 20 ];

[t, Y] = ode45(model, tspan, initial_y);

% calclulate the instantaneous rates
rates = zeros(length(t), length(rxns));

for i = 1:length(t)
    [temp, rates(i, :) ] = model(t(i), Y(i, :)); 
end

figure
plot(t, Y);
legend(identifiers);
title('Concentration With Arbitrary Rates')

figure
plot(t, rates);
legend(rxnidentifiers);
title('Elementary rxn rates (Arbitrary)')

% but the rates aint zero!
% let's fix that
%% Rates, the internally consistent way

backward(3) = forward(3)*backward(2)*forward(1)/(backward(1)*forward(2));

model = modelbuilder(identifiers, arrow, rxns, forward, backward);
    
initial_y = [ 1 0 0 ];

tspan = [ 0 20 ];

[t, Y] = ode45(model, tspan, initial_y);

% calclulate the instantaneous rates
rates = zeros(length(t), length(rxns));

for i = 1:length(t)
    [temp, rates(i, :) ] = model(t(i), Y(i, :)); 
end

figure
plot(t, Y);
legend(identifiers);
title('Consistent Rates');

figure
plot(t, rates);
legend(rxnidentifiers);
title('Elementary rxn rates (Consistent)');

