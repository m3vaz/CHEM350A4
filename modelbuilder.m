function [ modelhandle ]  ...
        = modelbuilder(identifiers, arrow, rxns, forward, backward)
%MODELBUILDER Build models to be used by ode solvers
%   Requires identifer list, arrow identifier, and a list of reactions
%   Returns a model for use
    
    if (iscellstr(rxns))
        rxns = {rxns};
    end
    
    % basic test of validity
    if length(rxns) ~= length(forward) || length(rxns) ~= length(backward)
        error('Input length mismatch');
    end
    
    % Assume that all identifiers are valid
    % if not, then errors going to be thrown
    % I'm not going to fix your data
    
    % this definition is done at runtime, so the values are bound
    % dynamically
    % this is the function that will be being passed back and is the
    % function to be passed to the ode solver
    function [ dydt,  rates ] = model(~, y)
        % equations are time invariant so t value is never used
        dydt = zeros(length(identifiers), 1);
        rates = zeros(1, length(rxns));
        
        % determine rxn rates
        for index = 1:length(rxns)
            % grab the rxn
            rxn = rxns{index};
            % grab the rates
            f = forward(index);
            b = backward(index);
            % grab the representation in the rxn
            forwardrep = zeros(length(identifiers), 1);
            backwardrep = zeros(length(identifiers), 1);
            
            % process the rxn
            pos = 1;
            forwarddiff = 1;
            badkwarddiff = 1;
            
            while true
                % run off the end of the array
                if (pos > length(rxn))
                    error(strcat('Arrow not found in rxn ', num2str(index)))
                end
                if (strcmp(rxn{pos}, arrow))
                    break;
                end
                
                identifier = rxn{pos};
                identpos = find(ismember(identifiers, identifier));
                if isempty(identpos)
                    error(strcat('Unexpected identifier ',identifier));
                end
                
                forwarddiff = forwarddiff*y(identpos);
                
                forwardrep(identpos) = forwardrep(identpos) + 1;
                
                pos = pos + 1;
            end
            
            % move past the arrow
            pos = pos+1;
            
            while pos <= length(rxn)
                
                identifier = rxn{pos};
                identpos = find(ismember(identifiers, identifier));
                if isempty(identpos)
                    error('Unexpected identifier');
                end
                
                badkwarddiff = badkwarddiff*y(identpos);
                
                backwardrep(identpos) = backwardrep(identpos) + 1;
                
                pos = pos + 1;
            end
            
            rate = f*forwarddiff - b*badkwarddiff;
            rates(index) = rate;
            dydt = dydt + rate*(backwardrep-forwardrep);
        end
        
        % dydt determined
        
    end
    
    % return them up
    modelhandle = @model;

end