function [ P, exitflag, output ] = LM_toy_nest( P0, data, options )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


[P, ~, ~, exitflag, output] = lsqnonlin(@nestedfun, P0,[],[],options);
foo = 1;

% Nested function
    function F = nestedfun(P0)
        
        [n_unique, n_S, ~] = size(data);
        
        % Find registration and translation matrices 
        S_R = cell(1,n_S);
        S_t = cell(1,n_S);
        ix = 1;
        for s = 2:n_S;
            rx = P0(ix);
            ry = P0(ix+1);
            rz = P0(ix+2);
            S_R{s} = compose_rotation(rx, ry, rz);
            S_t{s} = P0(ix+3:ix+5)';
            ix = ix + 6;
        end
        
        %data(:,1,:) = reshape(P0(ix:end),[(numel(P0)-ix+1)/3,3]);
        unique_1 = reshape(P0(ix:end),[(numel(P0)-ix+1)/3,3]);

        F = zeros(n_unique,n_S-1,3);
        for s = 2:n_S;
            unique_hat = (S_R{s}*squeeze(data(:,s,:))'+repmat(S_t{s},[1,n_unique]))';
           % error = unique - unique_hat;
            F(:,s-1,:) = unique_1 - unique_hat;
        end
        F = reshape(F,[],1);
        
    end
end
