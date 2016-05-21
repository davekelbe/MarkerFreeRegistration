function [ P,resnorm, residual, exitflag, output ] = LM_toy_nest4( P0, n_unique, n_S, dataix, options )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


[P, resnorm, residual, exitflag, output] = lsqnonlin(@nestedfun, P0,[],[],options);
foo3 = 1;

% Nested function
    function F = nestedfun(Pi)
        %disp(P0');
        %fprintf('\n iteration\n');
        %{
        fprintf('\n***************\n');
        fprintf('P(1) = %6.6f\t\n', Pi(1));
        fprintf('P(10) = %6.6f\t\n', Pi(10));
        fprintf('P(20) = %6.6f\t\n', Pi(20));
        fprintf('P(30) = %6.6f\t\n', Pi(30));
        fprintf('P(40) = %6.6f\t\n', Pi(40));
        fprintf('P(50) = %6.6f\t\n', Pi(50));
%}
        % Find registration and translation matrices
        S_R = cell(n_S,n_S);
        S_t = cell(n_S,n_S);
        ix = 1;
        for i = 1:n_S-1;
            for j = i+1:n_S;
                rx = Pi(ix);
                ry = Pi(ix+1);
                rz = Pi(ix+2);
                S_R{i,j} = compose_rotation(rx, ry, rz);
                S_t{i,j} = Pi(ix+3:ix+5)';
                ix = ix + 6;
            end
        end
        
        data = reshape(Pi(ix:end),[n_unique, n_S,3]);
        unique_1 = zeros(n_unique,3);
        for u = 1:n_unique;
            unique_1(u,:) = data(u,dataix(u),:);
        end
        %data(:,1,:) = reshape(P0(ix:end),[(numel(P0)-ix+1)/3,3]);
        
        F = nan(n_unique,n_S,n_S,3);
        for i = 1:n_S-1; % Reference
            for j = i+1:n_S;
                is_ucurrent = (dataix==i); % Points whose references are from camera i
                is_notnan = ~isnan(data(:,j,1)); % Points which camera j matches to i
                is_valid = is_ucurrent & is_notnan;
                n_valid = sum(is_valid);
                if n_valid== 0;
                    continue
                end
                foo = squeeze(data(is_valid,j,:))';
                if size(foo,1) ~= 3;
                    foo = reshape(foo,[3,1]);
                end
                unique_hat = (S_R{i,j}*foo+repmat(S_t{i,j},[1,n_valid]))';
                % Error between j and i
                F(is_valid,i,j,:) = unique_1(is_valid,:) - unique_hat;
                %{
                if i==1 && j==2;
                    cmax = .5;
                    cmin = 0;
                    figure;
                    legend_str{1} = sprintf('Reference points from %g', i);
                    legend_str{2} = sprintf('Matches of %g to %g',j,i);
                    plot3(unique_1(is_valid,1),unique_1(is_valid,2),unique_1(is_valid,3),'+k','markersize',10,...
                        'markerfacecolor',[0 0 0], 'linewidth',1.2);xlabel('x');ylabel('y');zlabel('z');
                    set(gca, 'Clim',[cmin cmax]);
                    view(0,90);
                    hold on
                    scatter3(unique_hat(:,1),unique_hat(:,2),unique_hat(:,3),40,...
                        sqrt(F(is_valid,i,j).^2),'filled');
                    legend(legend_str,'location','best');
                    axis equal
                    axis([-10 10 -10 10 -5 5]);
                    title('Points in WCS');
                    grid on
                    foo2 = 1;
                end
                %}
            end
        end
        %{
        for s = 2:n_S;
            unique_hat = (S_R{s}*squeeze(data(:,s,:))'+repmat(S_t{s},[1,n_unique]))';
           % error = unique - unique_hat;
            F(:,s-1,:) = unique_1 - unique_hat;
        end
        %}
        F = reshape(F,[],1);
        F = F(~isnan(F));
       % fprintf('\nRMSE = %6.6f\n', sqrt(mean(F.^2)));
        foo2 = 1; 
        % fprintf('\nF has numel %g\n',numel(F));
    end
end
