function [ P,exitflag, output ] = LM_toy_nest6( P0, data, dataix,R_empty,G_path, options )
%LM_TOY_NEST5 graphical adjustment
%   Detailed explanation goes here

G_path = G_path(~cellfun(@isempty,G_path));
n_g = numel(G_path);

% Check results before
%{
[~,data_WCSt, data_WCSj] =nestedfun(P0);
for g = 1:n_g;
    for i = 1:n_S;
        if  sum(sum(~isnan(squeeze(data_WCSt{1}(:,i,:,1)))))== 0 &&...
                sum(sum(~isnan(squeeze(data_WCSt{2}(:,i,:,1)))))  == 0;
            continue
        end
        figure;
        hold on
        plot3(reshape(data_WCSj{g}(:,i,:,1),[],1),...
            reshape(data_WCSj{g}(:,i,:,2),[],1),...
            reshape(data_WCSj{g}(:,i,:,3),[],1),'ok','markersize',5,...
            'markerfacecolor','r')
        plot3(reshape(data_WCSt{g}(:,i,:,1),[],1),...
            reshape(data_WCSt{g}(:,i,:,2),[],1),...
            reshape(data_WCSt{g}(:,i,:,3),[],1),'+','markersize',20,...
            'markerfacecolor',[0 0 0], 'linewidth',2);
        title(sprintf('Initial %g paths referenced to i=%g',g+1,i));
        grid on
        axis([-15 15 -15 15 -5 5]);
    end
end
%}

[P, resnorm, residual, exitflag, output] = lsqnonlin(@nestedfun, P0,[],[],options);
foo3 = 1;

% Check results after
%{
[~,data_WCSt0, data_WCSj0] =nestedfun(P0);
[~,data_WCSt, data_WCSj] =nestedfun(P);
for g = 1:n_g;
    for i = 1:n_S;
        if  sum(sum(~isnan(squeeze(data_WCSt{1}(:,i,:,1)))))== 0 &&...
                sum(sum(~isnan(squeeze(data_WCSt{2}(:,i,:,1)))))  == 0;
            continue
        end
        figure;
        hold on
        plot3(reshape(data_WCSj{g}(:,i,:,1),[],1),...
            reshape(data_WCSj{g}(:,i,:,2),[],1),...
            reshape(data_WCSj{g}(:,i,:,3),[],1),'ok','markersize',5,...
            'markerfacecolor','r')
        plot3(reshape(data_WCSj0{g}(:,i,:,1),[],1),...
            reshape(data_WCSj0{g}(:,i,:,2),[],1),...
            reshape(data_WCSj0{g}(:,i,:,3),[],1),'ok','markersize',5,...
            'markerfacecolor','b')
        plot3(reshape(data_WCSt{g}(:,i,:,1),[],1),...
            reshape(data_WCSt{g}(:,i,:,2),[],1),...
            reshape(data_WCSt{g}(:,i,:,3),[],1),'+r','markersize',20,...
            'markerfacecolor','r', 'linewidth',2);
        plot3(reshape(data_WCSt0{g}(:,i,:,1),[],1),...
            reshape(data_WCSt0{g}(:,i,:,2),[],1),...
            reshape(data_WCSt0{g}(:,i,:,3),[],1),'+b','markersize',20,...
            'markerfacecolor','b', 'linewidth',2);
        legend{1} = 'Points after';
        legend{2} = 'Points before';
        legend{3} = 'Truth after';
        legend{4} = 'Truth before';
        title(sprintf('Final %g paths referenced to i=%g',g+1,i));
        grid on
        axis([-15 15 -15 15 -5 5]);
    end
end
%}



% Nested function
    function [F, data_WCSt, data_WCSj] = nestedfun(P0)
        %disp(P0');
        %fprintf('\n iteration\n');
        [n_unique, n_S, ~] = size(data);
        color = jet(n_S);
        
        % Find registration and translation matrices
        S_R = cell(n_S,n_S);
        S_t = cell(n_S,n_S);
        ix = 1;
        for i = 1:n_S-1;
            for j = i+1:n_S;
                if ~R_empty(i,j);
                    rx = P0(ix);
                    ry = P0(ix+1);
                    rz = P0(ix+2);
                    S_R{i,j} = compose_rotation(rx, ry, rz);
                    S_t{i,j} = P0(ix+3:ix+5)';
                    ix = ix + 6;
                end
            end
        end
        
        % make S_R non-directed
        for i = 1:n_S;
            for j = i:n_S;
                if ~isempty(S_R{i,j});
                    S_R{j,i} = S_R{i,j}';
                    S_t{j,i} = -(S_R{i,j}')*S_t{i,j};
                end
            end
        end
        
        unique_1 = reshape(P0(ix:end),[(numel(P0)-ix+1)/3,3]);
        
        %F = nan(n_unique,n_S,n_S,3);
        data_WCSt = cell(n_g,1);
        data_WCSj = cell(n_g,1);
        for g = 1:n_g;
            data_WCSt{g} = nan(n_unique,n_S,n_S,3);
            data_WCSj{g} = nan(n_unique,n_S,n_S,3);
        end
        
        for g = 1:n_g;
            [n_perm, ~] = size(G_path{g});
            for p = 1:n_perm;
                path = G_path{g}(p,:);
                j = path(1);
                i = path(end);
                % Find effective R and t
                Reff = eye(3);
                teff = zeros(3,1);
                for k = 1:numel(path)-1;
                    if isempty(S_R{path(k+1),path(k)});
                        foo = 1;
                    end
                    Reff = S_R{path(k+1),path(k)}*Reff;
                    teff = S_R{path(k+1),path(k)}*(teff)+ S_t{path(k+1),path(k)};
                end
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
                % Mapped into l1 from l2
                unique_hat = (Reff*foo+repmat(teff,[1,n_valid]))';
                % Error between j and i
                %F(is_valid,i,j,:) = unique_1(is_valid,:) - unique_hat;
                data_WCSj{g}(is_valid,i,j,:) = unique_hat;
                data_WCSt{g}(is_valid,i,j,:) = unique_1(is_valid,:);
            end
            
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
        
        %
        %{
        for g = 1:n_g;
            for i = 1:n_S;
            if  sum(sum(~isnan(squeeze(data_WCSt{1}(:,i,:,1)))))== 0 &&...
               sum(sum(~isnan(squeeze(data_WCSt{2}(:,i,:,1)))))  == 0;
                continue
            end
            figure;
            hold on
            plot3(reshape(data_WCSj{g}(:,i,:,1),[],1),...
                reshape(data_WCSj{g}(:,i,:,2),[],1),...
                reshape(data_WCSj{g}(:,i,:,3),[],1),'ok','markersize',5,...
                'markerfacecolor','r')
            plot3(reshape(data_WCSt{g}(:,i,:,1),[],1),...
                reshape(data_WCSt{g}(:,i,:,2),[],1),...
                reshape(data_WCSt{g}(:,i,:,3),[],1),'+','markersize',20,...
                'markerfacecolor',[0 0 0], 'linewidth',2);
            title(sprintf('%g paths referenced to i=%g',g+1,i));
            grid on
            axis([-15 15 -15 15 -5 5]);
            end
        end
        %}
        
        F = [];
        for g  =1:n_g;
            diff = data_WCSj{g} - data_WCSt{g};
            diff = diff(:);
            diff = diff(~isnan(diff));
            F = [F ; diff];
        end
        
        %  F = reshape(F,[],1);
        %  F = F(~isnan(F));
        foo2 = 1;
    end
end
