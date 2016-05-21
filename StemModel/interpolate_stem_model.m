function [ pp_z, pp_y, pp_r1, pp_r2, pp_seg ] = interpolate_stem_model( reg_z, reg_y, reg_r1, reg_r2, reg_seg, t_z_resample, t_z_interp)
%INTERPOLATE_STEM_MODEL Summary of this function goes here
%   Detailed explanation goes here

n_seg = numel(unique(reg_seg));
model_z = cell(n_seg,1);
model_y = cell(n_seg,1);
model_r = cell(n_seg,1);

for s = 1:n_seg
    ix = (reg_seg == s);
    if sum(ix)>1;
        zy = [reg_z(ix,:);reg_y(ix,:)];
        [~,I] = sort(zy(:,3));
        zy = zy(I,:);
        minz = min(zy(:,3));
        maxz = max(zy(:,3));
        zstep = linspace(minz,maxz,(maxz-minz)/t_z_resample);
        zysmooth = zeros(numel(zstep),3);
        rr = [reg_r1(ix,:);reg_r2(ix,:)];
        rr = rr(I);
        if numel(zstep)>1
            linear_interp_z = (minz:t_z_interp:maxz)';
            linear_interp_x = interp1(zy(:,3),zy(:,1),linear_interp_z);
            linear_interp_y = interp1(zy(:,3),zy(:,2),linear_interp_z);
            zy_interp = [linear_interp_x linear_interp_y linear_interp_z];
            PP = splinefit(zy_interp(:,3),zy_interp(:,1:2)',zstep);
            %PP = splinefit(zy(:,3), zy(:,1:2)',zstep);
            zysmooth(:,1:2) = ppval(PP,zstep)';
            zysmooth(:,3) = zstep;
            P = polyfit( zy(:,3), rr, 1);
            rsmooth = P(2)+ P(1)*zstep;
            model_z{s} = zysmooth;
            model_r{s} = rsmooth';
            model_y{s} = circshift(model_z{s},-1);
            model_y{s}(end,:) = model_y{s}(end-1,:);
        else
            model_z{s} = reg_z(ix,:);
            tempy = circshift(reg_z(ix,:),-1);
            model_y{s} = tempy;
            tempend = reg_y(ix,:);
            model_y{s}(end,:) = tempend(end,:);
            model_r{s} = reg_r1(ix,:);
        end
    end
end


    pp_z = [];
    pp_y = [];
    pp_r1 = [];
    pp_r2 = [];
    pp_seg = [];
    
    for s = 1:n_seg
        if ~isempty(model_r{s})
            pp_z = [pp_z;model_z{s}];
            pp_y = [pp_y;model_y{s}];
            pp_r1 = [pp_r1;model_r{s}];
            tempr =  circshift(model_r{s},-1);
            tempr(end) = model_r{s}(end);
            pp_r2 = [pp_r2;tempr];
            pp_seg = [pp_seg;repmat(s,size(model_r{s}))];
        end
    end
    
    %pp.z = pp_z;
    %pp.y = pp_y;
    %pp.r1 = pp_r1; 
    %pp.r2 = pp_r2;
    %pp.seg = pp_seg;
    
    %pp.tf = true(size(pp_r1));
    %plot_cylinder_color3( pp, pp.r1, 'radius'  )
    %hold on
    %scatter3(data.x,data.y,data.z,5,[.5 .5 .5],'filled')
  
    

end

