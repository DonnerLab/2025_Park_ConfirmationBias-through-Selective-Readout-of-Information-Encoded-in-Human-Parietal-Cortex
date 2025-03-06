function null=shuffled_mat_NN(n_iter)

% load parcel-vertex matrix

load('vertex_parcel_L.mat');

% random rotation (0-360 degrees)
x = rand(n_iter,1)*(2*pi);
y = rand(n_iter,1)*(2*pi);
z = rand(n_iter,1)*(2*pi);

Xq=double(vertex_parcel(:, 1));
Yq=double(vertex_parcel(:, 2));
Zq=double(vertex_parcel(:, 3));

null=nan(180, n_iter);

for iter = 1:n_iter
    
    g_null = vertex_parcel;
    R = makehgtform('yrotate',x(iter),'xrotate',y(iter),'zrotate',z(iter)); %the rotation matrix
    R = R(1:3,1:3); %trim to 3D
    g_null = g_null * R; %rotate
    
    %rotated spherical vertices
    X = double(g_null(:, 1));
    Y = double(g_null(:, 2));
    Z = double(g_null(:, 3));
    
    %get rotated data points on our original surface points
    for i=randperm(180) % go through the parcel in random order
        
        dist_rot=sqrt((Xq-X(i)).^2+(Yq-Y(i)).^2+(Zq-Z(i)).^2); % distance between each rotated parcel to the original parcel-vertex
        [~, dist_rot_min_idx]=sort(dist_rot);
        
        null(i, iter)=dist_rot_min_idx(1);
        
        % in case this parcel was already assigned
        k=1;
        while sum(null(:, iter)==dist_rot_min_idx(k), 'omitnan')>1
            null(i, iter)=dist_rot_min_idx(k+1);
            k=k+1;
        end
        
    end
end



