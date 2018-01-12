% Dimensionality of the simulated problem (2 for 2D, 3 for 3D)
n_dims = 2;

% Starting frame for detections 
s_frame = 15;

% Number of detections per frame
n_dedFrame = length(XY_pts{s_frame});

% Plotting detection points


set(gca,'Ydir','reverse')
hold on
for i = 1:300
    str = num2str(i);
    for j = 1:size(XY_pts{i},1)
        pos = XY_pts{i}(j,:);
        plot(pos(2), pos(1), '.')
        
    end
end

max_linking_distance = 50;
max_gap_closing = Inf;
debug = true;
[ tracks adjacency_tracks ] = simpletracker(XY_pts,...
    'MaxLinkingDistance', max_linking_distance, ...
    'MaxGapClosing', max_gap_closing, ...
    'Debug', debug);

n_tracks = numel(tracks);
colors = hsv(n_tracks);
all_points = vertcat(XY_pts{:});

for i = 1:n_tracks
    
    track  = adjacency_tracks{i};
    track_pts = all_points(track,:);
    
    plot(track_pts(:,2), track_pts(:,1), 'Color', colors(i,:))
end

%https://www.mathworks.com/matlabcentral/fileexchange/34040-simple-tracker