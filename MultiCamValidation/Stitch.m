% stitch together two calibration results into a unified global coordinate system

clear all;

%pkg load statistics
%pkg load signal

% add necessary paths
addpath ../CommonCfgAndIO
addpath ./CoreFunctions
addpath ./InputOutputFunctions
addpath ../CalTechCal
addpath ../MultiCamSelfCal/CoreFunctions
addpath ../MultiCamSelfCal/OutputFunctions


arg_list = argv(); % Get command-line arguments
for i = 1:numel(arg_list)
  printf("Argument %d: %s\n", i, arg_list{i}); % Print each argument
end

found_cfg1 = 0;
found_cfg1 = 0;
filename = cell(1, 2); % Initialize a cell array to store two filenames.
shared_indexes = cell(1, 2); % Initialize a cell array to store two sets of shared indexes.
disable_plots = 0;

for i = 1:numel(arg_list)
  arg = arg_list{i};
  if length(arg) >= 10 && strcmp(arg(1:10), '--config1=')
    found_cfg1 = 1;
    filename{1} = arg(11:end); 
  elseif length(arg) >= 11 && strcmp(arg(1:11), '--indexes1=')
	found_indexes1 = 1;
    shared_indexes1 = arg(12:end)
    shared_indexes{1} = str2num(shared_indexes1); % Convert to numeric array
  elseif length(arg) >= 10 && strcmp(arg(1:10), '--config2=')
    found_cfg2 = 1;
    filename{2} = arg(11:end);
  elseif length(arg) >= 11 && strcmp(arg(1:11), '--indexes2=')
	found_indexes2 = 1;
    shared_indexes2 = arg(12:end)
    shared_indexes{2} = str2num(shared_indexes2); % Convert to numeric array
   % Check for --disable-plots
  elseif strcmp(arg, '--disable-plots')
    disable_plots = 1;
  end
end

disp(shared_indexes)

% since we are looking at the result of two calibrations at once, 
% we can;t use the normal config loading routine like in gocal/gorec and must set the fields manually
function config = set_points_and_pmat(filename, num_cameras)
  % Check if filename is provided
  if nargin < 1 || isempty(filename)
    error('Filename argument is required');
  end

  % Check if num_cameras is provided
  if nargin < 2 || isempty(num_cameras)
    error('Number of cameras argument is required');
  end

  % Initialize the config structure
  config = struct();
  config.files = struct();

  % Determine the directory name from the filename
  if isfolder(filename)
    config_dirname = filename;
  else
    config_dirname = fileparts(filename);
  end

  if (config_dirname(end) ~= '/')
    config_dirname = strcat(config_dirname, '/');
  end

  % Set config.paths.data
  config.paths.data = config_dirname;

  config.cal.cams2use = [1:num_cameras]; % hardcoded with 3 cameras as num_cameras for now

  % Generate the paths for each camera's Pmat.cal file
  config.files.CalPmat	= [config.paths.data,'camera%d.Pmat.cal'];

  config.files.Cst		= [config.paths.data,'Cst.dat'];

  % Generate the paths for each camera's points4cal.dat file
  config.files.points4cal = cell(1, num_cameras);
  for i = 1:num_cameras
    config.files.points4cal{i} = sprintf('%scam%d.points4cal.dat', config.paths.data, i);
  end
end

p = cell(1, 2);
x = cell(1, 2);
persistent_proj_matrices = []; % combined proj matrices from first calibration group- remain unchanged
combined_proj_matrices = [];
configArray = cell(1,2);

% we only want 3d points from each .dat where the 2d points are shared between both files
for i = 1:2
	configArray{i} = set_points_and_pmat(filename{i}, 3);
    p{i} = load(configArray{i}.files.points4cal{shared_indexes{i}(2)}); % Need to retrieve the points from one of the shared cameras from each trio- shouldnt matter which one
end

% Extract the x, y columns
xy1 = p{1}(:, 5:6); % x, y are columns 5 and 6 
xy2 = p{2}(:, 5:6);

% Find rows with matching x, y points
matches1 = ismember(xy1, xy2, 'rows'); % Check if each row of xy1 exists in xy2
matches2 = ismember(xy2, xy1, 'rows'); % Check if each row of xy2 exists in xy1

% Retain only rows with matching x, y points
p{1} = p{1}(matches1, :);
p{2} = p{2}(matches2, :);


for i = 1:2
    config = configArray{i}
    %p{i} = load(config.files.points4cal{shared_indexes{i}(2)}); % Need to retrieve the points from one of the shared cameras from each trio- doesn't matter which one
    p{i} = p{i}(:,1:3);
    %p{i} = p{i}';
    for j=1:3
        calPmatFile = sprintf(config.files.CalPmat, j);
        %fprintf('Processing CalPmat file: %s\n', calPmatFile);
        %calPmatFile = load(calPmatFile);
        %calPmatFile = config.files.CalPmat{j};
        Pmat = dlmread(calPmatFile, '', 5, 0); % Skip the first 5 lines, seems to be a problem with the metadata
        [c1Cam(j).K, c1Cam(j).R, c1Cam(j).t, c1Cam(j).C] = P2KRtC(Pmat);
        c1Cam(j).K
        if any(j == shared_indexes{i})
            p{i} = [p{i};c1Cam(j).C']; % append cam centers to point data- just shared ones? maybe unnecessary
        end
        if (i == 1)
            persistent_proj_matrices = [persistent_proj_matrices;Pmat]; % combine all pmatrices from the first calibration group 
		elseif (i == 2)
			combined_proj_matrices = [combined_proj_matrices;Pmat]; % combine all pmatrices from the second calibration group for alignment
        end
    end
 
end

if ~disable_plots
    figure(50),
    clf

    plot3(p{1}'(1,:),p{1}'(2,:),p{1}'(3,:),'*','Color', 'r');
    hold on;
    plot3(p{2}'(1,:),p{2}'(2,:),p{2}'(3,:),'*','Color', 'b');
    grid on
    pause
end



% get the similarity transform between point sets + shared camera centers
[align.simT.s, align.simT.R, align.simT.t] = estsimt(p{2}', p{1}');

align.simT

for i=1:2
    p{i} = p{i}';
    p{i}(4,:) = 1; % homogenize
    size(p{i})
end

% align 2 to 1 based on transform- apply to cameras and point set at the same time
[align.P, align.X] = align3d(combined_proj_matrices,p{2},align.simT);

align.P

% p2 is now aligned to p1
p{2} = align.X;


if ~disable_plots
    figure(100),
    clf
    plot3(p{1}(1,:),p{1}(2,:),p{1}(3,:),'*','Color', 'r');
    hold on;
    plot3(p{2}(1,:),p{2}(2,:),p{2}(3,:),'*','Color', 'b');
    grid on
    pause
end

disp(class(configArray{i}.files.CalPmat)); % Confirm the data type of config.files.CalPmat
disp(class(configArray{i}.cal.cams2use));
%[align.Cst,align.Rot] = savecalpar(align.P,configArray{2});

% save the calibration parameters for each subgroup (only one group gets modified, but we save to json as final step here)
%config = configArray{2};

for j=1:size(configArray,2)
	config = configArray{j};
   
    if (j == 1)
	    P = align.P;
    elseif (j== 2)
		P = persistent_proj_matrices;
	end
	idxused = config.cal.cams2use;
	CAMS = size(P,1)/3;


	Cst = zeros(CAMS,3);
	Pst = zeros(3*CAMS,3);
	Rot = [];
	for i=1:CAMS,
        
        outputFile = [config.paths.data,'Cam%d_extrinsics.json'];
        outputFile = sprintf(outputFile,idxused(i))
	    Pmat = P(i*3-2:i*3,:);
        % construct as 4x4, leave last row as 0,0,0,1
        Pmat = [Pmat;0,0,0,1];
        sc = norm(P(i*3,1:3));
        % first normalize the Projection matrices to get normalized pixel points
        P(i*3-2:i*3,:) = P(i*3-2:i*3,:)./sc;
        % decompose the matrix by using rq decomposition
        [K,R] = rq(P(i*3-2:i*3,1:3));
        tvec= inv(K)*P(i*3-2:i*3,4); % translation vector
        
        K

        fid = fopen(outputFile, 'w');

        jsonString = sprintf(['{\n', ...
            '  "m11": %.7f,\n', ...
            '  "m12": %.7f,\n', ...
            '  "m13": %.7f,\n', ...
            '  "m14": %.7f,\n', ...
            '  "m21": %.7f,\n', ...
            '  "m22": %.7f,\n', ...
            '  "m23": %.7f,\n', ...
            '  "m24": %.7f,\n', ...
            '  "m31": %.7f,\n', ...
            '  "m32": %.7f,\n', ...
            '  "m33": %.7f,\n', ...
            '  "m34": %.7f,\n', ...
            '  "m41": %.7f,\n', ...
            '  "m42": %.7f,\n', ...
            '  "m43": %.7f,\n', ...
            '  "m44": %.7f\n', ...
            '}'], ...
            R(1,1), R(1,2), R(1,3), tvec(1), ...
            R(2,1), R(2,2), R(2,3), tvec(2), ...
            R(3,1), R(3,2), R(3,3), tvec(3), ...
            0, 0, 0, 1);

        fprintf(fid, '%s', jsonString);


        % Close the file
        status = fclose(fid);
	end
end


%{
P = align.P;
idxused = config.cal.cams2use;
CAMS = size(P,1)/3;

Cst = zeros(CAMS,3);
Pst = zeros(3*CAMS,3);
Rot = [];
for i=1:CAMS,
  % do not save P matrices in separate files
  if 1
	Pmat = P(i*3-2:i*3,:);
    save(sprintf(config.files.CalPmat,idxused(i)),'Pmat');
  end

  sc = norm(P(i*3,1:3));
  % first normalize the Projection matrices to get normalized pixel points
  P(i*3-2:i*3,:) = P(i*3-2:i*3,:)./sc;
  % decompose the matrix by using rq decomposition
  [K,R] = rq(P(i*3-2:i*3,1:3));
  tvec= inv(K)*P(i*3-2:i*3,4);			% translation vector
  C	  = -R'*tvec;
  Psteph		   = R'*inv(K);
  Pst(i*3-2:i*3,:) = Psteph;
  Cst(i,:)		   = C';
  Rot	 = [Rot;R];
  %save(config.files.Pst,'Pst');
  save(config.files.Cst,'Cst');
end

%}

%{
jsonString = sprintf([

'{\n' ...
	'  "m11": %.7f,\n' ...
	'  "m12": %.7f,\n' ...
	'  "m13": %.7f,\n' ...
	'  "m14": %.7f,\n' ...
	'  "m21": %.7f,\n' ...
	'  "m22": %.7f,\n' ...
	'  "m23": %.7f,\n' ...
	'  "m24": %.7f,\n' ...
	'  "m31": %.7f,\n' ...
	'  "m32": %.7f,\n' ...
	'  "m33": %.7f,\n' ...
	'  "m34": %.7f,\n' ...
	'  "m41": %.7f,\n' ...
	'  "m42": %.7f,\n' ...
	'  "m43": %.7f,\n' ...
	'  "m44": %.7f\n' ...
	'}\n'], ...
	Pmat(1,1), Pmat(1,2), Pmat(1,3), Pmat(1,4), ...
	Pmat(2,1), Pmat(2,2), Pmat(2,3), Pmat(2,4), ...
	Pmat(3,1), Pmat(3,2), Pmat(3,3), Pmat(3,4), ...
	Pmat(4,1), Pmat(4,2), Pmat(4,3), Pmat(4,4));
        
	fprintf(fid, '%s', jsonString);
%}