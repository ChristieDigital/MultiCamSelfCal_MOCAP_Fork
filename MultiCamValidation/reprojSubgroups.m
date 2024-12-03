% do 3d reconstruction using selected camera subgroups

clear all;

% add necessary paths
addpath ../CommonCfgAndIO
addpath ./CoreFunctions
addpath ./InputOutputFunctions
addpath ../RansacM; % ./Ransac for mex functions (it is significantly faster for noisy data)

% preliminary parse just to grab the camera indexes shared between subgroups for reprojection

%{
arg_list = argv();
for i = 1:numel(arg_list)
  printf ("Argument %d: %s\n", i, arg_list{i});
endfor
%printf(arg_list{2})
%}

arg_list = argv(); % Get command-line arguments
for i = 1:numel(arg_list)
  printf("Argument %d: %s\n", i, arg_list{i}); % Print each argument
end

% Check for '--indexes=' argument
found_indexes = 0;
shared_indexes = []; % Initialize in case it isnt found
for i = 1:numel(arg_list)
  arg = arg_list{i};
  if length(arg) >= 10 && strcmp(arg(1:10), '--indexes=')
    found_indexes = 1;
    shared_indexes = arg(11:end); 
    shared_indexes = str2num(shared_indexes); % Convert to numeric array
    break;
  end
end

% Validate arguments
if ~found_indexes
  error('Missing --indexes=[shared camera trio indexes] command-line argument');
elseif isempty(shared_indexes) || length(shared_indexes) ~= 2 || any(isnan(shared_indexes))
  error('Invalid format for --indexes. Expected format: --indexes=x,y where x and y are numbers.');
end

disp('Shared indexes:');
disp(shared_indexes);

% Further checks after parsing
if isempty(shared_indexes)
  error('empty shared indexes')
end

disp(['Shared indexes: ', mat2str(shared_indexes)]);

% Read configuration from whatever is specified on command-line (via --config=FILENAME)
config = read_configuration();




% check if empty matrix
if isempty(shared_indexes)
  error('empty shared indexes');string.Join(",", shared_indexes)
end

UNDO_RADIAL = logical(config.cal.UNDO_RADIAL | config.cal.UNDO_HEIKK);

if UNDO_RADIAL
	% add functions dealing with radial distortion
	addpath ../RadialDistortions
end

% read the input data
loaded = loaddata(config);
linear = loaded;		% initalize the linear structure

CAMS = size(config.cal.cams2use,2);
FRAMES = size(loaded.IdMat,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See the README how to compute data
% for undoing of the radial distortion
if config.cal.UNDO_RADIAL
  for i=1:CAMS,
	[K,kc] = readradfile(sprintf(config.files.rad,config.cal.cams2use(i)));
	xn	   = undoradial(loaded.Ws(i*3-2:i*3,:),K,[kc,0]);
	linear.Ws(i*3-2:i*3,:) = xn;
  end
elseif config.cal.UNDO_HEIKK,
  for i=1:CAMS,
	heikkpar = load(sprintf(config.files.heikkrad,config.cal.cams2use(i)),'-ASCII');
	xn = undoheikk(heikkpar(1:4),heikkpar(5:end),loaded.Ws(i*3-2:i*3-1,:)');
	linear.Ws(i*3-2:i*3-1,:) = xn';
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detection of outliers
% RANSAC is pairwise applied
disp('RANSAC validation step running ...');

inl.IdMat = findinl(linear.Ws,linear.IdMat,config.cal.INL_TOL);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% fill cam(i) structures
for i=1:CAMS,
  cam(i).camId     = config.cal.cams2use(i);
  cam(i).ptsLoaded = find(loaded.IdMat(i,:)); % loaded structure
  cam(i).ptsInl	   = find(inl.IdMat(i,:));	% survived initial pairwise validation
  cam(i).xgt	   = loaded.Ws(3*i-2:3*i,cam(i).ptsLoaded);
  cam(i).xlin	   = linear.Ws(3*i-2:3*i,cam(i).ptsLoaded);
  cam(i).xgtin	   = loaded.Ws(3*i-2:3*i,cam(i).ptsInl);
  cam(i).P		   = loaded.Pmat{i};
  [cam(i).K, cam(i).R, cam(i).t, cam(i).C] = P2KRtC(cam(i).P);
end

disp('***********************************************************')
disp('Computing a robust 3D reconstruction via camera sampling ...')
% compute a reconstruction robustly

t1 = cputime;
reconstructed = estimateXSpecificCams(linear,inl.IdMat,cam,config,shared_indexes);
reconstructed.CamIds = config.cal.cams2use(reconstructed.CamIdx);
t2 = cputime;
disp(sprintf('Elapsed time for 3D computation: %d minutes %d seconds',floor((t2-t1)/60), round(mod((t2-t1),60))))

% compute reprojections
for i=1:CAMS,
  xe		= linear.Pmat{i}*reconstructed.X;
  cam(i).xe	= xe./repmat(xe(3,:),3,1);

  % these points were the input into Martinec and Pajdla filling
  mask.rec = zeros(1,FRAMES);	% mask of points that survived validation so far
  mask.vis = zeros(1,FRAMES); % maks of visible points
  mask.rec(reconstructed.ptsIdx)  = 1;
  mask.vis(cam(i).ptsLoaded) = 1;
  mask.both			   = mask.vis & mask.rec; % which points are visible and reconstructed for a particular camera
  unmask.rec			   = cumsum(mask.rec);
  unmask.vis			   = cumsum(mask.vis);
  cam(i).recandvis = unmask.rec(~xor(mask.rec,mask.both) & mask.rec);
  cam(i).visandrec = unmask.vis(~xor(mask.rec,mask.both) & mask.rec);
  cam(i).err2d	 = sqrt(sum([cam(i).xe(1:2,cam(i).recandvis) - cam(i).xlin(1:2,cam(i).visandrec)].^2));
  cam(i).mean2Derr = mean(cam(i).err2d);
  cam(i).std2Derr  = std(cam(i).err2d);
end

for i=1:CAMS,
  xe = loaded.Ws(i*3-2:i*3, reconstructed.ptsIdx(logical(loaded.IdMat(i,reconstructed.ptsIdx))));
  Xe = reconstructed.X(:, logical(loaded.IdMat(i,reconstructed.ptsIdx)));
  corresp = [Xe',xe'];
  save(sprintf(config.files.points4cal,config.cal.cams2use(i)),'corresp','-ascii');
end