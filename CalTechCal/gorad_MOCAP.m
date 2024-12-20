% main script to launch the estimation
% of the non-linear parameters by using the CalTech
% calibration toolbox and the output from the Svoboda's
% Multicamera self-calibration
%
% How to create the input data:
% 1) Run the MultiCamSelfCam
% 2) Run the MultiCamValidation
%
% $Id: gorad.m,v 2.0 2003/06/19 12:06:00 svoboda Exp $

clear all;

% TODO: isn't this unnecessary now?
addpath ../MultiCamSelfCal/Cfg
addpath ('../CommonCfgAndIO')

% Read configuration from whatever is specified on command-line (via --config=FILENAME)
config = read_configuration();

% if problem with desactivated images -> some problems with the estimation in general
desactivated_images = [];

idxcams = config.cal.cams2use;
selfcalib.goradproblem = 0;

for i = idxcams,
  [X_1,x_1] = preparedata(sprintf(config.files.points4cal,i));
  go_calib_optim_iter
  if any(isnan(param))
	  % when the iteration fails insert null distortion
	  % it is better than nonsense
	  kc(1:4) = [0,0,0,0];
	  selfcalib.goradproblem=1;
  else
	  %visualize_distortions
  end

  disp(sprintf('***** camera %d **********************************',i))
  %
  % save the intrinsic calibration data to a json file
  config.files.rad	= [config.paths.data,'Cam%d_intrinsics.json'];
  outputfile = sprintf(config.files.rad,i);

     fprintf(1,'\nExport of intrinsic calibration data to json file\n');
     % outputfile = input('File basename: ', 's');
     configfile = outputfile;
     disp(['Writing ' configfile]);

     fid = fopen(configfile, 'w');
        
     %{
        % Constructing the JSON structure as a string
        jsonString = sprintf([
            '{\n' ...
            '  "intrinsicMatrix": [\n' ...
            '    [%.16f, %.16f, %.16f],\n' ...
            '    [%.16f, %.16f, %.16f],\n' ...
            '    [%.16f, %.16f, %.16f]\n' ...
            '  ],\n' ...
            '  "distortionCoefficients": [\n' ...
            '    %.16f, %.16f, %.16f, %.16f\n' ...
            '  ]\n' ...
            '}\n'], ...
            KK(1,1), KK(1,2), KK(1,3), ...
            KK(2,1), KK(2,2), KK(2,3), ...
            KK(3,1), KK(3,2), KK(3,3), ...
            kc(1), kc(2), kc(3), kc(4));
     %}

        jsonString = sprintf([
            '{\n' ...
            '  "m11": %.7f,\n' ...
            '  "m12": %.7f,\n' ...
            '  "m13": %.7f,\n' ...
            '  "m21": %.7f,\n' ...
            '  "m22": %.7f,\n' ...
            '  "m23": %.7f,\n' ...
            '  "m31": %.7f,\n' ...
            '  "m32": %.7f,\n' ...
            '  "m33": %.7f,\n' ...
            '  "k1": %.7f,\n' ...
            '  "k2": %.7f,\n' ...
            '  "k3": %.7f,\n' ...
            '  "p1": %.7f,\n' ...
            '  "p2": %.7f\n' ...
            '}\n'], ...
            KK(1,1), KK(1,2), KK(1,3), ...
            KK(2,1), KK(2,2), KK(2,3), ...
            KK(3,1), KK(3,2), KK(3,3), ...
            kc(1), kc(2), 0.0, kc(3), kc(4)); %6th order radial distortion is not considered

        % Write the JSON string to the file
        fprintf(fid, '%s', jsonString);

        % Close the file
        status = fclose(fid);


%disp('Press any key to continue'),  pause

%%%
% clear already estimated parameters
clear fc kc alpha_c cc
end


