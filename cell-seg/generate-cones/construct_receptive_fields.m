function [RGC_fields, RGC_locs] = construct_receptive_fields(options)
    % inputs    
    arguments
        % neuron parameters
        options.ImageDimensions (1,1) {mustBeNumeric} = 200;  % size of input image being sampled (assumed to be square)
        options.RGC_sigma       (1,1) {mustBeNumeric} = 7;    % sigma of 2D gaussion of cones
        options.RGCHorSpacing   (1,1) {mustBeNumeric} = 20;   % spacing between RGCs in horizontal dim
        options.RGCTileRadius   (1,1) {mustBeNumeric} = 100;  % circular ring radius. Cells outside ring will not be simulated
        options.DoTranspose     (1,1) {mustBeNumericOrLogical} = false; % transpose the tiling array
        options.ProportionLtoM  (1,1) {mustBeNumeric} = 0.5;  % proportion of L to M cones
        options.ConeExcitiationFile   {mustBeText} = ''; 
        options.lambda          (1,1) {mustBeNumeric} = 605;
    end

    % returns 
    %  RGC_feilds - an [n_cells, ImageDimensions, ImageDimensions] matrix of
    %   that contains a gassuan in each image
    %  RGC_locs - a [n_cells, 3] matrix of the center of each cell, and its
    %   type (L == 1)
    %% construct RGC locations
    [X, Y] = meshgrid(0:options.RGCHorSpacing:(options.ImageDimensions+options.RGCHorSpacing+1));
    n = size(X,1);
       
    X = (sqrt(3) / 2) * X;
    Y = Y + repmat([0 round(options.RGCHorSpacing/2)],[n,n/2]);
    
    RGC_locs = zeros(numel(Y), 3);
    if options.DoTranspose
        RGC_locs(:, 1) = Y(:);
        RGC_locs(:, 2) = X(:);
    else
        RGC_locs(:, 1) = X(:);
        RGC_locs(:, 2) = Y(:);
    end

    temp_RGC_locs = [];
    for i=1:size(RGC_locs, 1)
        if sqrt((RGC_locs(i,1) - options.ImageDimensions/2)^2 ...
               +(RGC_locs(i,2) - options.ImageDimensions/2)^2) < options.RGCTileRadius
            temp_RGC_locs = [temp_RGC_locs; RGC_locs(i,:)];
        end
    end
    RGC_locs = temp_RGC_locs;
    n_cells = size(RGC_locs, 1);

    % define type
    RGC_locs(rand(n_cells,1)>options.ProportionLtoM,3) = 1;

    %% load activation 
    if numel(options.ConeExcitiationFile) ~= 0
        conesensitivity = readmatrix(options.ConeExcitiationFile);
        Lcone = interp1(conesensitivity(:,1),conesensitivity(:,2),options.lambda);
        Mcone = interp1(conesensitivity(:,1),conesensitivity(:,3),options.lambda);

    else
        warning  'No cone spectra provided, assuming same cone type'
        Lcone = 1;
        Mcone = 1;
    end

    %% construct RGC gaussians
    RGC_fields = zeros(n_cells, options.ImageDimensions, options.ImageDimensions);
    
    x=1:options.ImageDimensions; y=x;
    [X,Y] = meshgrid(x,y);
    for i=1:n_cells
        if RGC_locs(i,3) == 1
            RGC_fields(i, :,:) = Lcone*gauss2d(X,Y,RGC_locs(i,:), options.RGC_sigma);
        else
            RGC_fields(i, :,:) = Mcone*gauss2d(X,Y,RGC_locs(i,:), options.RGC_sigma);            
        end
    end

    
end
