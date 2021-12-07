settings.lambda             = 565;   % wavelength of light imageing at [nm]
settings.image_dimensions   = 240;   % size of image [pixels]
settings.reconstruction_dim = 240;   % size of S (neural "E" reconstruction)
settings.rgc_spacing        = 15;    % spacing between RGC
settings.proportionLM       = 0.5;   % proportion of L to M cones
settings.rgc_tile_radius    = 300; 
settings.RGC_sigma          = 1/(2*sqrt(2*log(2)))*settings.rgc_spacing;  % RGC receptive field sigma [pixels]
settings.conespectrafile    = 'linss2_10e_fine.csv'; 
settings.binary_colormap    = (linspace(1,0,11))'*ones(1,3); 
settings.threshold          = 3.1e-2;

[RGC_field, RGC_loc] = construct_receptive_fields( ...
    'ImageDimensions',      settings.image_dimensions, ...
    'RGC_sigma',            settings.RGC_sigma, ...
    'RGCHorSpacing',        settings.rgc_spacing, ...
    'RGCTileRadius',        settings.rgc_tile_radius, ...
    'ProportionLtoM',       settings.proportionLM, ...
    'ConeExcitiationFile',  settings.conespectrafile, ...
    'lambda',               settings.lambda ...
);

RGC_field(RGC_field<settings.threshold) = 0;

cones = 1e3*(1- squeeze(sum(RGC_field, 1)));
cones = (cones-min(cones(:)))/(max(cones(:))-min(cones(:)));

% plot and save receptive fields 
f = figure(1);
colormap(settings.binary_colormap);
imagesc(cones); 

axis equal off
xlim([0 , settings.image_dimensions])
ylim([0 , settings.image_dimensions])