clc; clear all; close all;
setenv('PATH', getenv('PATH')+":/usr/local/bin") % needed for finding ffmpeg to convert avi to gif

rootdir = '~/Desktop/cell-seg/data';
filelist = dir(fullfile(rootdir, '**/*1.avi'));

%% define constants
offsetx = -5;
offsety = 1;

%% load img
% for i = 1:size(filelist,1)
% if filelist(i).bytes == 0
%     continue
% end

fileidx = 46;% 36
mov=VideoReader([filelist(fileidx).folder '/' filelist(fileidx).name]);

frames = squeeze(mov.read());
img = mean(frames,3);

%% crop
img = img(192 + offsety:(192+127) + offsety, 128+offsetx:(128+255) + offsetx);
img = (img - min(img(:))) / (max(img(:)-min(img(:))));
I = img;

%% BEGIN SEGMENTATION
thresh = 0.0001e1;
cone_radius = 5; % radius of circular mask

fI = imgaussfilt(I,5);
p = FastPeakFind(imgaussfilt(fI,2), thresh);
loc = [p(1:2:end),p(2:2:end)];

segments = zeros(size(loc,1),size(fI,1),size(fI,2)); % scalar (nlocs,fI x, fI y)
for k = 1:size(loc,1)
    circlemask = createCirclesMask(fI,loc(k,:), cone_radius);
    conemask = circlemask;
    segments(k,:,:) = conemask(:,:);
end

figure(1)
colormap gray
subplot(2,1,1)
imagesc(fI); hold on
scatter(p(1:2:end),p(2:2:end),'r+')
hold off; 
axis off
subplot(2,1,2)
imagesc(squeeze(mean(segments,1)));
axis off
saveas(gcf,[filelist(fileidx).folder '/' filelist(fileidx).name '.eps'])


%% plot over time
% get all vid files 
vidfilelist = dir(fullfile(filelist(fileidx).folder, '*.avi'));

celldecay = zeros(size(segments,1),size(vidfilelist,1));
avemovie = zeros(size(vidfilelist,1),size(I,1),size(I,2));
for v = 1:size(vidfilelist,1)
    movitr=VideoReader([vidfilelist(v).folder '/' vidfilelist(v).name]);

    frames = squeeze(movitr.read());
    img = mean(frames,3);
    
    % crop
    img = img(192 + offsety:(192+127) + offsety, 128+offsetx:(128+255) + offsetx);
    %img = (img - min(img(:))) / (max(img(:)-min(img(:))));
    I_ = img;
    avemovie(v,:,:) = I_;

    celldecay(:,v) = segments(:,:)*I_(:);
end


%% plot decays
close 2
figure(2);
deltadecay = zeros(size(celldecay));
for i = 1:size(celldecay,1)
deltadecay(i,:) = celldecay(i,:)-celldecay(i,1);
plot(celldecay(i,1:(end-10))-celldecay(i,1),'LineWidth',2); hold on
%     plot(celldecay(i,1:(end-1)),'LineWidth',2); hold on
end
xlabel('iteration','FontSize',17)
ylabel('\Delta Intensity (a.u.)','FontSize',17)
set(gca,'FontSize',20)
set(gcf,'PaperType','A4')
% matlab.graphics.internal.setPrintPreferences('DefaultPaperPositionMode','manual')
% set(groot,'defaultFigurePaperPositionMode','manual')
saveas(gcf,[filelist(fileidx).folder '/' filelist(fileidx).name '_decay.eps'],'epsc')
% end

%% decision rule

timestep2classify = 8;
type1 = deltadecay(:,timestep2classify)>-1600;
type2 = deltadecay(:,timestep2classify)<-1600;

figure(3)
for i = 1:size(celldecay,1)
    if type1(i)
        plot(celldecay(i,1:(end-10))-celldecay(i,1),'g','LineWidth',2); hold on
    else
        plot(celldecay(i,1:(end-10))-celldecay(i,1),'r','LineWidth',2); hold on
    end
end
xlabel('iteration','FontSize',17)
ylabel('\Delta Intensity (a.u.)','FontSize',17)
set(gca,'FontSize',20)
set(gcf,'PaperType','A4')
% matlab.graphics.internal.setPrintPreferences('DefaultPaperPositionMode','manual')
% set(groot,'defaultFigurePaperPositionMode','manual')
saveas(gcf,[filelist(fileidx).folder '/' filelist(fileidx).name '_decay_c.eps'],'epsc')

figure(4)
imagesc(squeeze(avemovie(1,:,:))); hold on; colormap gray; axis image; axis off
scatter(loc(type1,1),loc(type1,2),'g.','AlphaData',0.5,SizeData=400,LineWidth=30)
scatter(loc(type2,1),loc(type2,2),'r.','AlphaData',0.5,SizeData=400,LineWidth=30)
saveas(gcf,[filelist(fileidx).folder '/' filelist(fileidx).name '_class.eps'],'epsc')


%% save gif
savegif(avemovie,['figures/' filelist(fileidx).name '_ave_movie.gif'])


function savegif(movie, fname)    
    % convert movie to gif and save to fname
    % parameters:
    % movie - scalar (n frames, vert. pix, horiz.pix)
    % fname - string

    mmax = max(movie, [], 'all');
    mmin = min(movie, [], 'all');
    nmovie = uint8(256*(movie - mmin)/(mmax-mmin));
    
    fps = 5;

    v = VideoWriter(fname);
    v.FrameRate = fps;
    open(v);
    for i = 1:size(nmovie,1)
       writeVideo(v,squeeze(nmovie(i,:,:)));
    end
    close(v);
    % relies on ffmpeg, install with brew or something if you dont have it
    cmd = sprintf('ffmpeg -y -i %s.avi -filter_complex "fps=%d,scale=320:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 %s',fname,fps,fname);
    system(cmd)
    system(sprintf('rm %s.avi',fname))
end