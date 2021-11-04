folder = "D:\Video_Folder\Alex2\11_3_2021_15_0_4\";
vidmats = {};

offsetx = -5;
offsety = 1;

indy = (196 + offsety):((192+80) + offsety);
indx = (135+offsetx):((128+80) + offsetx);
for i = 1:11
    fname = sprintf('%sAlex2_%003d.avi',folder,i);
    hvid = VideoReader(fname);
    vidmat = squeeze(hvid.read());
    vidmat_crop = vidmat(indy, indx, :);
    vidmats{i} = vidmat_crop(:,:,2:end);
end

%%
mm = [];
ss = [];
for i = 1:11
    mm = [mm, mean(double(vidmats{i}(:)))];
    ss = [ss, std(double(vidmats{i}(:)))];
end
%%
ff = @(a,b,c,x) a ./ (1+exp(-b*(x + c))); 
L = @(params,x,y) norm(ff(params(1), params(2), params(3),x) - y)^2;
loss_factory = @(x,y) @(params) L(params,x,y);
fun = loss_factory(xx,mm);
%%
% D:\Video_Folder\Alex2\11_3_2021_15_0_4\
% resulted in x = [135.476847369744,9.38609473656521,-0.295548692871276]
x0 = [1,1,0];
x = fminsearch(fun,x0)

figure; plot(xx, ff(x(1), x(2), x(3),xx))
hold on;
plot(xx,mm)

%% 
inverse_map = @(params,y) (-log((params(1) ./ y) -1) ./ params(2)) - params(3);