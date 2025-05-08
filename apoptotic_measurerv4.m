clc;
clear variables;
close all;

currentdir = pwd;
addpath(pwd);
main = uigetdir();
[filepath, name] = fileparts(main);

pouch_dir = [main, '/mask_pouch'];
cherry_dir = [main, '/mask_cherry'];
signal_dir = [main, '/avg_basal'];
if exist([main, '/dcp1_binary'],'dir') == 0
    mkdir(main, '/dcp1_binary');
end
binary_dir = [main, '/dcp1_binary'];

cd(signal_dir);
files = dir('*.tif');
numberfiles= numel(files);
summary =  zeros(numberfiles, 13);
resolution = 1.76;
conversion = (1 /(resolution^2));


for i=1:numberfiles
    
    % read this https://stackoverflow.com/questions/3008326/background-subtracting-in-matlab
    
    currentfile= [num2str(i),'.tif'];
    cd(pouch_dir);
    % se= strel("diamond",10);
    mask_pouch= imbinarize(imread(currentfile));
    % mask_pouch=imerode(mask_pouch,se);
    pouch_region = regionprops(mask_pouch);
    cd(cherry_dir);
    mask_cherry= imbinarize(imread(currentfile));
    % mask_cherry=imerode(mask_cherry,se);
    cherry_region = regionprops(mask_cherry);
    cd(signal_dir);
    disc = bfopen(currentfile);
    disc_image= disc{1,1};
    dcp1_signal= disc_image{3};
    
    dcp1_signal2 = imgaussfilt(dcp1_signal,1);
    dcp1_signal3 = imtophat(dcp1_signal2, offsetstrel('ball', 5,5));
    % test= imbinarize(dcp1_signal3);
    % background= imsubstract(dcp1_signal3,test);
    background = imopen(dcp1_signal3, strel('ball', 10, 10));
    % dcp1_signal3 = imtophat(dcp1_signal2, strel('disk',25));
    dcp1_signal3 = imsubtract(dcp1_signal3, 0.3*background);
    dcp1_signal3(mask_pouch==0)=0;
    % dcp1_signal3 = imadjust(dcp1_signal3);
    th= graythresh(dcp1_signal3);
    dcp1bw= imbinarize(dcp1_signal3, th);
    dcp1bw= imdilate(dcp1bw,strel("disk",1));
    % dcp1bw = bwareaopen(dcp1bw, 3);

    imwrite(dcp1bw,"matlab.tif")
    

    image = figure;
    imshowpair(dcp1_signal, dcp1bw);


    % th = graythresh(dcp1_signal2(mask_pouch==1));
    % th = graythresh(dcp1_signal3);
    % dcp1bw2= imbinarize(dcp1_signal3,th);
    % dcp1bw2= imdilate(dcp1bw2,strel("disk",1));
    % dcp1bw2 = bwareaopen(dcp1bw2, 3);
    % image = figure;
    % imshowpair(dcp1_signal, dcp1bw2);

    hold on;
    cd(binary_dir);
    image_name = [num2str(i),'_dcp1_signal.tif'];
    print(image, '-dtiff', '-r150', image_name);
    close all

    % dcp1_signal2 = imgaussfilt(dcp1_signal,1);
    % 
    % background = imopen(dcp1_signal2, offsetstrel('ball', 25, 25));
    % dcp1_signal2 = imsubtract(dcp1_signal2, background);
    % dcp1bw= imbinarize(dcp1_signal2);
    % dcp1bw = bwareaopen(dcp1bw, 5);
    

    summary(i,1)= i;
    summary(i,2)= conversion * pouch_region.Area;
    summary(i,4) = conversion * cherry_region.Area;
    summary(i,3) =  summary(i,2) - summary(i,4);

    pouchbw =  dcp1bw;
    pouchbw(mask_pouch==0) = 0;
    bwant = pouchbw;
    bwant(mask_cherry==1) = 0;
    dcp1_ant = regionprops(bwant);
    summary(i, 5) = size(dcp1_ant, 1); % calculate total number of clusters
    summary(i, 7) = conversion * sum([dcp1_ant.Area]); % calculate total apoptotic area

    bwpos = pouchbw;
    bwpos(mask_cherry==0) = 0;
    dcp1_pos = regionprops(bwpos);
    summary(i, 6) = size(dcp1_pos, 1); % calculate total number of clusters
    summary(i, 8) = conversion * sum([dcp1_pos.Area]); % calculate total apoptotic area

    summary(i,9) =  summary(i,5) / summary(i,3);
    summary(i,10) =  summary(i,6) / summary(i,4);
    summary(i,11) =  (summary(i,7) / summary(i,3)) * 100; % calculate % of anterior death area
    summary(i,12) =  (summary(i,8) / summary(i,4)) *100; % calculate % of posterior death area
    summary(i,13) =  summary(i,12) - summary(i,11); % calculate pouch size ratio
    % summary(i,14) = (max(0.025,summary(i,12))) / (max(0.025,summary(i,11)));  % calculate death ratio
    cd(signal_dir);

end

cd(main);
results = array2table(summary);
results.Properties.VariableNames = {'Disc', 'Area pouch', 'Anterior Pouch' 'Posterior Pouch', 'Anterior Cells',...
    'Posterior Cells', 'Anterior death area', 'Posterior death area', 'Anterior index', 'Posterior Index',...
    '% Anterior apoptotic area', '% Posterior apoptotic area', 'death difference'};
writetable(results,[name,'_results.csv']);



close all;
clear variables;
clc;