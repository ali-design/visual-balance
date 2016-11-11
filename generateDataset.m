function generateDataset(dataPath,N)

dataPath = [dataPath,'/'];
dirDataPath = dir([dataPath,'*.png']);

if exist('dataset.mat', 'file')
    dataset = load('dataset.mat');
    dataset = dataset.dataset;
else
    dataset = {};
end

a = size(dataset,2);

%I = size(dirDataPath,1);
I = a + N;
a = a + 1;
tic
parfor c = a:I
    %[pt,name,ext] = fileparts(dirDataPath(c).name);
    salImg = imread(dirDataPath(c).name);
    salImg = salImg(:,:,1);
    [Ys,Xs] = size(salImg);
    salQV = reshape(salImg,Ys*Xs,1);
    salQV = uint16(floor(double(salQV)/32));
    %salQV = uint8(salQV);
    
    [X,Y] = meshgrid(1:Xs,1:Ys);
    pos = uint16(([X(:),Y(:)]));
    points = [pos,salQV];
    idx = salQV~=0;
    points = points(idx,:);
    img = {};
    img.salValue = uint8(points(:,3));
    img.pixels = points(:,1:2);
    img.w = uint16(Xs);
    img.h = uint16(Ys);
    dataset{c} = img;
end

fprintf('\n Saving dataset to dataset.mat ... \n');
save('dataset.mat','dataset');
toc
% plot dataset
% plotDatasetScatterplot(dataset);
