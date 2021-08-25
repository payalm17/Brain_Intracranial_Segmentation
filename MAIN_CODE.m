% Brain Intracranial Bleed Segmentation 
% Main Code

%% Reading Volume Files
tic
[img,imginfo]=load_3d_image_dicom('*.*', ...
'C:\Users\payal\Desktop\Research_PhD\Brain_Intracranial_Seg\IPH\cSDH 24\SDH 24\1');
rsinfo=dicominfo('C:\Users\payal\Desktop\Research_PhD\Brain_Intracranial_Seg\IPH\cSDH 24\SDH 24\DICOMDIR');
view3dgui(img)

%% Intensity Rescaling
[row ,column, n]=size(img);
img=im2double(img);
I3=0*img;
for i=1:n
I=img(:,:,i);
ma=max(I(:));
I=I/(ma);
I2=imadjustn(I);
I3(:,:,i)=I2;
end
%view3dgui(I3)

%% Skull Stripping/Removal

% Convert Image to grayscale 0-1
[row ,column, n]=size(I3);
get_image=0*I3;
binaryImage=0*I3;
binaryImage3=0*I3;
binaryImage4=0*I3;
binaryImage5=0*I3;
skull_strippedImage=0*I3;
for i=1:n

% Threshold to create a binary image
binaryImage1 = I3(:,:,i) > .5;
binaryImage(:,:,i)=binaryImage1;

% Get rid of small specks of noise
binaryImage2 = bwareaopen(binaryImage(:,:,i), 20);
binaryImage3(:,:,i)=binaryImage2;

%Fill the image
binaryImage4(:,:,i) = imfill(binaryImage3(:,:,i), 'holes');

% Erode away 15 layers of pixels.
se = strel('disk', 30, 0);
binaryImage5(:,:,i) = imerode(binaryImage4(:,:,i), se);

% Mask the gray image
finalImage1 = I3(:,:,i); % Initialize.
binaryImage6=binaryImage5(:,:,i);
finalImage1(~binaryImage6) = 0;
skull_strippedImage(:,:,i)=finalImage1;
end
% view3dgui(I3);
% view3dgui(binaryImage5);
% view3dgui(skull_strippedImage);

%% Contrast-limited adaptive histogram equalization...
[row,column,slice]=size(skull_strippedImage);
pout_adapthisteq=zeros(row,column,slice);
for k=1:slice
        	pout_adapthisteq(:,:,k) = adapthisteq(skull_strippedImage(:,:,k));
 end
%view3dgui(pout_adapthisteq)

%% Contrast-limited adaptive histogram equalization...
[row,column,slice]=size(pout_adapthisteq);
pout_adapthisteq1=zeros(row,column,slice);
for k=1:slice
        	pout_adapthisteq1(:,:,k) = adapthisteq(pout_adapthisteq(:,:,k));
end
view3dgui(pout_adapthisteq1)

pause(30.0);
%% Selecting Start and End Slice with bleed for segmentation 
answer = inputdlg('Enter Starting and Ending Slice Number for Segmentation:',...
             'IPH_Slices', [5 80]);
user_val = str2num(answer{1});
st=user_val(1);
en=user_val(2);
pout_adapthisteq1=pout_adapthisteq1(:,:,st:en);

%% Dehazing
[row,column,slice]=size(pout_adapthisteq1);
dehazed_image=zeros(row,column,slice);
for k=1:slice
        	dehazed_image(:,:,k) = imreducehaze(pout_adapthisteq1(:,:,k));
 end
%view3dgui(dehazed_image)

%% Gaussian Smoothing
[row,column,slice]=size(dehazed_image);
smooth_image=zeros(row,column,slice);
for k=1:slice
        	smooth_image(:,:,k) = imgaussfilt(dehazed_image(:,:,k),1);
 end
%view3dgui(smooth_image)

%% Deep Denoising
% Get the dimentional size
[~,~,m]=size(smooth_image);

% Generate a new EMPTY output (3D)
deep_denoisedImage=double(0*smooth_image);

for i=1:m
	% Display the current slice
	fprintf('Slice No. : %d\n',i');  
	get_image=smooth_image(:,:,i);
	% 2D Denosing Function
	cd('C:\Users\payal\Desktop\Research_PhD\Brain_Intracranial_Seg\FFDNet-master\FFDNet-master');
	[output] = DeepDenosing(get_image);
	deep_denoisedImage(:,:,i)=(output);
	clc;
end
%view3dgui(deep_denoisedImage)

%% Active Contour Segmentation

% mask creation
imshow(deep_denoisedImage(:,:,1));
r = drawrectangle; % Create customizable rectangular ROI
mask = createMask(r);

imshow(deep_denoisedImage(:,:,m));
r1 = drawrectangle; % Create customizable rectangular ROI
mask1 = createMask(r1);

% segmentation
final_Image=double(0*deep_denoisedImage);
for i=1:m
p1=deep_denoisedImage(:,:,i);
if i<m
BW = activecontour(p1,mask,500,'Chan-Vese','ContractionBias',0.73); % Segment image
else 
BW = activecontour(p1,mask1,500,'Chan-Vese','ContractionBias',0.73); % Segment image 
end
final_Image(:,:,i)=BW;
end
view3dgui(final_Image)

%% Overlaying image of original contrast enhanced and segmented image 
BW4=bwperim(final_Image(:,:,:));
I4=pout_adapthisteq1(:,:,:);
C1=I4+BW4;
overlayed_image=C1(:,:,:);
view3dgui(overlayed_image);

%% Calculating Volume of Total Brain Tissue
[row,column,slice]=size(skull_strippedImage);
non_boon_image1 = imgaussfilt3(skull_strippedImage,4);
count=nnz(non_boon_image1);
x=imginfo.Voxel_Size(1);
y=imginfo.Voxel_Size(2);
z=imginfo.Voxel_Size(3);
count1=count*x*y*z;

%% Calculating Volume of Blood Clot
bleed=sum(final_Image(:));
bleed1=bleed*x*y*z;

%% msg box for displaying calculated volume 
msgbox(sprintf('The total volume of Brain tissue region in mm^3 is = %d\n The total volume of Hemorrhage region in mm^3 is = %d',count1,bleed1),'Volume');
toc

%% Manual Ground Truth Contouring  

% lidx=find(strcmpi({handles.contours.Label},'test')==1);
% test_mask=bitget(handles.mask_from_contours,lidx);
% view3dgui(test_mask)
% 
% 
%% Evaluation qualitative analysis 

% bigbinaryImage1 = double(test_mask(:,:,5));
% binaryImage2 = final_Image(:,:,5);
% andImage = bigbinaryImage1 & binaryImage2;
% orImage = bigbinaryImage1 | binaryImage2;
% diceCoeff = 1 * (sum(andImage) / sum(orImage));
% binaryImage3 = abs(bigbinaryImage1 - binaryImage2);
% subplot(1,3,1);
% imshow(bigbinaryImage1);
% title('Original Mask (Ground Truth)')
% subplot(1,3,2);
% imshow(binaryImage2);
% title('Constructed Mask')
% subplot(1,3,3);
% imshow(binaryImage3);
% title('Mask Difference')
