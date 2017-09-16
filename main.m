%  This Matlab code demonstrates an edge-based active contour model as an application of 
%  the Distance Regularized Level Set Evolution (DRLSE) formulation in the following paper:
%
%  C. Li, C. Xu, C. Gui, M. D. Fox, "Distance Regularized Level Set Evolution and Its Application to Image Segmentation", 
%     IEEE Trans. Image Processing, vol. 19 (12), pp. 3243-3254, 2010.
%
% Author: Chunming Li, all rights reserved
% E-mail: lchunming@gmail.com   
%         li_chunming@hotmail.com 
% URL:  http://www.imagecomputing.org/~cmli/

clear all;
close all;
for i = 78:1:500
    imgs(i,:) = char([num2str(i, '%03d'), '.bmp']);
    %dicoms(271 - i, :) = char([num2str(i, '%03d'), '.dcm']);
end
%imgs = ['CT2\336.jpg'; 'CT2\335.jpg'];
c0=2;
initialLSF=c0*ones([512, 512]);
% generate the initial region R0 as a rectangle
% maybe the region can be initialize in other shape
% from 300
%initialLSF(300:346, 235:290)=-c0;
%772 -0.8 50
%initialLSF(253:273, 280:335)=-c0;
initialLSF(245:258, 280:295)=-c0;

finalLSF = initialLSF;
initC = contour(initialLSF, [0 0], 'visible', 'off');

warning('off', 'Images:initSize:adjustingMag');
tic;
for i = 78:1:length(imgs)

fprintf('%s\n', imgs(i, :));
Img = imread(['CT4\', imgs(i, :)]); % real miscroscope image of cells
Img=double(Img(:,:,1));  % from integer to double, for computing facility
%Idcm = dicomread(['C:\Users\wwjwh\Desktop\ThirdAttempt\00044775\', dicoms(i,:)]);
%Idcm = double(Idcm);
%parameter setting
timestep= 5;  % time step
mu=0.2/timestep;  % coefficient of the distance regularization term R(phi)
iter_inner=10;
iter_outer=5;
lambda = 5; % coefficient of the weighted length term L(phi)
alfa = -1.5 * 4;  % coefficient of the weighted area term A(phi)
epsilon = 1.5; % papramater that specifies the width of the DiracDelta function

sigma= 1.5;     % scale parameter in Gaussian kernel
%Img_smooth = medfilt2(Img, [3 3]);
G=fspecial('gaussian', 5, sigma);
Img_smooth=conv2(Img, G, 'same');  % smooth image by Gaussiin convolution
Img_smooth = imfill(Img_smooth);
%[Ix,Iy]=gradient(Img_smooth);

Ix = DGradient(Img_smooth, 1, 2);
Iy = DGradient(Img_smooth, 1, 1);

f=Ix.^2+Iy.^2;
g=1./(1+f);  % edge indicator function.

% initialize LSF as binary step function
%phi=initialLSF;
phi = finalLSF;
% figure(1);
% mesh(-phi);   % for a better view, the LSF is displayed upside down
% hold on;  contour(phi, [0,0], 'r','LineWidth',2);
% title('Initial level set function');
% view([-80 35]);

potential=2;  
if potential ==1
    potentialFunction = 'single-well';  % use single well potential p1(s)=0.5*(s-1)^2, which is good for region-based model 
elseif potential == 2
    potentialFunction = 'double-well';  % use double-well potential in Eq. (16), which is good for both edge and region based models
else
    potentialFunction = 'double-well';  % default choice of potential function
end


%% start level set evolution

for n=1:iter_outer
    if mod(n - 1, 2) == 0
        figure(1);
        imshow(Img_smooth); axis off; axis equal; colormap(gray); hold on;  
        C = contour(phi, [0,0], 'r');
        fprintf('%d\n', size(initC, 2));
        rate = double(abs(size(initC, 2) - size(C, 2)));
        initC = C;
        fprintf('%d\n', size(C, 2));
        title(num2str(i));
        pause(0.5);
    end
    phi = drlse2d(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction);
    ind1 = sub2ind(size(phi), find(phi < 0));
    ind2 = sub2ind(size(Img_smooth), find(Img_smooth < 50));
    ind12 = intersect(ind1, ind2);
    phi(ind12) = c0;
    fprintf('iter_outer is %d\n', n);
    if rate <= 1
        break;
    end
end

%%
% refine the zero level contour by further level set evolution with alfa=0
% alfa=0;
% iter_refine = 10;
% phi = drlse_edge(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction);

figure('visible', 'off');
imshow(Img); axis off; axis equal;
colormap(gray); hold on;  contour(phi, [0,0], 'g', 'LineWidth', 3);
pos = strfind(imgs(i, :), '.');
saveas(gcf, ['C:\Users\wwjwh\Desktop\DRLSE_v0\DRLSE_v0\result3\',imgs(i, 1:pos), 'bmp']);
finalLSF=phi;

end
toc;