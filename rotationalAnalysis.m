%  File and function name : rotationalAnalysis
%  This program tracks the bright field images of beads in three
%  dimensions. In x and y, 2D Gaussian fitting is used and in z axis,
%  intensity of the beads is used. 
%
%  Date of completion     : 1 January 2014
%  
%  Written by    :   Sinan Can    sinan.can@berkeley.edu
%    
%  Inputs        :
%
%           file= Analized file name
%           pixelsize= Pixel size of recorded video.
%
%  Outputs       :   
%           x,y,z= x,y,z positions of tracked beads. 
%             
%
%  To Run >>    - Video file with bead image
%               - Spotlist file
%                       spotlist: is the name of the spotlist file that contains
% a list of spots to be tracked as well as header information. This file is required
% to run the 2D Gaussian tracking. You should define rough estimate of x position (Sx),
% y position (Sy), start (Start) and end (End) frames for tracking as well
% as the FileName of video. 
%               - Defaults file
%                       defaults.mat has the necessary parameters for
% Gaussian tracking.
             

function [x y z intensity]=rotationalAnalysis(file,pixelsize)

name= strtok(file, '.');
filecorrect=strcat(name,'corrected.tif');

% gets the information of video (size, length)
info=imfinfo(file);
w=info(1).Width;
h=info(1).Height;


% Takes the magnitude of the image for black white correcting and saves it into a new file. 
% Mean intensity is taken as zero level.
for k=1:length(info);
    a=imread(file,'Index',k);
    avg=mean2(a);
          
    for i=1:h         
        for j=1:w 
            if a(i,j)< avg
                a(i,j)=a(i,j)+2*(avg-a(i,j));
            end
        end
    end
    
    if k==1
       mkdir('temp')
       imwrite(a,strcat('temp/',filecorrect));
    else
        imwrite(a,strcat('temp/',filecorrect),'WriteMode','append');
    end
end
movefile(strcat('temp/',filecorrect));
rmdir('temp');
%WHTTrackHighRes_PCfix finds the peak positions by using 2D Gaussian
%fitting. This functions returns x and y values of bead images and saves it
%into a text file.

WHTrackHighRes_PCfix; % This function is written by the members of Yildiz Lab

%Loads the output of WHTrackHighRes_PCfix function for analysis

data=load(strcat(name,'corrected_spot1_1_xy_only.txt'));

% WHTrackHighRes_PCfix returns NaN for some regions where it cannot fit the
% image. NaN's are cleared out for further analysis

dataNew= data(0== sum(isnan(data), 2), :);

x_pos=dataNew(:,1);
y_pos=dataNew(:,2);

% Rotates the x and y to find the on and off-axis parts of motion. Function
% line_rotation rotates the given data. First I fit it to a line and rotate the
% data based on the angle and center point. 
Mol=[x_pos y_pos];
    if mod(length(x_pos),2)==0
        CenterMol=[x_pos((length(Mol))/2)  y_pos((length(Mol))/2)];
    end

    if mod((length(x_pos)-1),2)==0
        CenterMol=[x_pos((length(Mol)-1)/2)  y_pos((length(Mol)-1)/2)];
    end

pmol=polyfit(Mol(:,1),Mol(:,2),1);
mmol=atan(pmol(1));
% 57.2957795 is forchange from radian to degrees
molrotated=rotation(Mol,CenterMol,-mmol*57.2957795);  
x=molrotated(:,1)-min(molrotated(:,1));
y=molrotated(:,2)-mean(molrotated(:,2));
% plots and saves the x-y traces of of tracked beads. 
h=figure;
plot(x,y);
axis equal;
saveas(h,name,'fig');


%% Z position Tracking
% Intensities corresponding to every x-y point is found by using original
% video file. intensity array has the intensity of each x-y point.
%pixelsize=160;
l=0;
n=0;
intensity=0;

for k=1:length(info)
         
    a=imread(file,'Index',k);
     
    if isfinite(data(k,:))==1
        n=k-l;
        intensity(n)=a(round(dataNew(n,2)/pixelsize),round(dataNew(n,1)/pixelsize)); 
    end
    if isfinite(data(k,:))==0
        k=k+1;
        l=l+1;
    end
    intensity=cast(intensity,'double');   
end

%% Z information is found by using the intensity data. 
% Intensities are calibrated by recording the intensities of beads at 
% different heights with respect to objective.Based on this calibration, the
% intensities are fitted to polynomial p given below. Since the intensity
% varies based on the bright field light source power, they are scaled to
% fit the calibration data intensities (-58 to 89). And z positions are
% evaluated by using fit polynomial and scaledintensity of images. 

p=[0.0002 -0.0277 2.9805 47.363];
scaledintensity=scaledata(intensity,-58,89);
z=polyval(p,scaledintensity);

%% Plot the 3D figure and save with the variables. 
h=figure;
plot3(smooth(x,5),smooth(y,5),smooth(z,5));
axis equal
saveas(h,strcat(name,'ThreeD'),'fig');

save variables.mat




