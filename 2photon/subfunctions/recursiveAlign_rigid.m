function [displace,m] = recursiveAlign_rigid(stack,idx,optimizer,metric,fixedRefObj)

%Function for recursive alignement of 2-photon images in a t-stack.
%aligned image(m).
%
%INPUTS:
%   stack: YXT image(only one channel). 
%   idx: image index for recursive branches. Input a vector from 1:T for initial function call. 
%   optimizer and metric:  are obtained using imregconfig('monomodal')
%   fixedRefObj: obtained using imref2d on the first frame of the channel.
%
%OUTPUTs:
%   displace: cell array of tform objects
%   m: average of the fixed and moving image of each recursion
%
%
% Based on the method published on Dr Dario Ringach's blog
% (https://scanbox.org/2016/06/30/non-rigid-image-alignment-in-twenty-lines-of-matlab/)
%
% Adapted and modified by Eric Martineau - Universite de Montreal

%%
if length(idx)==1 % just one frame, return same image
    A = squeeze(stack(:,:,idx));
    m = A;
    displace = {affine2d(diag([1,1,1]))};
elseif length(idx) == 2 % if two frames only, align them together
    A = squeeze(stack(:,:,idx(1)));
    B = squeeze(stack(:,:,idx(2)));    
    
    [tform] = imregtform(A,B,'rigid', optimizer, metric); %caculate transform

    RegIm = imwarp(A,tform,'Outputview',fixedRefObj); %apply to A
    m = (RegIm/2 + B/2); %output average image after registeration
    displace = {tform affine2d(diag([1,1,1]))}; %output transformation in cell array
else
    idx0 = idx(1:floor(end/2)); %split data in two
    idx1 = idx(floor(end/2)+1 : end);

    [D0,m0]= recursiveAlign_rigid(stack,idx0,optimizer,metric, fixedRefObj);
    [D1,m1] = recursiveAlign_rigid(stack,idx1,optimizer,metric, fixedRefObj);
    
    [tform] = imregtform(m0,m1,'rigid', optimizer, metric);
    RegIm = imwarp(m0,tform,'Outputview',fixedRefObj);
    m = (RegIm/2+m1/2);
    
    D0 = cellfun(@(x) (affine2d(x.T*tform.T)),D0,'UniformOutput',false);
    displace = [D0, D1];
end
   