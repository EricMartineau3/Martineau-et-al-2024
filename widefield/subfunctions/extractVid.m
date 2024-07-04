function vid = extractVid(Trial, TrialType)

 for i = 1:length(Trial)
    vid(:,:,:,i) = Trial(i).vid;
 end
 
 vid = vid(:,:,:,TrialType);
 vid = mean(vid,4);
 vid = mat2gray(vid); 