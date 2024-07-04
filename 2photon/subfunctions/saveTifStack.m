function saveTifStack(file,filename)

filename = strcat(filename, ".tif");

if isfile(filename) > 0
    delete(filename);
end

for i = 1:size(file,3)
    imwrite(file(:,:,i), filename ,'WriteMode', 'append');
end