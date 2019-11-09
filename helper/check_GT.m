GTDir = 'C:\Users\xzh87\OneDrive\Reference\Lab\code\HARG-SSVM_v0.9\ground_truth';

class = 'Duck';

listOfGTFile = dir(fullfile(GTDir, class, '*.mat'));

for i = 1:length(listOfGTFile)
   GTFile = fullfile(GTDir, class, listOfGTFile(i).name);
   load(GTFile, 'groundTruth');
   
   if rank(groundTruth.assign) ~= size(groundTruth.assign, 1)
       fprintf("Wrong in %s\n", GTFile); 
       groundTruth.assign = eye(10, 10);
       
       save(GTFile, 'groundTruth');
   end
end