function preprocess
p = 'C:\Users\xzh87\OneDrive\Reference\Lab\code\WILLOW-ObjectClass_dataset\feature';
Class = ["Car", "Duck", "Face", "MotorBike", "Winebottle"];
for cls = Class
    listOfFile = dir(fullfile(p, cls, "*.mat"));
    fprintf("processing %s...\n", cls);
    for ii = 1:length(listOfFile)
        pp = fullfile(listOfFile(ii).folder, listOfFile(ii).name);
        viewInfo = load(pp);
        view = viewInfo.cdata.view;
        save(pp, 'view');
    end
end
end