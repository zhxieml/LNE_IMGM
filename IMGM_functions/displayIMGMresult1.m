clear all;
fileList = dir('imgm*.mat');
dirLength = length(fileList);
acc = zeros(3,10);
con = zeros(3,10);
scr = zeros(3,10);
tim = zeros(3,10);

for i=1:length(fileList)
    load(fileList(i).name);
    acc=acc+accResult;
    con=con+conResult;
    scr=scr+scrResult;
    tim=tim+timeResult;
    clear accResult; clear conResult; clear scrResult; clear timeResult;
end

acc = acc/dirLength;
con = con/dirLength;
scr = scr/dirLength;
tim = tim/dirLength;

figure, title('Accuracy'); hold on;
acc1 = plot(acc(1,:),'Color','r');
acc2 = plot(acc(2,:),'Color','g');
acc3 = plot(acc(3,:),'Color','b');
legend([acc1, acc2, acc3],'Raw CAO','Incre CAO','IMGM');
hold off;

figure, title('Time'); hold on;
tim1 = plot(tim(1,:),'Color','r');
tim2 = plot(tim(2,:),'Color','g');
tim3 = plot(tim(3,:),'Color','b');
legend([tim1, tim2, tim3],'Raw CAO','Incre CAO','IMGM');
hold off;