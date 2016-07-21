function class=svm(C)
load Trainfeature;
load Truelabel;
figure; title('Train');
svmStruct = svmtrain(round(Trainfeature),Truelabel,'showplot',true);
class=svmclassify(svmStruct,C,'showplot',true);
hold on;
plot(C(:,1),C(:,2),'ro','MarkerSize',12);
title('Train');
hold off