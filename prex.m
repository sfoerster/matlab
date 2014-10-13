g=7;
childrenPerParent=5;
p=0.3;

figure(23); h=gca; 
figure(24); h2=gca;
figure(25); h3=gca;

[dist,shis,ex] = geneticDist([0 1],p,childrenPerParent,g);
axes(h); cla; plot(h,shis(:,1),'b-'); hold on; grid on; xlabel('Generation'); ylabel('Probability of extinction'); title(['Children per parent: ' num2str(childrenPerParent) '; Each child has a ' num2str(p) ' probability of reproducing'])
axes(h2); cla; plot(h2,ex,'b-'); hold on; grid on; xlabel('Generation'); ylabel('Expected population'); title(['Children per parent: ' num2str(childrenPerParent) '; Each child has a ' num2str(p) ' probability of reproducing'])
axes(h3); cla; plot(h3,diff(shis(:,1)),'b-'); hold on; grid on; xlabel('Generation'); ylabel('Probability of extinction in the next generation'); title(['Children per parent: ' num2str(childrenPerParent) '; Each child has a ' num2str(p) ' probability of reproducing'])

[dist,shis,ex] = geneticDist([0 0 1],p,childrenPerParent,g);
plot(h,shis(:,1),'r-');
plot(h2,ex,'r-');
plot(h3,diff(shis(:,1)),'r-');

[dist,shis,ex] = geneticDist([0 0 0 1],p,childrenPerParent,g);
plot(h,shis(:,1),'k--');
plot(h2,ex,'k--');
plot(h3,diff(shis(:,1)),'k--');

[dist,shis,ex] = geneticDist([0 0 0 0 1],p,childrenPerParent,g);
plot(h,shis(:,1),'m-.');
plot(h2,ex,'m-.');
plot(h3,diff(shis(:,1)),'m-.');

[dist,shis,ex] = geneticDist([0 0 0 0 0 1],p,childrenPerParent,g);
plot(h,shis(:,1),'b-.');
plot(h2,ex,'b-.');
plot(h3,diff(shis(:,1)),'b-.');

legend(h,{'1 initial parent','2 initial parents','3 initial parents','4 initial parents','5 initial parents'},'Location','Northwest')
legend(h2,{'1 initial parent','2 initial parents','3 initial parents','4 initial parents','5 initial parents'},'Location','Northwest')
legend(h3,{'1 initial parent','2 initial parents','3 initial parents','4 initial parents','5 initial parents'})