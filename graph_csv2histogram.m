data = readmatrix('stat_azimuthal_vel_l1.csv');              % �t�@�C�������݂̃t�H���_�[�܂��� MATLAB �p�X��̃t�H���_�[�ɂȂ��ꍇ�́Afilename �ɐ�΃p�X���܂��͑��΃p�X�����w�肵�ĂˁB

vp=data(:,1);

f1 = figure;                                % initial position plot

figure(f1);
nbins = 20;
histogram(vp*1e2,nbins)
hold on
% XX = -1e-3:1e-5:1e-3;
% mu = 5;
% sigma = 1/3*1e-3;
% f = 1000/sqrt(2*pi*sigma^2)*exp(-XX.^2./(2*sigma^2));
% plot(XX*1e3,f,'LineWidth',2);
xlabel(' azimuthal velocity unification [%] '); ylabel(' Number (SUM = 100) ');
ax = gca;
ax.FontSize = 13;

