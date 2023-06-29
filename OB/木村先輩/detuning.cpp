data = readmatrix('fileName.csv');              % ファイルが現在のフォルダーまたは MATLAB パス上のフォルダーにない場合は、filename に絶対パス名または相対パス名を指定します。

x=data(1,:);
y=data(2,:);
z=data(3,:);

f1 = figure;                                %traj plot
f2 = figure;                                %traj plot in xy


%traj plot

figure(f1);
plot3(x*1e6,y*1e6,z*1e2);
hold on
xlabel(' x[um] '); ylabel(' y[um] '); zlabel(' z[cm] ');
plot3(x(1)*1e6,y(1)*1e6,0, '.k')
ax = gca;
ax.FontSize = 13;


%traj plot in xy

figure(f2);
plot(x*1e6,y*1e6);
hold on
xlabel(' x[um] '); ylabel(' y[um] ');
plot(x(1)*1e6,y(1)*1e6, 'ko')
plot(x(10000)*1e6,y(10000)*1e6, 'ro')
ax = gca;
ax.FontSize = 13;
% text( 20,31, "detune:"+ int16(detune/(2*pi*1000000))+"MHz")
% text( 9.5, -9, "confinement radius :"+int16(r_conf*1e6)+"um",'Fontsize',12)
xticks([-800 -400 0 400 800])
yticks([-1000 -500 0 500 1000])
% t = linspace(0,2*pi,100);
% plot( w0*1e6 *sin(t),w0*1e6 *cos(t),'-r');
