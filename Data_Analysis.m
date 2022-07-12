
% written by Fabrice Gielen, 01/2022

% 'A parallel droplet generation and absorbance detection platform for the
% simultaneous screening of enzymatic reactions in nanoliter droplets'

% This software generates figures shown in Figure 3 of the paper. It identifies absorbance values for all droplets, match substrate concentration and extract Michaelis-Menten parameters 

%%requires Curve Fitting Toolbox, getsubs.m 

%% Manual input required:

    %1-set time for end of gradients manually in 'grad'
    
    %2-set duration of gradients in 'grad_time' variable  - based on Figure 1 fits
    
    %3- set signal for 1 mM PNP calibration in 'sig_1mM' variable (signal difference with oil)

    %4 change number of droplets/gradient based on first gradient

    %5 change initial enzyme and substrate concentration

close all
clear all

Path            =     'D:\Fabrice\Collaborations\Stefanie Neun\Matlab analysis\b-xylose_example-analysis\Example data analysis\'; % change folder accordingly
File            =     '1uM_243_400mM_b-Xyl_1.txt'; % raw data file

A=load([Path File]);
time=A(:,1);

% display raw data
figure(1); set(figure(1),'Color','w');
font_size=18;

k=12; % select line to analyse 
            
sig=A(:,k+1); %load next channel
oil=mode(sig); %oil level should be most represented value in signal
plot(time,sig); 
hold on;

xlabel('Time (s)');
ylabel('Signal (A.U.)');
set(gca,'FontSize',font_size,'FontWeight','bold','XColor', [0,0,0], 'YColor', [0,0,0])
set(gcf,'color','w');

%%%%%%%%%%%%%%% Change time for end of gradients below %%%%%%%%%%%%%%%%%%

%gradient times 
grad_1=[114;304.4;528.2;706.2;933.2;1112.3]; %channel 1
grad_2=[134.5;283.6;549;685.8;953.8;1091.9]; %channel 2
grad_3=[158;258.5;574;660.8;978.9;1067]; %channel 3
grad_4=[183;234.1;598.1;636.5;1003;1042.7]; %channel 4
grad_5=[95;320.1;508.5;721.7;913.7;1127.4]; %channel 5
grad_6=[120;296.5;533.2;698.2;938.5;1103.8]; %channel 6
grad_7=[146;267;562.9;668.7;968.2;1074.3]; %channel 7
grad_8=[169.5;243.1;586.5;644.8;991.8;1050.4]; %channel 8
grad_9=[64.6;352.2;481.6;753.9;886.6;1159.7]; %channel 9
grad_10=[84.3;333.2;500.9;735.2;906;1140.9]; %channel 10
grad_11=[109.6;307.4;526.7;709;932.2;1115]; %channel 11
grad_12=[126.8;290.5;543.4;692.3;948.6;1098]; %channel 12

% Change apparent duration of gradients per channel and 1 mM PNP
% calibration signal
switch k
    
    case 1
grad=grad_1;channel=1;sig_1mM=242;grad_time=31;
    case 2
grad=grad_2;channel=2;sig_1mM=198.5;grad_time=30.5;
    case 3
grad=grad_3;channel=3;sig_1mM=237.5;grad_time=30.5;
    case 4
grad=grad_4;channel=4;sig_1mM=221;grad_time=30.7;
    case 5
grad=grad_5;channel=5;sig_1mM=179;grad_time=30.2;
    case 6
grad=grad_6;channel=6;sig_1mM=183;grad_time=29; 
    case 7
grad=grad_7;channel=7;sig_1mM=110;grad_time=29.5; 
    case 8
grad=grad_8;channel=8;sig_1mM=170;grad_time=29.8; 
    case 9
grad=grad_9;channel=9;sig_1mM=136;grad_time=29.5;   
    case 10
grad=grad_10;channel=10;sig_1mM=186;grad_time=29.8;    
    case 11
grad=grad_11;channel=11;sig_1mM=100;grad_time=30.8;    
    case 12
grad=grad_12;channel=12;sig_1mM=96;grad_time=29.8; 
end

N=60; % change number of droplets/gradient based on first gradient
subs_i=400; %initial substrate concentration in mM
enz_i=1; %initial enzyme concentration in uM

filt_length=5; %length of vector for STD, start with 5
filt_std=1; %stringency of filter, start with 1
%%%%%%%%%%%%%%%%%%%%%%%%%%

B=A(:,channel+1); % open correct line
oil_baseline=mode(B); % oil level 
sampling=200; % sampling rate =200Hz
t_res=1/sampling;
time_end=grad(end)+grad_time+10; % end time for analyzing trace after last gradient

time2=time(1:time_end*sampling);
B3=B(1:time_end*sampling);
nb_grad=length(grad);

%%%%filter signal based on local standard deviation, eliminates high
%%%%scatter

h=waitbar(0,'Filtering Signal..');

for i=1:nb_grad %loop through the gradients
    
   waitbar(i/nb_grad)
   
     if mod(i,2)==1 %load in right direction 
        
         for j=(grad(i)-grad_time)*sampling:grad(i)*sampling % process only gradients, not all signal !
             
             
       if (std(B(j:j+filt_length))>filt_std || std(B(j-filt_length:j)) >filt_std) 
            
                B3(uint32(j))=NaN;
       end
         end
            
     else
           
         for j=(grad(i)*sampling):(grad(i)+grad_time)*sampling % process only gradients, not all signal 

               

       if (std(B(j:j+filt_length))>filt_std || std(B(j-filt_length:j)) >filt_std) 
            
                B3(uint32(j))=NaN;
       end
         end
           
     end 
    
end

close(h)

drops=zeros(N, nb_grad);
time_drops=zeros(N, nb_grad);
reac=zeros(grad_time*sampling+1, nb_grad);


for i=1:nb_grad % loop through the gradients to find fits

    
    if mod(i,2)==1 %load in right direction 
        grad_pos=uint32(grad(i)*sampling-grad_time*sampling:grad(i)*sampling);
        B2=B3(grad_pos); %load next gradient
        time2=time(grad_pos); %load time of gradient
        
            for j=round(length(B2)/2):1:length(B2) % this loop eliminates oil wrongly identified as drop   
            
                if abs(B2(j)-oil_baseline)<10
                    B2(j)=NaN;
                end
     end 
        
    else
        grad_pos=uint32(grad(i)*sampling:(grad(i)*sampling+grad_time*sampling));
        B2=B3(grad_pos); %load next gradient
        time2=time(grad_pos); %load time of gradient
        
            for j=1:round(length(B2)/2) % this loop eliminates oil wrongly identified as drop   
                 if abs(B2(j)-oil_baseline)<10
                    B2(j)=NaN;
                 end
            end 
    end
    
     time3=time2; %used for fitting
     time2(isnan(B2))=[];
     B2(isnan(B2))=[]; 
 
    %%%%%%%%%%%%%%%%%%%%%%
    figure(2); % plot individual droplet signals 
    set(figure(2),'Color','w');
    
    B_drops=(B2-oil_baseline);
    plot(time2,B_drops,'.','MarkerEdgeColor',[0 0 1]);%hold on;
    
    ylim([min(B_drops) 0]);
    
    xlabel('Time (s)');
    ylabel('Droplet signal (A.U.)');
    set(gca,'FontSize',font_size,'FontWeight','bold','XColor', [0,0,0], 'YColor', [0,0,0])
    set(gcf,'color','w');
    %%%%%%%%%%%%%%%%%%%%%%
  
    if mod(i,2)==1 %find exponential fit
        
         ft = fittype('a*(exp(b*(x+c)))+d');
            options = fitoptions(ft);
            options.StartPoint=[8,-0.17,-grad(i),B2(1)];
            options.Lower = [-10 -10 -10000 0];
            options.upper = [20 0 0 1e5];

            p=fit(time2,B2,ft,options);
            hold on;
            
            legend off;
            b=p.a*(exp(p.b.*(time3+p.c)))+p.d;
    
       
     
    else
       
         ft = fittype('a*(exp(b*(x+c)))+d');
            options = fitoptions(ft);
            options.StartPoint=[8,0.17,-grad(i),B2(end)];
            options.Lower = [-10 0 -10000 0];
            options.upper = [20 10 0 1e5];

            p=fit(time2,B2,ft,options);
            hold on;
           
            legend off;
            b=p.a*(exp(p.b.*(time3+p.c)))+p.d;
    end
    
    figure(1);hold on;
    plot(time3,b); % overlay the fits onto the raw data
    reac(:,i)=max(b)-b; 
      
    %%now convert to uM PNP 
    reac(:,i)=reac(:,i)*1000/sig_1mM; 

    %%%% now segment reactions into N separate equally spaced droplets
    for j=1:N 
        bins=size(reac,1)/N;
        drops(j,i)=mean(reac((j-1)*bins+1:j*bins,i));
        time_drops(j,i)=mean(time3((j-1)*bins+1:j*bins));
    end

    %%%%%%%%%%%%%%%%%%%%%%
    figure(3) % plot all droplet absorbance minus oil baseline values
    plot(time3,reac(:,i));
    hold on
    plot(time_drops,drops,'ro');
   
    xlabel('Time (s)');
    ylabel('Signal difference  (A.U.)');
    set(gca,'FontSize',font_size,'FontWeight','bold','XColor', [0,0,0], 'YColor', [0,0,0])
    set(gcf,'color','w');
    %%%%%%%%%%%%%%%%%%%%%%

end

    %flip every second column of drops /time_drops to extract slopes

    for i=1:nb_grad
        if mod(i,2)==0
            drops(:,i)=flipud(drops(:,i));
            time_drops(:,i)=flipud(time_drops(:,i));
        end    
    end

%%%%%%%%%%%%%%%%%%%%%%
figure(4) % Signal difference time courses
slope=zeros(N,1);
intercept=zeros(N,1);

for i=1:N
    plot(time_drops(i,:), drops(i,:));
    hold on
    
    p = polyfit(time_drops(i,:), drops(i,:), 1);
    slope(i,1) = p(1);
    intercept(i,1) = p(2);
    
end

xlabel('Time (s)');
ylabel('Signal difference time courses  (A.U.)');
set(gca,'FontSize',font_size,'FontWeight','bold','XColor', [0,0,0], 'YColor', [0,0,0])
set(gcf,'color','w');
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%Convert time_drops to substrate concentration
[subs,enz]=getsubs(time_drops(:,1)-time_drops(1,1),subs_i);
drops2=drops(:,1)-drops(1,1);

% fit slopes with Michaelis-Menten equation

fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0],...
               'Upper',[Inf,20],...
               'StartPoint',[1 5]);
ft = fittype('a*x/(x+b)','options',fo);

[curve2,gof2] = fit(subs,slope,ft);

%%%%%%%%%%%%%%%%%%%%%% 
figure(5) % plot slopes versus substrate concentration
plot(subs,slope);
hold on
plot(curve2,subs,slope);

%extract kcat and Km
kcatKm= coeffvalues(curve2);
kcat=kcatKm(1);
Km=kcatKm(2);

%correct kcat fro enzyme concentration
kcat=kcat/enz_i; % change according to enzyme concentration
text(0.6,0.5,strcat('k_{cat}= ',num2str(round(kcat,2)),' s^{-1}'),'Units','normalized','FontSize',16, 'FontWeight', 'bold'); % kcat
text(0.6,0.35,strcat('K_{M}= ',num2str(round(Km,2)),' mM'),'Units','normalized','FontSize',16, 'FontWeight', 'bold'); % KM
legend('off'); 

xlabel('Substrate concentration (mM)');
ylabel('Slope (uM/s)');
set(gca,'FontSize',font_size,'FontWeight','bold','XColor', [0,0,0], 'YColor', [0,0,0])
set(gcf,'color','w');
%%%%%%%%%%%%%%%%%%%%%%

