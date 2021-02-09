%% Phase II starts at line 96
%% Phase III starts at line 198

function [am, bm, cm, dm, al, bl, cl, dl, l] = e190661_aliyev(txt_path, lbr_path)
Xb=-1;
Bb=-1;
    %% Open the file and get the first line
    fileID = fopen(txt_path,'r');

    tline = fgetl(fileID);
    %% assign d2XX values to -1 
    d2AB=-1;
    d2AC=-1;
    d2BC=-1;
    %% get the library
    lbr = table2cell(readtable(lbr_path, 'PreserveVariableNames',true));
    %% some constants and variables
    eps=8.854*10^-12 % permitivity of free space
    f=50; % frequency of the system,(Hz)
    
    %% exporting variables from the file
    while ~strcmp(tline,'-999')
        if (strcmp(tline,'Number of circuits'))
            k=str2double(fgetl(fileID));     % Number of circuits
        end

        if (strcmp(tline,'Number of bundle conductors per phase'))
            n=str2double(fgetl(fileID));     % Number of bundle conductors per phase
        end

        if (strcmp(tline,'Bundle distance (m)'))
            d=str2double(fgetl(fileID));     % Bundle distance (m)
        end

        if (strcmp(tline,'Length of the line (km)'))
            l=str2double(fgetl(fileID));     % Length of the line (km)
        end

        if (strcmp(tline,'ACSR conductor name'))
            name=fgetl(fileID);              % ACSR conductor name
            
            while (~any(strcmp(lbr(:,1),name)))  % give warning in case of wrong conductor name
                warning('WARNING!!! Conductor name is not listed in the library. Please try another conductor name.');
                name=input('Enter anoter name: ', 's');
            end
            %% get the parameters of the conductor
            names=lbr(:,1); % names of the conductors in double format
            IndexC = strfind(names,name); % find row of the conductor
            Index = find(not(cellfun('isempty',IndexC))); % row of the conductor in the library
            sAL=lbr{Index,2}*5.067e-10;   % Aliminium area (m^2) 
            stranding=lbr{Index,3};       % Stranding 
            layers=lbr{Index,4};          % Layers of Aliminium
            Do=lbr{Index,5}*0.0254;       % Outside diameter (m) 
            ro = Do/2;                    % Outside radius (m)
            
            Rac=lbr{Index,7}/1609.344;    % AC 50 Hz Resistance 20?C (ohm/m)
            GMR=lbr{Index,8}*0.3048;      % GMR (m)
        end

        %% get the coordinates of each phase and distances between centers of each bundles
        if (strcmp(tline,'C1 Phase A (centre)')) %% C1 Phase A (centre)
            x1=str2double(fgetl(fileID));  
            y1=str2double(fgetl(fileID));
        end

        if (strcmp(tline,'C1 Phase B (centre)')) %% C1 Phase B (centre)
            x2=str2double(fgetl(fileID));  
            y2=str2double(fgetl(fileID));
        end

        if (strcmp(tline,'C1 Phase C (centre)')) %% C1 Phase C (centre)
            x3=str2double(fgetl(fileID));
            y3=str2double(fgetl(fileID));
        end

        if (strcmp(tline,'C2 Phase A (centre)')) %% C2 Phase A (centre)
            x4=str2double(fgetl(fileID));
            y4=str2double(fgetl(fileID));
        end

        if (strcmp(tline,'C2 Phase B (centre)')) %% C2 Phase B (centre)
            x5=str2double(fgetl(fileID));
            y5=str2double(fgetl(fileID));
        end

        if (strcmp(tline,'C2 Phase C (centre)')) %% C2 Phase C (centre)
            x6=str2double(fgetl(fileID));
            y6=str2double(fgetl(fileID));
        end
        tline = fgetl(fileID); %% get the next line (corresponds to text lines)
    end
    fclose('all');
    
    
    %% Phase II    
    % distances between phases
    d1AB=sqrt((x2-x1)^2+(y2-y1)^2);
    d1AC=sqrt((x3-x1)^2+(y3-y1)^2);
    d1BC=sqrt((x2-x3)^2+(y2-y3)^2); 
    GMD = nthroot(d1AB*d1AC*d1BC,3);    
    if k==2  
        d2AB=sqrt((x5-x4)^2+(y5-y4)^2);
        d2AC=sqrt((x6-x4)^2+(y6-y4)^2);
        d2BC=sqrt((x6-x5)^2+(y6-y5)^2);
        dA1A2=sqrt((x6-x1)^2+(y6-y1)^2);
        dA1B2=sqrt((x5-x1)^2+(y5-y1)^2);
        dA1C2=sqrt((x4-x1)^2+(y4-y1)^2);    
        dB1A2=sqrt((x6-x2)^2+(y6-y2)^2);
        dB1B2=sqrt((x5-x2)^2+(y5-y2)^2);
        dB1C2=sqrt((x4-x2)^2+(y4-y2)^2);
        dC1A2=sqrt((x6-x3)^2+(y6-y3)^2);
        dC1B2=sqrt((x5-x3)^2+(y5-y3)^2);
        dC1C2=sqrt((x4-x3)^2+(y4-y3)^2); 
        DAC=nthroot(dA1C2*dC1A2*d1AC*d2AC,4);
        DAB=nthroot(d1AB*d2BC*dB1A2*dA1B2,4)
        DBC=nthroot(dB1C2*dC1B2*d1BC*d2AB,4)
        GMD=nthroot(DAB*DAC*DBC,3);
    end      
 
 %
    if n==1
        GMRL=GMR;  
        GMRC=ro;
    elseif n==2
        GMRL=sqrt(GMR*d);
        GMRC=sqrt(ro*d);
    elseif n==3
        GMRL=nthroot(GMR*d^2,3);
        GMRC=nthroot(ro*d^2,3);
    elseif n==4
        GMRL=nthroot(GMR*d^3*sqrt(2),4);
        GMRC=nthroot(ro*d^3*sqrt(2),4);
    elseif n==5
        GMRL=nthroot(GMR*d^4*(sin(72*pi/180)/cos(54*pi/180))^2,5);
        GMRC=nthroot(ro*d^4*(sin(72*pi/180)/cos(54*pi/180))^2,5);
    elseif n==6
        GMRL = nthroot(GMR*2*d^5*3,6);
        GMRC= nthroot(ro*2*d^5*sqrt(3)^2,6);
    elseif n==7
        GMRL=nthroot(GMR*d^2*((2*sin((450/7)pi/180)*d)^2)((d*sin((540/7)*pi/180)/cos((450/7)*pi/180))^2),7);
        GMRC=nthroot(ro*d^2*((2*sin((450/7)pi/180)*d)^2)((d*sin((540/7)*pi/180)/cos((450/7)*pi/180))^2),7);
    elseif n==8
        GMRL=nthroot(GMR*d^2*((2*dbundle*sin((135/2)*pi/180))^2)*((d*(sqrt(2)+1))^2)*(d/cos((135/2)*pi/180)),8);
        GMRC=nthroot(ro*d^2*((2*d*sin((135/2)*pi/180))^2)*((d*(sqrt(2)+1))^2)*(d/cos((135/2)*pi/180)),8);
    end
    
    if k==2
        % GMR of each phase for inductor
        GMRa=sqrt(dA1A2*GMRL);
        GMRb=sqrt(dB1B2*GMRL);
        GMRc=sqrt(dC1C2*GMRL);
        %
        GMRCa=sqrt(dA1A2*GMRC);
        GMRCb=sqrt(dB1B2*GMRC);
        GMRCc=sqrt(dC1C2*GMRC);
        GMRL=nthroot(GMRa*GMRb*GMRc,3); 
        GMRC=nthroot(GMRCa*GMRCb*GMRCc,3); 
    end
    
%% Earth effect capacitance    
if k == 1
    H1 = sqrt((x2 - x1)^2 + (y1 + y2)^2)*sqrt((x2 - x3)^2 + (y2 + y3)^2)*sqrt((x1 - x3)^2 + (y1 + y3)^2); 
    Cearth = (H1/(2*y1*2*y2*2*y3))^(1/3);
elseif k == 2
    H2= (sqrt((x2-x1)^2 + (y1+y1)^2)*...
            sqrt((x4-x1)^2 + (y5+y1)^2)*...
            sqrt((x2-x6)^2 + (y1+y6)^2)*...
            sqrt((x6-x4)^2 + (y6+y5)^2)*...
            sqrt((x3-x2)^2 + (y3+y2)^2)*...
            sqrt((x4-x2)^2 + (y4+y2)^2)*...
            sqrt((x3-x4)^2 + (y3+y5)^2)*...
            sqrt((x4-x4)^2 + (y5+y4)^2)*...
            sqrt((x3-x1)^2 + (y3+y1)^2)*...
            sqrt((x4-x1)^2 + (y4+y1)^2)*...
            sqrt((x3-x6)^2 + (y3+y6)^2)*...
            sqrt((x4-x6)^2 + (y4+y6)^2))^(1/12);
    H3 = (64*y1*y3*y2*y4*y6*y5*...
               ((x6-x1)^2 + (y1+y6)^2)*...
               ((x3-x4)^2 + (y3+y4)^2)*...
               ((x2-x4)^2 + (y1+y5)^2))^(1/12);
    Cearth = (H2/H3);
end
 
    L=2*10^-4*log(GMD/GMRL);
    C=2*pi*eps/(log(GMD/GMRC)-log(Cearth));
    if k==2 
        C=2*pi*eps/(log(GMD/GMRC));
    end
    Rf=lbr{Index,7}/((1.60934)*n*k);          % DC Resistance 20?C (ohm/km)
    Xf=2*pi*f*L;
    Bf=2*pi*10^3*f*C;    
    
    
    
%% Phase III

R = Rf*l/1000;
X = Xf*l/1000;
Z = R+j*X;
Y = j*Bf*l/1000;
y = Y/l;
z = Z/l;
Zo = sqrt(z/y);
g = sqrt(z*y);
l = l*1000;
% if (l<241)
    am = (Z*Y/2)*10^6+1;
    bm = Z;
    cm = Y*(1+Z*Y/4);
    dm = am;

% else
    al = cosh(g*l);
    bl = Zo*sinh(g*l);
    cl = sinh(g*l)/Zo;
    dl = al;
% end
end