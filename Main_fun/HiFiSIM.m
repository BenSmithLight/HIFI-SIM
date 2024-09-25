function varargout = HiFiSIM(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HiFiSIM_OpeningFcn, ...
                   'gui_OutputFcn',  @HiFiSIM_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before HiFiSIM is made visible.
function HiFiSIM_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
%% 
handles.x=[];
handles.y=[];
handles.fig=[];
h = handles.figure1;

movegui(h, 'northwest');
set(0,'Units','normalized'); 

set(handles.checkbox1,'Value',0);
set(handles.checkbox2,'Value',0);
set(handles.checkbox6,'Value',0);
set(handles.edit1,'String',' ');
set(handles.edit2,'String',' ');
set(handles.edit3,'String',' ');
set(handles.edit4,'String',' ');
set(handles.edit5,'String',' ');
set(handles.edit6,'String',' ');
set(handles.edit8,'String',' ');

set(handles.text22,'Visible','off');
set(handles.edit14,'Visible','off');

set(handles.text12,'Visible','off');
set(handles.edit11,'Visible','off');
set(handles.text14,'Visible','off');

set(handles.text13,'Visible','off');
set(handles.edit10,'Visible','off');
set(handles.text15,'Visible','off');

set(handles.text10,'Visible','off');
set(handles.edit9,'Visible','off');
set(handles.text11,'Visible','off');

clc;            

FlagDataType=0;
param.FlagDataType=FlagDataType;
param.FlagApoFWHM=0;
handles.param=param;
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = HiFiSIM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.checkbox1,'Value')
    FlagDataType=1;
    set(handles.checkbox2,'Value',0);
    set(handles.checkbox6,'Value',0);
    
    set(handles.edit1,'String',3);
    set(handles.edit2,'String',3);
    set(handles.edit3,'String',' ');
    set(handles.edit4,'String',525);
    set(handles.edit5,'String',78.6);
    set(handles.edit6,'String',1.42);
    set(handles.edit8,'String',0.90);
    set(handles.edit9,'String',1.2);
    set(handles.edit10,'String',1.0);
    set(handles.edit11,'String',0.5);
else
    FlagDataType=0;
    set(handles.edit1,'String',' ');
    set(handles.edit2,'String',' ');
    set(handles.edit3,'String',' ');
    set(handles.edit4,'String',' ');
    set(handles.edit5,'String',' ');
    set(handles.edit6,'String',' ');
    set(handles.edit8,'String',' ');
    set(handles.edit9,'String',' ');
    set(handles.edit10,'String',' ');
    set(handles.edit11,'String',' ');
end

set(handles.edit9,'Enable','off');
set(handles.edit10,'Enable','off');
set(handles.edit11,'Enable','off');

param.phaOff=0;   
param.fac=ones(1,2);    
param.nrBands=2;
param.FlagDataType=FlagDataType;
param.FlagApoFWHM=0;
handles.param=param;
guidata(hObject,handles);

% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.checkbox2,'Value')
    set(handles.checkbox1,'Value',0);
    set(handles.checkbox6,'Value',0);
    
    set(handles.edit1,'String',3);
    set(handles.edit2,'String',5);
    set(handles.edit3,'String',' ');
    set(handles.edit4,'String',525);
    set(handles.edit5,'String',79.5);
    set(handles.edit6,'String',1.40);
    set(handles.edit8,'String',0.90);
    set(handles.edit9,'String',1.2);
    set(handles.edit10,'String',1.0);
    set(handles.edit11,'String',0.5);
    FlagDataType=1;
else
    FlagDataType=0;
    set(handles.edit1,'String',' ');
    set(handles.edit2,'String',' ');
    set(handles.edit3,'String',' ');
    set(handles.edit4,'String',' ');
    set(handles.edit5,'String',' ');
    set(handles.edit6,'String',' ');
    set(handles.edit8,'String',' ');
    set(handles.edit9,'String',' ');
    set(handles.edit10,'String',' ');
    set(handles.edit11,'String',' ');
end
set(handles.edit9,'Enable','off');
set(handles.edit10,'Enable','off');
set(handles.edit11,'Enable','off');

param.phaOff=0;   
param.fac=ones(1,3);    
param.nrBands=3;
param.FlagDataType=FlagDataType;
param.FlagApoFWHM=0;
handles.param=param;
guidata(hObject,handles);

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
param=handles.param;
if param.FlagDataType==0
    warndlg('Please set the raw  data type!','!!warn!!','modal'); 
    return;
end

%% Load  raw SIM data
FlagLoad=1;
param=handles.param;
param.nrDirs=str2double(get(handles.edit1,'String')); 
param.nrPhases=str2double(get(handles.edit2,'String'));
N=param.nrDirs*param.nrPhases;
[filename, pathname] = uigetfile({'*.tiff;*.tif'},'Select the raw data sequence','MultiSelect','off');  
if isequal(filename,0)
    disp('Cancel data reading！');
    return;
else
    info = imfinfo(fullfile(pathname, filename));
    if numel(info)==N
        NPixel=max(info(1).Width,info(1).Height);
        Iraw0=zeros(info(1).Height,info(1).Width,N);
        Iraw=zeros(NPixel,NPixel,N); 
        for j=1:N
            Iraw0(:,:,j)=double(imread(fullfile(pathname, filename),j));
            if info(1).Width == info(1).Height || info(1).Width > info(1).Height
                Iraw(1:info(1).Height,:,j)=Iraw0(1:info(1).Height,:,j);
            else
                Iraw(:,1:info(1).Width,j)=Iraw0(:,1:info(1).Width,j);
            end 
%             MIJ.createImage(Iraw(:,:,j));
        end
    else
        warndlg('Raw data frames is not correct!','!!warn!!','modal');
        return;
    end
end

set(handles.edit3,'String',pathname);
Format=info.Format;
param.Format=Format;
param.filename=filename;
param.imgSize=NPixel;
param.Size1 = info(1).Height; 
param.Size2 = info(1).Width;   
param.Iraw = Iraw;  
FlagParameter=0;                    
param.FlagParameter=FlagParameter;
param.FlagLoad=FlagLoad;
handles.param=param;
guidata(hObject,handles);

function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
param=handles.param;
if isfield (param,'FlagLoad')==0
    warndlg('Please load the raw  data!','!!warn!!','modal'); 
    return;
end

FlagParameter=param.FlagParameter;
if FlagParameter==0
    Progressbar = waitbar(1/4, 'Preprocessing...','Name','SI Reconstruction');
   
    %% Generate approximate OTF/PSF
    Iraw=param.Iraw;
    NPixel=size(Iraw,1);  
    param.micronsPerPixel=str2double(get(handles.edit5,'String'))*10^(-3);    % 单位：微米
    param.cyclesPerMicron=1/(NPixel*param.micronsPerPixel);
    param.NA=str2double(get(handles.edit6,'String'));
    param.lambda=str2double(get(handles.edit4,'String'));
    
    param.attStrength=0;
    param.OtfProvider=SimOtfProvider(param,param.NA,param.lambda,1);
    
    psf=abs(otf2psf((param.OtfProvider.otf)));
    handles.param.OtfProvider=param.OtfProvider;
    handles.param=param;
    
    %% Preprocessing
    N=param.nrDirs*param.nrPhases;
    Temp=importImages(Iraw);  
    IIraw=deconvlucy(Temp,psf,5);
    for I=1:N
        IIrawFFT(:,:,I)=FFT2D(IIraw(:,:,I),false);
    end;
    
    WF=zeros(NPixel,NPixel);
    Tdeconv=zeros(NPixel,NPixel);
    WFdeconv=zeros(NPixel,NPixel,param.nrDirs);
    WFdeconvFFT=zeros(NPixel,NPixel,param.nrDirs);
    for i=1:param.nrDirs
        for j=1:param.nrPhases
            Tdeconv=Tdeconv+IIraw(:,:,(i-1)*param.nrDirs+j);
            WF(:,:)=WF(:,:)+Iraw(:,:,(i-1)*param.nrDirs+j);
        end
        WFdeconv(:,:,i)=Tdeconv/param.nrPhases;
        WFdeconvFFT(:,:,i)=FFT2D(WFdeconv(:,:,i),false);
    end
    WF=WF/N;
    WF=importImages(WF); 
    
    %% Wide-field
    Size1=2*param.Size1;
    Size2=2*param.Size2;
    WF2=zeros(Size1,Size2);
    Temp=zeros(2*size(WF,1),2*size(WF,2));
    fftWF=fftshift(fft2(WF));
    Temp(size(WF,1)/2+1:size(WF,1)/2+size(WF,1),size(WF,2)/2+1:size(WF,2)/2+size(WF,2))=fftWF;
    Temp=abs(ifft2(Temp));
    Temp=255*Temp/max(max(Temp));
    WF2(1:Size1,1:Size2)=Temp(1:Size1,1:Size2);
    WF2=importImages(WF2);
    handles.WF2=WF2;
    
    %% Parameter Estimation
    waitbar(2/4,Progressbar, 'Parameter Estimation...');
    
    % Kc region MASK
    cnt=[NPixel/2+1,NPixel/2+1];
    param.cutoff=1000/(0.5*param.lambda/param.NA);
    [x,y]=meshgrid(1:NPixel,1:NPixel);
    rad=sqrt((y-cnt(1)).^2+(x-cnt(2)).^2);
    Mask=double(rad<=1.0*(param.cutoff/param.cyclesPerMicron+1));
    NotchFilter0=getotfAtt(NPixel,param.OtfProvider.cyclesPerMicron,0.5*param.cutoff,0,0);
    NotchFilter=NotchFilter0.*Mask;
    Mask2=double(rad<=1.10*(param.cutoff/param.cyclesPerMicron+1));
    NotchFilter2=NotchFilter0.*Mask2;
    
    CrossCorrelation=zeros(size(Mask2,1),size(Mask2,2),param.nrDirs);
    k0=zeros(1,param.nrDirs);
    
    for I=1:param.nrDirs
        lb=2;
        if param.nrBands==2
            hb=2;
            fb=lb;
        elseif param.nrBands==3
            hb=4;
            fb=hb;
        end
         
        param.phaOff=0;    
        param.fac=ones(1,param.nrBands); 
        separateII=separateBands(IIrawFFT(:,:,(I-1)*param.nrPhases+1:I*param.nrPhases),param.phaOff,param.nrBands,param.fac);
        
        SeparateII{1,I}=separateII;
        
        if param.nrBands==2
            c0=separateII(:,:,1);
            c2=separateII(:,:,lb);
           
            c0=(c0./(max(max(abs(c0)))));
            c2=(c2./(max(max(abs(c2)))));
          
            c0=c0.*NotchFilter;
            c2=c2.*NotchFilter;
            
            c0=FFT2D(c0,false);
            c2=FFT2D(c2,false);
            c2=c2.*conj(c0); 
            c2=c2./max(max(c2));              

            vec=fftshift(FFT2D(c2,true));
        elseif param.nrBands==3
            c0=separateII(:,:,1);
            c3=separateII(:,:,hb);
            
            c0=c0./(max(max(abs(separateII(:,:,1))))).*NotchFilter;
            c3=c3./(max(max(abs(separateII(:,:,hb))))).*NotchFilter;
 
            c0=FFT2D(c0,false);
            c3=FFT2D(c3,false);
            c3=c3.*conj(c0);
            c3=c3./max(max(c3)); 

            vec=fftshift(FFT2D(c3,true));
        end
        CrossCorrelation(:,:,I)=vec;
        %% 
%         temp=vec.*NotchFilter; 
        temp=vec.*NotchFilter2;
        temp=log(1+abs(temp));
        temp=temp./max(max(temp));
%         MIJ.createImage(temp);
    
        [yPos,xPos]=find(temp==max(max(temp)));
        peak.xPos=xPos(1);
        peak.yPos=yPos(1);
        k0(I)=sqrt((peak.xPos-cnt(1))^2+(peak.yPos-cnt(2))^2);
    end
    
    Flag=0;
    if param.nrDirs>2           % For very few special cases
        if max(k0)-min(k0)>8
            Flag=1;
            Kobject=min(k0); 
%             Kobject=208; 
            Mask1=rad>=(Kobject+1);
            Mask2=rad<=(Kobject-1);
        end
    end
    
    for I=1:param.nrDirs
        vec=CrossCorrelation(:,:,I);
        if Flag==1
            vec(Mask1)=0;
            vec(Mask2)=0;
        end
        temp=vec.*NotchFilter2;
        temp=log(1+abs(temp));
        temp=temp./max(max(temp));
%         MIJ.createImage(temp);
    
        [yPos,xPos]=find(temp==max(max(temp)));
        peak.xPos=xPos(1);
        peak.yPos=yPos(1);

        cntrl=zeros(10,30);
        overlap=0.15;
        step=2.5;
        bn1=(param.nrBands-1)*2;
        kx=(peak.xPos-cnt(2));
        ky=(peak.yPos-cnt(1));
        
        separateII=SeparateII{1,I};
        [peak,cntrl]=fitPeak(separateII(:,:,1)/(max(max(abs(separateII(:,:,1))))),separateII(:,:,fb)/(max(max(abs(separateII(:,:,fb))))),1,bn1,param.OtfProvider,-kx,-ky,overlap,step,cntrl);
        
        if lb~=hb
            if param.nrBands==2
                peak.kx=peak.kx*2;
                peak.ky=peak.ky*2;
            end
            
            p1=getPeak(separateII(:,:,1),separateII(:,:,lb),0,1,param.OtfProvider,peak.kx/2,peak.ky/2,overlap);
            p2=getPeak(separateII(:,:,1)/(max(max(abs(separateII(:,:,1))))),separateII(:,:,hb)/(max(max(abs(separateII(:,:,1))))),0,2,param.OtfProvider,peak.kx,peak.ky,overlap);
            
            param.Dir(I).px=-peak.kx/2;    
            param.Dir(I).py=-peak.ky/2;
            param.Dir(I).phaOff=-phase(p1);
            Temp_m1=abs(p1);
            Temp_m2=abs(p2);
            
            Temp_m1(Temp_m1>1.0)=1;
            Temp_m2(Temp_m2>1.0)=1.0;
            param.Dir(I).modul(1)=Temp_m1;
            param.Dir(I).modul(2)=Temp_m2;
        end
        if lb==hb
            p1=getPeak(separateII(:,:,1),separateII(:,:,lb),1,lb,param.OtfProvider,peak.kx,peak.ky,overlap);
            param.Dir(I).px=-peak.kx;
            param.Dir(I).py=-peak.ky;
            param.Dir(I).phaOff=-phase(p1);
            Temp_m=abs(p1);
            Temp_m(Temp_m>1.0)=1.0;
            param.Dir(I).modul=Temp_m;
        end
        K0(I)=sqrt((param.Dir(I).px)^2+(param.Dir(I).py)^2);
        %% 
%             fittedPeak=vec;
%             for x=1:30;
%                 for y=1:10
%                     fittedPeak(y,x)=cntrl(y,x);
%                 end
%             end
            
%             figure,imshow(temp,[]);
%             colormap('hot');
%             hold on;
%             
%             if lb~=hb
%                 f=(param.nrBands-1)/2;
%             else
%                 f=1;
%             end
%             t=0:0.1:2.1*pi;
%             x=cnt(2)-peak.kx+10*sin(t);
%             y=cnt(1)-peak.ky+10*cos(t);
%             plot(x,y,'-w','LineWidth',2);
% %             title(['Orientation',num2str(I)]);
%             %         plotbrowser('on');
% %         end
%         %     plotbrowser('off');
    end
  
    %%
    SIMparam=zeros(3,5);
    if param.nrPhases==3
        for i=1:param.nrDirs
            SIMparam(i,1)=atan(param.Dir(i).py/param.Dir(i).px)*180/pi;
            SIMparam(i,2)=sqrt((param.Dir(i).px)^2+(param.Dir(i).py)^2);
            SIMparam(i,3)=param.Dir(i).phaOff;
            SIMparam(i,4)=param.Dir(i).modul;
            if param.Dir(i).modul<0.35
                param.Dir(i).modul=0.7;
            end
        end
    elseif param.nrPhases==5
        for i=1:param.nrDirs
            SIMparam(i,1)=atan(param.Dir(i).py/param.Dir(i).px)*180/pi;
            SIMparam(i,2)=sqrt((param.Dir(i).px)^2+(param.Dir(i).py)^2)*2;
            SIMparam(i,3)=param.Dir(i).phaOff;
            SIMparam(i,4)=param.Dir(i).modul(1);
            SIMparam(i,5)=param.Dir(i).modul(2);
            if param.Dir(i).modul(1)<0.35
                param.Dir(i).modul(1)=0.7;
            end
            if param.Dir(i).modul(2)<0.35
                param.Dir(i).modul(2)=0.7;
            end
        end
    end 
    set(handles.uitable1,'Data',SIMparam);
    FlagParameter=1;
    handles.K0=K0;
    handles.IIrawFFT=IIrawFFT;
    handles.param.Dir=param.Dir;
    param.WF2=WF2;
    handles.param=param;
    guidata(hObject,handles);

    waitbar(3/4,Progressbar, 'Reconstruction...');
end

%% Reconstruction
if param.FlagParameter==1
    Progressbar = waitbar(1/2, 'Reconstruction...','Name','SI Reconstruction');
end

param=handles.param;
IIrawFFT=handles.IIrawFFT;
Iraw=param.Iraw;
siz=size(Iraw(:,:,1));
w=siz(2);
h=siz(1);

param.attStrength=str2double(get(handles.edit8,'String'));
param.a=str2double(get(handles.edit10,'String'));           %  阻尼因子“damping” factor：β
param.attFWHM=1.0;
param.OtfProvider=SimOtfProvider(param,param.NA,param.lambda,param.a);
handles.param=param;
guidata(hObject,handles);
       
%% HiFi-SIM：Spectrum optimization
fftDirectlyCombined=zeros(h*2,w*2);
for I=1:param.nrDirs
    par=param.Dir(I);
    param.fac(2:param.nrBands)=param.Dir(I).modul(1:param.nrBands-1);   
    param.fac(2:param.nrBands)=param.Dir(I).modul(1:param.nrBands-1);
    separate=separateBands(IIrawFFT(:,:,(I-1)*param.nrPhases+1:I*param.nrPhases),par.phaOff,param.nrBands,param.fac);
    
    shifted=zeros(2*h,2*w,param.nrPhases);
    shifted(:,:,1)=placeFreq(separate(:,:,1));
    
    for b=2:param.nrBands
        pos=b*2-2;
        neg=b*2-1;
        shifted(:,:,pos)=placeFreq(separate(:,:,pos));
        shifted(:,:,neg)=placeFreq(separate(:,:,neg));
        
        shifted(:,:,pos)=NfourierShift(shifted(:,:,pos),-(b-1)*par.px,-(b-1)*par.py);
        shifted(:,:,neg)=NfourierShift(shifted(:,:,neg),(b-1)*par.px,(b-1)*par.py);
    end
        shifted(:,:,1)=applyOtf(shifted(:,:,1),param.OtfProvider,1,0,0,1,0);
        for b=2:param.nrBands
            pos=b*2-2;
            neg=b*2-1;
            shifted(:,:,pos)=applyOtf(shifted(:,:,pos),param.OtfProvider,b,-(b-1)*par.px,-(b-1)*par.py,1,0);
            shifted(:,:,neg)=applyOtf(shifted(:,:,neg),param.OtfProvider,b,(b-1)*par.px,(b-1)*par.py,1,0);
        end
    for J=1:param.nrBands*2-1
        fftDirectlyCombined=fftDirectlyCombined+shifted(:,:,J);
    end
end
% Temp1=real(ifft2(fftshift((fftDirectlyCombined))));
% Temp1(Temp1<0)=0;
% MIJ.createImage(Temp1);

w1=str2double(get(handles.edit9,'String'));  % Initial optimization Wiener constant：[0.9-2.5]
w2=0.1;
%% 
param.cutoff=1000/(0.5*param.lambda/param.NA);                        
param.sampleLateral=ceil(param.cutoff/param.cyclesPerMicron)+1;  
K0=handles.K0;
K=max([ceil(K0)]);
if param.nrBands==2
    cutoff=floor(1*K)/param.sampleLateral+1.0;
    R=K;
elseif	param.nrBands==3
    cutoff=floor(2*K)/param.sampleLateral+1.0;
    R=2*K;
end

otfHiFi=zeros(2*h,2*w);
otfHiFi=writeApoVector(otfHiFi,param.OtfProvider,cutoff);      % Ideal OTF
Mask=zeros(2*h,2*w);
Mask(otfHiFi~=0)=1;

%% Traditional Wiener-SIM
% if size(Iraw,3)==9
%     wFilter0=WienerFilterWiener_3D(param);
% else
%     wFilter0=WienerFilterWiener_3D(param);
% end
% Wk0=otfHiFi./(wFilter0.wDenom+w2^2);
% fftWiener=real(ifft2(fftshift((fftDirectlyCombined.*Wk0.*Mask))));
% Temp=fftWiener;
% Temp(Temp<0)=0;
% Wiener=255*Temp/max(max(Temp));
% MIJ.createImage(Wiener);                  % The reconstruction results are displayed in the imageJ window
% % figure, imshow(Wiener,[]);    % The reconstruction results are displayed in the matlab window
% % colormap('hot');


%% HiFi-SIM
% Step 1
% if size(Iraw,3)==9
if param.nrBands==2
    wFilter1=WienerFilterW1_2D(param);
else
    wFilter1=WienerFilterW1_3D(param);
end

Wk1=otfHiFi./(wFilter1.wDenom+w1^2);

fftInitialHiFi=fftDirectlyCombined.*Wk1.*Mask;
% Temp3=real(ifft2(fftshift((fftInitialHiFi))));
% Temp3(Temp3<0)=0;
% MIJ.createImage(Temp3);

% Step 2
if size(Iraw,3)==9
    wFilter2=WienerFilterW2_2D(param);
else
    wFilter2=WienerFilterW2_3D(param);
end

% 
if get(handles.checkbox6,'Value')
    ApoFWHM=str2double(get(handles.edit11,'String'));
else
    ApoFWHM=0.5*(cutoff-1);
    ApoFWHM=min(0.5,round(ApoFWHM*100)/100);
    set(handles.edit11,'String',ApoFWHM);
end
apo= apodize_gauss([2*h,2*w], struct('rad',ApoFWHM));  
Wk2=apo./(wFilter2.wDenom+w2^2);
fftHiFi=real(ifft2(fftshift((fftInitialHiFi.*Wk2.*Mask))));

% Results of wide field
if get(handles.checkbox4,'Value')==1
    WF=handles.WF2;
    if size(WF,1)~=size(WF,2)
        WF=importImages2(WF);
    end    
    MIJ.createImage(WF);      % The reconstruction results are displayed in the imageJ window
%     figure, imshow(WF,[]);    % The reconstruction results are displayed in the matlab window
%     colormap('hot');
    param.WF2=WF;
end
% Results of HiFi-SIM
Size1=2*param.Size1;
Size2=2*param.Size2;
HiFi=zeros(Size1,Size2);
Temp=fftHiFi;
Temp(Temp<0)=0;
Temp=255*Temp/max(max(Temp));
HiFi(1:Size1,1:Size2)=Temp(1:Size1,1:Size2);
HiFi=importImages2(HiFi);
if get(handles.checkbox5,'Value')==1
    MIJ.createImage(HiFi);  % The reconstruction results are displayed in the imageJ window
%     figure,imshow(HiFi,[]);   % The reconstruction results are displayed in the matlab window
    colormap('hot');
end

if FlagParameter==0
    waitbar(4/4,Progressbar, 'Reconstruction done!');
else
    waitbar(2/2,Progressbar, 'Reconstruction done!');
end

close(Progressbar);
FlagReconstruction=1;

param.HiFi=HiFi;
param.FlagParameter=FlagParameter;
param.FlagReconstruction=FlagReconstruction;
handles.param=param;
guidata(hObject,handles);
% save param param 

function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double

param=handles.param;
FlagApoFWHM=1;
param.FlagApoFWHM=FlagApoFWHM;
handles.param=param;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in radiobutton2.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2

if get(handles.checkbox6,'Value')
    
    set(handles.text22,'Visible','on');
    set(handles.edit14,'Visible','on');
    
    set(handles.text12,'Visible','on');
    set(handles.edit11,'Visible','on');
    set(handles.text14,'Visible','on');
    
    set(handles.text13,'Visible','on');
    set(handles.edit10,'Visible','on');
    set(handles.text15,'Visible','on');
    
    set(handles.text10,'Visible','on');
    set(handles.edit9,'Visible','on');
    set(handles.text11,'Visible','on');
    
    set(handles.edit9,'Enable','on');
    set(handles.edit10,'Enable','on');
    set(handles.edit11,'Enable','on');
else
    set(handles.text22,'Visible','off');
    set(handles.edit14,'Visible','off');
    
    set(handles.text12,'Visible','off');
    set(handles.edit11,'Visible','off');
    set(handles.text14,'Visible','off');
    
    set(handles.text13,'Visible','off');
    set(handles.edit10,'Visible','off');
    set(handles.text15,'Visible','off');
    
    set(handles.text10,'Visible','off');
    set(handles.edit9,'Visible','off');
    set(handles.text11,'Visible','off');
    
    
    set(handles.edit9,'String',1.2);
    set(handles.edit10,'String',1.0);
    set(handles.edit11,'String',0.5);

    set(handles.edit9,'Enable','off');
    set(handles.edit10,'Enable','off');
    set(handles.edit10,'Enable','off');
end
    
% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% % --- Executes on button press in checkbox6.
% function radiobutton2_Callback(hObject, eventdata, handles)
% % hObject    handle to checkbox6 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
param=handles.param;

if isfield (param,'FlagLoad')==0
    warndlg('Please load the raw  data!','!!warn!!','modal'); 
    return;
end


Format=param.Format;
filename=param.filename;
Filename=strrep(filename,['.',Format],'');

if isfield (param,'FlagReconstruction')==0
    warndlg('No results to save!','!!warn!!','modal'); 
    return;
end

param=handles.param;
WF=param.WF2;
HiFi=param.HiFi;

DataDir=get(handles.edit3,'String');
mkdir([DataDir,'\Results']);
ResultDir=[DataDir,'\Results\'];

NmaeWideField=[Filename,'_WF.tif'];
NmaeSR=[Filename,'_SIM.tif'];

imgsave32(WF,[ResultDir NmaeWideField]);     
imgsave32(HiFi,[ResultDir NmaeSR]);     

msgbox('Saved successfully！','modal');

% --- Executes during object creation, after setting all properties.
function uitable1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject, 'Data', cell(3));
