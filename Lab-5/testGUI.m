function varargout = testGUI(varargin)
% TESTGUI MATLAB code for testGUI.fig
%      TESTGUI, by itself, creates a new TESTGUI or raises the existing
%      singleton*.
%
%      H = TESTGUI returns the handle to a new TESTGUI or the handle to
%      the existing singleton*.
%
%      TESTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TESTGUI.M with the given input arguments.
%
%      TESTGUI('Property','Value',...) creates a new TESTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before testGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to testGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help testGUI

% Last Modified by GUIDE v2.5 04-Aug-2015 13:40:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @testGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @testGUI_OutputFcn, ...
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
end

% --- Executes just before testGUI is made visible.
function testGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to testGUI (see VARARGIN)

% Choose default command line output for testGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes testGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

axes(handles.axes3);
A=imread('problem.png','png');
imshow(A)
end

% --- Outputs from this function are returned to the command line.
function varargout = testGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end



function [ temp ] = initialCondition( Xj, a, b, c )
           temp = sin(Xj);
       end
        
       function [ temp ] = fi1( tk, a, b, c )
           temp = exp((c-a)*tk)*(cos(b*tk)+sin(b*tk));
       end
       
       function [ temp ] = fi2( tk, a, b, c )
           temp = (-1)*exp((c-a)*tk)*(cos(b*tk)+sin(b*tk));
       end
       
       function [ temp ] = analiticAnswer( Xj, tk, a, b, c )
           temp = exp((c-a)*tk)*sin(Xj + b*tk);
       end
       
       
       function [  ] = solve( scheme, approximation, handles, T_frame, K_frame, l_frame, N_frame, a_frame, b_frame, c_frame, alpha_frame, beta_frame, gamma_frame, delta_frame, tk_frame )
           % задание параметров
           T = T_frame;
           K = K_frame;
           l = l_frame;
           N = N_frame;
           tau = T / (K - 1);
           h = l / (N - 1);
           
           a = a_frame;
           b = b_frame;
           c = c_frame;
           
           alpha = alpha_frame;
           beta = beta_frame;
           gamma = gamma_frame;
           delta = delta_frame;
           
           sigma = a*tau/(h*h)
           if( sigma > 0.5)
               disp('Явная конечно-разностная схема неустойчива, задайте другой шаг');
           end

           U = zeros(N, K);
           % начальные условия
           f = 1;
           for Xj = 0:h:l
               U(f, 1) = initialCondition(Xj, a, b, c);  
               f = f + 1;
           end
           
           %% явная конечно-разностная схема
           if(scheme == 1)
               % двухточечная аппроксимация с первым порядком
               if(approximation == 1)
                   for k = 1:K - 1
                       for j = 2:N - 1
                           U(j, k+1) = tau*(U(j+1, k)*(a/(h*h) + b/(2*h)) + U(j, k)*((-2*a)/(h*h) + c + 1/tau) + U(j-1, k)*(a/(h*h) - b/(2*h)));
                       end
                       % граничные условия
                       tk = (k+1)*tau;
                       U(1, k+1) = (fi1(tk, a, b, c) - alpha*U(2,k+1)/h)/(beta - alpha/h);
                       U(N, k+1) = (fi2(tk, a, b, c) + gamma*U(N-1,k+1)/h)/(delta + gamma/h);
                   end
               end 
               
               % трехточечная аппроксимация со вторым порядком
               if(approximation == 2)
                   for k = 1:K - 1
                       for j = 2:N - 1
                           U(j, k+1) = tau*(U(j+1, k)*(a/(h*h) + b/(2*h)) + U(j, k)*((-2*a)/(h*h) + c + 1/tau) + U(j-1, k)*(a/(h*h) - b/(2*h)));
                       end
                       % граничные условия
                       tk = (k+1)*tau;
                       U(1, k+1) = (fi1(tk, a, b, c)*2*h + alpha*(U(3, k+1) - 4*U(2, k+1)))/(2*beta*h - 3*alpha);
                       U(N, k+1) = (fi2(tk, a, b, c)*2*h + gamma*(4*U(N-1, k+1) - U(N-2, k+1)))/(2*delta*h + 3*gamma);
                   end
               end 
               
               % двухточечная аппроксимация со вторым порядком
               if(approximation == 3)
                   for k = 1:K - 1
                       for j = 2:N - 1
                           U(j, k+1) = tau*(U(j+1, k)*(a/(h*h) + b/(2*h)) + U(j, k)*((-2*a)/(h*h) + c + 1/tau) + U(j-1, k)*(a/(h*h) - b/(2*h)));
                       end
                       % граничные условия
                       tk = (k+1)*tau;
                       U(1, k+1) = (h/tau*U(1,k) - fi1(tk, a, b, c)*(2*a*a-b*h)/alpha + 2*a*a/h*U(2,k+1)) / (2*a*a/h + h/tau - c*h - beta/alpha*(2*a*a-b*h));
                       U(N, k+1) = (h/tau*U(N,k) + fi2(tk, a, b, c)*(2*a*a+b*h)/gamma + 2*a*a/h*U(N-1,k+1)) / (2*a*a/h + h/tau - c*h + delta/gamma*(2*a*a+b*h)); 
                   end
               end                     
           end

           
           %% неявная конечно-разностная схема
           if(scheme == 2)
               % двухточечная аппроксимация с первым порядком
               if(approximation == 1)
                   A = zeros(N, N);
                   d = zeros(N, 1);
                   for k = 1:K - 1
                       tk = (k+1)*tau;                      
                       
                       A(1, 1) = 1;
                       A(1, 2) = alpha/(beta*h - alpha);
                       A(N, N-1) = -gamma/(delta*h + gamma);
                       A(N, N) = 1;
                       d(1) = fi1(tk, a, b, c)/(beta - alpha/h);
                       d(N) = fi2(tk, a, b, c)/(delta + gamma/h);
                                             
                       q = 2;
                       for j = 1:N-2
                           A(q, j) = b/(2*h) - a/(h*h);
                           A(q, j+1) = 2*a/(h*h) + 1/tau - c;
                           A(q, j+2) = -(a/(h*h) + b/(2*h));
                           d(j+1) = 1/tau*U(j+1, k);
                           q = q + 1;
                       end
                                              
                       X = linsolve(A, d);
                       for j = 1:N
                          U(j, k+1) = X(j); 
                       end                     
                   end
               end 
               
               % трехточечная аппроксимация со вторым порядком
               if(approximation == 2)
                   A = zeros(N, N);
                   d = zeros(N, 1);
                   for k = 1:K - 1
                       tk = (k+1)*tau;                                    
                       
                         A(1, 1) = 1;
                         A(1, 2) = 4*alpha/(2*beta*h - 3*alpha);
                         A(1, 3) = -alpha/(2*beta*h - 3*alpha);
                         A(N, N-2) = gamma/(2*delta*h + 3*gamma);
                         A(N, N-1) = -4*gamma/(2*delta*h + 3*gamma);
                         A(N, N) = 1;
                         d(1) = 2*h*fi1(tk, a, b, c)/(2*beta*h - 3*alpha);
                         d(N) = 2*h*fi2(tk, a, b, c)/(2*delta*h + 3*gamma);
                                             
                       q = 2;
                       for j = 1:N-2
                           A(q, j) = b/(2*h) - a/(h*h);
                           A(q, j+1) = 2*a/(h*h) + 1/tau - c;
                           A(q, j+2) = -(a/(h*h) + b/(2*h));
                           d(j+1) = 1/tau*U(j+1, k);
                           q = q + 1;
                       end
                                              
                       X = linsolve(A, d);
                       for j = 1:N
                          U(j, k+1) = X(j); 
                       end                     
                   end                   
               end 
               
               % двухточечная аппроксимация со вторым порядком
               if(approximation == 3)
                   A = zeros(N, N);
                   d = zeros(N, 1);
                   for k = 1:K - 1
                       tk = (k+1)*tau;                      
                       
                       A(1, 1) = 2*a*a/h + h/tau -c*h - beta/alpha*(2*a*a - b*h);
                       A(1, 2) = -2*a*a/h;
                       A(N, N-1) = -2*a*a/h;
                       A(N, N) = 2*a*a/h + h/tau -c*h + delta/gamma*(2*a*a + b*h);
                       d(1) = h/tau*U(1, k) - fi1(tk, a, b, c)*(2*a*a - b*h)/alpha;
                       d(N) = h/tau*U(N, k) + fi2(tk, a, b, c)*(2*a*a + b*h)/gamma;
                                             
                       q = 2;
                       for j = 1:N-2
                           A(q, j) = b/(2*h) - a/(h*h);
                           A(q, j+1) = 2*a/(h*h) + 1/tau - c;
                           A(q, j+2) = -(a/(h*h) + b/(2*h));
                           d(j+1) = 1/tau*U(j+1, k);
                           q = q + 1;
                       end
                                              
                       X = linsolve(A, d);
                       for j = 1:N
                          U(j, k+1) = X(j); 
                       end                     
                   end
               end                     
           end
           
           %% конечно-разностная схема Кранка-Николсона
           if(scheme == 3)
               % двухточечная аппроксимация с первым порядком
               if(approximation == 1)
                   A = zeros(N, N);
                   d = zeros(N, 1);
                   for k = 1:K - 1
                       tk = (k+1)*tau;                      
                       
                       A(1, 1) = 1;
                       A(1, 2) = alpha/(beta*h - alpha);
                       A(N, N-1) = -gamma/(delta*h + gamma);
                       A(N, N) = 1;
                       d(1) = fi1(tk, a, b, c)/(beta - alpha/h);
                       d(N) = fi2(tk, a, b, c)/(delta + gamma/h);
                                             
                       q = 2;
                       for j = 1:N-2
                           A(q, j) = b/(4*h) - a/(2*h*h);
                           A(q, j+1) = a/(h*h) + 1/tau - c/2;
                           A(q, j+2) = -(a/(2*h*h) + b/(4*h));
                           d(j+1) = 1/tau*U(j+1, k) + a/(2*h*h)*(U(j+2, k) - 2*U(j+1, k) + U(j, k)) + b/(4*h)*(U(j+2, k) - U(j, k)) + c/2 * U(j+1, k);
                           q = q + 1;
                       end
                                              
                       X = linsolve(A, d);
                       for j = 1:N
                          U(j, k+1) = X(j); 
                       end                     
                   end
               end 
               
               % трехточечная аппроксимация со вторым порядком
               if(approximation == 2)
                   A = zeros(N, N);
                   d = zeros(N, 1);
                   for k = 1:K - 1
                       tk = (k+1)*tau;                                    
                       
                         A(1, 1) = 1;
                         A(1, 2) = 4*alpha/(2*beta*h - 3*alpha);
                         A(1, 3) = -alpha/(2*beta*h - 3*alpha);
                         A(N, N-2) = gamma/(2*delta*h + 3*gamma);
                         A(N, N-1) = -4*gamma/(2*delta*h + 3*gamma);
                         A(N, N) = 1;
                         d(1) = 2*h*fi1(tk, a, b, c)/(2*beta*h - 3*alpha);
                         d(N) = 2*h*fi2(tk, a, b, c)/(2*delta*h + 3*gamma);
                                             
                       q = 2;
                       for j = 1:N-2
                           A(q, j) = b/(4*h) - a/(2*h*h);
                           A(q, j+1) = a/(h*h) + 1/tau - c/2;
                           A(q, j+2) = -(a/(2*h*h) + b/(4*h));
                           d(j+1) = 1/tau*U(j+1, k) + a/(2*h*h)*(U(j+2, k) - 2*U(j+1, k) + U(j, k)) + b/(4*h)*(U(j+2, k) - U(j, k)) + c/2 * U(j+1, k);
                           q = q + 1;
                       end
                                              
                       X = linsolve(A, d);
                       for j = 1:N
                          U(j, k+1) = X(j); 
                       end                     
                   end                   
               end 
               
               % двухточечная аппроксимация со вторым порядком
               if(approximation == 3)
                   A = zeros(N, N);
                   d = zeros(N, 1);
                   for k = 1:K - 1
                       tk = (k+1)*tau;                      
                       
                       A(1, 1) = 2*a*a/h + h/tau -c*h - beta/alpha*(2*a*a - b*h);
                       A(1, 2) = -2*a*a/h;
                       A(N, N-1) = -2*a*a/h;
                       A(N, N) = 2*a*a/h + h/tau -c*h + delta/gamma*(2*a*a + b*h);
                       d(1) = h/tau*U(1, k) - fi1(tk, a, b, c)*(2*a*a - b*h)/alpha;
                       d(N) = h/tau*U(N, k) + fi2(tk, a, b, c)*(2*a*a + b*h)/gamma;
                                             
                       q = 2;
                       for j = 1:N-2
                           A(q, j) = b/(4*h) - a/(2*h*h);
                           A(q, j+1) = a/(h*h) + 1/tau - c/2;
                           A(q, j+2) = -(a/(2*h*h) + b/(4*h));
                           d(j+1) = 1/tau*U(j+1, k) + a/(2*h*h)*(U(j+2, k) - 2*U(j+1, k) + U(j, k)) + b/(4*h)*(U(j+2, k) - U(j, k)) + c/2 * U(j+1, k);
                           q = q + 1;
                       end
                                              
                       X = linsolve(A, d);
                       for j = 1:N
                          U(j, k+1) = X(j); 
                       end                     
                   end
               end                     
           end
        
           
           %% построение слоя, соответствующего tk
           analiticAnsw = zeros(1, N);
           UAnsw = zeros(1, N);
           x = zeros(1, N);

           k = tk_frame;
           tk = k*tau;
          

           f = 1;
           for Xj = 0:h:l
               analiticAnsw(f) = analiticAnswer( Xj, tk, a, b, c );
               UAnsw(f) = U(f, k);
               x(f) = Xj;
               f = f + 1;
           end

           %figure;
           axes(handles.axes1);
           plot(x, UAnsw, 'r', x, analiticAnsw, 'k');
           %hold on;
           %plot(x, analiticAnsw, 'k');
            
           % зависимость погрешности от tk
           for n = 0:N - 1   
               tk = (n * floor((K-1)/(N-1)) + 1)* tau;
               f = 1;
               for Xj = 0:h:l
                   analiticAnsw(f) = analiticAnswer( Xj, tk, a, b, c );
                   UAnsw(f) = U(f, n*floor((K-1)/(N-1)) + 1);
                   x(f) = Xj;
                   f = f + 1;
               end
               for i = 1:N
                   maxt(i) = abs(UAnsw(i) - analiticAnsw(i));
               end
               error(n + 1) = max(maxt);
           end

           %figure;
           axes(handles.axes2);
           plot(x, error, 'r');

           

       end





function SchemeChoice_Callback(hObject, eventdata, handles)
% hObject    handle to SchemeChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SchemeChoice as text
%        str2double(get(hObject,'String')) returns contents of SchemeChoice as a double
end

% --- Executes during object creation, after setting all properties.
function SchemeChoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SchemeChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function ApproximationChoice_Callback(hObject, eventdata, handles)
% hObject    handle to ApproximationChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ApproximationChoice as text
%        str2double(get(hObject,'String')) returns contents of ApproximationChoice as a double
end

% --- Executes during object creation, after setting all properties.
function ApproximationChoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ApproximationChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%scheme = str2double(get(handles.SchemeChoice, 'String'));
%approximation = str2double(get(handles.ApproximationChoice, 'String'));
T_frame = str2double(get(handles.TChoice, 'String'));
K_frame = str2double(get(handles.KChoice, 'String'));
l_frame = str2double(get(handles.lChoice, 'String'));
N_frame = str2double(get(handles.NChoice, 'String'));
a_frame = str2double(get(handles.aChoice, 'String'));
b_frame = str2double(get(handles.bChoice, 'String'));
c_frame = str2double(get(handles.cChoice, 'String'));
alpha_frame = str2double(get(handles.AlphaChoice, 'String'));
beta_frame = str2double(get(handles.BetaChoice, 'String'));
gamma_frame = str2double(get(handles.GammaChoice, 'String'));
delta_frame = str2double(get(handles.DeltaChoice, 'String'));
tk_frame = str2double(get(handles.tkChoice, 'String'));

if(get(handles.radiobutton2, 'Value') == 1)
    scheme = 1
end
if(get(handles.radiobutton3, 'Value') == 1)
    scheme = 2
end
if(get(handles.radiobutton4, 'Value') == 1)
    scheme = 3
end

if(get(handles.radiobutton5, 'Value') == 1)
    approximation = 1
end
if(get(handles.radiobutton6, 'Value') == 1)
    approximation = 2
end
if(get(handles.radiobutton7, 'Value') == 1)
    approximation = 3
end

solve( scheme, approximation, handles, T_frame, K_frame, l_frame, N_frame, a_frame, b_frame, c_frame, alpha_frame, beta_frame, gamma_frame, delta_frame, tk_frame )
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles(GUI_EVENT_CLOSE)
end




           
            

       



function TChoice_Callback(hObject, eventdata, handles)
% hObject    handle to TChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TChoice as text
%        str2double(get(hObject,'String')) returns contents of TChoice as a double
end

% --- Executes during object creation, after setting all properties.
function TChoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function KChoice_Callback(hObject, eventdata, handles)
% hObject    handle to KChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of KChoice as text
%        str2double(get(hObject,'String')) returns contents of KChoice as a double
end

% --- Executes during object creation, after setting all properties.
function KChoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function lChoice_Callback(hObject, eventdata, handles)
% hObject    handle to lChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lChoice as text
%        str2double(get(hObject,'String')) returns contents of lChoice as a double
end

% --- Executes during object creation, after setting all properties.
function lChoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function NChoice_Callback(hObject, eventdata, handles)
% hObject    handle to NChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NChoice as text
%        str2double(get(hObject,'String')) returns contents of NChoice as a double
end

% --- Executes during object creation, after setting all properties.
function NChoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function AlphaChoice_Callback(hObject, eventdata, handles)
% hObject    handle to AlphaChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AlphaChoice as text
%        str2double(get(hObject,'String')) returns contents of AlphaChoice as a double
end

% --- Executes during object creation, after setting all properties.
function AlphaChoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AlphaChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function BetaChoice_Callback(hObject, eventdata, handles)
% hObject    handle to BetaChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BetaChoice as text
%        str2double(get(hObject,'String')) returns contents of BetaChoice as a double
end

% --- Executes during object creation, after setting all properties.
function BetaChoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BetaChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function GammaChoice_Callback(hObject, eventdata, handles)
% hObject    handle to GammaChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GammaChoice as text
%        str2double(get(hObject,'String')) returns contents of GammaChoice as a double
end

% --- Executes during object creation, after setting all properties.
function GammaChoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GammaChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function aChoice_Callback(hObject, eventdata, handles)
% hObject    handle to aChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of aChoice as text
%        str2double(get(hObject,'String')) returns contents of aChoice as a double
end

% --- Executes during object creation, after setting all properties.
function aChoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to aChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function bChoice_Callback(hObject, eventdata, handles)
% hObject    handle to bChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bChoice as text
%        str2double(get(hObject,'String')) returns contents of bChoice as a double
end

% --- Executes during object creation, after setting all properties.
function bChoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function cChoice_Callback(hObject, eventdata, handles)
% hObject    handle to cChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cChoice as text
%        str2double(get(hObject,'String')) returns contents of cChoice as a double
end

% --- Executes during object creation, after setting all properties.
function cChoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function DeltaChoice_Callback(hObject, eventdata, handles)
% hObject    handle to DeltaChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DeltaChoice as text
%        str2double(get(hObject,'String')) returns contents of DeltaChoice as a double
end

% --- Executes during object creation, after setting all properties.
function DeltaChoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DeltaChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function tkChoice_Callback(hObject, eventdata, handles)
% hObject    handle to tkChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tkChoice as text
%        str2double(get(hObject,'String')) returns contents of tkChoice as a double
end

% --- Executes during object creation, after setting all properties.
function tkChoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tkChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
