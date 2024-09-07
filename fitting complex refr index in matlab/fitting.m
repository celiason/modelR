% This program makes fitting of dielectric function by arbitrary number 
% of Debye, Drude, Lorentz and modified Lorentz terms

% Using modified Lorentz terms, in which dielectric polarization depends 
% both on the electric field and its first time derivative,
% makes fitting more flexible and allows to fit complicated dispersion profiles
% For details of its realization in FDTD, see
% A. Deinega and S. John, "Effective optical response of silicon to sunlight 
% in the finite-difference time-domain method," Opt. Lett. 37, 112-114 (2012) 

% there are five different types of terms
% 1 - just eps0
% 2 - Debye, deps/(1-1.i*gamma*w)
% 3 - Drude, deps*w0^2/(-2.i*gamma*w-w^2)
% 4 - Lorentz, deps*w0^2/(w0^2-2.i*gamma*w-w^2)
% 5 - modified Lorentz, deps*(w0^2-1.i*gamma_num*w)/(w0^2-2.i*gamma*w-w^2)

% name of each type
terms_name={'' 'debye' 'drude' 'lorentz' 'm-lorentz'};
% number of coefficients in each type
terms_sz=[1 2 3 3 4];
% name of coefficients in each type
terms_coef_name={{'eps0'} {'deps' 'gamma'} {'deps' 'w0' 'gamma'} {'deps' 'w0' 'gamma'} {'deps' 'w0' 'gamma' 'gamma'''}};

%%%%%%%%%%%%%%%% THIS IS VARIABLES ARE MODIFIED BY USER %%%%%%%%%%%%%%%%%%%
% start values of coefficients
start_val={{1} {1 1} {1 1 1} {1 1 1} {1 1 1 1}};
% minimal and maximal values for each coefficient
min_val={{1} {0 0} {-.2 0 0} {-.2 0 0} {-10 0 0 0}};
max_val={{1.00001} {10 2.5} {10 5 2.5} {10 5 2.5} {10 5 2.5 25}};
% number of terms of each type
% first type (eps0) is only one, therefore first
% element of this array should be 1
terms_num=[1 0 0 20 0];
% number of iterations
% at each iteration matlab function fminsearch is called
iterations_number=30;

% data to be fitted is stored in file which can be found by path and file
Path='/Users/Chad/Desktop/';
FileName='melanin.nk';
% file should contain table with 3 columns
% first column coresponds to wavelength (if arg="lambda') or frequency (otherwise)
arg='lambda';
% second and third columns correspond to real and imaginary part of epsilon (if nk=0)
% or n and k (otherwise)
nk=1;
% fitting is made only between vmin and vmax values in the first column
arg_min=300;
arg_max=1000;
% values in the first column are in increasing (inv=0) 
% or decreasing (inv=1) order
inv=0;
% if table_step>1 than table from file will be thinned out with this step
table_step=1;
% you can alse modify residue.m to control fitting process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load table from file
Data=load([Path FileName]);
% calculate size of the table
[Data_nrow,Data_ncol]=size(Data);

% if n and k are given, rewrite 2nd and 3rd slices of Data with real and imaginary part of epsilon     
if (nk==1)
    for i=1:Data_nrow-1
        eps_re=Data(i,2)*Data(i,2)-Data(i,3)*Data(i,3);
        eps_im=2*Data(i,2)*Data(i,3);
        Data(i,2)=eps_re;
        Data(i,3)=eps_im;
    end
end

% rearange table in decreasing arguments order case
if (inv==1)
    Data_inv = Data;
    for i=1:Data_nrow
        for j=1:Data_ncol
            Data_inv(i,j)=Data(Data_nrow+1-i,j);
        end
    end
    Data = Data_inv;
end

% Data2 takes values from Data in [arg_min,arg_max] arguments range 
% and thin them out if table_step>1
j=-1;
Data2=0;
for i=1:table_step:(Data_nrow-1)
    if ((j<1)&(Data(i,1)>=arg_min))
        j=1;
    end
    if ((Data(i,1)>arg_max))
        break;
    end
    if (j>0)
        Data2(j,1)=Data(i,1);
        Data2(j,2)=Data(i,2);
        Data2(j,3)=Data(i,3) ;       
        j=j+1;
    end
end
[Data2_nrow,Data2_ncol]=size(Data2);

% total number of coefficients
coef_size=0;
for i=1:5
    coef_size=coef_size+terms_num(i)*terms_sz(i);
end

% x - values of all coefficients
x=ones(coef_size,1);
% lower and upper limit values for coefficients 
lb=ones(coef_size,1); 
ub=ones(coef_size,1);

ii=1;
for i=1:5
    for k=1:terms_num(i)    
        for j=1:terms_sz(i)
            lb(ii)=min_val{i}{j};
            ub(ii)=max_val{i}{j};
            x(ii)=start_val{i}{j};
            ii=ii+1;
        end
    end
end


for it=0:iterations_number % fitting cycle
    
    options=optimset('TolX',1E-15);
%   options=optimset('MaxFunEvals', 10000);

    if (it>0)
        [x,f] = fminsearch(@residue,x,options,Data2,terms_num,terms_sz,arg,lb,ub,it);
    end

    % table with epsilon calculated by current fiting
    Result=zeros(Data2_nrow,3);
    for i=1:Data2_nrow
        eps=calculate_eps(x,Data2(i),terms_num,terms_sz,arg,lb,ub);
        Result(i,2)=real(eps);
        Result(i,3)=imag(eps);
    end

    plot_title=sprintf('iteration number: %d',it);
    plot_eps(Data2(:,1),Data2(:,2),Data2(:,3),Result(:,2),Result(:,3),arg,plot_title);
    
    % if x is outside interval [lb,ub], move it back,
    % assuming that calculate_eps function is periodic with one period [lb,ub]
    xx = move_inside_period(x,lb,ub);
    
    format long;
    
    message=sprintf('iteration number: %d\n',it);
    ii=1;
    for i=1:5
        for k=1:terms_num(i)    
            if(terms_num(i)>=1 && length(terms_name{i})>0)
                message=sprintf('%s%d %s: ', message, k, terms_name{i});
            end
            for j=1:terms_sz(i)
                message=sprintf('%s%s=%f ', message, terms_coef_name{i}{j}, xx(ii));
                ii=ii+1;
            end
            message=sprintf('%s\n', message);
        end
    end
    message
    pause(1);
    
% sc=1;
% xls=zeros(5,1+2+5);
% xls(1,1)=xx(1);
% j=2;
% for i=1:terms_num(3)
%     xls(2,i+1)=xx(j);
%     xls(3,i+1)=xx(j+1)*sc;
%     xls(4,i+1)=2*xx(j+2)*sc;
%     xls(5,i+1)=0;    
%     j=j+terms_sz(3);
% end
% for i=1:terms_num(4)
%     xls(2,i+1+2)=xx(j);
%     xls(3,i+1+2)=xx(j+1)*sc;
%     xls(4,i+1+2)=2*xx(j+2)*sc;
%     xls(5,i+1)=0;    
%     j=j+terms_sz(4);
% end
% for i=1:terms_num(5)
%     xls(2,i+1+2+terms_num(4))=xx(j);
%     xls(3,i+1+2+terms_num(4))=xx(j+1)*sc;
%     xls(4,i+1+2+terms_num(4))=2*xx(j+2)*sc;
%     xls(5,i+1+2+terms_num(4))=xx(j+3)*sc;
%     j=j+terms_sz(5);
% end
% fid = fopen('xls.txt','w');
% [xls_nrow,xls_ncol]=size(xls);
% for i=1:xls_nrow
%     for j=1:xls_ncol
%         fprintf(fid,'%f',xls(i,j));
%         if j<xls_ncol
%             fprintf(fid,'%c',9);
%         end
%     end
%     if i<xls_nrow
%         fprintf(fid,'\n');
%     end
% end
% fclose(fid);
    
end