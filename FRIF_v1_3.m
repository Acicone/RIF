function [IMF,stats] = FRIF_v1_3(f,L,options,PickPower)


%
%  function IMF = FRIF_v1_3(f,options,L)
%
% It generates the decomposition of the signal f :
%
%  f = IMF(1,:) + IMF(2,:) + ... + IMF(K, :)
%
% where the last row in the matrix IMF is the trend and the other rows
% are actual IMFs
%
%                                Inputs
%
%   f         Signal to be decomposed
%
%   L         Mask length values for each Inner Loop
%
%   options    Structure, generated using function Settings_FRIF_v1, containing
%              all the parameters needed in the various algorithms
%
%                               Output
%
%   IMF       Matrices containg in row i the i-th IMF. The last row
%              contains the remainder-trend.
%
%   stats     Statistics regarding the IMFs
%               logG     Resempling function used for each IMF
%               inStepN  Inner loop steps number for each IMF
%
%   See also Settings_FRIF_v1, GET_MASK_V1_1, MAXMINS_v3_8, PLOT_IMF_V10.
%
%  Please cite:
%
%  A. Cicone. 'Iterative Filtering as a direct method for the decomposition
%  of nonstationary signals'. Numerical Algorithms, Volume 373, 2020,  112248.
%  doi: 10.1007/s11075-019-00838-z
%  arXiv http://arxiv.org/abs/1811.03536
%
%  G. Barbarino, A. Cicone. 'Stabilization and Variations to the Adaptive 
%  Local Iterative Filtering Algorithm: the Fast Resampled Iterative Filtering Method'
%  Submitted 2021
%  arXiv http://arxiv.org/abs/2111.02764


%% we deal with the input

tol=10^-14;
%load('TriangularFilter')
load('prefixed_double_filter','MM');

if nargin < 1,  help FRIF_v1_3; return; end
if nargin < 2
    disp('FRIF needs as input a mask lenght matrix')
    IMF=[];
    stats=[];
    return
end
if nargin < 3, options = Settings_FRIF_v1; end

N = length(f);
if size(f,1)>size(f,2)
    f = f.';
end
if size(f,1)>1
    disp('Wrong dataset, the signal must be a single row vector')
    disp('If you have a multivariate signal you can try the MvFIF code')
    disp('If, instead, is a multidimensional data set, plese try the FIF2 code')
    IMF=[];
    stats=[];
    return
end

sL=size(L);
if sL(1)==N
    L=L.';
    sL=size(L);
elseif sL(2)==N
    if options.NIMFs==sL(1)
        options.NIMFs=sL(1);
    end
else
    disp('Wrong size for the mask lenght matrix')
    disp('L should contain as many rows as IMFs required to extract and as')
    disp('many columns as the length of the input signal')
    IMF=[];
    stats=[];
    return
end

if nargin <4, 
    PickPower = -10*ones(1,sL(1)); 
end

Norm1f=norm(f,1); % to avoid dealing with way too small values
f=f/Norm1f;

% we initialize the matrix containing the IMFs
IMF =zeros(sL(1)+1,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Risampled Iterative Filtering 2         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:sL(1)
    % we construct the G resampling function starting from the mask length
    % L(x) provided by the user
    G=zeros(1,N);
    counter=1;
    G(counter)=1;
    while G(counter)<N
        counter=counter+1;
        % we devide by 2 L values to increase the number of sample points
        % captured by the mask (we clearly double the mask length)
        G(counter)=G(counter-1)+interp1(1:N,L(i,:),G(counter-1),'linear')/options.UpSampling;
        % ATT! we interpolate here because, otherwise, we might be problems if L values are less than 4
    end
    % we remove the extra entries in G and the final one if bigger than N
    if G(counter)>N
        G(counter:end)=[];
    else
        G(counter+1:end)=[];
    end
    % we save the function G
    stats(i).logG=G;
    SD=1;
    h=interp1(1:N,f,G,'pchip');
    Nh = length(h);
    if options.plots>0.5
    figure
    plot(h)
    title(['Signal - IMF ' num2str(i)])
    end
    h_pp=h;
    h_pp(abs(h)<=tol)=[];
    if isempty(h_pp)
        disp('Signal too small')
        IMF=[];
        stats=[];
        return
    end
    maxmins_pp=Maxmins_v3_8(h_pp,tol);
    diffMaxmins_pp=diff(maxmins_pp);
    %     N_pp=length(h_pp);
    %     k_pp = length(maxmins_pp);
    Num_iter=2;
%     if i==1
%         Num_iter=2;
%     else
%         Num_iter=1;
%     end
    for jj=1:Num_iter
        length_fixed_Mask = round(2*options.Xi*prctile(diffMaxmins_pp,options.alpha));
        
        inStepN=0;
        if options.verbose>0
            fprintf('\n IMF # %1.0d\n',i)
            fprintf('\n  step #            SD\n\n')
        end
        a = get_mask_v1_1(MM,length_fixed_Mask,options.verbose,tol);
        ExtendSig=1==0;
        
        if Nh < length(a) % we need to extend the signal
            ExtendSig=1==1;
            Nxs=ceil(length(a)/Nh);
            N_old=Nh;
            if rem(Nxs,2)==0
                Nxs=Nxs+1;
            end
            h_n=[];
            for ii=1:Nxs
                h_n=[h_n h];
            end
            h=h_n;
            Nh=Nxs*Nh;
        end
        
        Nza=Nh-length(a);
        if rem(Nza,2)==0
            a = [zeros(1,Nza/2) a zeros(1,Nza/2)];
            fftA=real(fft([a((length(a)-1)/2+1:end) a(1:(length(a)-1)/2)]));
            % figure,plot(circshift(a,(length(a)-1)/2+1)-ifft(real(fft(circshift(a,(length(a)-1)/2+1)))),'r')
        else
            a = [zeros(1,(Nza-1)/2) a zeros(1,(Nza-1)/2+1)];
            %csA=circshift(a,(length(a))/2+1);
            fftA=real(fft([a((length(a))/2:end) a(1:(length(a))/2-1)]));
            % figure,plot(circshift(a,(length(a))/2+1)-ifft(real(fft(circshift(a,(length(a))/2+1)))),'r')
        end
        
        fftH=fft(h);
        fft_h_new=fftH;
        
        
        if jj==1 && Num_iter==2
            % We fix Xi in order to allign the filter minimum with the signal
            % maximum in the Fourier domain
            temp_s=log(abs(fftH(1:(end-1)/2)).^2);
            Sm = find(islocalmax(temp_s));
            Pos_s = find(temp_s(Sm)>PickPower(i));
            temp_w=log(abs(fftA(1:(end-1)/2)).^2);
            Fm = find(islocalmin(temp_w));                       
            [val_w,pos_w]=mink(abs(Fm-Sm(Pos_s(end))),3); % we look for the first 
            % three closest minima in the filter fourier transform to the 
            % highest frequency maximun in the signal fourier transform with value bigger than 10^-10
            %[val2_w,pos2_w]=min(temp_w(Fm(pos_w)));
            [val2_w,pos2_w]=min(pos_w);
            Xi_new = Fm(pos_w(pos2_w))/Sm(Pos_s(end));
            if options.plots>0.5
            figure
            plot(temp_s)
            hold on
            plot(temp_w)
            plot(Fm(pos_w(pos2_w)),temp_w(Fm(pos_w(pos2_w))),'rx')
            plot(Sm(Pos_s(end)),temp_s(Sm(Pos_s(end))),'rx')
            end
            
            options.Xi=options.Xi*Xi_new;
        end
    end
    if options.plots>0.5 %&& rem(inStepN,5)==0
        if gcf > 30
            close all
        end
        figN=figure;
        set(figN,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    end
    if options.plots>=1
        figMask=figure;
        figRem=figure;
        set(figMask,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    end
    
    % Plotting of the FFT of the filter and signal
    if options.plots>=0.5
        figure
        semilogy(abs(fftA(1:(end-1)/2)).^2,'r')
        hold on
        semilogy(abs(fftH(1:(end-1)/2)).^2,'k')
        legend('Filter FFT','Signal FFT')
    end
    while SD>options.delta && inStepN < options.MaxInner
        inStepN=inStepN+options.NumSteps;
        
        fft_h_old=(1-fftA).^(inStepN-1).*fftH;
        fft_h_new=(1-fftA).^inStepN.*fftH;
        
        %%%%%%%%%%%%%%%% Updating stopping criterium %%%%%%%%%%%%%%%%%
        
        SD=norm(fft_h_new-fft_h_old)^2/norm(fft_h_old)^2;
        
        %%%%%%%%%%%%%%%%%% generating f_n %%%%%%%%%%%%%%%%%%
        
        if options.verbose>0
            fprintf('    %2.0d      %1.40f\n',inStepN,SD)
        end
        
        if options.plots>=1  && rem(inStepN,2)==0
            figure(figMask)
            title(['IMF ' num2str(i) ' step # ' num2str(inStepN) ])
            plot(ifft(fft_h_new)*Norm1f,'linewidth',2)
            figure(figRem)
            title(['Remainder after IMF ' num2str(i) ' step # ' num2str(inStepN) ])
            plot((f-ifft(fft_h_new))*Norm1f,'linewidth',2)
            
            pause(0.01)
        end
        
    end
    
    h=ifft(fft_h_new);
    
    if ExtendSig % we reduce the signal
        Nh=N_old;
        h=h(Nh*(Nxs-1)/2+1:Nh*((Nxs-1)/2+1));
    end
    if inStepN >= options.MaxInner && options.verbose>0
        disp('Max # of inner steps reached')
        %return
    end
    stats(i).inStepN=inStepN;
    
    % we resample the IMF back to the original sampling rate
    h_res=interp1(G,h,1:N,'pchip');
    
    IMF(i,:) = h_res;
    f=f-h_res;
    
    if options.saveInter==1
        save([nameFile '_intermediate_FRIF_v1_3.mat'],'IMF','f','stats','-v7.3');
    end
    
end


IMF = [IMF(1:i,:); f];

IMF=IMF*Norm1f; % we scale back to the original values

if options.plots>=1
    if gcf > 30
        close all
    end
    figN=plot_imf_v10(IMF,1:N);
    for ii=1:length(figN)
        set(figN(ii),'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        if options.saveplots>0
            saveas(figN(ii),[nameFile '_IMFs'], 'fig')
            saveas(figN(ii),[nameFile '_IMFs'], 'epsc')
            saveas(figN(ii),[nameFile '_IMFs'], 'png')
        end
        
    end
end

if options.saveEnd == 1
    save([ 'Final_' nameFile '_FRIF_v1_3.mat'],'IMF','stats','-v7.3');
end

end


%% Auxiliar functions


function a=get_mask_v1_1(y,k,verbose,tol)
%
% Rescale the mask y so that its length becomes 2*k+1.
% k could be an integer or not an integer.
% y is the area under the curve for each bar

n=length(y);
m=(n-1)/2;

if k<=m % The prefixed filter contains enough points
    
    if mod(k,1)==0     % if the mask_length is an integer
        
        a=zeros(1,2*k+1);
        
        for i=1:2*k+1
            s=(i-1)*(2*m+1)/(2*k+1)+1;
            t=i*(2*m+1)/(2*k+1);
            
            %s1=s-floor(s);
            s2=ceil(s)-s;
            
            t1=t-floor(t);
            %t2=ceil(t)-t;
            
            if floor(t)<1
                disp('Ops')
            end
            a(i)=sum(y(ceil(s):floor(t)))+s2*y(ceil(s))+t1*y(floor(t));
        end
        
    else   % if the mask length is not an integer
        new_k=floor(k);
        extra = k-new_k;
        c=(2*m+1)/(2*new_k+1+2*extra);
        
        a=zeros(1,2*new_k+3);
        
        t=extra*c+1;
        t1=t-floor(t);
        %t2=ceil(t)-t;
        if k<0
            disp('Ops')
            a=[];
            return
        end
        a(1)=sum(y(1:floor(t)))+t1*y(floor(t));
        
        for i=2:2*new_k+2
            s=extra*c+(i-2)*c+1;
            t=extra*c+(i-1)*c;
            %s1=s-floor(s);
            s2=ceil(s)-s;
            
            t1=t-floor(t);
            
            
            a(i)=sum(y(ceil(s):floor(t)))+s2*y(ceil(s))+t1*y(floor(t));
        end
        t2=ceil(t)-t;
        
        a(2*new_k+3)=sum(y(ceil(t):n))+t2*y(ceil(t));
    end
else % We need a filter with more points than MM, we use interpolation
    dx=0.01;
    % we assume that MM has a dx = 0.01, if m = 6200 it correspond to a
    % filter of length 62*2 in the physical space
    f=y/dx; % function we need to interpolate
    dy=m*dx/k;
    b=interp1(0:m,f(m+1:2*m+1),0:m/k:m);
    if size(b,1)>size(b,2)
        b=b.';
    end
    if size(b,1)>1
        fprintf('\n\nError!')
        disp('The provided mask is not a vector!!')
        a=[];
        return
    end
    a=[fliplr(b(2:end)) b]*dy;
    if abs(norm(a,1)-1)>tol
        if verbose>0
            fprintf('\n\n Warning!\n\n')
            fprintf(' Area under the mask equals %2.20f\n',norm(a,1))
            fprintf(' it should be equal to 1\n We rescale it using its norm 1\n\n')
        end
        a=a/norm(a,1);
    end
end

end