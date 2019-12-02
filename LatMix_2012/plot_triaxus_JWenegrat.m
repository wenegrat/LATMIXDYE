% Example code to show the dye data I'm sending to Jacob Wenegrat from LatMix 2012
%
% Written by Miles A. Sundermeyer, 6/7/19, for Wenegrat et al., paper

transectflag = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load all four injections in inj_data_2012
load inj_data_2012;		% variable name is InjData{1-4}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identify the times we want for a particular survey, and get the files
% encompassing those times
exptName{1} = 'ITE';
exptName{2} = 'SI1';
exptName{3} = 'SI2';
exptName{4} = 'Fil';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables in my_triaxus_**.mat files are:
% Time series versions
% S_ts			(PSU) Salinity
% T_ts			(deg C) Temperature
% depth_ts		(m) depth
% fluor_ts		(ppb) un-corrected fluorescein
% fluorClean_ts		(ppb) corrected fluorescein
% fluorPPB_ts		(ppb) corrected + calibrated fluorescein
% heading_ts		(degrees N) ships heading
% jday_ts		(days) decimal day, with Jan 1, noon = 0.5
% lat_ts		(deg N) Latitude
% lon_ts		(deg E) Longitude
% pdens_ts		(kg/m^3) potential density (-1000)
% rhd_ts		(ppb) un-corrected rhodamine
% shiplog_ts		(km) ship log (integrated distance traveled)
% pr			misc units - Triaxus data with various degrees of corrections, binned in pressure
% ts			misc units - Triaxus time series data with various degrees of corrections
%
% Binned in pressure, with successive casts in successive columns of the arrays
% DEPTH			(m) depth
% FLUOR			(ppb) un-corrected fluorescein
% FLUORCLEAN		(ppb) corrected fluorescein
% FLUORPPB		(ppb) corrected + calibrated fluorescein
% HEADING		(degrees N) ships heading
% JDAY			(days) decimal day, with Jan 1, noon = 0.5
% LAT			(deg N) Latitude
% LON			(deg E) Longitude
% PDENS			(kg/m^3) potential density (-1000)
% RHD			(ppb) un-corrected rhodamine
% SALIN			(PSU) Salinity
% SHIPLOG		(km) ship log (integrated distance traveled)
% TEMP			(deg C) Temperature
%
% ii			index of time series for successive "transects" across GS
% II			index of pressure binned arrays (upper case variables) for successive "transects" across GS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nF_list = [1 2 3 4];	% run all realizations/drifts
nF_list = [2 3];	% just run SI drifts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nF = nF_list	% or a range based on above lists
  if nF==1 		% Intra-thermocline eddy
    start_t = 56.916;		% Tue, Mar 13, first thing a.m. GMT
    end_t = 60.016;		% Thur, Mar 15, end of day GMT
    end_t2 = 60.016;
    %adcp_depth = ; 		% Set depth bin to use for ADCP

    clims_temp = [14 21];
    clims_salin = [34.5 36.5];
    clims_salin = [35.5 36.5];
    clims_pdens = [24.0 27.0];
    clims_press = [0 200];
    clims_fluor_log = [-2 0.1];

    ylims_temp = [14 15.5];
    ylims_salin = [36 36.5];
    ylims_salin = [35.4 35.9];
    ylims_pdens = [25.8 26.8];
    ylims_depth = [0 200];

    Injnum = 1;
  elseif nF==2 		% SI-1
    start_t = 60.0;		% Thur, Mar 1, first thing a.m. GMT
    start_t = 60.7196;		% Thur, Mar 1, first thing a.m. GMT
    end_t = 62.99999;		% Sat, Mar 3, end of day GMT
    end_t = 63.0714;		% Sat, Mar 3, end of day GMT
    end_t2 = 62.8;
    %adcp_depth = 27; 		% Set depth bin to use for ADCP

    clims_temp = [14 21];
    clims_salin = [35.5 36.5];
    clims_pdens = [24.0 27.0];
    clims_press = [0 200];
    clims_fluor_log = [-2 0.1];

    ylims_temp = [18 22];
    ylims_salin = [36 36.5];
    ylims_pdens = [25.5 26];
    ylims_depth = [0 100];

    Injnum = 2;
  elseif nF==3 		% SI-2
    start_t = 64.0;		% Mon, Mar 5, first thing a.m. GMT
    end_t = 68.99999;		% Fri, Mar 9, end of day GMT
    end_t2 = 67.7;
    %adcp_depth = 27; 		% Set depth bin to use for ADCP

    clims_temp = [14 21];
    clims_salin = [35.5 36.5];
    clims_pdens = [24.0 27.0];
    clims_press = [0 200];
    clims_fluor_log = [-2.3 -1];

    ylims_temp = [18 20.75];
    ylims_salin = [35.8 36.5];
    ylims_pdens = [25.5 26];
    ylims_depth = [0 100];

    Injnum = 3;
  elseif nF==4 		% Filament
    start_t = 72.5;		% Tue, Mar 13, first thing a.m. GMT
    end_t = 74.99999;		% Thur, Mar 15, end of day GMT
    end_t2 = 74.3;
    %adcp_depth = 55; 		% Set depth bin to use for ADCP

    clims_temp = [14 21];
    clims_salin = [34.5 36.5];
    clims_salin = [35.5 36.5];
    clims_pdens = [24.0 27.0];
    clims_press = [0 200];

    ylims_temp = [12 18.5];
    ylims_salin = [34.5 36];
    ylims_pdens = [25.8 26.5];
    ylims_depth = [0 100];

    Injnum = 4;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Load the triaxus data for this experiment
  load(['my_triaxus_',exptName{nF}]);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % zero out all negative fluorPPB_ts values
  toss = find(fluorPPB_ts<0);
  fluorPPB_ts(toss) = 0;
  % and all negative FLUORPPB values
  toss = find(FLUORPPB<0);
  FLUORPPB(toss) = 0;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  titstr = {'Triaxus Survey';[datestr(jday_ts(1)+datenum(2012,1,1)),' - ',...
 			datestr(jday_ts(end)+datenum(2012,1,1))]};

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % pick some interesting transect times
  contflag = 1;
  counter = 0;
  if exptName{nF}=='ITE'
    translist = [1:length(ii)];
    translist = [34 41:44 52:54];	% first post-injection is 34, many @100-200 m
    translist = [34:54];	% first post-injection is 34, many @100-200 m
    translist = [34 37 40 41 42 44 45 46 48 50 52 54];  % first post-injection transect is 34, many transects 100-200 m
  elseif exptName{nF}=='SI1'
    translist = [10 43];
    translist = [10 40];
    translist = [10 12 14:22 24:40 43:44];
  elseif exptName{nF}=='SI2'
    translist = [1:length(ii)];
    translist = [5 41];
    translist = [3 5 7:9 11:14 16:18 20:28 30 31 33 34 40:43 45 47 49 51];
  elseif exptName{nF}=='Fil'
    translist = [1:length(ii)];
    translist = [9 22];		% for Klymak et al paper
    translist = [6 9 10 12:15 17 19:27 28:36];
  else
    translist = [1:length(ii)];
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Make some overview plots
  ss = 10;
  sss = 100;
  dyethresh = 0.001;				% (ppb)
  dyethresh = 0.005;				% (ppb)
  keep = find(fluorPPB_ts>dyethresh);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % plot pretty view of fluor data (ppb) as function of depth and time with no dye in gray, dye in color
  figure(nF*100+6)
  clf

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(3,1,1)
  cla
  plot(jday_ts(1:ss:end),depth_ts(1:ss:end),'.','color',0.8*[1 1 1])
  hold on
  %scatter(jday_ts(1:ss:end),depth_ts(1:ss:end),100*fluorPPB_ts(1:ss:end)+1e-10, ...
  %  						fluorPPB_ts(1:ss:end)+1e-10)
  
  if(1)		% Note: here "keep" are the indices where fluor is above a certain threshold (see above)
    thisfluorPPB = log10(fluorPPB_ts(keep(1:ss:end)));
    thisjday = jday_ts(keep(1:ss:end));
    thisdepth = depth_ts(keep(1:ss:end));
    thispdens = pdens_ts(keep(1:ss:end));
    [Y,I] = sort(thisfluorPPB);
    scatter(thisjday(flipud(I)),thisdepth(flipud(I)),50,thisfluorPPB(flipud(I))+1e-10,'filled')
  else
    scatter(jday_ts(keep(1:ss:end)),depth_ts(keep(1:ss:end)),50,log10(fluorPPB_ts(keep(1:ss:end))+1e-10),'filled')
  end

  plot(InjData{Injnum}.jday(InjData{Injnum}.injind(1:ss:end)), ...
  		InjData{Injnum}.press(InjData{Injnum}.injind(1:ss:end)),'r.')

  caxis(clims_fluor_log)
  xlim([InjData{nF}.jday(InjData{nF}.injind(1)),end_t2])
  ylim(ylims_depth)
  set(gca,'ydir','rev')
  %xlabel('Yearday')
  ylabel('Depth (m)')
  title(['Experiment: ',exptName{nF},'Fluorescein (log_{10}(ppb)) vs. time'])

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % similar to previous plot, but this time as function of density rather than depth
  subplot(3,1,2)
  cla
  %plot(jday_ts(1:ss:end),pdens_ts(1:ss:end),'.','color',0.8*[1 1 1])
  plot(jday_ts(1:ss:end),pdens_ts(1:ss:end),'.','color',0.8*[1 1 1])
  hold on
  %scatter(jday_ts(1:ss:end),pdens_ts(1:ss:end),100*fluorPPB_ts(1:ss:end)+1e-10, ...
  %  						fluorPPB_ts(1:ss:end)+1e-10)
  if(1)
    scatter(thisjday(flipud(I)),thispdens(flipud(I)),50,thisfluorPPB(flipud(I))+1e-10,'filled')
  else
    scatter(jday_ts(keep(1:ss:end)),pdens_ts(keep(1:ss:end)),20,log10(fluorPPB_ts(keep(1:ss:end))+1e-10),'filled')
  end

  plot(InjData{Injnum}.jday(InjData{Injnum}.injind(1:ss:end)), ...
  		InjData{Injnum}.sigma1(InjData{Injnum}.injind(1:ss:end)),'r.');

  caxis(clims_fluor_log)
  xlim([InjData{nF}.jday(InjData{nF}.injind(1)),end_t2])
  ylim(ylims_pdens)
  set(gca,'ydir','rev')
  xlabel('Yearday')
  ylabel('Potential Density (kg/m^3)')
  ylabel('\sigma_{\theta} (kg/m^3)')

  clear S10 T10
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure(nF*100+2)
  clf
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % plot a TS plot with all the data, plus injection
  %subplot(3,2,5)
  if nF==1
    tsdiagram([34.7 37.2],[12 22],0,0)
  elseif nF==2
    tsdiagram([34.5 37.0],[13 23],0,0)
  elseif nF==3
    tsdiagram([34.7 37.0],[12.5 22],0,0)
  elseif nF==4
    tsdiagram([34.5 37.0],[12 22],0,0)
  end
  hold on

  if(1)			% do as fancy meshplot
    % use Matlab's sparse matrix formulation to create a 2-d histogram of the T-S data
    % create sparse matrix, using rounded T & S as ordinate values, then convert to full matrix
    ss=1;
    T10 = real(T_ts*10);
    S10 = real(S_ts*20);
    keep = find((S10 > 0.5) & (T10 > 0.5));
    M = full(sparse(round(S10(keep(1:ss:end))), round(T10(keep(1:ss:end))), 1));
    Tordinate = [1:max(round(T10(keep)))]/10;
    Sordinate = [1:max(round(S10(keep)))]/20;

    pcolor(Sordinate,Tordinate,log10(M'));
    colormap(flipud(gray))
    colormap(flipud(hot))
    shading flat;

    caxis([0 5])
  end
  ss=10;
  plot(InjData{Injnum}.sal1(InjData{Injnum}.injind(1:ss:end)), ...
  		InjData{Injnum}.temp1(InjData{Injnum}.injind(1:ss:end)),'g.')

  title(['Experiment: ',exptName{nF}])

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(transectflag)		% TOGGLE THIS PART OF CODE
  % Plot transects of individual GS crossings within an experiment
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(999)
    clf

    for thistrans=translist
      keepII = [II(thistrans,1):II(thistrans,2)];

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subplot(2,2,1)
      cla
      pcolor(SHIPLOG(:,keepII),DEPTH(:,keepII),SALIN(:,keepII));

      xlabel('Ship Log (km)')
      ylabel('Depth (m)')
      title('Salinity (PSU)')
      shading flat
      set(gca,'ydir','rev')

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subplot(2,2,2)
      cla
      pcolor(SHIPLOG(:,keepII),DEPTH(:,keepII),FLUORPPB(:,keepII));

      xlabel('Ship Log (km)')
      ylabel('Depth (m)')
      title('Fluorescein (ppb)')
      shading flat
      set(gca,'ydir','rev')

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subplot(2,2,3)
      cla
      plot(FLUORPPB(:,keepII),DEPTH(:,keepII))
      xlabel('Fluor (ppb)')
      ylabel('Depth (m)')
      set(gca,'ydir','rev')

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subplot(2,2,4)
      cla
      plot(LON(:,:),LAT(:,:),'k-')
      hold on
      plot(LON(:,keepII),LAT(:,keepII),'r-','linewidth',2)

      xlabel('Lon (^oE)')
      ylabel('Lat (^oN)')
      title(['Experiment: ',exptName{nF}])

      set(gca,'DataAspectRatioMode','Manual');
      ax = axis;
      axrat = 1/cos(mean(ax(3:4))*pi/180);
      set(gca,'DataAspectRatio',[1 1/axrat 1]);

      disp('*** PAUSE (hit any key to continue) ***')
      pause

    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end 		% TOGGLE THIS PART OF CODE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end 		% end looping through dye experiments
