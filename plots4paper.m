close all
clear all

%% produces plots for "Causal inference in Environmental sound recognition", James Traer, Sam V. Norman-Haignere & Josh H McDermott
%% -- July, 2020

%% first specify information for each experiment to plot. Store information (experiment_number, subplot, plot_type, data_subset, and data_path) in a structure "exp". Plotter will scroll through elements of "exp".
%% -- plot_type options are: 'level','level-louder','reverb','reverb-louder','distance','naturalism','reverb-In-vs-Out','level-with-noise'
%% -- data_subset options are: 'difficulty_matched', 'distance_matched', 'excitation_matched'

e=0;
%%% == Experiment 1: Source recognition against level ==
%%% --- all data
%e=e+1; exp(e).exp=1; exp(e).subplot='2B'; exp(e).plot_type='level'; exp(e).data_subset='all_data'; exp(e).data_path='exp01_fig2B_2C_2F_3D_6D.csv';
%%% --- difficulty-matched subset of sources
%e=e+1; exp(e).exp=1; exp(e).subplot='2C'; exp(e).plot_type='level'; exp(e).data_subset='difficulty_matched'; exp(e).data_path='exp01_fig2B_2C_2F_3D_6D.csv';
%%% --- distance-matched subset of sources
%e=e+1; exp(e).exp=1; exp(e).subplot='3D'; exp(e).plot_type='level'; exp(e).data_subset='distance_matched'; exp(e).data_path='exp01_fig2B_2C_2F_3D_6D.csv';
%%% --- excitation-matched subset of sources
%e=e+1; exp(e).exp=1; exp(e).subplot='6D'; exp(e).plot_type='level'; exp(e).data_subset='excitation_matched'; exp(e).data_path='exp01_fig2B_2C_2F_3D_6D.csv';
%
%%% == Experiment 2: Source recognition against level (multiple choice) ==
%e=e+1; exp(e).exp=2; exp(e).subplot='2D'; exp(e).plot_type='level'; exp(e).data_subset='all_data'; exp(e).data_path='exp02_fig2D.csv';
%
%%% == Experiment 3: Fraction of Exp 1 results judged louder ==
%e=e+1; exp(e).exp=3; exp(e).subplot='2F'; exp(e).plot_type='level-louder'; exp(e).data_subset='all_data'; exp(e).data_path='exp01_fig2B_2C_2F_3D_6D.csv';
%
%%% == Experiment 4: Source recognition in reverberation ==
%%% --- all data
%e=e+1; exp(e).exp=4; exp(e).subplot='2J'; exp(e).plot_type='reverb'; exp(e).data_subset='all_data'; exp(e).data_path='exp04_fig2J_2K_2L_3E_5D_6E.csv';
%%% --- difficulty-matched subset of sources
%e=e+1; exp(e).exp=4; exp(e).subplot='2K'; exp(e).plot_type='reverb'; exp(e).data_subset='difficulty_matched'; exp(e).data_path='exp04_fig2J_2K_2L_3E_5D_6E.csv';
%%% --- distance-matched subset of sources
%e=e+1; exp(e).exp=4; exp(e).subplot='3E'; exp(e).plot_type='reverb'; exp(e).data_subset='distance_matched'; exp(e).data_path='exp04_fig2J_2K_2L_3E_5D_6E.csv';
%%% --- excitation-matched subset of sources
%e=e+1; exp(e).exp=4; exp(e).subplot='6E'; exp(e).plot_type='reverb'; exp(e).data_subset='excitation_matched'; exp(e).data_path='exp04_fig2J_2K_2L_3E_5D_6E.csv';
%%% --- Comparison of indoor vs outdoor sounds: with distance-matched subset of sources
%e=e+1; exp(e).exp=4; exp(e).subplot='5D'; exp(e).plot_type='reverb-In-vs-Out'; exp(e).data_subset='distance_matched'; exp(e).data_path='exp04_fig2J_2K_2L_3E_5D_6E.csv';
%
%%%% == Experiment 5: Source loudness judgments ==
%e=e+1; exp(e).exp=5; exp(e).subplot='2L'; exp(e).plot_type='reverb-louder'; exp(e).data_subset='all_data'; exp(e).data_path='exp04_fig2J_2K_2L_3E_5D_6E.csv';
%
%%% == Experiments 6+7: Source distance judgments ==
%e=e+1; exp(e).exp=6; exp(e).subplot='3B'; exp(e).plot_type='distance'; exp(e).data_subset='all_data'; exp(e).data_path='exp06_fig3B.csv';
%e=e+1; exp(e).exp=7; exp(e).subplot='4C'; exp(e).plot_type='distance'; exp(e).data_subset='all_data'; exp(e).data_path='exp07_fig4C.csv';
%
%%%% == Experiment 8: Source recognition (studio-recordings) ==
%%% --- all data
%e=e+1; exp(e).exp=8; exp(e).subplot='4E'; exp(e).plot_type='reverb'; exp(e).data_subset='all_data'; exp(e).data_path='exp08_fig4E_4F_5D.csv';
%%% --- distance matched
%e=e+1; exp(e).exp=8; exp(e).subplot='4E'; exp(e).plot_type='reverb'; exp(e).data_subset='distance_matched'; exp(e).data_path='exp08_fig4E_4F_5D.csv';
%%% --- distance and difficulty matched
%e=e+1; exp(e).exp=8; exp(e).subplot='4E'; exp(e).plot_type='reverb'; exp(e).data_subset='distance_and_difficulty_matched'; exp(e).data_path='exp08_fig4E_4F_5D.csv';
%%% --- Comparison of indoor vs outdoor sounds: with distance-matched subset of sources
%e=e+1; exp(e).exp=8; exp(e).subplot='5D'; exp(e).plot_type='reverb-In-vs-Out'; exp(e).data_subset='distance_matched'; exp(e).data_path='exp08_fig4E_4F_5D.csv';
%
%%%% == Experiment 9: Source loudness judgments (distance-matched studio recordings) ==
%e=e+1; exp(e).exp=8; exp(e).subplot='4F'; exp(e).plot_type='reverb-louder'; exp(e).data_subset='distance_matched'; exp(e).data_path='exp08_fig4E_4F_5D.csv';
%
%%% == Experiments 10+11: Source naturalism judgments ==
%e=e+1; exp(e).exp=10; exp(e).subplot='5C'; exp(e).plot_type='naturalism'; exp(e).data_subset='distance_matched'; exp(e).data_path='exp10_fig5C.csv';
%e=e+1; exp(e).exp=11; exp(e).subplot='5C'; exp(e).plot_type='naturalism'; exp(e).data_subset='distance_matched'; exp(e).data_path='exp11_fig5C.csv';

%% == Experiment 12: Source recognition against level and noise ==
e=e+1; exp(e).exp=12; exp(e).subplot='7B'; exp(e).plot_type='level-with-noise'; exp(e).data_subset='all_data'; exp(e).data_path='exp12_fig7B.csv';

%% scroll through experiments and plot
for jj=1:length(exp);
  %% load data
  fprintf('Loading data for Exp %d: subplot %s (%s)\n',exp(jj).exp,exp(jj).subplot,exp(jj).data_subset);
  D=readtable(exp(jj).data_path);

  %% extract specified data subsets
  if strcmp(exp(jj).data_subset,'difficulty_matched')
    index=find(D.difficulty_matched_subset==1);
  elseif strcmp(exp(jj).data_subset,'distance_matched')
    if strcmp(exp(jj).plot_type,'reverb-In-vs-Out')||strcmp(exp(jj).plot_type,'naturalism')
      index=find(rem(D.distance_matched_subset,1)==0.5);
    else
      index=find(D.distance_matched_subset>=1);
    end
  elseif strcmp(exp(jj).data_subset,'distance_and_difficulty_matched')
    index=find(D.distance_matched_subset>=1);
    index2=find(D.difficulty_matched_subset==1);
    index=intersect(index,index2);
  elseif strcmp(exp(jj).data_subset,'excitation_matched')
    index=find(D.excitation_matched_subset==1);
  else
    index=1:size(D,1);
  end
  D=D(index,:);
  Dsz=size(D);
  N=max([D.listener]);

  %% plot data against presentation level (intensity) (if appropriate): Figs 2B,2C,2D,2F,3D,6D
  if strcmp(exp(jj).plot_type,'level')||strcmp(exp(jj).plot_type,'level-louder') 
    level=[1:7];
    %% pre-allocate data for High- (H) and Low-intensity (L) sources 
    H=zeros(N,length(level)); Hcnt=H; 
    L=zeros(N,length(level)); Lcnt=L;
    % scroll through data and group into conditions
    last_participant_index=0;
    last_trial=0;
    for jv=1:Dsz(1);
      level_index=find(level==D(jv,:).level_condition);
      participant_index=D(jv,:).listener;
      if participant_index~=last_participant_index||D(jv,:).listener_trial~=last_trial
        if strcmp(D(jv,:).typical_source_intensity{1},'high-intensity') || (strcmp(exp(jj).data_subset,'distance_matched')&&strcmp(D(jv,:).typical_source_intensity{1},'3'));
          if strcmp(exp(jj).plot_type,'level')
            H(participant_index,level_index)=H(participant_index,level_index)+D(jv,:).proportion_correct;
          elseif strcmp(exp(jj).plot_type,'level-louder')
            H(participant_index,level_index)=H(participant_index,level_index)+D(jv,:).proportion_judged_louder;
          end
          Hcnt(participant_index,level_index)=Hcnt(participant_index,level_index)+1;
        elseif strcmp(D(jv,:).typical_source_intensity{1},'low-intensity') || (strcmp(exp(jj).data_subset,'distance_matched')&&strcmp(D(jv,:).typical_source_intensity{1},'2'));
          if strcmp(exp(jj).plot_type,'level')
            L(participant_index,level_index)=L(participant_index,level_index)+D(jv,:).proportion_correct;
          elseif strcmp(exp(jj).plot_type,'level-louder')
            L(participant_index,level_index)=L(participant_index,level_index)+D(jv,:).proportion_judged_louder;
          end
          Lcnt(participant_index,level_index)=Lcnt(participant_index,level_index)+1;
        end
      end
      last_participant_index=participant_index;
      last_trial=D(jv,:).listener_trial;
    end
    H=H./Hcnt;
    L=L./Lcnt;
    % plot
    figure;
    if max(isnan([H(:); L(:)]))==0;
      errorbar(level,mean(H,1),std(H,1)/sqrt(N),'r'); hold on;
      errorbar(level+0.2,mean(L,1),std(L,1)/sqrt(N),'b'); hold on;
    else
      fprintf('WARNING: not all particpants heard all conditions\n')
      errorbar(level,nanmean(H,1),nanstd(H,1)/sqrt(N),'r'); hold on;
      errorbar(level+0.2,nanmean(L,1),nanstd(L,1)/sqrt(N),'b'); hold on;
    end
    set(gca,'ylim',[0.3 1])
    text(35,0.35,1,sprintf('N=%d',N))
    title(sprintf('figure %s: Experiment %d (%s)',exp(jj).subplot,exp(jj).exp,exp(jj).data_subset));
    set(gca,'xtick',level,'xticklabel',{'30';'40';'50';'60';'70';'80';'90'});
    xlabel('Presentation Intensity (dB)');
    if strcmp(exp(jj).plot_type,'level')
      ylabel('Proportion correct');
    elseif strcmp(exp(jj).plot_type,'level-louder')
      ylabel('Proportion judged louder');
    end
    hold off; drawnow
    if strcmp(exp(jj).data_subset,'distance_matched');
      %% block levels into groups of 3 to ensure we don't have empty sets
      H2(:,1)=nanmean(H(:,1:3),2);
      H2(:,2)=nanmean(H(:,3:5),2);
      H2(:,3)=nanmean(H(:,5:7),2);
      L2(:,1)=nanmean(L(:,1:3),2);
      L2(:,2)=nanmean(L(:,3:5),2);
      L2(:,3)=nanmean(L(:,5:7),2);
      %% write data for stats
      save_stats(3,H2,L2,exp(jj))
    else
      %% write data for stats
      save_stats(7,H,L,exp(jj))
    end

  %% == plot data against reverb condition (if appropriate): Figs 2J,2K,3E,4E,4F,6E ==
  elseif strcmp(exp(jj).plot_type,'reverb')||strcmp(exp(jj).plot_type,'reverb-louder')
    rvrb={'unaltered', 'reverb'};
    H=zeros(N,length(rvrb)); Hcnt=H; 
    L=zeros(N,length(rvrb)); Lcnt=L;
    last_sbj=0;
    last_trial=0;
    for jv=1:Dsz(1);
      rvrb_ndx=find(strcmp(rvrb,D(jv,:).reverb_condition));
      sbj=D(jv,:).listener;
      if sbj~=last_sbj||D(jv,:).listener_trial~=last_trial
        if strcmp(D(jv,:).typical_source_intensity,'high-intensity');
          if strcmp(exp(jj).plot_type,'reverb')
            H(sbj,rvrb_ndx)=H(sbj,rvrb_ndx)+D(jv,:).proportion_correct;
          elseif strcmp(exp(jj).plot_type,'reverb-louder')
            H(sbj,rvrb_ndx)=H(sbj,rvrb_ndx)+D(jv,:).proportion_judged_louder;
          end
          Hcnt(sbj,rvrb_ndx)=Hcnt(sbj,rvrb_ndx)+1;
        elseif strcmp(D(jv,:).typical_source_intensity,'low-intensity');
          if strcmp(exp(jj).plot_type,'reverb')
            L(sbj,rvrb_ndx)=L(sbj,rvrb_ndx)+D(jv,:).proportion_correct;
          elseif strcmp(exp(jj).plot_type,'reverb-louder')
            L(sbj,rvrb_ndx)=L(sbj,rvrb_ndx)+D(jv,:).proportion_judged_louder;
          end
          Lcnt(sbj,rvrb_ndx)=Lcnt(sbj,rvrb_ndx)+1;
        end
        last_sbj=sbj;
        last_trial=D(jv,:).listener_trial;
      end
    end
    H=H./Hcnt;
    L=L./Lcnt;
    figure;
    errorbar([0 1],mean(H,1),std(H,1)/sqrt(N),'r'); hold on;
    errorbar([0 1]+0.02,mean(L,1),std(L,1)/sqrt(N),'b'); hold on;
    set(gca,'ylim',[0.0 1])
    text(0.15,0.35,1,sprintf('N=%d',N))
    title(sprintf('figure %s: Experiment %d (%s)',exp(jj).subplot,exp(jj).exp,exp(jj).data_subset));
    set(gca,'xtick',[0 1],'xticklabel',{'unaltered';'Added reverberation'});
    if strcmp(exp(jj).plot_type,'reverb')
      ylabel('Proportion correct');
    elseif strcmp(exp(jj).plot_type,'reverb-louder')
      ylabel('Proportion judged louder');
    end
    hold off; drawnow
    %% write data for stats
    save_stats(2,H,L,exp(jj)); 

  %% == plot distance ratings: fig 3B, 4C ==
  elseif strcmp(exp(jj).plot_type,'distance')
    rvrb={'unaltered', 'reverb'};
    H=zeros(N,length(rvrb)); Hcnt=H; 
    L=zeros(N,length(rvrb)); Lcnt=L;
    last_sbj=0;
    last_trial=0;
    for jv=1:Dsz(1);
      rvrb_ndx=find(strcmp(rvrb,D(jv,:).reverb_condition));
      sbj=D(jv,:).listener;
      if sbj~=last_sbj||D(jv,:).listener_trial~=last_trial
        if strcmp(D(jv,:).typical_source_intensity,'high-intensity');
          H(sbj,rvrb_ndx)=H(sbj,rvrb_ndx)+D(jv,:).distance_rating;
          Hcnt(sbj,rvrb_ndx)=Hcnt(sbj,rvrb_ndx)+1;
        elseif strcmp(D(jv,:).typical_source_intensity,'low-intensity');
          L(sbj,rvrb_ndx)=L(sbj,rvrb_ndx)+D(jv,:).distance_rating;
          Lcnt(sbj,rvrb_ndx)=Lcnt(sbj,rvrb_ndx)+1;
        end
      end
      last_sbj=sbj;
      last_trial=D(jv,:).listener_trial;
    end
    H=H./Hcnt;
    L=L./Lcnt;
    figure;
    errorbar([0 1],mean(H,1),std(H,1)/sqrt(N),'r'); hold on;
    errorbar([0 1]+0.02,mean(L,1),std(L,1)/sqrt(N),'b'); hold on;
    set(gca,'ylim',[1 7],'ytick',[1:7],'yticklabel',{'10cm';'30cm';'1m';'3m';'10m';'30m';'100m'});
    text(0.15,1.35,1,sprintf('N=%d',N));
    title(sprintf('figure %s: Experiment %d (%s)',exp(jj).subplot,exp(jj).exp,exp(jj).data_subset));
    set(gca,'xtick',[0 1],'xticklabel',{'unaltered';'Added reverberation'});
    ylabel('Distance ratings');
    hold off; drawnow
    %% write data for stats
    save_stats(2,H,L,exp(jj))

  %% == plot naturalism ratings: fig 5C ==
  elseif strcmp(exp(jj).plot_type,'naturalism')
    rvrb={'unaltered', 'reverb'};
    HI=zeros(N,length(rvrb)); HIcnt=HI;
    HO=zeros(N,length(rvrb)); HOcnt=HO;
    LI=zeros(N,length(rvrb)); LIcnt=LI;
    LO=zeros(N,length(rvrb)); LOcnt=LO;
    last_sbj=0;
    last_trial=0;
    for jv=1:Dsz(1);
      rvrb_ndx=find(strcmp(rvrb,D(jv,:).reverb_condition));
      sbj=D(jv,:).listener;
      if sbj~=last_sbj||D(jv,:).listener_trial~=last_trial
        if strcmp(D(jv,:).typical_source_intensity{1}(1:5),'high-')&&strcmp(D(jv,:).typical_source_location,'indoor');
          HI(sbj,rvrb_ndx)=HI(sbj,rvrb_ndx)+D(jv,:).natural_rating;
          HIcnt(sbj,rvrb_ndx)=HIcnt(sbj,rvrb_ndx)+1;
        elseif strcmp(D(jv,:).typical_source_intensity{1}(1:5),'high-')&&strcmp(D(jv,:).typical_source_location,'outdoor');
          HO(sbj,rvrb_ndx)=HO(sbj,rvrb_ndx)+D(jv,:).natural_rating;
          HOcnt(sbj,rvrb_ndx)=HOcnt(sbj,rvrb_ndx)+1;
        elseif strcmp(D(jv,:).typical_source_intensity{1}(1:4),'low-')&&strcmp(D(jv,:).typical_source_location,'indoor');
          LI(sbj,rvrb_ndx)=LI(sbj,rvrb_ndx)+D(jv,:).natural_rating;
          LIcnt(sbj,rvrb_ndx)=LIcnt(sbj,rvrb_ndx)+1;
        elseif strcmp(D(jv,:).typical_source_intensity{1}(1:4),'low-')&&strcmp(D(jv,:).typical_source_location,'outdoor');
          LO(sbj,rvrb_ndx)=LO(sbj,rvrb_ndx)+D(jv,:).natural_rating;
          LOcnt(sbj,rvrb_ndx)=LOcnt(sbj,rvrb_ndx)+1;
        end
      end
      last_sbj=sbj;
      last_trial=D(jv,:).listener_trial;
    end
    HI=HI./HIcnt;
    HO=HO./HOcnt;
    LI=LI./LIcnt;
    LO=LO./LOcnt;
    figure;
    errorbar([0 1],nanmean(HO,1),nanstd(HO,1)/sqrt(N),'g'); hold on;
    errorbar([0 1]+0.02,nanmean(HI,1),nanstd(HI,1)/sqrt(N),'m'); hold on;
    errorbar([2 3],nanmean(LO,1),nanstd(LO,1)/sqrt(N),'g'); hold on;
    errorbar([2 3]+0.02,nanmean(LI,1),nanstd(LI,1)/sqrt(N),'m'); hold on;
    set(gca,'ylim',[1 5])
    text(0.15,1.35,1,sprintf('N=%d',N))
    title(sprintf('figure %s: Experiment %d',exp(jj).subplot,exp(jj).exp));
    set(gca,'xtick',[0 1 2 3],'xticklabel',{'unaltered';'Added reverberation';'unaltered';'Added reverberation'});
    ylabel('Naturalism ratings (1-5)');
    hold off; drawnow
    %% write data for stats
    save_stats(2,HI,HO,exp(jj),'_High-Intensity')
    save_stats(2,LI,LO,exp(jj),'_Low-Intensity')

  %% == plot data against reverb condition for indoor-vs-outdoor sounds (if appropriate): Fig 5D ==
  elseif strcmp(exp(jj).plot_type,'reverb-In-vs-Out')
    rvrb={'unaltered', 'reverb'};
    HI=zeros(N,length(rvrb)); HIcnt=HI;
    HO=zeros(N,length(rvrb)); HOcnt=HO;
    LI=zeros(N,length(rvrb)); LIcnt=LI;
    LO=zeros(N,length(rvrb)); LOcnt=LO;
    last_sbj=0;
    last_trial=0;
    for jv=1:Dsz(1);
      rvrb_ndx=find(strcmp(rvrb,D(jv,:).reverb_condition));
      sbj=D(jv,:).listener;
      if sbj~=last_sbj||D(jv,:).listener_trial~=last_trial
        if strcmp(D(jv,:).typical_source_intensity{1}(1:5),'high-')&&strcmp(D(jv,:).typical_source_location,'indoor');
          HI(sbj,rvrb_ndx)=HI(sbj,rvrb_ndx)+D(jv,:).proportion_correct;
          HIcnt(sbj,rvrb_ndx)=HIcnt(sbj,rvrb_ndx)+1;
        elseif strcmp(D(jv,:).typical_source_intensity{1}(1:5),'high-')&&strcmp(D(jv,:).typical_source_location,'outdoor');
          HO(sbj,rvrb_ndx)=HO(sbj,rvrb_ndx)+D(jv,:).proportion_correct;
          HOcnt(sbj,rvrb_ndx)=HOcnt(sbj,rvrb_ndx)+1;
        elseif strcmp(D(jv,:).typical_source_intensity{1}(1:4),'low-')&&strcmp(D(jv,:).typical_source_location,'indoor');
          LI(sbj,rvrb_ndx)=LI(sbj,rvrb_ndx)+D(jv,:).proportion_correct;
          LIcnt(sbj,rvrb_ndx)=LIcnt(sbj,rvrb_ndx)+1;
        elseif strcmp(D(jv,:).typical_source_intensity{1}(1:4),'low-')&&strcmp(D(jv,:).typical_source_location,'outdoor');
          LO(sbj,rvrb_ndx)=LO(sbj,rvrb_ndx)+D(jv,:).proportion_correct;
          LOcnt(sbj,rvrb_ndx)=LOcnt(sbj,rvrb_ndx)+1;
        end
      end
      last_sbj=sbj;
      last_trial=D(jv,:).listener_trial;
    end
    HI=HI./HIcnt;
    HO=HO./HOcnt;
    LI=LI./LIcnt;
    LO=LO./LOcnt;
    figure;
    errorbar([0 1],nanmean(HO,1),nanstd(HO,1)/sqrt(N),'g'); hold on;
    errorbar([0 1]+0.02,nanmean(HI,1),nanstd(HI,1)/sqrt(N),'m'); hold on;
    errorbar([2 3],nanmean(LO,1),nanstd(LO,1)/sqrt(N),'g'); hold on;
    errorbar([2 3]+0.02,nanmean(LI,1),nanstd(LI,1)/sqrt(N),'m'); hold on;
    set(gca,'ylim',[0.0 1])
    text(0.15,0.35,1,sprintf('N=%d',N))
    title(sprintf('figure %s: Experiment %d (%s)',exp(jj).subplot,exp(jj).exp,exp(jj).data_subset));
    set(gca,'xtick',[0 1 2 3],'xticklabel',{'unaltered';'Added reverberation';'unaltered';'Added reverberation'});
    ylabel('Proportion correct');
    hold off; drawnow
    %% write data for stats
    save_stats(2,HI,HO,exp(jj),'_High-Intensity')
    save_stats(2,LI,LO,exp(jj),'_Low-Intensity')

  elseif strcmp(exp(jj).plot_type,'level-with-noise')
    level=[2:7];
    %% pre-allocate data for High- (H) and Low-intensity (L) sources 
    H=zeros(N,length(level)); Hcnt=H; 
    L=zeros(N,length(level)); Lcnt=L;
    Hn=zeros(N,length(level)); Hncnt=H; 
    Ln=zeros(N,length(level)); Lncnt=L;
    noise_index=find(D.noise_condition==1);
    no_noise_index=find(D.noise_condition==0);
    % scroll through data and group into conditions
    last_participant_index=0;
    last_trial=0;
    for jv=no_noise_index.';
      level_index=find(level==D(jv,:).level_condition);
      participant_index=D(jv,:).listener;
      if participant_index~=last_participant_index||D(jv,:).listener_trial~=last_trial
        if strcmp(D(jv,:).typical_source_intensity{1},'high-intensity')
          H(participant_index,level_index)=H(participant_index,level_index)+D(jv,:).proportion_correct;
          Hcnt(participant_index,level_index)=Hcnt(participant_index,level_index)+1;
        elseif strcmp(D(jv,:).typical_source_intensity{1},'low-intensity') 
          L(participant_index,level_index)=L(participant_index,level_index)+D(jv,:).proportion_correct;
          Lcnt(participant_index,level_index)=Lcnt(participant_index,level_index)+1;
        end
      end
      last_participant_index=participant_index;
      last_trial=D(jv,:).listener_trial;
    end
    H=H./Hcnt;
    L=L./Lcnt;
    % scroll through trials with noise
    last_participant_index=0;
    last_trial=0;
    for jv=noise_index.';
      level_index=find(level==D(jv,:).level_condition);
      participant_index=D(jv,:).listener;
      if participant_index~=last_participant_index||D(jv,:).listener_trial~=last_trial
        if strcmp(D(jv,:).typical_source_intensity{1},'high-intensity')
          Hn(participant_index,level_index)=Hn(participant_index,level_index)+D(jv,:).proportion_correct;
          Hncnt(participant_index,level_index)=Hncnt(participant_index,level_index)+1;
        elseif strcmp(D(jv,:).typical_source_intensity{1},'low-intensity') 
          Ln(participant_index,level_index)=Ln(participant_index,level_index)+D(jv,:).proportion_correct;
          Lncnt(participant_index,level_index)=Lncnt(participant_index,level_index)+1;
        end
      end
      last_participant_index=participant_index;
      last_trial=D(jv,:).listener_trial;
    end
    Hn=Hn./Hncnt;
    Ln=Ln./Lncnt;
    % plot
    figure;
    errorbar(level,nanmean(H,1),nanstd(H,1)/sqrt(N),'r'); hold on;
    errorbar(level+0.2,nanmean(L,1),nanstd(L,1)/sqrt(N),'b'); hold on;
    errorbar(level,nanmean(Hn,1),nanstd(Hn,1)/sqrt(N),'r--'); hold on;
    errorbar(level+0.2,nanmean(Ln,1),nanstd(Ln,1)/sqrt(N),'b--'); hold on;
    set(gca,'ylim',[0.3 1])
    text(35,0.35,1,sprintf('N=%d',N))
    title(sprintf('figure %s: Experiment %d (%s)',exp(jj).subplot,exp(jj).exp,exp(jj).data_subset));
    set(gca,'xtick',level,'xticklabel',{'30';'40';'50';'60';'70';'80';'90'});
    xlabel('Presentation Intensity (dB)');
    if strcmp(exp(jj).plot_type,'level')
      ylabel('Proportion correct');
    elseif strcmp(exp(jj).plot_type,'level-louder')
      ylabel('Proportion judged louder');
    end
    hold off; drawnow
    %% write data for stats
    save_stats(6,H,L,exp(jj))

  end
end

function save_stats(No_presentation_conditions,source1,source2,exp_info,modifier)
if nargin<5;
modifier='';
end

fid=fopen('tmp.csv','w');
fprintf(fid,'source_condition,presentation_condition,fraction_correct,subject\n')
for jSrc=[1:2];
for jPrs=[1:No_presentation_conditions];
  for jSbj=1:size(source1,1);
    if jSrc==1;
      fprintf(fid,'%d,%d,%0.8f,%d\n',jSrc,jPrs,source1(jSbj,jPrs),jSbj);
    else
      fprintf(fid,'%d,%d,%0.8f,%d\n',jSrc,jPrs,source2(jSbj,jPrs),jSbj);
    end
  end
end
end
fclose(fid);
if ~isempty(modifier)
fprintf('===== %s =====\n',modifier)
end
unix('python3 stts.py');
%unix(sprintf('mv tmpStts.csv exp%02d_fig%s_stats%s.csv',exp_info.exp,exp_info.plots,modifier))
end %EOF
