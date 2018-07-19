function [] = adi_avg_MEG(path, dataName, condition)

  load ([path, dataName])
  cfg = [];  
  cfg.vartrllength = 2;      
  switch condition
      case 'dislike'
          tdata = ft_timelockanalysis(cfg, dislike_allRuns)
      case 'dontcare'
          tdata = ft_timelockanalysis(cfg, dontcare_allRuns)
      case 'like'
          tdata = ft_timelockanalysis(cfg, like_allRuns)
  end
  
  figure
  plot(tdata.time, mean(abs(tdata.avg(1:248,:))))
  axis tight
  xlabel ('time')
  ylabel ('mean absolute activation')
  title(dataName(1:end-4))
  savefig([path, 'avg_abs_', dataName(1:end-4),'.fig'])
  fig = ([path, 'avg_abs_', dataName(1:end-4)]);
  print('-dpng', fig); 
  close all

end