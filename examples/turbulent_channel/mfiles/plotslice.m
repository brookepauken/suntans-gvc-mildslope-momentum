function [x,z,d,data] = plotslice(PLOT,datadir,n)

  EMPTY=999999;         % Empty cells are defined by this
  precision='float64';  % Precision for reading in data
  
  % cellcentered data contains the voronoi points and the depths
  % at those points.
  dz = load([datadir,'/vertspace.dat']);
  dv_data = load([datadir,'/depth.dat-voro']);
  xv = dv_data(:,1);
  yv = dv_data(:,2);
  dv = dv_data(:,3);
  
  % Total number of cells in the horizontal Nc and vertical Nk
  Nc = length(xv);
  Nk = length(dz);
  
  % Length and depth of domain
  L = max(xv);
  dmax = max(dv);

  % Set up the Cartesian grid
  z = -cumsum(dz)/sum(dz);
  [xs,is]=sort(xv);
  [x,z]=meshgrid(xs,z);
  zeta = get_zeta(datadir,n);
  d = -zeta(:,end);
  
  z = 0.5*(zeta(:,1:end-1)+zeta(:,2:end))';

  % Read free surface to create sigma grid
  %  hdata = fread(fopen([datadir,'/fs.dat'],'rb'),'float64');
  
  % Open up file descriptors for binary files
  switch(PLOT)
   case 'q'
    file = [datadir,'/q.dat'];
   case 's'
    file = [datadir,'/s.dat'];    
   case 'T'
    file = [datadir,'/T.dat'];    
   case {'u','v','w','U'}
    file = [datadir,'/u.dat'];    
   case 's0'
    file = [datadir,'/s0.dat'];    
   case 'h'
    file = [datadir,'/fs.dat'];    
   case 'nut'
    file = [datadir,'/nut.dat'];    
   case 'kappat'
    file = [datadir,'/kappat.dat'];    
   otherwise
    fprintf('Unrecognized plot variable.\n');
    fprintf('Use one of ''q'',''s'',''u'',''w'',''s0'',''h'',''nut'',''kappat'', ''U''.\n');
    data=zeros(Nk,Nc);
    return;
  end

  fid = fopen(file,'rb');

  switch(PLOT)
   case {'u','v','w','U'}
    data = reshape(getcdata(fid,Nk*3*Nc,n,precision),Nc,Nk,3);
    if(PLOT=='u')
      data = squeeze(data(:,:,1));
    elseif(PLOT=='v')
      data = squeeze(data(:,:,2));
    elseif(PLOT=='w')
      data = squeeze(data(:,:,3));
    end
    if(PLOT~='U')
        data = data(is,:)';
    else
        data = data(is,:,:);
    end
    data(find(data==EMPTY))=nan;
   case 'h'
    data = getcdata(fid,Nc,n,precision);
    data = data(is);
   otherwise
    data = reshape(getcdata(fid,Nc*Nk,n,precision),Nc,Nk);
    data = data(is,:)';
    data(find(data==EMPTY))=nan;
  end

  fclose(fid);
