function zeta = get_zeta(dirname,n)

    addpath("~/suntans/suntans-master/mfiles/");

    suntansfile = [dirname,'/suntans.dat'];

    depthdata = load([dirname,'/depth.dat-voro'],'-ascii');

    dv = depthdata(:,3);

    Nkmax = getvalue(suntansfile,'Nkmax');
    Nc = length(dv);

    fid = fopen([dirname,'/vert_dzz.dat'],'rb');
    fseek(fid,8*Nkmax*Nc*(n-1),'bof');
    dzz = reshape(fread(fid,Nkmax*Nc,'float64'),Nc,Nkmax);
    fclose(fid);
    
    fid = fopen([dirname,'/fs.dat'],'rb');
    fseek(fid,8*Nc*(n-1),'bof');
    h = fread(fid,Nc,'float64');
    fclose(fid);
    
    zeta = zeros(Nc,Nkmax+1);
    zeta(:,1) = h;
    %zeta(:, 1) = h*0;
    for k=2:Nkmax+1
        zeta(:,k) = zeta(:,k-1)-dzz(:,k-1);
    end
    
    


    
