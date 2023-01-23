function plotter_nd(X,Y,Z,V,name,ratio,vari_size,caxismin,caxismax,percentile_min,percentile_max,symflag,cut)

if isempty(V), delete(gca), return, end
if numel(V) ~= prod(vari_size), delete(gca), return, end
if prod(vari_size)==1, delete(gca), return, end
V = reshape(V,[vari_size 1]);

% assign default value to percentile_min if not provided
if nargin < 10  || isempty(percentile_min)
  percentile_min = 2;
end
% assign default value to percentile_max if not provided
if nargin < 11 || isempty(percentile_max)
  percentile_max = 98;
end
% assign percentile_min value to caxismin if not provided
if nargin < 8 || isempty(caxismin)
  caxismin  = my_prctile(V(:),percentile_min);
end
% assign percentile_max value to caxismax if not provided
if nargin < 9 || isempty(caxismax)
  caxismax  = my_prctile(V(:),percentile_max);
end
% exchange caxismin and caxismax if not a valid interval
if caxismax < caxismin
  caxismin2 = caxismin;
  caxismin  = caxismax;
  caxismax  = caxismin2;
end
% apply unit interval if interval is too small
if caxismax-eps < caxismin
  caxismin = mean([caxismin,caxismax])-0.5; caxismax = caxismin + 1;
end
% force limits symmetric about zero if symflag is set
if nargin > 11 && symflag==1
  caxismin = mean([caxismin,-caxismax]);
  caxismax = -caxismin;
end

% check cutting parameters
if nargin < 13 || isempty(cut)
  cut = [0 0 0];
end
n1           = 1:length(X);
n2           = 1:length(Y);
n3           = 1:length(Z);
if cut(1) == true
  n1         = floor(length(X)/2):length(X);
end
if cut(2) == true
  n2         = floor(length(Y)/2):length(Y);
end
if cut(3) == true
  n3         = 1:floor(length(Z)/2);
end

if length(vari_size)==1
  plot(X,V)
  ylim([caxismin caxismax])
  title(name)
  xlabel('x [m]')
elseif length(vari_size)==2
  pcolor(X,Y,V)
  shading interp
  caxis([caxismin caxismax])
  colorbar
  daspect(ratio)
  set(gcf,'renderer','zbuffer')
  xlabel('x [m]')
  ylabel('y [m]')
  title(name)
elseif length(vari_size)==3
  caxisrng   = abs(caxismax-caxismin);
  levels     = caxismin:caxisrng*0.1:caxismax;
  colors     = round(max(min(1+(levels-caxismin)/caxisrng*63,64),1));

  % handling caxis in case of constant values to be plotted
  if var(V(:)) > eps*eps
    s=slice(X,Y,Z,permute(V,[2 1 3]),X(end),Y(end),Z(1));
    shading interp
    set(s,'facealpha',0.5)
  end
 
  % plotting isosurfaces of slices depending on noisiness of plots
  noise      = sum(sum(sum(abs(del2(min(max(V,caxismin),caxismax))))))/numel(V)/caxisrng;
  if noise < 0.3 && var(V(:)) > eps*eps
    contourslice(X(n1),Y(n2),Z(n3),permute(V(n1,n2,n3),[2 1 3]),[X(1) X(end)],[Y(1) Y(end)],[Z(1) Z(end)],levels)
    map=colormap('jet');
    for ii=1:length(levels)
      p=patch(isosurface(X(n1),Y(n2),Z(n3),permute(V(n1,n2,n3),[2 1 3]),levels(ii)));
      isonormals(X(n1),Y(n2),Z(n3),permute(V(n1,n2,n3),[2 1 3]),p)
      set(p,'facecolor',map(colors(ii),:),'edgecolor','none','facealpha',0.5);
    end
  elseif var(V(:)) > eps*eps
    slice(X,Y,Z,permute(V,[2 1 3]),[X(end)/2 X(end)],[Y(end)/2 Y(end)],[0 Z(end)/2])
    hold on
    for i=0:0.5:1
      for j=0:0.5:1
        for k=0:0.5:1
          plot3([0 1]*X(end),[j j]*Y(end),[k k]*Z(end),'k')
          plot3([i i]*X(end),[0 1]*Y(end),[k k]*Z(end),'k')
          plot3([i i]*X(end),[j j]*Y(end),[0 1]*Z(end),'k')
        end
      end
    end
    hold off
    shading interp
  end

  view(-30,25)
  camproj perspective
  daspect(ratio);
  xlims = (X([1 end]));
  ylims = (Y([1 end]));
  zlims = (Z([1 end]));
  xlim(xlims)
  ylim(ylims)
  zlim(zlims)

  a=get(gca,'xtick');
  set(gca,'xtick',a(1:end-1));
  
  camlight left
  camlight head
  lighting phong
  material([0,1,0])

  title(name,'position',[0 50 40],'horizontalalignment','left')
  box on, grid off
  caxis([caxismin caxismax])
  colorbar
  set(gcf,'renderer','opengl')
%   xlabel('x [m]','position',[         (xlims(2)-xlims(1))*0.5  , ylims(1)-(ylims(2)-ylims(1))*0.2 , zlims(1)               ])
%   ylabel('y [m]','position',[xlims(1)-(xlims(2)-xlims(1))*0.15 ,          (ylims(1)+ylims(2))*0.5 , zlims(1)               ])
%   zlabel('z [m]','position',[xlims(1)-(xlims(2)-xlims(1))*0.15 , ylims(2)                         , (zlims(1)+zlims(2))*0.5])
end
