function add_abc_tofig(abcd,xoff,yoff,fnsize)
%function add_abc_tofig(abcd,xoff,yoff,fnsize)
%
% DGM 3/7/18
if nargin<4,fnsize=8;end
if nargin<3,yoff=-0.01;end
if nargin<2,xoff=-0.015;end
zz=get(gca,'OuterPosition');
textboxloc=[zz(1)+xoff zz(2)+zz(4)+yoff .025 .025];
hi=annotation('textbox',textboxloc,'String',['(' char('a'+abcd-1), ')']);
set(hi,'LineStyle','none');set(hi,'FontSize',fnsize);
