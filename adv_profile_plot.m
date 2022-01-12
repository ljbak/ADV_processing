function lh = adv_profile_plot(q_prof,rng,c,s,ms,ls)
% plot mean or rms profile q_prof along the ADV profile range using color
% c, symbol s, linestyle l, and markersize ms. returns handle to line lh

if nargin==2
    lh = plot(q_prof,rng,'.','linewidth',1,'markersize',10);
elseif nargin==3
    lh = plot(q_prof,rng,'color',c,'symbol','+','linewidth',1);
elseif nargin==4
    lh = plot(q_prof,rng,'color',c,'symbol',s,'linewidth',1);
elseif nargin==5
    lh = plot(q_prof,rng,'color',c,'symbol',s,'markersize',ms,'linewidth',1);
elseif nargin==6
    lh = plot(q_prof,rng,'color',c,'symbol',s,'markersize',ms,'linestyle',ls,'linewidth',1);
else
    error('Incorrect number of input arguments')
end

ylabel('z [m]')