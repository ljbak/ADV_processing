function adv_profile_plot(q_prof,rng,c,s,ms,ls)
% plot mean or rms profile q_prof along the ADV profile range using color
% c, symbol s, linestyle l, and markersize ms

if nargin==2
    plot(q_prof,rng,'+','linewidth',1);
elseif nargin==3
    plot(q_prof,rng,'color',c,'symbol','+','linewidth',1);
elseif nargin==4
    plot(q_prof,rng,'color',c,'symbol',s,'linewidth',1);
elseif nargin==5
    plot(q_prof,rng,'color',c,'symbol',s,'markersize',ms,'linewidth',1);
elseif nargin==6
    plot(q_prof,rng,'color',c,'symbol',s,'markersize',ms,'linestyle',ls,'linewidth',1);
else
    error('Incorrect number of input arguments')
end

ylabel('z [m]')