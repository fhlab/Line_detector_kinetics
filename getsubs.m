function [cths,cthe]=getsubs(time,subs_i)
%%%calculate substrate concentration based on flow parameters


qw=2; % loading withdrawal rate in uL/min (ratio water/oil in drip drop)
qi=10; % infuse rate in uL/min 
Vi=40; % initial volume in uL
%qr=4; %reading flow rate in ul/min

enz_i=1; %initial enzyme concentration - used later to correct slopes

cthe=enz_i*(Vi./(Vi+(qi-qw).*time./60));
cths=subs_i.*(qi.*time./60./(Vi+(qi-qw).*time./60));

end