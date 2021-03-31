function init = get_init(p, r, i, s, gps, prm)

init = zeros(1,i.nx);

% --- Allocate initial immunity -------------------------------------------

mat = zeros(length(gps.age),length(gps.risk), length(gps.vac));
mat(:,1,1) = prm.N.*(1-p.risk).*(1-p.vacc1);
mat(:,1,2) = prm.N.*(1-p.risk).*p.vacc1;
mat(:,2,1) = prm.N.*p.risk.*(1-p.vacc2);
mat(:,2,2) = prm.N.*p.risk.*p.vacc2;

% keyboard;
% mat(:,:,1) has two column vector; 1st: non-vaccinated/non risk; 2nd: non-vaccinated/risk; 
% mat(:,:,2) has two column vector; 1st: vaccinated/non risk; 2nd: vaccinated/risk; 

tmp = mat.*repmat((1-p.seropos'),[1,length(gps.risk), length(gps.vac)]);
init(s.U) = tmp(:)';

tmp = mat.*repmat(p.seropos',[1,length(gps.risk), length(gps.vac)]);
init(s.R) = tmp(:)';


% --- Allocate seeded infection -------------------------------------------
seed = 10;
init(i.MS.v0.r0.a2) = seed; 
init(i.U.v0.r0.a2)  = init(i.U.v0.r0.a2) - seed;
% init(i.SS.el.v0) = seed; init(i.U.el.v0) = init(i.U.el.v0) - seed;
% init(i.SS.hcw.v0) = seed; init(i.U.hcw.v0) = init(i.U.hcw.v0) - seed;

% keyboard;

end














