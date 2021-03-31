function [p,r] = alloc_parameters(x, xi, p, r, R0, i, s, gps, prm)

% --- Allocate parameters -------------------------------------------------
r.incub           = 1/x(xi.incub_period);                                                
p.sympto          = x(xi.p_sympto);                                                           
p.c               = x(xi.p_c);                                                           
r.eta             = 1/x(xi.presym_period);                                               
r.gamma           = 1/x(xi.infec_period);
%r.vcm             = x(xi.vcm);
%r.vmu             = x(xi.vmu);
r.mu      = (1+x(xi.vmu)).*p.cfr*r.gamma./(1-p.cfr);
%prm.contact = (1+x(xi.vcm)).*[[1.3710,1.4317,0.0500];[2.5175,2.9007,0.0993];[0.2789,0.3350,0.0160]];% %Age group(0-24,25-60,60+) 
mat = 1+reshape(x(xi.vcm),3,3);
prm.contact = [[1.3710,1.4317,0.0500];[2.5175,2.9007,0.0993];[0.2789,0.3350,0.0160]].*mat;

% % --- Find beta so that we get given value of R0 --------------------------

r.beta = 1;
p1 = p; p1.vacc1 = p.vacc1*0; p1.vacc2 = p.vacc2*0;

M = make_model2(p1, r, i, s, gps, prm);

inds = [s.E, s.infectious];

% Construct the new entries to E due to infection
tmp = repmat(M.lam,length(gps.risk)*length(gps.vac),1);
tmp((size(tmp,1)/2+1):end,:) = tmp((size(tmp,1)/2+1):end,:)*p.c1;

m = zeros(i.nstates);
m(s.E,:) = tmp;
% Account for the initial number in each compartment
mat = zeros(length(gps.age), length(gps.risk), length(gps.vac));
tmp = prm.N.*(1-p.seropos);
mat(:,1,1) = tmp.*(1-p.risk).*(1-p.vacc1);
mat(:,1,2) = tmp.*(1-p.risk).*p.vacc1;
mat(:,2,1) = tmp.*p.risk.*(1-p.vacc2);
mat(:,2,2) = tmp.*p.risk.*p.vacc2;
% Collapse them into vector form that lines up with elements of s.U
m(s.E,:) = m(s.E,:).*repmat(mat(:),1,size(m,2));
F = m(inds,inds);
% Construct V
tmp = -(M.lin - diag(M.mortvec));
V   = tmp(inds,inds);
% Bring them together
r.beta = R0/max(eigs(V\F));

