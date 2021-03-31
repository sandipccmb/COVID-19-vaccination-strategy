function out = goveqs_basis3(t, in, M, i, s, p, r, agg, sel, prm)

invec = in(1:i.nstates);

% Get the force of infections
lam = M.lam*invec;

%allmat = M.lin + lam(1)*M.nlin.urb.a1 + lam(2)*M.nlin.urb.a2 + lam(3)*M.nlin.urb.a3 + ...
%               + lam(4)*M.nlin.rur.a1 + lam(5)*M.nlin.rur.a2 + lam(6)*M.nlin.rur.a3;

allmat = M.lin + lam(1)*M.nlin.a1 + lam(2)*M.nlin.a2 + lam(3)*M.nlin.a3;

out = allmat*invec;

% Implement deaths
morts = M.mortvec.*invec;
out(1:i.nstates) = out(1:i.nstates) - morts;

% Get the auxiliaries
out(i.aux.inc)   = agg.inc*(sel.inc.*allmat)*invec;
% out(i.aux.hosp)  = agg.hosp*(sel.hosp.*allmat)*invec;
out(i.aux.mort)  = agg.mor*morts; 
% out(i.aux.hosp2) = sum((sel.hosp2.*Mlin)*invec);

