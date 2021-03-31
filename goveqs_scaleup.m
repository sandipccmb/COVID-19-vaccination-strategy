function out = goveqs_scaleup(t, in, M0, M1, i, s, p, r, agg, sel, interval)

scale = min(max((t-interval(1))/(interval(2)-interval(1)),0),1);

M1.lin = M0.lin + scale*(M1.lin - M0.lin);
M1.lam = M0.lam + scale*(M1.lam - M0.lam);

out = goveqs_basis3(t, in, M1, i, s, p, r, agg, sel);