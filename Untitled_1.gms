

*$Offtext
sets
        f  factors    /labor, capital/
        s  sectors    /mfrs,  nonmfrs/
        h  households /rich,  poor/;

alias (h,k), (s,ss), (f,ff);
*
*       demand function parameters.
*
parameter sigmac(h)
       / rich    1.5 ,  poor    0.75/;

table alpha(s,h)
                rich    poor
        mfrs    0.5     0.3
        nonmfrs 0.5     0.7;

table e(f,h)
                rich    poor
        labor             60
        capital   25  ;
*
*       production function parameters.
*
parameter phi(s)
        / mfrs 1.5,  nonmfrs 2.0 /;

table delta(f,s)
                        mfrs    nonmfrs
        labor           0.6     0.7
        capital         0.4     0.3;

parameter sigma(s)
      /  mfrs 2.0,   nonmfrs 0.5/;

parameter       tshr(h) share of tax revenue /rich 0.4, poor 0.6/,
                t(f,s)  ad-valorem tax rates;

t(f,s) = 0;

positive
variables
        w(f)            factor price,
        p(s)            commodity price,
        y(s)            production level,
        i(h)            income;

equations
        fmkt(f)         factor market,
        cmkt(s)         commodity market,
        profit(s)       zero profit,
        income(h)       income equation;

fmkt(f)..       sum(h, e(f,h)) =g=
                sum(s, y(s) * phi(s)**(sigma(s)-1) *
                (delta(f,s) * (sum(ff, delta(ff,s)**sigma(s) *
                (w(ff)*(1 + t(ff,s)))**(1 - sigma(s)))
                **(1/(1-sigma(s)))/phi(s))
                / (w(f) * (1 + t(f,s))))**sigma(s));

cmkt(s)..       y(s) =g= sum(h,
                (i(h)/sum(ss, alpha(ss,h) * p(ss)**(1-sigmac(h)))) *
                alpha(s,h) * (1 /p(s))**sigmac(h));

profit(s)..     sum(f, delta(f,s)**sigma(s) *
                (w(f)*(1 + t(f,s)))**(1 - sigma(s)))**(1/(1-sigma(s)))/phi(s)
                =g= p(s);

income(h)..     i(h) =g= sum(f, e(f,h) * w(f)) + tshr(h) *
                sum((s,f), t(f,s) * w(f) * y(s) * phi(s)**(sigma(s)-1) *
                (delta(f,s) * (sum(ff, delta(ff,s)**sigma(s) *
                (w(ff)*(1 + t(ff,s)))**(1 - sigma(s)))
                **(1/(1-sigma(s)))/phi(s))/(w(f) * (1 + t(f,s))))**sigma(s));

model jel / fmkt.w, cmkt.p, profit.y, income.i/;

*       compute solution for this dimension problem:

w.lo(f) = 0.0001;
p.lo(s) = 0.0001;

w.l(f) = 1;
p.l(s) = 1;
y.l(s) = 10;
i.l(h) = sum(f, w.l(f) * e(f,h));

*       solve the reference case:

i.fx(h) = i.l(h);
solve jel using mcp;
