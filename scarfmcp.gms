$title  Scarf's Activity Analysis Example (SCARFMCP,SEQ=132)
$Ontext

Scarf's Activity Analysis Example


Scarf, H, and Hansen, T, The Computation of Economic Equilibria.
Yale University Press, 1973.

$Offtext

sets
        c       commodities

                /capeop, nondurbl, durable, capbop, skillab, unsklab/

        h       consumers

                /agent1,agent2,agent3,agent4,agent5/

        s       sectors

                /d1, d2, n1, n2, n3, cd, c1, c2/;

alias (c,cc);

table e(c,h)  commodity endowments

             agent1     agent2     agent3     agent4     agent5
capbop         3         0.1         2          1           6
skillab        5         0.1         6          0.1         0.1
unsklab        0.1       7           0.1        8           0.5
durable        1         2           1.5        1           2

table d(c,h) reference demands

             agent1     agent2     agent3     agent4     agent5
capeop         4         0.4          2         5          3
skillab        0.2                    0.5
unsklab                  0.6                    0.2        0.2
nondurbl       2         4           2          5          4
durable        3.2       1           1.5        4.5        2

parameter  esub(h)  elasticities in demand

/       agent1          1.2,
        agent2          1.6,
        agent3          0.8,
        agent4          0.5,
        agent5          0.6  /;

table data(*,c,s)  activity analysis matrix

                         d1          d2          n1          n2          n3

output.nondurbl                                 6.0         8.0         7.0
output.durable          4.0         3.5
output.capeop           4.0         4.0         1.6         1.6         1.6
input .capbop           5.3         5.0         2.0         2.0         2.0
input .skillab          2.0         1.0         2.0         4.0         1.0
input .unsklab          1.0         6.0         3.0         1.0         8.0

              +          cd          c1          c2

output.capeop           0.9         7.0         8.0
input .capbop           1.0         4.0         5.0
input .skillab                      3.0         2.0
input .unsklab                      1.0         8.0;


parameter       alpha(c,h)      demand function share parameter;
alpha(c,h) = d(c,h) / sum(cc, d(cc,h));

parameter  a(c,s)  activity analysis matrix;
a(c,s) = data("output",c,s) - data("input",c,s);

positive
variables

        p(c)    commodity price,
        y(s)    production,
        i(h)    income;

equations
        mkt(c)          commodity market,
        profit(s)       zero profit,
        income(h)       income index;

*       distinguish ces and cobb-douglas demand functions:

mkt(c)..        sum(s, a(c,s) * y(s)) + sum(h, e(c,h)) =g=

                sum(h$(esub(h) ne 1),
                (i(h)/sum(cc, alpha(cc,h) * p(cc)**(1-esub(h)))) *
                alpha(c,h) * (1/p(c))**esub(h)) +

                sum(h$(esub(h) eq 1),
                i(h) * alpha(c,h) / p(c));

profit(s)..     -sum(c, a(c,s) * p(c)) =g= 0;

income(h)..     i(h) =g= sum(c, p(c) * e(c,h));

model scarf / mkt.p, profit.y, income.i/;

p.l(c) = 1;
y.l(s) = 1;
i.l(h) = sum(c, p.l(c) * e(c,h));

p.lo(c)  = 0.00001$(smax(h, alpha(c,h)) gt eps);

*       fix the price of numeraire commodity:

i.fx(h)$(ord(h) eq 1) = i.l(h);

solve scarf using mcp;


