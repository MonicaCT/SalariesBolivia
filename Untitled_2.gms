  SETS Factor Basic factors of production /Labor, Capital/
Sector Producing industries /Food, NonFood/
households Households Types /Farm, NonFarm/ ;
ALIAS (households,households1), (Sector,OtherSector), (Factor,Factor1);
PARAMETER SigmaC(households) Elasticity of substitution in household CES
/NonFarm 1.5
Farm 1.00/ ;
TABLE Alpha(Sector,households) Consumption share in household CES
NonFarm Farm
Food 0.5 0.4
NonFood 0.5 0.6 ;
TABLE Endowment(Factor,households) Household endowments of factors
NonFarm Farm
Labor 60 5
Capital 25 2 ;
Table Intermediateuse(OtherSector,sector) Use in othersector when producing in sector
Nonfood Food
Nonfood 0.00 0.0
Food 0.00 0.0;
PARAMETER Phi(Sector) scale parameter in CES production function
/Food 2.2
NonFood 2.0/ ;
TABLE Delta(Factor,Sector) distribution parameter in CES production
Food NonFood
Labor 0.6 0.7
Capital 0.4 0.3 ;
PARAMETER Sigma(Sector) Elasticity of production factor substitution
/Food 2.0
NonFood 0.5/
Parameter TaxShare(households) Household share of tax disbursements
/NonFarm 0.7
Farm 0.3/ ;
Parameter Governmentpurch(sector) Government goods purchase dependence on revenue
/Food 0.00
NonFood 0.00 /
Scalar Incometax Household tax level /0.001/;
Parameter CorporateTax(Sector) Tax on corporate earnings
/food 0, nonfood 0/;
Parameter TaxDeduction(Households) Household Tax Deductions
/farm 0, NonFarm 0/;
POSITIVE VARIABLES
FACTORPRICE(Factor) Prices for the factors labor and capital
FACTORQUAN(Factor,Sector) Factor use by a producing sector
COMPRICE(Sector) Prices of commodities
DemCommod(Households,Sector) Commodity Demand by Households
PRODUCTION(Sector) Production Quantity
HHIncome(households) Household income
TAXREVENUE Total government tax revenue ;
EQUATIONS
FactorMkt(Factor) Factor market balances
FactorDem(Factor,Sector) Factor demand by a sector
CommodMkt(Sector) Commodity market balance
CommodDem(Households,Sector) Commodity Demand by Households
Profit(Sector) Zero profit condition
Income(households) Household budget constraint
GovBal Government budget constraint ;
FactorMkt(Factor)..
sum(households,Endowment(Factor,households)) =G=
sum(Sector,FactorQuan(Factor,Sector)) ;
FactorDem(Factor,Sector)..
FACTORQUAN(Factor,Sector) =g=
Production(Sector)*Phi(Sector)**(sigma(Sector)-1)
*(Delta(Factor,Sector)
*( sum(Factor1,Delta(Factor1,Sector)**sigma(Sector)
*(FactorPrice(Factor1))**(1 - sigma(Sector)))
**(1/(1-sigma(Sector)))/Phi(Sector))/
FactorPrice(Factor))**sigma(Sector) ;
CommodDem(Households,Sector)..
DemCommod(Households,Sector) =E=
(HHIncome(households)
/sum(otherSector,alpha(otherSector,households)
*ComPrice(OtherSector)**(1-SigmaC(households)))
)*Alpha(Sector,households) * (1 /ComPrice(Sector))**sigmaC(households);
CommodMkt(Sector)..
Production(Sector) =G=
sum(households,DemCommod(Households,Sector))
+ Governmentpurch(sector)*TaxRevenue/COMPRICE(Sector)
+ sum(othersector,intermediateuse(sector,othersector)*production(othersector));
Profit(Sector)..
sum(Factor,FactorPrice(Factor)*FactorQuan(Factor,Sector))
+corporatetax(sector)*(
ComPrice(Sector)* Production(Sector)
-sum(othersector,ComPrice(OtherSector)*intermediateuse(othersector,sector)
*production(sector))
-sum(Factor,FactorPrice(Factor)*FactorQuan(Factor,Sector)))
=G= ComPrice(Sector)* Production(Sector) ;
Income(households)..
HHIncome(households) =G=
(1-incometax)*sum(Factor,Endowment(Factor,households) * FactorPrice(Factor))
+ incometax * TaxDeduction(Households)
+ TaxShare(households) *TaxRevenue ;
GovBal..
TaxRevenue =G=
sum(sector,corporatetax(sector)*(
ComPrice(Sector)* Production(Sector)
-sum(Othersector,ComPrice(OtherSector)*intermediateuse(Othersector,sector)
*production(sector))
-sum(Factor,FactorPrice(Factor)*FactorQuan(Factor,Sector))))
+sum(households,incometax
*(sum(Factor,Endowment(Factor,households) * FactorPrice(Factor))
-TaxDeduction(Households)));
MODEL CGEModel /FactorMkt.FACTORPRICE, FactorDem.FACTORQUAN,
CommodMkt.COMPRICE,commoddem.DemCommod
Profit.PRODUCTION, Income.Hhincome, GovBal.TAXREVENUE/;
* lower bounds
FACTORPRICE.LO(Factor) = 0.01;
COMPRICE.LO(Sector) = 0.01;
HHincome.Lo(households) = 0.01;
*starting point
FACTorPRICE.L(Factor) = 1 ;
FACTorQuan.L(Factor,sector) = 1 ;
DemCommod.l(Households,Sector)=1;
COMPRICE.L(Sector) = 1 ;
PRODUCTION.L(Sector) = 10;
HHincome.L(households) =
sum(Factor, FACTorPRICE.l(Factor) * Endowment(Factor,households));
HHincome.fx(households)$(ord(households) =1) = HHincome.L(households);
OPTION MCP = PATH;
SOLVE CGEModel USING MCP;