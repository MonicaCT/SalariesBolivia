*Illustrate use of spreadsheet interfaces XLIMPORT and XLDUMP

 option limrow=0;option limcol=0;
 set places
          /newyork,chicago,topeka,totalsupply/
    sources
          /seattle,sandiego,totalneed/
*$libinclude xlimport places myspread.xls input!a1:e1 cdim=1
*$libinclude xlimport places myspread.xls input!a1:e1 Row
display places

*$libinclude xlimport sources myspread.xls input!a2:a4 Col
display sources
set destinaton(places)
    source(sources);

    destinaton(places)=yes;
    destinaton("totalsupply")=no;
    source(sources)=yes;
    source("totalneed")=no;
 parameter trandata (sources,places) transport data from spreadsheet
           Supply(Sources)   Supply at each source plant in cases
           Need(places)      Amount neeeded at each market destination in cases;
$libinclude xlimport trandata myspread.xls input!a1:e4

  supply(source) = trandata(source,"totalsupply");
  Need(destinaton) = trandata("totalneed",destinaton);
  Scalar
     prmilecst   freight cost in $ per case per 1000 miles /90/
     loadcost    freight loading cost in $ per case /25/        ;
  Parameter trancost(Sources,places)  transport cost in dollars per case ;
     trancost(Source,Destinaton) =
           loadcost + prmilecst * trandata(Source,Destinaton) ;
  Positive Variable
      transport(Sources,places) shipment quantities in cases;
  Variable
      totalcost total transportation costs in dollars ;
  Equations
         Costsum                 total transport cost -- objective function
         supplybal(sources)      supply at sources
         Demandbal(places)       demand at destinations ;

  Costsum ..  totalcost  =e=  sum((Source,Destinaton),
           trancost(Source,Destinaton)*transport(Source,Destinaton));

  Supplybal(Source) ..
      sum(destinaton, transport(Source,Destinaton))
              =l=  supply(Source) ;

    demandbal(Destinaton) ..
     sum(Source, transport(Source,Destinaton))
             =g=  need(Destinaton) ;

    Model tranport /all/ ;
    Solve tranport using lp minimizing totalcost ;

*$libinclude xlexport transport.l myspread.xls output3!a1..d4
*$libinclude xlexport transport.m myspread.xls output3!f1..i4
*$libinclude xlexport transport.l myspread.xls output3!a6:d8  /m
*$libinclude xlexport transport.l myspread.xls output3!f6:i9  /m
*$libinclude xlexport transport.l myspread.xls output3!j1
$libinclude xldump transport.l myspread.xls output2!a1
$libinclude xldump transport.l myspread.xls output!a1..D4
$libinclude xldump transport.l myspread.xls output4!a1..D4

$ontext
#user model library stuff
Main topic Linking to other programs
Featured item 1 Spreadsheet
Featured item 2 XLDUMP.gms
Featured item 3 XLIMPORT.gms
Featured item 4
include myspread.xls
Description
Illustrate use of spreadsheet interfaces XLIMPORT and XLDUMP

$offtext