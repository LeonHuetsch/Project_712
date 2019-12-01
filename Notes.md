Notes

1. Use equidistant grid for assets first, exp(ln) makes it quite weird.
2. For some reason savings do not cross 45 degree line how do I get this?


3. If i go from sigma=1 to sigma=1.5 running time decreases drastically. Why?

4. If maxGridCap is very high it takes way longer. Probably because convergence in these areas takes a long time?

5. If policy function for assets does not cross 45 degree line, need to increase max grid value. The necessary max grid value can also be reduced by lowering precautionary savings motive, for example lower variance of income process (risk), lower risk aversion (sigma), lower persistence of income process.


6. Problem when delta=0 (persistence of shock): transition matrix by itself looks find, it is iid. vgridShock does not change much compared to very low persistence and looks find. However stationary distribution (eigenvector of transition matrix) sums up to almost zero, so normalizing makes the values huge. Thus, the product of the grid vector and stationary distribution, which is the average income and/or the aggregate income is huge. Then normalizing the grid with that number yields only grid value basically equal to zero.
   Solution: Took wrong eigenvector, need to take the one with eigenvalue 1, which for delta=0 is the 2nd not first one!!!


7. When persistence delta is low (tried 0 and 0.3) it cannot find the equilibrium r, the difference does no go down anymore ad expected asset holdings seem to jump a lot. WHY?


8. Problems when maxAsset is higher: If 15 it all works and finds GE interest rate for all parameter specifications. When large, e.g. 35, neither fsolve nor fzero can find them for some parameter specifications: particularly for small delta (=0.3). Possible reason: less saving desire means interested rate has to be higher to get people to save and there is a problem with that. Is it becasue r is larger then rho? or initial guess way too small?


9. Any interest rate here is net of depreciation: so r = r1 - depreciation, where r1 is the real net interest rate. Write in PDF!!!
