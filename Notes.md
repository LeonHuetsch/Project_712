Notes

1. Use equidistant grid for assets first, exp(ln) makes it quite weird.
2. For some reason savings do not cross 45 degree line how do I get this?


3. If i go from sigma=1 to sigma=1.5 running time decreases drastically. Why?

4. If maxGridCap is very high it takes way longer. Probably because convergence in these areas takes a long time?

5. If policy function for assets does not cross 45 degree line, need to increase max grid value. The necessary max grid value can also be reduced by lowering precautionary savings motive, for example lower variance of income process (risk), lower risk aversion (sigma), lower persistence of income process.


6. Problem when delta=0 (persistence of shock): transition matrix by itself looks find, it is iid. vgridShock does not change much compared to very low persistence and looks find. However stationary distribution (eigenvector of transition matrix) sums up to almost zero, so normalizing makes the values huge. Thus, the product of the grid vector and stationary distribution, which is the average income and/or the aggregate income is huge. Then normalizing the grid with that number yields only grid value basically equal to zero.
   Solution: Took wrong eigenvector, need to take the one with eigenvalue 1, which for delta=0 is the 2nd not first one!!!


7. When persistence delta is low (tried 0 and 0.3) it cannot find the equilibrium r, the difference does no go down anymore ad expected asset holdings seem to jump a lot. WHY?   
