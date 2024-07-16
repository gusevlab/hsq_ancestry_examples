Re-analysis of Ritchie, Bates, Deary et al.
================
2024-07-15

## Load in the data

``` r
input_mat = '
1                                           
0.42    1                                       
0.43    0.31    1                                   
0.44    0.3 0.31    1                               
0.46    0.31    0.33    0.37    1                           
0.46    0.31    0.28    0.39    0.57    1                       
0.32    0.22    0.48    0.24    0.31    0.27    1                   
0.47    0.27    0.33    0.62    0.45    0.48    0.22    1               
0.45    0.25    0.4 0.41    0.44    0.4 0.3 0.45    1           
0.42    0.21    0.3 0.3 0.4 0.34    0.27    0.34    0.54    1       
0.28    0.14    0.24    0.31    0.38    0.4 0.16    0.41    0.42    0.32    1   
0.67    0.39    0.46    0.49    0.58    0.51    0.35    0.53    0.51    0.4 0.36    1'
input_names = c("Age11_MHT","Edu","LM","DSS","MR","BD","VPA","SS","LNS","DSB","SpSp","Age70_MHT")
cor_mat = getCov(input_mat, names = input_names, lower=T, diag=T  )
ntotal = 1028
```

## Fit and confirm the published models

“*Using modification indexes calculated in Mplus, we found five residual
covariances that were significant in the baseline model. Four of these
described clear content overlap in the tests (between Matrix Reasoning
and Block Design, Logical Memory and Verbal Paired Associates,
Digit-Symbol and Symbol Search, Digit Span Backwards and Letter-Number
Sequencing) and one was unexpectedly negative (between the MHT and
Spatial Span).*”

``` r
# ---
EduOnlyG.model <- '
g_adult =~ Age70_MHT + LM + DSS + MR + BD + VPA + SS + LNS + DSB + SpSp
g_adult + Edu ~ Age11_MHT
g_adult ~ Edu
Age70_MHT ~~ SpSp
MR ~~ BD
LM ~~ VPA
DSS ~~ SS
DSB ~~ LNS'

EduGandS.model <- '
g_adult =~ Age70_MHT + LM + DSS + MR + BD + VPA + SS + LNS + DSB + SpSp
g_adult + Edu ~ Age11_MHT
g_adult + DSS + LM ~ Edu
Age70_MHT ~~ SpSp
MR ~~ BD
LM ~~ VPA
DSS ~~ SS
DSB ~~ LNS'

EduOnlyS.model <- '
g_adult =~ Age70_MHT + LM + DSS + MR + BD + VPA + SS + LNS + DSB + SpSp
g_adult + Edu ~ Age11_MHT
Age70_MHT + LM + DSS + MR + BD + VPA + SS ~ Edu
MR ~~ BD
LM ~~ VPA
DSS ~~ SS
DSB ~~ LNS
Age70_MHT ~~ SpSp'

EduOnlyG.fit <- cfa(EduOnlyG.model, sample.cov = cor_mat , sample.nobs = ntotal , std.lv = T )
EduOnlyS.fit <- cfa(EduOnlyS.model, sample.cov = cor_mat , sample.nobs = ntotal , std.lv = T )
EduGandS.fit <- cfa(EduGandS.model, sample.cov = cor_mat , sample.nobs = ntotal , std.lv = T )

cat("Model 1" , fitMeasures(EduOnlyG.fit,c("rmsea","aic")),'\n')
```

    ## Model 1 0.05806203 27644.45

``` r
cat("Model 2" , fitMeasures(EduGandS.fit,c("rmsea","aic")),'\n')
```

    ## Model 2 0.05712551 27634.42

``` r
cat("Model 3" , fitMeasures(EduOnlyS.fit,c("rmsea","aic")),'\n')
```

    ## Model 3 0.05735391 27626.13

``` r
cat ("dAIC Only S vs. G and S" , fitMeasures(EduOnlyS.fit,"aic") - fitMeasures(EduGandS.fit,"aic") , '\n' )
```

    ## dAIC Only S vs. G and S -8.289836

``` r
cat ("dAIC Only S vs. Only G" , fitMeasures(EduOnlyS.fit,"aic") - fitMeasures(EduOnlyG.fit,"aic") , '\n' )
```

    ## dAIC Only S vs. Only G -18.32223

We confirm the lower (better) AIC for a model with paths only on S
compared to a model with paths to g and S, and much lower AIC for a
model only on g.

“*We then compared Model C to the previous models. It had significantly
better fit than both Model A, ΔAIC = 19.08, χ2(6) = 31.08, p \< .001;
and Model B, ΔAIC = 9.90, χ2(4) = 17.90, p = .001.*”

## Plot Model 1 : Only on g

``` r
lavaanPlot( model = EduOnlyG.fit , coefs=T , stand=T )
```

![](ritchie_reanalysis_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->
\## Plot Model 2 : On g and s

``` r
lavaanPlot( model = EduGandS.fit , coefs=T , stand=T )
```

![](ritchie_reanalysis_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->
\## Plot Model 3 : Only on s

``` r
lavaanPlot( model = EduOnlyS.fit , coefs=T , stand=T )
```

![](ritchie_reanalysis_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## Do model selection for Model 2

``` r
# start without the effect of Edu on S

EduGandS.baseline <- '
g_adult =~ Age70_MHT + LM + DSS + MR + BD + VPA + SS + LNS + DSB + SpSp
g_adult + Edu ~ Age11_MHT
Age70_MHT ~~ SpSp
MR ~~ BD
LM ~~ VPA
DSS ~~ SS
DSB ~~ LNS'

all_vals = c("LM","DSS","MR","BD","VPA","SS","LNS","DSB","SpSp","Age70_MHT")
min_p = NA
included_vals = "g_adult"
# Iterate until no more significant paths
while( is.na(min_p) | min_p < 0.05 ) {
  # Iterate through all specific tests
  min_p = 1
  for ( i in 1:length(all_vals) ) {
      new_vals = c(included_vals,all_vals[i])
      eqn = paste(paste(new_vals,collapse=" + ")," ~ Edu",sep='')
      new_model = paste(EduGandS.baseline,'\n',eqn,sep='')
      EduGandS.fit <- cfa(new_model, sample.cov = cor_mat , sample.nobs = ntotal , std.lv = T)
      coefs = standardizedSolution(EduGandS.fit)
      cur_pv = coefs[coefs$lhs == all_vals[i] & coefs$op == "~" & coefs$rhs == "Edu","pvalue"]
      if ( cur_pv < min_p ) {
        min_val = all_vals[i]
        min_p = cur_pv
      }
  }
  
  # if there's a significant path, add it
  cat(min_val,min_p,'\n')
  if ( min_p < 0.05 ) {
    included_vals = c(included_vals,min_val)
    all_vals = setdiff(all_vals,min_val)
  }
  cat('\n')
}
```

    ## SpSp 3.793093e-05 
    ## 
    ## LM 0.01059601 
    ## 
    ## DSS 0.02503372 
    ## 
    ## LNS 0.08537644

``` r
# print final model
cat( "Final specific paths:" , included_vals , '\n' )
```

    ## Final specific paths: g_adult SpSp LM DSS

Interestingly, SpSp was the first selected path and highly significant.
LM and DSS were included as paths in Ritchie et al. but SpSp was not
(reasons unstated).

Now let’s test the relative model fit.

``` r
EduGandS_new.model <- '
g_adult =~ Age70_MHT + LM + DSS + MR + BD + VPA + SS + LNS + DSB + SpSp
g_adult + Edu ~ Age11_MHT
g_adult + DSS + LM + SpSp ~ Edu
Age70_MHT ~~ SpSp
MR ~~ BD
LM ~~ VPA
DSS ~~ SS
DSB ~~ LNS'

EduGandS_new.fit <- cfa(new_model, sample.cov = cor_mat , sample.nobs = ntotal , std.lv = T)
cat ("dAIC Only S vs. G and S (updated)" , fitMeasures(EduOnlyS.fit,"aic") - fitMeasures(EduGandS_new.fit,"aic") , '\n' )
```

    ## dAIC Only S vs. G and S (updated) 3.484087

The updated model fits better than the model only on S. Let’s visualize
it:

``` r
lavaanPlot( model = EduGandS_new.fit , coefs=T , stand=T )
```

![](ritchie_reanalysis_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->
There’s a significant positive path from Edu to g, and a negative path
from Edu to SpSp. Let’s look at the stanardized coefficients:

``` r
standardizedSolution(EduGandS.fit)
```

    ##          lhs op       rhs est.std    se      z pvalue ci.lower ci.upper
    ## 1    g_adult =~ Age70_MHT   0.825 0.019 44.304  0.000    0.788    0.861
    ## 2    g_adult =~        LM   0.497 0.028 17.594  0.000    0.441    0.552
    ## 3    g_adult =~       DSS   0.557 0.026 21.243  0.000    0.506    0.609
    ## 4    g_adult =~        MR   0.679 0.019 36.216  0.000    0.642    0.715
    ## 5    g_adult =~        BD   0.640 0.020 31.479  0.000    0.600    0.679
    ## 6    g_adult =~       VPA   0.415 0.028 15.066  0.000    0.361    0.469
    ## 7    g_adult =~        SS   0.661 0.019 34.178  0.000    0.623    0.699
    ## 8    g_adult =~       LNS   0.649 0.020 32.717  0.000    0.610    0.687
    ## 9    g_adult =~       DSB   0.530 0.024 21.902  0.000    0.482    0.577
    ## 10   g_adult =~      SpSp   0.590 0.029 20.113  0.000    0.532    0.647
    ## 11   g_adult  ~ Age11_MHT   0.667 0.021 31.980  0.000    0.626    0.708
    ## 12       Edu  ~ Age11_MHT   0.420 0.025 17.123  0.000    0.372    0.468
    ## 13 Age70_MHT ~~      SpSp  -0.195 0.042 -4.610  0.000   -0.278   -0.112
    ## 14        MR ~~        BD   0.241 0.032  7.430  0.000    0.177    0.304
    ## 15        LM ~~       VPA   0.334 0.028 11.725  0.000    0.278    0.390
    ## 16       DSS ~~        SS   0.386 0.028 13.573  0.000    0.330    0.441
    ## 17       LNS ~~       DSB   0.304 0.030 10.095  0.000    0.245    0.363
    ## 18   g_adult  ~       Edu   0.155 0.030  5.216  0.000    0.097    0.214
    ## 19      SpSp  ~       Edu  -0.117 0.032 -3.647  0.000   -0.180   -0.054
    ## 20        LM  ~       Edu   0.082 0.029  2.835  0.005    0.025    0.139
    ## 21       DSS  ~       Edu   0.065 0.027  2.405  0.016    0.012    0.117
    ## 22 Age70_MHT  ~       Edu   0.031 0.026  1.201  0.230   -0.020    0.082
    ## 23 Age70_MHT ~~ Age70_MHT   0.297 0.021 13.999  0.000    0.255    0.338
    ## 24        LM ~~        LM   0.711 0.025 28.125  0.000    0.662    0.761
    ## 25       DSS ~~       DSS   0.654 0.026 25.314  0.000    0.603    0.705
    ## 26        MR ~~        MR   0.539 0.025 21.208  0.000    0.490    0.589
    ## 27        BD ~~        BD   0.591 0.026 22.736  0.000    0.540    0.642
    ## 28       VPA ~~       VPA   0.828 0.023 36.229  0.000    0.783    0.873
    ## 29        SS ~~        SS   0.563 0.026 22.034  0.000    0.513    0.613
    ## 30       LNS ~~       LNS   0.579 0.026 22.538  0.000    0.529    0.630
    ## 31       DSB ~~       DSB   0.719 0.026 28.049  0.000    0.669    0.769
    ## 32      SpSp ~~      SpSp   0.698 0.028 25.036  0.000    0.644    0.753
    ## 33       Edu ~~       Edu   0.824 0.021 39.973  0.000    0.783    0.864
    ## 34   g_adult ~~   g_adult   0.444 0.023 19.699  0.000    0.400    0.489
    ## 35 Age11_MHT ~~ Age11_MHT   1.000 0.000     NA     NA    1.000    1.000

## Sensitivity analysis: Remove SpSp entirely from the analysis

``` r
EduOnlyG.model <- '
g_adult =~ Age70_MHT + LM + DSS + MR + BD + VPA + SS + LNS + DSB
g_adult + Edu ~ Age11_MHT
g_adult ~ Edu
MR ~~ BD
LM ~~ VPA
DSS ~~ SS
DSB ~~ LNS'

EduGandS.model <- '
g_adult =~ Age70_MHT + LM + DSS + MR + BD + VPA + SS + LNS + DSB
g_adult + Edu ~ Age11_MHT
g_adult + DSS + LM ~ Edu
MR ~~ BD
LM ~~ VPA
DSS ~~ SS
DSB ~~ LNS'

EduOnlyS.model <- '
g_adult =~ Age70_MHT + LM + DSS + MR + BD + VPA + SS + LNS + DSB
g_adult + Edu ~ Age11_MHT
Age70_MHT + LM + DSS + MR + BD + VPA + SS ~ Edu
MR ~~ BD
LM ~~ VPA
DSS ~~ SS
DSB ~~ LNS'

EduOnlyG.fit <- cfa(EduOnlyG.model, sample.cov = cor_mat , sample.nobs = ntotal , std.lv = T )
EduOnlyS.fit <- cfa(EduOnlyS.model, sample.cov = cor_mat , sample.nobs = ntotal , std.lv = T )
EduGandS.fit <- cfa(EduGandS.model, sample.cov = cor_mat , sample.nobs = ntotal , std.lv = T )

cat("Model 1" , fitMeasures(EduOnlyG.fit,c("rmsea","aic")),'\n')
```

    ## Model 1 0.04781269 24981.4

``` r
cat("Model 2" , fitMeasures(EduGandS.fit,c("rmsea","aic")),'\n')
```

    ## Model 2 0.04685773 24975.26

``` r
cat("Model 3" , fitMeasures(EduOnlyS.fit,c("rmsea","aic")),'\n')
```

    ## Model 3 0.04997851 24980.49

``` r
fitMeasures(EduOnlyS.fit,"aic") - fitMeasures(EduGandS.fit,"aic")
```

    ##   aic 
    ## 5.223

``` r
fitMeasures(EduOnlyS.fit,"aic") - fitMeasures(EduOnlyG.fit,"aic")
```

    ##    aic 
    ## -0.915

The model with Edu paths to g and S still fits better than a model with
only paths to S. Moreover, the model with paths only to g now has
approximately the same fit as the model with paths only to S. All of the
primary results were apparently driven by the inclusion of SpSp.
