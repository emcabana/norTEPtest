# TEPnormality: Normality Tests Based on Transformed Empirical Processes

This R package implements univariate normality tests using Transformed Empirical Processes (TEPs). 

## Scientific Credits

The methodology behind these tests is based on a general family of consistent and focused tests for normality using isometries that make the empirical process insensitive to unknown location and scale parameters. 

## Key Features

*  Focused Power: The tests are designed to be specifically sensitive to departures in skewness (using the third Hermite polynomial $H_3$) or kurtosis (using the fourth Hermite polynomial $H_4$).
*   Consistency:Despite being focused, these tests remain consistent against all non-normal alternatives.
*   Invariance: The statistics are invariant to unknown location and scale parameters.

## Implemented Statistics

### Quadratic Cramer-von Mises-Watson Type 
Based on quadratic functionals (integrals) of the process:
*   `TEPnorm(X,type,...)`
** 'type=skew': Focused on skewness.
*   `type=kurt`: Focused on kurtosis.

## References

*   Caba&ntilde;a, A., & Caba&ntilde;a, E. M. (2003).Tests of Normality based on transformed empirical processes. Methodology and Computing in applied probability 5, 309-335.
