# dea_extensions

This is a collection of R functions for Data Envelopment Analysis (DEA) models that I have developed either for my publications or for further data exploration in the background.

## dea_ar_type1_primal()
Primal of DEA model with Assurance Regions (AR) of Type I.

### Arguments

| Argument  | Description |
| ------------- | ------------- |
| base          | A data frame with *n* rows and *s* + *m* columns, where *n* is the number of Decision-Making Units (DMUs), *s* is the number of outputs and *m* is the number of inputs.  |
| noutput       | Numeric. The number of outputs produced by the DMUs.  |
| rts           | Returns to scale specification. 1 for constant returns to scale and 2 for variable returns to scale.  |
| orientation   | Orientation of the DEA model. 1 for input orientation and 2 for output orientation.  |
| bounds        |   |
| duals        |   |


## dea_gle()

Global efficiencies (GLE) model for ranking Decision-Making Units (DMUs) with a common set of input and output weights (Despotis, 2002).

### Arguments

| Argument  | Description |
| ------------- | ------------- |
| base          | A data frame with *n* rows and *s* + *m* columns, where *n* is the number of Decision-Making Units (DMUs), *s* is the number of outputs and *m* is the number of inputs.  |
| noutput       | Numeric. The number of outputs produced by the DMUs.  |
| eff           | A numeric vector of the efficiency scores.  |
| t_values      | A numeric vector of the values of *t* for which GLE should be run. Should be increasing incrementally from 0 to 1, e.g. by 0.01 or otherwise. Defaults to `seq(0, 1, 0.01)`. |

### Details
The function calculates, for each DMU, the GLE scores and sets of intput and output weights for the different values of *t*, the average GLE score, as well the ranking factor for each DMU. In his original paper, Despotis (2002) proposes to consider the unique sets of weights when these are duplicated for different values of *t*. However, this does not work well in R, because seemingly identical sets of weights for different values of *t* are often not considered equal. Consider the classic example `sqrt(2) ^ 2 == 2` which returns `FALSE`. So if, for example, in a single-input single-output case the set of weights was `c(2, 3)` for *t* = 0.02 and `c(sqrt(2) ^ 2, sqrt(3) ^ 2)` for *t* = 0.03, then `c(2, 3) == c(sqrt(2) ^ 2, sqrt(3) ^ 2)` would return `c(FALSE, FALSE)`. It is therefore much safer to average over *all* GLE scores rather than remove the duplicated ones.

### Value
A list with elements `GLE` and `weights`.

| Element  | Description |
| ------------- | ------------- |
| GLE          | A data frame with *n* rows and 1 + `length(t_values)` + 1 + 1 + 1 columns with the DEA efficiency scores, the GLE scores for the different values of *t*, the average GLE score (AVGLE), variable *q* indicating the number of times a GLE score was equal to 1, and the ranking factor *q* + AVGLE.  |
| weights       | a data frame with `length(t_values)` rows and *s* + *m* columns, containing the sets of common weights for each value of *t*. |

### Note
For *t* = 0 the non-linear model is solved by bisection search (see Despotis, 2002, p. 317). The implementation of bisection search in `dea_gle()` sometimes results in GLE scores that are greater than 1. When that happens, these scores are set to 1. A future improvement of `dea_gle()` should ensure that the bisection search returns scores that are at most 1.

## References
Despotis DK (2002). Improving the discriminating power of DEA: focus on globally efficient units. *Journal of the Operational Research Society* 53: 314–323. DOI: 10.1057=palgrave=jors=2601253.
