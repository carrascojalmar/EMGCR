# EMGLMLasso
## Mixture cure models with flexible survival and link functions: Parameter estimation via the EM algorithm
### Authors
Chaeyeon Yoo,
Dipak K. Dey,
Víctor H. Lachos,
Jalmar M. F. Carrasco.
### Abstract
The mixture cure model is a special type of survival model, and it assumes that the studied population is a mixture of susceptible individuals who may experience the event of interest and
cure/non-susceptible individuals who will never experience the event. The Weibull model is typically specified as a latency component for the event time and symmetric logistic regression as an incidence component for the probability of uncured. Analytical challenges arise when incorporating alternative survival distributions and link functions. To address this, we propose a flexible mixture cure rate model, where parameter estimation is performed using the Expectation-Maximization algorithm. This approach is easily implementable in statistical software that supports weighted survival regression models along with the efficient routine glm(), which allows the definition of custom link functions. The proposed method delivers notable gains in computational efficiency, with standard errors as a byproduct. In simulations under varied scenarios and in an analysis of real-world liver cancer data, we demonstrate the efficiency and flexibility of our proposed model.
