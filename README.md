# GibbsPy
This is a Python function that calculates several thermodynamic parameters like gas phase or aqueous phase gibbs energy of formation.

## Functions
I suggest using only the following 4 fuctions: gibbs, gibbsAQ, henrys_coef, and fugcoef. Gibbs energies calculated with these functions use the "normal" convention for Gibbs energy of formation as used on [this wikipedia page](https://en.wikipedia.org/wiki/Standard_Gibbs_free_energy_of_formation).

**gibbs(species, T)** - Calculates the Gibbs energy of formation (in units J/mol) of a GAS PHASE species at temperature T (kelvin). This function uses the [NASA thermodynamic database](https://publications.anl.gov/anlpubs/2005/07/53802.pdf).

**gibbsAQ(species, T, P)** - Calculates the Gibbs energy of formation (in units J/mol) of a AQUEOUS PHASE species at temperature T (kelvin), and pressure P (bar). This function uses the [SUPCRT thermodynamic database](https://www.sciencedirect.com/science/article/pii/009830049290029Q).

**henrys_coef(species, T, P)** - Calculates the Henry law coefficient (mol L<sup>-1</sup> bar<sup>-1</sup>) of a species at temperature T (kelvin), and pressure P (bar).

**fugcoef(temperature, pressure, names, n)** - This function uses the Soave Equation (see Equation 27 and 28 in [Krissansen-Totton et al. 2016](https://www.liebertpub.com/doi/full/10.1089/ast.2015.1327?casa_token=LVXXfDQznwEAAAAA%3AiKsmU_1GdJO4ifGxVA30LzbPhFSMKexWc0gaJrDNxCf7D_i0ae9ym2ylWijp8dkpZOwp4ZC_Ivcr)) to calculate fugacity coefficients. Inputs are temperature (kelvin), pressure (bar), names (species names), and mol fraction of each species (n).

NOTE - I suggest not using the function gibbsG, because this calculates gibbs energy of formation of a gas phase species in a DIFFERENT convention than the aqueuous phase function (gibbsAQ). Therefore, you would get the wrong answer if you used gibbsG and gibbsAQ to calculate equilibrium constant of gas-aqueous reaction. You get the right answer when you use the functions gibbs and gibbsAQ.

