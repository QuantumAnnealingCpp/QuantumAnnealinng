## Simulation of quantum annealing and its application for traveling salesman problem

Suzuki-Trotter formalism maps the quantum problem to a classical one. It is essentially a method to transform a d-dimensional
quantum Hamiltonian into a (d+1)-dimensional effective classical Hamiltonian. To obtain the effective hamiltonian one starts with Transverse Ising
Hamiltonian:

![equation](http://latex.codecogs.com/gif.latex?H%20%3D%20-%5CGamma%20%5Csum%5EN_%7Bi%3D1%7D%20%5Csigma%5Ex_i%20-%20%5Csum_%7B%28i%2Cj%29%7D%20J_%7Bij%7D%20%5Csigma%5Ez_i%20%5Csigma%5Ez_j)

After applying the Suzuki-Trotter decomposition one obtains the following effective hamitonian. It is used as a goal function which is optimized using the Metropolis algorithm:

![equation](http://latex.codecogs.com/gif.latex?H_%7Beff%7D%28%5Csigma%29%20%3D%20%5Csum%5EN_%7B%28i%2Cj%29%7D%20%5Csum%5EM_%7Bk%3D1%7D%20%5Cbigg%5B-%5Cfrac%7BJ_%7Bij%7D%7D%7BM%7D%5Csigma_%7Bik%7D%5Csigma_%7Bjk%7D%20-%20%5Cfrac%7B%5Cdelta_%7Bij%7D%7D%7B2%5Cbeta%7D%20%5Cln%20%5Ccoth%5CBig%28%5Cfrac%7B%5Cbeta%20%5CGamma%20%7D%7BM%7D%5CBig%29%5Csigma_%7Bik%7D%5Csigma_%7Bik&plus;1%7D%20%5Cbigg%5D)

Setting the following parameters one can obtain a phase transition for the modulus of the magnetization of the system:
```
double T		= 0.1; //thermodynamical temperature
double gamma	    = 1.5;//transverse field strenght factor
double kB		= 1.0; //Boltzman constan
int N = 8; //N quantum spins
int M			= 32; //additional dimension of spins
int NT			= 1000000; //liczba kroków czasowych, w których losowany jest jeden spin
int snapNT		= NT/1000; //a value of magnetization is saved every 1000 steps of the simulation
int Gsteps		= 20; //number of evaluation points of gamma for the | <s> | = f(gamma) plot
```

The point at which the average magnetization |M| > 0 starts approximate at gamma/J = 1.2 which is in agreement with theory found in the references listed below.

![alt text](https://github.com/QuantumAnnealingCpp/QuantumAnnealinng/blob/develop/Graphs/M16DOBRYWYKRES.png?raw=true)




After rewritng the problem in a way which is useful for solving the traveling salesman problem we obtain the goal funtion which can be interpreted as a length of path. However, one gets an
additional term which takes into account "quantum effects":

![equation](http://latex.codecogs.com/gif.latex?L%20%3D%20%5Csum%5EN_%7B%28ij%29%7D%5Csum%5EM_%7Bk%3D1%7D%20%5Cbigg%5B%20d_%7Bij%7Dn_%7Bi%28t%29j%7D%20n_%7Bj%28t&plus;1%29k%7D%20&plus;%20%5Cfrac%7B%5Cdelta_%7Bi%28t%29j%7D%7D%7B2%5Cbeta%7D%20%5Cln%5Ccoth%5CBig%28%5Cfrac%7B%5Cbeta%20%5CGamma%7D%7BM%7D%5CBig%29%20n_%7Bi%28t%29k%29n_%7Bi%28t%29k&plus;1%7D%7D%20%5Cbigg%5D)

Setting the following parameters one can obtain the following results
```
double T		= 0.05; //thermodynamical temperature
double gamma	    = 1.5;//transverse field strenght factor
double kB		= 1.0; //Boltzman constan
int N = 16;//N = 32 //N quantum spins
int M			= 16; //additional dimension of spins
int NT			= 1000000; //liczba kroków czasowych, w których losowany jest jeden spin
int Gsteps		= 20; //number of evaluation points of gamma for the | <s> | = f(gamma) plot
```

![alt text](https://github.com/QuantumAnnealingCpp/QuantumAnnealinng/blob/develop/Graphs/M16N16bestpath.png?raw=true)

![alt text](https://github.com/QuantumAnnealingCpp/QuantumAnnealinng/blob/develop/Graphs/M16N32bestpath.png?raw=true)

One does not have certainty whether it is the shortes possible path. However it is a local mnimum of the goal function.

Points are chosen at rendom (in the unit square), hence your output may be different.

### Prerequisites

To open the solution VS 2017 edition is required. No special software is required. 


## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Authors

* **Paweł Kukliński**  - [pawelkuk](https://github.com/pawelkuk)
* **Jakub Sobolewski**  - [WelcomeToMyVirtualHome](https://github.com/WelcomeToMyVirtualHome)

