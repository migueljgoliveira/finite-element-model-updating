# ParamID-2.0
 PARAMeter IDentification 2.0

## Input

### INFO

`*INFO` defines the type of procedure and database of tests. 

#### Procedure
Three keywords are used to define the type of procedure. Only one keyword is required.

`*CALIBRATION` sets the procedure to calibration.

`*SIMULATION` sets the procedure to simulation. A single iteration is performed using the initial set of parameters.

`*SENSITIVITY` sets the procedure to sensitivity. A sensitivity study of variables is performed.

#### Database

`*TESTS` defines the database of tests. 

Data lines (required) to define the database of tests (separated by commas):

  1. Number of tests
  2. Name of tests
  3. Weight of tests

##### *Example*

  ```
  *INFO
  *CALIBRATION
  *TESTS
   2
   Test1, Test2
   0.5, 0.5
  *END INFO
  ```

### METHOD

`*METHOD` defines the type of inverse methods used in the procedure. 

##### Required parameters:

 * `method`
    * `method=<string>` defines the inverse method(s) to use.
      * `method=vfm`: virtual fields method.
      * `method=femu`: finite element model updating method.
      * `method=hybrid`: finite element model updating method and virtual fields method.

* `cost`
    * `cost=<string>` defines the type of cost function. 
      * `cost=1`: difference.
      * `cost=2`: weighted difference.
      * `cost=3`: quadratic weighted difference.

#### Finite Element Model Updating

`*FEMU` defines settings of the finite element model updating method.

###### Optional parameters:
  * `norm`
    * `norm=<bool>` defines the type of normalisation used in the objective function.
      * `norm=false`: difference between experimental and numerical strain components and force not normalised.
      * `norm=true`: difference between experimental and numerical strain components normalised by maximum strain of all time increments, points, and components; and difference between experimental and numerical force normalised by maximum force of all time increments (default).

  * `residual`
    * `residual=<integer>` defines the type of residuals used in the objective function.
      *  `residual=1`: residuals for each time increment of strain plus force.
      *  `residual=2`: residuals for each time increment of strain and force.
      *  `residual=3`: residuals for each time increment, point and component of strain and for each time increment of force (default).

  * `cpus`
    * `cpus=<integer>` defines the number of CPUs used in the finite element analysis (default is 1).

###### Optional parameter when `method=hybrid`:
  * `weight`
      * `weight=<float>` defines the weight of the finite element model updating. This value must range between 0 and 1. Default is 0.5.

#### Virtual Fields Method

`*VFM` defines settings of the virtual fields method.

##### Optional parameters:
  * `norm`
    * `norm=<bool>` defines the type of normalisation used in the objective function.
      * `norm=false`: difference between internal and external virtual work not normalised (default).
      * `norm=true`: difference between internal and external virtual work normalised by external virtual work

  * `virtualfields`
    * `virtualfields=<string>` defines the type of virtual fields.
      * `virtualfields=manual`: manual virtual fields (default.
      * `virtualfields=auto`: sensivity-based virtual fields (not yet available).


##### Optional parameter when `method=hybrid`:
  * `weight`
      * `weight=<float>` defines the weight of the virtual fields method. This value must range between 0 and 1. Default is 0.5.

##### *Example*

  ```
  *METHOD, method=hybrid, objective=3
  *FEMU, norm=true, residual=3, cpus=2, weight=0.6
  *VFM, norm=false, weight=0.4
  *END METHOD
  ```

### OPTIMISATION

`*OPTIMISATION` defines settings of the optimisation procedure.

  * `parallel`
      * `parallel=<bool>` defines the type of mutation strategy.
        * `parallel=false`: evaluation of tests is performed in sequence (default).
        * `parallel=true` evaluation of tests is performed in parallel.
#### Differential Evolution

`*DE` defines settings of the differential evolution algorithm.

###### Optional parameters:

  * `strategy`
    * `strategy=<string>` defines the type of mutation strategy.
      * `strategy=best1bin`: (default)
      * `strategy=best1exp`:
      * `strategy=rand1exp`:
      * `strategy=randtobest1exp`:
      * `strategy=currenttobest1exp`:
      * `strategy=best2exp`:
      * `strategy=rand2exp`:
      * `strategy=randtobest1bin`:
      * `strategy=currenttobest1bin`:
      * `strategy=best2bin`:
      * `strategy=rand2bin`:
      * `strategy=rand1bin`:

  * `polish`
    * `polish=<bool>` defines if the L-BFGS-B method is used to polish the best population member at the end of the optimisation.
      * `polish=false`: best population member is not polished at the end of the optimisation (default).
      * `polish=true`: best population member is polished at the end of the optimisation.

  * `init`
    * `init=<string>` defines the type of population initialization.
      * `init=latinhypercube`: population is generated in a near-random way, maximising the search space coverage (default).
      * `init=random`: population is randomly generated.

  * `updating`
    * `updating=<string>` defines the type of scheme used to update the best solution vector.
      * `updating=immediate`: best solution vector is continuously updated within a single generation; it can lead to faster convergence, but it is not compatible with parallelization (default).
      * `updating=deferred`: best solution vector is updated once per generation; it is compatible with parallelization.

#### Dual Annealing

`*DA` defines settings of the dual-annealing algorithm.

#### Least-Squares

`*LS` defines settings of the least-squares algorithms.

#### Minimize

`*MINI` defines settings of the minimize algorithms.

##### *Example*

  ```
  *OPTIMISATION, algorithm=de, parallel=true
  *DE, strategy=best1bin, polish=false, init=latinhypercube, updating=immediate
  *END OPTIMISATION
  ```

### MATERIAL
