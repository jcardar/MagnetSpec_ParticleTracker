# Magnetic Spectrometer Particle Tracker
This project is a fully 3D, relativistic simulation of a magnetic spectrometer often used in Laser Wakefield Acceleration experiments. It simulates user-defined particle beams traversing through user defined magnetic dipole geometries, and landing on imaging screens as set up by the user. A genetic algorithm is implimented to manipulate input parameters of the genetic algorithm and converge on an local optimal setup for a given figure of merit.
UNDER DEVELOPMENT.
Particle tracker (.cpp) and visualization (.ipynb) tools for magnetic spectrometer, as well as genetic algorithm (.py) to run and interpret magnetic spectrometer data.

To make input deck for particle stepper code, from the repository base directory navigate to './data/analysis/' and open the Jupyter Notebook "Magnetic Spectrometer Interface". If you desire to use the genetic algorithm, the parameters for that will be set up in this Jupyter Notebook, too.

To do a single run of the magnetic spectrometer code run in terminal by changing directory to 'source_code' directory, executing the 'make' command, and then typing '.\run' in terminal. This can also currently be done in the first line of the developer Jupyter Notebook './data/analysis/analsis_mag_spec.ipynb'.

To use the genetic algorithm (under development), change directory to './source_code/', and run in terminal 'python geneticalgorithm_run.py'. Prompts will appear to guide the remaining setup of the genetic algorithm process. Data will be saved in './souce_code/save/', where the original input deck and the best input deck according to the figure of merit chosen will also be saved.

File organization:

./source_code/ - directory of .cpp and .py codes required to run the magnetic particle tracker and the genetic algorithm. 

./data/ - directory where data from the last run of the magnetic particle tracker code is stored, as well as the analysis subdirectory.

./analysis/ - folder where the Jupyter notebooks for setting up an input deck for the code is stored, as well as (in-development) visualization tools for checking the results of running the spectrometer code. Also in this directory is the input deck is stored for each run, as well as the best input deck of the last genetic algorithm run and .


Data output into 'data' directory.

Inside 'data' directory is analysis tool 'analysis_mag_spec.ipynb', this can run the code and import data. (Currently quite messy, am cleaning up).
<br>
<br>
NORMALIZATION:

Time           -> Average Classical Cyclotron Frequency (\omega_{c,0}) = 1

Velocity       -> Speed of light (c)   = 1

Position       -> (c / \omega_{c,0})    = 1

Charge         -> Particle charge (qe) = 1

Mass           -> Particle mass (me)   = 1

Magnetic field -> Average magnetic field (B_0) = 1
<br>
<br>
TO USE: (Will change!)<br>
<br>
Currently values in particle stepper must be hard-coded in.<br>
Some effort has been made to make sure user-inputted values are near the top of the main loop in the file "mag_spec_tracker.cpp".<br>
In particular, users may be interested in changing:<br>
<br>
  Time Properties:<br>
  del_time (line 32)        -time step for particle stepping algorithm<br>
  time_step_test (line 33)  -whether to half the time-step for each subsequent particle (=true) or not (=false)<br>
  particle_time_limit (line 150) -user-defined value for a time value to have the particle stepper stop (useful if particle is stuck in magnetic region doing gyroradius circles forever.)<br><br>
  
  beam properties (lines 47 - 55) -All beam properties for simulation (example given below)<br>
  <br>
  example beam properties:<br>
    int num_par          {1};<br>
    double charge        {-1.0};<br>
    double mass          {1.0};<br>
    double energy0       {10.0};  //Normalized Energy = gamma<br>
    double energy_spread {0.0};<br>
    ThreeVec initial_position(0.0, 0.0, 0.0);<br>
    ThreeVec initial_position_spread(0.0, 0.0, 0.0);<br>
    ThreeVec initial_angular_direction(0.0, M_PI/2.000000, M_PI/2.000000);<br>
    ThreeVec initial_angular_spread(0.0, 0.0, 0.0);<br>
    <br>
      
  magnet properties (lines 61 - 90)<br>
    - Properties of the magnets in the system. First define number of magnets, then change their properties in the switch loop. This will soon be changed to a for loop when code is made to read values from user-generated input deck.<br>
    <br>
  screen properties (lines 95 - 119)<br>
    - Properties of the screens in the system.<br>

