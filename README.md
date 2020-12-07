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

./analysis/ - folder where the Jupyter notebooks for setting up an input deck for the code is stored, as well as (in-development) visualization tools for checking the results of running the spectrometer code. Also in this directory is the input deck is stored for each run, as well as the best input deck of the last genetic algorithm run and the original input deck for the last genetic algorithm run.

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

In particular, users may be interested in changing some values in the mag_spec_tracker.cpp code:<br>
<br>
  Time Properties:<br>
  del_time            - normalized time step for particle stepping algorithm<br>
  time_step_test      - whether to half the time-step for each subsequent particle (=true) or not (=false)<br>
  particle_time_limit - user-defined value for a time value to have the particle stepper stop (useful if particle could be stuck in magnetic region doing gyroradius circles forever.)<br><br>
  

