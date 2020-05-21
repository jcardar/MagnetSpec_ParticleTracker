# Magnetic Spectrometer Particle Tracker
UNDER DEVELOPMENT.
The particle tracker (.cpp) and visualization (.ipynb) tools for magnetic spectrometer. 

Most up-to-date version is "continuous field" version.
Code run by changing directory to 'source_code' directory and executing the 'make' command.
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

