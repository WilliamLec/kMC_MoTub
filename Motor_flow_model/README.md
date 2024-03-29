## Motor Flow Model

Here you will find the script used to simulate with a kMC algorithm a flow of processive motors along a lattice.

The motor is considered to be two-headed and walking with a hand-over-hand fashion. 
The lattice is composed of m lines of n sites and with periodic conditions (ie. when the motor exits the lattice, it re-enters on the other side).
In the following clip (real time), you can see a flow of motor with a velocity of 10 sites/s. 
The color code is the following: blue = motor's front head, lemon = motor's rear head, ligth green = lattice and dark green = lattice vacancy.

![](../gifs/dynein_flow_with_defect.gif)

The kMC algorithm was highly optimized thanks to two tricks: 
1. after each reaction the lattice update is restricted to the neighborhood of the reaction site (ie. a 7x3 rectangle) 
2. the use of mirror lists to create redundance gives us access at any time of the position of any active reaction 

## Motor Reactions

This model is inspired from the one developped by [Rank *and al.*](https://www.sciencedirect.com/science/article/pii/S0006349518308269) 
1. motor attach, active when 2 adjacent sites are vacant, rate constant = kp
2. walking, active when the motor in front of the motor is vacant, rate constant = kw 
3. motor detachment, active for any motor present on the lattice, rate constant = km 
4. modified motor detachment, active when the motor has only one head bound to the lattice (ie. facing a lattice vacancy), rate constant = kM

## Motor Characteristics

To tune motor characteristics, you can play with: 
- kw = mean motor velocity in site/s in very low density phase
- km = inverse of mean motor binding time 
- kw/km = motor processivity = mean distance travelled by a single motor before detachment of the lattice
- kM/km = detachment coefficient = impact how lattice vacancy disturbs the motor flux / fixed to 1 you get a traffic jam but to 100 you don't 

