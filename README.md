# 2D ball balancing

Work carried out as part of the Mechatronics Lab course on a ball balancing
system with a tiltable table. In particular, the goal was to design, develop and test some suitable control
strategies to stabilize a small ball in a certain position or along a certain trajectory on a moving
flat surface, by exploiting the feedback on the position of the ball and actuators to tilt the table.  

The project was implemented:
1. Through Simulink simulations, for comparison with experimental results
2. On a provided experimental bench, using Simulink+PoliArd (Arduino-based architecture) and touchpad for position feedback
3. On a newly constructed bench, using Python (for image processing and data handling) and Simulink (for controller deployment) along with TI Launchpad microcontroller and camera for position feedback
