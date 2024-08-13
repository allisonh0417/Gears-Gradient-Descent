Gear's Method (1st to 5th Order) for Gradient Descent
Equations Used: https://www.chemecomp.com/Gear.pdf

Nesterov's Accelerated Gradient Descent:
https://cseweb.ucsd.edu/classes/sp24/cse291-e/papers/Nesterov/DifferentialML.pdf

* To initialize and generate matrices with random values, use 'init_func.m' before proceeding

Function List:
* NAG_Adaptive: Nesterov's with Adaptive Step Size
* NAG_Const: Nesterov's with Constant Step Size
* Gears_5th_old: Full Run of 5th Order Gear's with Adaptive Step Size (less stable)
* Gears_5th_new: Full Run of 5th Order Gear's with Adaptive Step Size w/ always decreasing error (more stable)
* Gears_5th_iter: Full Run of 5th Order Gear's with Adaptive Step Size w/ testing new step size every X iterations
* BDF1_Adaptive: 1st Order Gear's with Adaptive Step Size
