### What this is

An event-driven simulation of an ideal gas in a box interacting with a movable piston. The piston is under the influence of a gravity-like force, $ F = - M g $. Gravity does not affect the particles.

### How to

If you don't have zig installed, do  `brew install zig`.  If you don't have brew installed, install it. 

Then do:

`zig build`

This will run the simulation but won't output anything. It's just a skeleton for producing outputs and can be modified as you need.

`./zig-out/bin/ergsim <number-of-particles> <time>`

will run the simulation.

I wrote some python code (not included) to make some histograms. These are with particle gravity turned off:

https://github.com/alecstein/phys/assets/16236421/cc232826-662c-405e-bda4-f017c42fab71

https://github.com/alecstein/phys/assets/16236421/e56aa705-617d-4697-a1a9-e7163f423df2


### Diary

##### Dec. 24

* Fixed two bugs related to particle collisions 
* Removed discriminant positivity check (should always be positive, panic otherwise)
* Stopped using `eval` to parse floats in python, this was reading in the numbers wrong and creating weird graphs
* Fixed more bugs in the zig code that got there via ChatGPT
* Rearranged some things for clarity

##### Dec. 25

"I am not impressed by [Zig's] performance" (Georges St. Pierre voice).

* Tried to use a priority queue for speed gainz but the problem is you have to update the entire queue every time there's a collision with the piston (this potentially reorganizes all of the events.) There's probably a math insight you can use to speed this up but I couldn't think of it.
* Tried to use `@Vector`s but this made everything hella slow instead of fast
* Putting things into structs for organizational clarity slows things down quite a bit

##### Dec. 30-31

* Successfully (I think) added gravity to the simulation
* Worked out some of the theory
* Updated the code to be a little cleaner -- removed file writing for now to make it easy for @josh to check for bugs
