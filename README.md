### What this is

An event-driven simulation of an ideal gas in a box interacting with a movable piston. The piston is under the influence of a gravity-like force, $ F = - M g $. Gravity does not affect the particles.

### How to

In a nutshell:

1. Install homebrew
2. Install zig
3. Build the binary
4. Run it

If you don't have [homebrew](https://brew.sh) installed, install it by going to that page and running the install script:

```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Then clone this repo and `cd` into the directory. In one stroke:

```
brew install zig
git clone https://github.com/alecstein/ergsim.git
cd ergsim
zig build -Drelease
```

Note: if you get an error about the has upon building, this is a bug with the zig package manager. They're fixing it but here's a workaround.

#### Workaround 1
1. Delete the entire repo
2. `rm -rf ~/.cache/zig`
3. Follow the steps above

#### Workaround 2
1. `rm -rf ~/.cache/zig`
1. Open up `build.zig.zon`
3. Comment out the entire line with the hash in it
4. Run `zig build -Drelease`. This will output a (different) hash
5. Copy and paste that new hash back where the old hash was in `build.zig.zon`
6. Uncomment that line and run `zig build -Drelease`

Then to run the simulation,

`./zig-out/bin/ergsim <number-of-particles>`

https://github.com/alecstein/ergsim/assets/16236421/431e808b-3be6-4a03-b814-9a7a5ddd9a1d

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

* Successfully (I think) added gravity to the particles in the simulation
* Worked out some of the theory. Double factorials got me like !!
* Updated the code to be a little cleaner -- removed file writing for now to make it easy for @josh to check for bugs

##### Jan. 1

* Worked out some more optimizations
* A good way of unit testing to write code that you're sure is correct, count the collisions (for fixed initial conditions), then do the same thing for your optimized code. If the counts are the same, then you're probably ok. The counts will not be exactly the same due to floating point errors -- sometimes, collision times can vary in the 10th decimal place or so between optimizations. This shouldn't matter too much

##### Jan. 3

* Rewrote using zig-webui. Kinda weird and annoying to use, but we have liftoff. Total freedom from python and matplotlib
* Running this is as simple as doing `vis <particle number>` and letting 'er rip
* Archived old visuals [here](https://github.com/alecstein/phys/assets/16236421/cc232826-662c-405e-bda4-f017c42fab71) and [here](https://github.com/alecstein/phys/assets/16236421/e56aa705-617d-4697-a1a9-e7163f423df2)


#### Jan. 10

* Figured out D3.js and rewrote the visualization to be a little nicer. Lots of potential there. We can move stuff around now
* Rewrote the simulation in terms of a single parameter, mu
* Removed some code
* Cache old [visuals](https://github.com/alecstein/ergsim/assets/16236421/8d250379-dfd1-4092-82fc-6ace1afd737f)
