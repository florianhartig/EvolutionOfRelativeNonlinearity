;==================================================================
; GLOBAL VARIABLES
; global variables are either defined here, or in the GUI.

globals
[
  b-datapoints             
  disturbance-strength
  init-rate 
  population
  population1
  population2
  temporal-disturbance
  time-per-1000-steps
  stratcounter
  spreadcounter
  agecounter
  bcounter
  reproductioncounter
  invasioncounter
  break
  lastvalue                  
  individuals-die
]

;==================================================================
; AGENT DEFINITIONS
; (defines the state variables of the model)

turtles-own [
  individuals    ; number of individuals represented by the turtle
  b-strategy     ; density compensation b-strategy
  rate           ; growth rate
  mutants        ; temporary variable 
  interaction    ; niche overlap with other species
  age
]

; derives "sub-species"
breed [species1]
breed [species2]
breed [invadors]

; invadors have additional property "growth"
invadors-own[
  growth         ; for summing up average growth rates
]


;==================================================================
; SETUP PROCEDURES
; (define the initialization of the model, different options are avaiable 
; that are selected by the global variables init-b and species)

to setup 
  clear-turtles 
  clear-patches
  clear-all-plots 
  reset-ticks  
  reset-timer
  set stratcounter 0
  set spreadcounter 0
  set agecounter 0
  set bcounter 0
  set reproductioncounter 0
  set invasioncounter 0
  set break False 
  set individuals-die True
  
  ;======================================
  ; set the initial size of the population
  ; pop is divided evenly withing species
  let popsize capacity / species
  ;; times chosen initial popsize
  ;; random value of +/- 0.05% is added to avoid artefacts e.g. from setting the population right at capacity, or withing the cycles
  set popsize ceiling(initial-population * popsize * (0.95 + random-float 0.1)) 
    
  ;======================================
  ; set initial population sizes
  if init-b = "flat" 
  [
    create-species1 popsize
    [
      set b-strategy random-float (15 - cutoff) + cutoff
      set individuals 1
      set interaction interaction12
      set rate initial-rate-1
    ]
    if species > 1 [
      create-species2 popsize
      [
        set b-strategy random-float (15 - cutoff) + cutoff
        set individuals 1
        set interaction interaction21
        set rate initial-rate-2
      ]
    ]
  ]
  if init-b = "single" [
    create-species1 1
    [
      set b-strategy initial-strategy-1
      set individuals popsize
      set interaction interaction12
      set rate initial-rate-1
    ]
    if species > 1 [
      create-species2 1
      [
        set b-strategy initial-strategy-2
        set individuals popsize
        set interaction interaction21
        set rate initial-rate-2
      ]
    ]
  ]
  if init-b = "invasion"
  [
    ifelse species > 1 [
      user-message("only one species allowed with invasion")
    ][ 
      create-species1 1
      [
        set b-strategy initial-strategy-1
        set individuals popsize
        set rate initial-rate-1
        set interaction interaction12
      ]
    ] 
  ]    
  if init-b = "double"
  [
    create-species1 1
    [
      set b-strategy initial-strategy-1
      set individuals round(popsize / 2)
      set rate initial-rate-1
      set interaction interaction12
    ]
    create-species2 1
    [
      set b-strategy initial-strategy-2
      set individuals round(popsize / 2)
      set rate initial-rate-2
      set interaction interaction21
    ]   
  ]
  if visualization [visualize]
end


;==================================================================
; GO PROCEDURES 
; (define the subprocedures that are called at each time step;
; different procedures are available for different model experiments
; standard is the procedure that is called "go")



; go is the "normal" model scheduling
to go
  update-population
  tick
  if visualization [visualize]
end

; for invasion calculations
to invasion-fitness-go
  update-population
  tick
  ifelse(demographic-stochasticity)[
    ;This part is hard-coded. Invasions appear every 600 years, with 100 years time inbetween for the resident population to stabilize
    ;
    ;
    if (((ticks + 600) mod 600) = 100)[  
        create-species2 1
        [
          set b-strategy initial-strategy-2
          set individuals 3
          set interaction interaction21
          set rate initial-rate-2
        ]
    ]
    if ((((ticks + 600) mod 600) = 0) and (ticks > 100))[
        ifelse (count species2 = 0)
          [set invasioncounter (invasioncounter + 1) ]
          [ ask species2 [die]         ]
        if (count species1 = 0)
        [
          create-species1 1
          [
            set b-strategy initial-strategy-1
            set individuals initial-population * capacity
            set interaction interaction12
            set rate initial-rate-1
          ]  
        ]
    ]     
  ]
  [
    let invasion-time 500
    if (ticks = invasion-time)[
      create-invadors 1
      [
        set b-strategy initial-strategy-2
        set individuals 1
        set growth 0
        set interaction interaction21
        set rate initial-rate-2
      ]
    ]
    if (ticks > invasion-time)[
      ask invadors[
        set growth (growth + ln(individuals))
        set individuals 1 
        set invasioncounter invasioncounter + 1
      ]
    ]  
  ]
  if visualization [visualize]
end

; for IBM ESS calculations
to ess-record
  go
  ask turtles [set age (age + 1)]
  set population sum [individuals] of turtles
  if (ticks > 100000)[
    set stratcounter (stratcounter + number-of-strategies) 
    set spreadcounter (spreadcounter + average-b-spread)
    set agecounter (agecounter + sum [age] of turtles)
    set bcounter (bcounter + average-b-strategy)
  ]  
end



;==================================================================
; SUB PROCEDURES
; subprocedures are called by the different go procedures



; update population is the main procedure that describes the biological 
; processes that act on the population at each time step
to update-population
  set population sum [individuals] of (turtle-set species1 species2)
  ask turtles
  [
    calculate-reproduction
    if (disturbance > 0)[
      set individuals (( 1 - disturbance) * individuals ) 
    ]
    if demographic-stochasticity [set individuals random-poisson(individuals)]
    if (mutate-b or mutate-rate or mutate-alpha) [mutate]
    ifelse individuals-die 
      [if (individuals <= 0)[die]]
      [if (individuals <= 1)[set individuals 1]]
  ]
end

;    set population1 sum [individuals] of species1
;    set population2 sum [individuals] of species2

;    ask species1[ set individuals random-poisson ( rate * individuals / (1 + ( rate - 1) * ( (population1 + interaction * population2) / capacity) ^ b-strategy )) ]
;    ask species2[ set individuals random-poisson ( rate * individuals / (1 + ( rate - 1) * ( (population2 + interaction * population1) / capacity) ^ b-strategy ))]




; calculate-reproduction calculates the growth rate for different growth models
to calculate-reproduction
  if (function = "MSS") [ set individuals (  rate * individuals / (1 + ( rate - 1) * (population / capacity) ^ b-strategy))]
  if (function = "Hassel") [ set individuals ( rate * individuals / (1 + ( rate ^ ( 1 /  b-strategy) - 1) * (population / capacity)) ^ b-strategy )]
  if (function = "Ricker") [ set individuals  ( individuals * exp ( rate * (1 - (population / capacity)^ b-strategy ))) ]
  if (function = "Gompertz") [ set individuals  (  individuals * exp ( rate * (1 - (ln (population) / ln(capacity)))^ b-strategy )) ]
  if (function = "Experimental") [ set individuals  ( (rate - 2 * exp(-(( b-strategy - 2.3) / 0.5) ^ 2 )) * individuals / (1 + ( rate - 1) * (population / capacity) ^ b-strategy ))  ]
  if (function = "MSS-Mod") [ 

    let lower 1.0
    let upper 5.5
        
    let scaling 4
    
    ifelse( b-strategy < lower or b-strategy > upper) 
      [set individuals (  rate * individuals / (1 + ( rate - 1) * (population / capacity) ^ b-strategy))]
      [ifelse((population / capacity) < 1)
        [ set individuals (  rate * individuals / (1 + ( rate - 1) * (population / capacity) ^ (lower + (b-strategy - lower)^(1 / scaling)/(upper - lower)^(1 / (scaling - 1)))))]  
        [ set individuals (  rate * individuals / (1 + ( rate - 1) * (population / capacity) ^ (lower + (b-strategy - lower)^ scaling /(upper - lower)^(scaling - 1))))]
      ]
  ]
end



; mutate implements the mutation
to mutate
  if (individuals > 0)
      [ 
        set mutants random-poisson(individuals * mutation-probability / capacity)
        if mutants > individuals [set mutants individuals]
        set individuals (individuals - mutants)
        hatch mutants 
        [ 
          set individuals 1
          set age 0
          if mutate-b 
          [ 
            set b-strategy (b-strategy + random-normal 0 mutation) 
            if b-strategy < cutoff [ set b-strategy  cutoff]
          ]
          if mutate-rate
          [ 
            set rate (rate + random-normal 0 mutation) 
            if rate < 1 [ set rate 1]
          ] 
          if mutate-alpha 
          [ 
            set interaction (interaction + random-normal 0 mutation) 
            if interaction < 0 [ set interaction 0]
            if interaction > 1 [ set interaction 1]
          ] 
        ] 
      ]   
end 

; visualize implements the visualization of the model
to visualize
  set population sum [individuals] of turtles
  if population > 0 [
   set population1 sum [individuals] of species1
   set population2 sum [individuals] of species2
  
   set-current-plot "Population size"
   if ticks > 200 [set-plot-x-range (ticks - 150)  (ticks + 50) ]
   set-current-plot-pen "species1"
   plot population1
   if (species > 1 )[
     set-current-plot-pen "species2"
     plot population2
   ]
   set-current-plot-pen "capacity"
   plot capacity
   
   
   
   if true [
   set-current-plot "b-histogram population"
   clear-plot 
   ;set-plot-y-range 0 2 * capacity
   set-current-plot-pen "frequency"
   set-histogram-num-bars 300
   
   let pop [] let a 0 let b 0
   (foreach ([b-strategy] of turtles) ([individuals] of turtles) [ set a ?1 set b ?2 (set pop sentence pop (n-values b [a]) )])
   histogram pop
   ]
   

    set-current-plot "b-development"
    if ticks > 500 [ set-plot-x-range (ticks - 500)  (ticks + 50) ]
    set-current-plot-pen "species1"
    plot average-b-strategy1
    if (species = 2) [
      set-current-plot-pen "species2"
      plot average-b-strategy2
    ]
    
    set-current-plot "r-development"
    if ticks > 500 [ set-plot-x-range (ticks - 500)  (ticks + 50) ]
    set-current-plot-pen "species1"
    plot average-rate-1
    if (species = 2) [
      set-current-plot-pen "species2"
      plot average-rate-2
    ]
    

    set-current-plot "genetic-spread"
    if ticks > 500 [ set-plot-x-range (ticks - 500)  (ticks + 50) ]
    set-current-plot-pen "b"
    plot average-b-spread
    set-current-plot-pen "rate"
    plot average-rate-spread


;    set-current-plot "invasion-fitness"
;    set-current-plot-pen "fitness"
;    plot (sum [growth] of invadors)

   
   ; visualize rate histogram
   if true [
     set-current-plot "r-histogram population"
     clear-plot 
     set-plot-y-range 0 2 * capacity
     set-plot-x-range 0 ceiling (([rate] of max-one-of turtles [rate]) / 20 ) * 20
     set-current-plot-pen "rate"
     set-histogram-num-bars 300
     let pop [] let a 0 let b 0
     (foreach ([rate] of turtles) ([individuals] of turtles) [ set a ?1 set b ?2 (set pop sentence pop (n-values b [a]) )])
     histogram pop
   ]
   
   
;   if true [
;    set-current-plot "alpha"
;    if ticks > 5000 [ set-plot-x-range ticks - 5000  ticks + 50 ]
;    set-current-plot-pen "species1"
;    plot average-interaction1
;    if (species > 1) [
;      set-current-plot-pen "species2"
;      plot average-interaction2
;    ]
;  ]
 ]
end


;==================================================================
; Reporters 
; (in NL language, reporters are functions that calculate (report) some 
; aggregated variable from the state variables of the model)


to-report average-strategies
  report (stratcounter / (ticks - 100000))
end

to-report average-spread
  report (spreadcounter / (ticks - 100000))
end

to-report average-age
  report (agecounter / (ticks - 100000))
end

to-report average-b
  report (bcounter / (ticks - 100000))
end

to-report average-b-strategy
  ifelse (population > 0) 
    [ report (sum ( map [?1 * ?2] [b-strategy] of turtles  [individuals] of turtles )) / population ]
    [ report -1 ] 
end

to-report average-b-strategy1
  ifelse (population1 > 0) 
    [
      report (sum ( map [?1 * ?2] [b-strategy] of species1  [individuals] of species1 )) / population1 ]
    [ report -1 ] 
end

to-report average-b-strategy2
  ifelse (population2 > 0) 
    [ report (sum ( map [?1 * ?2] [b-strategy] of species2  [individuals] of species2 )) / population2 ]
    [ report -1 ] 
end


to-report average-rate
  ifelse (population > 0) 
    [ report (sum ( map [?1 * ?2] [rate] of turtles  [individuals] of turtles )) / population ]
    [ report -1 ] 
end

to-report average-rate-1
  ifelse (population1 > 0) 
    [ report (sum ( map [?1 * ?2] [rate] of species1  [individuals] of species1 )) / population1 ]
    [ report -1 ] 
end

to-report average-rate-2
  ifelse (population2 > 0) 
    [ report (sum ( map [?1 * ?2] [rate] of species2  [individuals] of species2 )) / population2 ]
    [ report -1 ] 
end

to-report average-b-spread
  ; requires population to be updated
  ifelse (population > 0) [ 
      let average average-b-strategy
      report (sum ( map [abs(?1 - average) * ?2] [b-strategy] of turtles  [individuals] of turtles )) / population 
    ][ 
      report -1 
    ] 
end

to-report average-rate-spread
  ifelse (count turtles > 0) [ 
      let average average-rate
      report (sum ( map [abs(?1 - average) * ?2] [rate] of turtles  [individuals] of turtles )) / population 
    ][ 
      report -1 
    ] 
end

to-report number-of-strategies
  report count turtles
end

to-report average-interaction1
  ifelse (count species1 > 0) 
    [
      report (sum ( map [?1 * ?2] [interaction] of species1  [individuals] of species1 )) / (sum [individuals] of species1) ]
    [ report -1 ] 
end

to-report average-interaction2
  ifelse (count species2 > 0) 
    [ report (sum ( map [?1 * ?2] [interaction] of species2  [individuals] of species2 )) / (sum [individuals] of species2) ]
    [ report -1 ] 
end


to-report critical-b
  ifelse break [report (mean [b-strategy] of turtles)][report 0]
end
  
  
; reporter for branching  
to-report occurrences [x the-list]
  report reduce [ifelse-value (abs (?2 - x) < resolution / 2) [?1 + 1] [?1]] (fput 0 the-list)
end






;==================================================================
; OLD MAIN PROCEDURES
; not used for current model results

to strategy-count-go
  update-population
  tick
  if ((ticks > 10000) and ((ticks mod 20) = 0)) [
     set stratcounter (stratcounter + number-of-strategies)
     set invasioncounter (invasioncounter + 1)
  ]
end

to critical-go
  update-population
  tick
  if ((ticks mod 50) = 48) [set lastvalue  sum [individuals] of turtles]
  if ((ticks mod 50) = 49) [
    ifelse (abs (lastvalue - (sum [individuals] of turtles)) > 5)[set break True][ask turtles [set b-strategy  (b-strategy + b-strategy / 50)]]
  ]
  if visualization [visualize]
end




to branching-setup
 reset-timer
 setup
 set b-datapoints ceiling(range / resolution)
 file-open (word "output/rate-" initial-rate-1 "_disturbance-" disturbance "startingpoint-" initial-strategy-1 ".dat " )
 file-print (word "# File: output/rate-" initial-rate-1 "_disturbance-" disturbance "startingpoint-" initial-strategy-1 ".dat " )
 file-print (word "# Time steps between data recording: " data-distance)
 file-print (word "# Histogramm per time step: Range 0-" range " datapoints = " b-datapoints )
 let raster n-values b-datapoints [(? * resolution) + resolution / 2 ]
 file-print raster
end

to branching-go
  update-population
  if (ticks = 30000)[set mutation-probability 0.3] 
  tick
  if ((ticks mod data-distance) = 0) [
    ;visualize
    ;export-plot "b-b-strategy" (word "output/rate-" initial-rate "_disturbance-" disturbance "startingpoint-" initial-b-strategy-1 "_plot-" (floor ticks / data-distance) )
    let a 0 let b 0
    
    if mutate-b[    
      let pop [] set a 0 set b 0
      (foreach ([b-strategy] of turtles) ([individuals] of turtles) [ set a ?1 set b ?2 (set pop sentence pop (n-values b [a]) )])
      let raster n-values b-datapoints [(? * resolution) + resolution / 2 ]
      let pop-histogram (map [occurrences ? pop] raster)
      file-print pop-histogram
    ]
    if mutate-rate [
      let pop-rate [] set a 0 set b 0
      (foreach ([rate] of turtles) ([individuals] of turtles) [ set a ?1 set b ?2 (set pop-rate sentence pop-rate (n-values b [a]) )])
      let raster n-values b-datapoints [(? * resolution) + resolution / 2 ]
      let pop-rate-histogram (map [occurrences ? pop-rate] raster)
      file-print pop-rate-histogram
    ]
    visualize
  ]
end

to branching-end
  file-print (word "# Experiment ended after: " timer)
  file-close-all 
end
@#$#@#$#@
GRAPHICS-WINDOW
1480
410
1725
527
1
1
28.7
1
10
1
1
1
0
1
1
1
-1
1
-1
1
0
0
0
ticks
30.0

BUTTON
6
10
111
55
setup
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
5
110
110
160
Go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
385
420
545
540
b-histogram population
b histogram
Frequency
0.0
20.0
0.0
10.0
true
false
"" ""
PENS
"frequency" 1.0 1 -16777216 true "" ""

SLIDER
5
260
190
293
initial-population
initial-population
0
1
1
0.01
1
NIL
HORIZONTAL

PLOT
385
10
875
260
Population size
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"species1" 1.0 0 -16777216 true "" ""
"capacity" 1.0 0 -2674135 true "" ""
"species2" 1.0 0 -10899396 true "" ""
"niche1" 1.0 2 -16777216 true "" ""
"niche2" 1.0 2 -10899396 true "" ""

SLIDER
225
355
365
388
mutation
mutation
0
1
0.5
0.001
1
NIL
HORIZONTAL

SLIDER
185
470
366
503
disturbance
disturbance
0
2
0.05
0.01
1
NIL
HORIZONTAL

CHOOSER
5
210
190
255
init-b
init-b
"flat" "invasion" "single" "double"
2

SLIDER
225
390
365
423
mutation-probability
mutation-probability
0
1
0.3
0.01
1
NIL
HORIZONTAL

PLOT
550
420
710
540
r-histogram population
rate histogram
frequency
0.0
200.0
0.0
10.0
false
false
"" ""
PENS
"rate" 1.0 1 -2674135 true "" ""

SLIDER
225
320
365
353
cutoff
cutoff
-3
3
0.17
0.01
1
NIL
HORIZONTAL

SWITCH
225
250
365
283
mutate-rate
mutate-rate
1
1
-1000

INPUTBOX
5
300
95
360
initial-strategy-1
4
1
0
Number

INPUTBOX
100
300
190
360
initial-strategy-2
6.0496475
1
0
Number

INPUTBOX
100
365
190
425
initial-rate-2
5
1
0
Number

INPUTBOX
5
365
95
425
initial-rate-1
5
1
0
Number

INPUTBOX
115
110
207
170
capacity
1000
1
0
Number

BUTTON
5
60
110
105
Step
go\n
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
225
215
365
248
mutate-b
mutate-b
0
1
-1000

BUTTON
1480
275
1575
308
NIL
branching-go\n
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1480
240
1575
273
NIL
branching-setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
1480
40
1630
71
For data acquisation and testing only\n
13
0.0
1

INPUTBOX
1480
90
1565
150
range
7
1
0
Number

INPUTBOX
1570
90
1655
150
resolution
0.05
1
0
Number

INPUTBOX
1480
155
1565
215
data-distance
100
1
0
Number

MONITOR
315
10
378
55
NIL
ticks
17
1
11

PLOT
385
270
620
415
b-development
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"species1" 1.0 0 -16777216 true "" ""
"species2" 1.0 0 -2674135 true "" ""

PLOT
715
420
875
540
genetic-spread
NIL
NIL
0.0
4.0
0.0
4.0
true
true
"" ""
PENS
"b" 1.0 0 -16777216 true "" ""
"rate" 1.0 0 -2674135 true "" ""

CHOOSER
115
10
207
55
species
species
1 2
0

SLIDER
5
470
131
503
interaction12
interaction12
-1
1
1
0.01
1
NIL
HORIZONTAL

SLIDER
5
509
131
542
interaction21
interaction21
-1
1
1
0.01
1
NIL
HORIZONTAL

SWITCH
225
285
365
318
mutate-alpha
mutate-alpha
1
1
-1000

PLOT
630
270
875
415
r-development
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"species1" 1.0 0 -16777216 true "" ""
"species2" 1.0 0 -2674135 true "" ""

SWITCH
270
60
380
93
visualization
visualization
1
1
-1000

CHOOSER
115
60
207
105
function
function
"MSS" "Hassel" "Ricker" "Gompertz" "Experimental" "MSS-Mod"
5

BUTTON
1645
235
1750
268
NIL
ess-record\n
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1645
270
1750
303
NIL
invasion-fitness-go\n
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
184
507
367
540
Demographic-Stochasticity
Demographic-Stochasticity
0
1
-1000

BUTTON
1645
305
1750
338
NIL
critical-go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1480
310
1575
343
NIL
branching-end
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
10
185
160
203
Initialization parameters
13
0.0
1

TEXTBOX
225
185
375
203
Mutation parameters
13
0.0
1

TEXTBOX
10
445
160
463
Limited niche-overlap
13
0.0
1

TEXTBOX
185
445
375
476
Disturbance and demogr. stoch.
13
0.0
1

@#$#@#$#@
## WHAT IS IT?

This section could give a general understanding of what the model is trying to show or explain.

## HOW IT WORKS

This section could explain what rules the agents use to create the overall behavior of the model.

## HOW TO USE IT

This section could explain how to use the model, including a description of each of the items in the interface tab.

## THINGS TO NOTICE

This section could give some ideas of things for the user to notice while running the model.

## THINGS TO TRY

This section could give some ideas of things for the user to try to do (move sliders, switches, etc.) with the model.

## EXTENDING THE MODEL

This section could give some ideas of things to add or change in the procedures tab to make the model more complicated, detailed, accurate, etc.

## NETLOGO FEATURES

This section could point out any especially interesting or unusual features of NetLogo that the model makes use of, particularly in the Procedures tab.  It might also point out places where workarounds were needed because of missing features.

## RELATED MODELS

This section could give the names of models in the NetLogo Models Library or elsewhere which are of related interest.

## CREDITS AND REFERENCES

This section could contain a reference to the model's URL on the web if it has one, as well as any other necessary credits or references.
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 5.0.3
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="competitive-exclusion-zerodisturbance" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <exitCondition>number-of-strategies &lt; 2</exitCondition>
    <metric>(count species1) &gt; 0</metric>
    <metric>(count species2) &gt; 0</metric>
    <metric>ticks</metric>
    <enumeratedValueSet variable="Demographic-Stochasticity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-rate">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="species">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-b">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="data-distance">
      <value value="5200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resolution">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="function">
      <value value="&quot;MSS&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disturbance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction12">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation-probability">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="range">
      <value value="28.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cutoff">
      <value value="0.17"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="capacity">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-population">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-1">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-alpha">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualization">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction21">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-1">
      <value value="0.1353353"/>
      <value value="0.1380692"/>
      <value value="0.1408584"/>
      <value value="0.1437039"/>
      <value value="0.146607"/>
      <value value="0.1495686"/>
      <value value="0.1525901"/>
      <value value="0.1556726"/>
      <value value="0.1588174"/>
      <value value="0.1620258"/>
      <value value="0.1652989"/>
      <value value="0.1686381"/>
      <value value="0.1720449"/>
      <value value="0.1755204"/>
      <value value="0.1790661"/>
      <value value="0.1826835"/>
      <value value="0.186374"/>
      <value value="0.190139"/>
      <value value="0.19398"/>
      <value value="0.1978987"/>
      <value value="0.2018965"/>
      <value value="0.2059751"/>
      <value value="0.2101361"/>
      <value value="0.2143811"/>
      <value value="0.2187119"/>
      <value value="0.2231302"/>
      <value value="0.2276377"/>
      <value value="0.2322363"/>
      <value value="0.2369278"/>
      <value value="0.241714"/>
      <value value="0.246597"/>
      <value value="0.2515786"/>
      <value value="0.2566608"/>
      <value value="0.2618457"/>
      <value value="0.2671353"/>
      <value value="0.2725318"/>
      <value value="0.2780373"/>
      <value value="0.283654"/>
      <value value="0.2893842"/>
      <value value="0.2952302"/>
      <value value="0.3011942"/>
      <value value="0.3072787"/>
      <value value="0.3134862"/>
      <value value="0.319819"/>
      <value value="0.3262798"/>
      <value value="0.3328711"/>
      <value value="0.3395955"/>
      <value value="0.3464558"/>
      <value value="0.3534547"/>
      <value value="0.3605949"/>
      <value value="0.3678794"/>
      <value value="0.3753111"/>
      <value value="0.3828929"/>
      <value value="0.3906278"/>
      <value value="0.398519"/>
      <value value="0.4065697"/>
      <value value="0.4147829"/>
      <value value="0.4231621"/>
      <value value="0.4317105"/>
      <value value="0.4404317"/>
      <value value="0.449329"/>
      <value value="0.458406"/>
      <value value="0.4676664"/>
      <value value="0.4771139"/>
      <value value="0.4867523"/>
      <value value="0.4965853"/>
      <value value="0.506617"/>
      <value value="0.5168513"/>
      <value value="0.5272924"/>
      <value value="0.5379444"/>
      <value value="0.5488116"/>
      <value value="0.5598984"/>
      <value value="0.5712091"/>
      <value value="0.5827483"/>
      <value value="0.5945205"/>
      <value value="0.6065307"/>
      <value value="0.6187834"/>
      <value value="0.6312836"/>
      <value value="0.6440364"/>
      <value value="0.6570468"/>
      <value value="0.67032"/>
      <value value="0.6838614"/>
      <value value="0.6976763"/>
      <value value="0.7117703"/>
      <value value="0.726149"/>
      <value value="0.7408182"/>
      <value value="0.7557837"/>
      <value value="0.7710516"/>
      <value value="0.7866279"/>
      <value value="0.8025188"/>
      <value value="0.8187308"/>
      <value value="0.8352702"/>
      <value value="0.8521438"/>
      <value value="0.8693582"/>
      <value value="0.8869204"/>
      <value value="0.9048374"/>
      <value value="0.9231163"/>
      <value value="0.9417645"/>
      <value value="0.9607894"/>
      <value value="0.9801987"/>
      <value value="1"/>
      <value value="1.0202013"/>
      <value value="1.0408108"/>
      <value value="1.0618365"/>
      <value value="1.0832871"/>
      <value value="1.1051709"/>
      <value value="1.1274969"/>
      <value value="1.1502738"/>
      <value value="1.1735109"/>
      <value value="1.1972174"/>
      <value value="1.2214028"/>
      <value value="1.2460767"/>
      <value value="1.2712492"/>
      <value value="1.2969301"/>
      <value value="1.3231298"/>
      <value value="1.3498588"/>
      <value value="1.3771278"/>
      <value value="1.4049476"/>
      <value value="1.4333294"/>
      <value value="1.4622846"/>
      <value value="1.4918247"/>
      <value value="1.5219616"/>
      <value value="1.5527072"/>
      <value value="1.584074"/>
      <value value="1.6160744"/>
      <value value="1.6487213"/>
      <value value="1.6820276"/>
      <value value="1.7160069"/>
      <value value="1.7506725"/>
      <value value="1.7860384"/>
      <value value="1.8221188"/>
      <value value="1.858928"/>
      <value value="1.8964809"/>
      <value value="1.9347923"/>
      <value value="1.9738777"/>
      <value value="2.0137527"/>
      <value value="2.0544332"/>
      <value value="2.0959355"/>
      <value value="2.1382762"/>
      <value value="2.1814723"/>
      <value value="2.2255409"/>
      <value value="2.2704998"/>
      <value value="2.316367"/>
      <value value="2.3631607"/>
      <value value="2.4108997"/>
      <value value="2.4596031"/>
      <value value="2.5092904"/>
      <value value="2.5599814"/>
      <value value="2.6116965"/>
      <value value="2.6644562"/>
      <value value="2.7182818"/>
      <value value="2.7731948"/>
      <value value="2.829217"/>
      <value value="2.886371"/>
      <value value="2.9446796"/>
      <value value="3.004166"/>
      <value value="3.0648542"/>
      <value value="3.1267684"/>
      <value value="3.1899333"/>
      <value value="3.2543742"/>
      <value value="3.3201169"/>
      <value value="3.3871877"/>
      <value value="3.4556135"/>
      <value value="3.5254215"/>
      <value value="3.5966397"/>
      <value value="3.6692967"/>
      <value value="3.7434214"/>
      <value value="3.8190435"/>
      <value value="3.8961933"/>
      <value value="3.9749016"/>
      <value value="4.0552"/>
      <value value="4.1371204"/>
      <value value="4.2206958"/>
      <value value="4.3059595"/>
      <value value="4.3929457"/>
      <value value="4.4816891"/>
      <value value="4.5722252"/>
      <value value="4.6645903"/>
      <value value="4.7588212"/>
      <value value="4.8549558"/>
      <value value="4.9530324"/>
      <value value="5.0530903"/>
      <value value="5.1551695"/>
      <value value="5.2593108"/>
      <value value="5.365556"/>
      <value value="5.4739474"/>
      <value value="5.5845285"/>
      <value value="5.6973434"/>
      <value value="5.8124374"/>
      <value value="5.9298564"/>
      <value value="6.0496475"/>
      <value value="6.1718584"/>
      <value value="6.2965383"/>
      <value value="6.4237368"/>
      <value value="6.5535049"/>
      <value value="6.6858944"/>
      <value value="6.8209585"/>
      <value value="6.958751"/>
      <value value="7.0993271"/>
      <value value="7.242743"/>
      <value value="7.3890561"/>
      <value value="7.5383249"/>
      <value value="7.6906092"/>
      <value value="7.8459698"/>
      <value value="8.0044689"/>
      <value value="8.1661699"/>
      <value value="8.3311375"/>
      <value value="8.4994376"/>
      <value value="8.6711377"/>
      <value value="8.8463063"/>
      <value value="9.0250135"/>
      <value value="9.2073309"/>
      <value value="9.3933313"/>
      <value value="9.5830892"/>
      <value value="9.7766804"/>
      <value value="9.9741825"/>
      <value value="10.1756743"/>
      <value value="10.3812366"/>
      <value value="10.5909515"/>
      <value value="10.8049029"/>
      <value value="11.0231764"/>
      <value value="11.2458593"/>
      <value value="11.4730407"/>
      <value value="11.7048115"/>
      <value value="11.9412644"/>
      <value value="12.182494"/>
      <value value="12.4285967"/>
      <value value="12.679671"/>
      <value value="12.9358173"/>
      <value value="13.1971382"/>
      <value value="13.463738"/>
      <value value="13.7357236"/>
      <value value="14.0132036"/>
      <value value="14.2962891"/>
      <value value="14.5850933"/>
      <value value="14.8797317"/>
      <value value="15.1803222"/>
      <value value="15.4869851"/>
      <value value="15.7998429"/>
      <value value="16.1190209"/>
      <value value="16.4446468"/>
      <value value="16.7768507"/>
      <value value="17.1157655"/>
      <value value="17.4615269"/>
      <value value="17.8142732"/>
      <value value="18.1741454"/>
      <value value="18.5412875"/>
      <value value="18.9158463"/>
      <value value="19.2979718"/>
      <value value="19.6878166"/>
      <value value="20.0855369"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-2">
      <value value="0.1353353"/>
      <value value="0.1380692"/>
      <value value="0.1408584"/>
      <value value="0.1437039"/>
      <value value="0.146607"/>
      <value value="0.1495686"/>
      <value value="0.1525901"/>
      <value value="0.1556726"/>
      <value value="0.1588174"/>
      <value value="0.1620258"/>
      <value value="0.1652989"/>
      <value value="0.1686381"/>
      <value value="0.1720449"/>
      <value value="0.1755204"/>
      <value value="0.1790661"/>
      <value value="0.1826835"/>
      <value value="0.186374"/>
      <value value="0.190139"/>
      <value value="0.19398"/>
      <value value="0.1978987"/>
      <value value="0.2018965"/>
      <value value="0.2059751"/>
      <value value="0.2101361"/>
      <value value="0.2143811"/>
      <value value="0.2187119"/>
      <value value="0.2231302"/>
      <value value="0.2276377"/>
      <value value="0.2322363"/>
      <value value="0.2369278"/>
      <value value="0.241714"/>
      <value value="0.246597"/>
      <value value="0.2515786"/>
      <value value="0.2566608"/>
      <value value="0.2618457"/>
      <value value="0.2671353"/>
      <value value="0.2725318"/>
      <value value="0.2780373"/>
      <value value="0.283654"/>
      <value value="0.2893842"/>
      <value value="0.2952302"/>
      <value value="0.3011942"/>
      <value value="0.3072787"/>
      <value value="0.3134862"/>
      <value value="0.319819"/>
      <value value="0.3262798"/>
      <value value="0.3328711"/>
      <value value="0.3395955"/>
      <value value="0.3464558"/>
      <value value="0.3534547"/>
      <value value="0.3605949"/>
      <value value="0.3678794"/>
      <value value="0.3753111"/>
      <value value="0.3828929"/>
      <value value="0.3906278"/>
      <value value="0.398519"/>
      <value value="0.4065697"/>
      <value value="0.4147829"/>
      <value value="0.4231621"/>
      <value value="0.4317105"/>
      <value value="0.4404317"/>
      <value value="0.449329"/>
      <value value="0.458406"/>
      <value value="0.4676664"/>
      <value value="0.4771139"/>
      <value value="0.4867523"/>
      <value value="0.4965853"/>
      <value value="0.506617"/>
      <value value="0.5168513"/>
      <value value="0.5272924"/>
      <value value="0.5379444"/>
      <value value="0.5488116"/>
      <value value="0.5598984"/>
      <value value="0.5712091"/>
      <value value="0.5827483"/>
      <value value="0.5945205"/>
      <value value="0.6065307"/>
      <value value="0.6187834"/>
      <value value="0.6312836"/>
      <value value="0.6440364"/>
      <value value="0.6570468"/>
      <value value="0.67032"/>
      <value value="0.6838614"/>
      <value value="0.6976763"/>
      <value value="0.7117703"/>
      <value value="0.726149"/>
      <value value="0.7408182"/>
      <value value="0.7557837"/>
      <value value="0.7710516"/>
      <value value="0.7866279"/>
      <value value="0.8025188"/>
      <value value="0.8187308"/>
      <value value="0.8352702"/>
      <value value="0.8521438"/>
      <value value="0.8693582"/>
      <value value="0.8869204"/>
      <value value="0.9048374"/>
      <value value="0.9231163"/>
      <value value="0.9417645"/>
      <value value="0.9607894"/>
      <value value="0.9801987"/>
      <value value="1"/>
      <value value="1.0202013"/>
      <value value="1.0408108"/>
      <value value="1.0618365"/>
      <value value="1.0832871"/>
      <value value="1.1051709"/>
      <value value="1.1274969"/>
      <value value="1.1502738"/>
      <value value="1.1735109"/>
      <value value="1.1972174"/>
      <value value="1.2214028"/>
      <value value="1.2460767"/>
      <value value="1.2712492"/>
      <value value="1.2969301"/>
      <value value="1.3231298"/>
      <value value="1.3498588"/>
      <value value="1.3771278"/>
      <value value="1.4049476"/>
      <value value="1.4333294"/>
      <value value="1.4622846"/>
      <value value="1.4918247"/>
      <value value="1.5219616"/>
      <value value="1.5527072"/>
      <value value="1.584074"/>
      <value value="1.6160744"/>
      <value value="1.6487213"/>
      <value value="1.6820276"/>
      <value value="1.7160069"/>
      <value value="1.7506725"/>
      <value value="1.7860384"/>
      <value value="1.8221188"/>
      <value value="1.858928"/>
      <value value="1.8964809"/>
      <value value="1.9347923"/>
      <value value="1.9738777"/>
      <value value="2.0137527"/>
      <value value="2.0544332"/>
      <value value="2.0959355"/>
      <value value="2.1382762"/>
      <value value="2.1814723"/>
      <value value="2.2255409"/>
      <value value="2.2704998"/>
      <value value="2.316367"/>
      <value value="2.3631607"/>
      <value value="2.4108997"/>
      <value value="2.4596031"/>
      <value value="2.5092904"/>
      <value value="2.5599814"/>
      <value value="2.6116965"/>
      <value value="2.6644562"/>
      <value value="2.7182818"/>
      <value value="2.7731948"/>
      <value value="2.829217"/>
      <value value="2.886371"/>
      <value value="2.9446796"/>
      <value value="3.004166"/>
      <value value="3.0648542"/>
      <value value="3.1267684"/>
      <value value="3.1899333"/>
      <value value="3.2543742"/>
      <value value="3.3201169"/>
      <value value="3.3871877"/>
      <value value="3.4556135"/>
      <value value="3.5254215"/>
      <value value="3.5966397"/>
      <value value="3.6692967"/>
      <value value="3.7434214"/>
      <value value="3.8190435"/>
      <value value="3.8961933"/>
      <value value="3.9749016"/>
      <value value="4.0552"/>
      <value value="4.1371204"/>
      <value value="4.2206958"/>
      <value value="4.3059595"/>
      <value value="4.3929457"/>
      <value value="4.4816891"/>
      <value value="4.5722252"/>
      <value value="4.6645903"/>
      <value value="4.7588212"/>
      <value value="4.8549558"/>
      <value value="4.9530324"/>
      <value value="5.0530903"/>
      <value value="5.1551695"/>
      <value value="5.2593108"/>
      <value value="5.365556"/>
      <value value="5.4739474"/>
      <value value="5.5845285"/>
      <value value="5.6973434"/>
      <value value="5.8124374"/>
      <value value="5.9298564"/>
      <value value="6.0496475"/>
      <value value="6.1718584"/>
      <value value="6.2965383"/>
      <value value="6.4237368"/>
      <value value="6.5535049"/>
      <value value="6.6858944"/>
      <value value="6.8209585"/>
      <value value="6.958751"/>
      <value value="7.0993271"/>
      <value value="7.242743"/>
      <value value="7.3890561"/>
      <value value="7.5383249"/>
      <value value="7.6906092"/>
      <value value="7.8459698"/>
      <value value="8.0044689"/>
      <value value="8.1661699"/>
      <value value="8.3311375"/>
      <value value="8.4994376"/>
      <value value="8.6711377"/>
      <value value="8.8463063"/>
      <value value="9.0250135"/>
      <value value="9.2073309"/>
      <value value="9.3933313"/>
      <value value="9.5830892"/>
      <value value="9.7766804"/>
      <value value="9.9741825"/>
      <value value="10.1756743"/>
      <value value="10.3812366"/>
      <value value="10.5909515"/>
      <value value="10.8049029"/>
      <value value="11.0231764"/>
      <value value="11.2458593"/>
      <value value="11.4730407"/>
      <value value="11.7048115"/>
      <value value="11.9412644"/>
      <value value="12.182494"/>
      <value value="12.4285967"/>
      <value value="12.679671"/>
      <value value="12.9358173"/>
      <value value="13.1971382"/>
      <value value="13.463738"/>
      <value value="13.7357236"/>
      <value value="14.0132036"/>
      <value value="14.2962891"/>
      <value value="14.5850933"/>
      <value value="14.8797317"/>
      <value value="15.1803222"/>
      <value value="15.4869851"/>
      <value value="15.7998429"/>
      <value value="16.1190209"/>
      <value value="16.4446468"/>
      <value value="16.7768507"/>
      <value value="17.1157655"/>
      <value value="17.4615269"/>
      <value value="17.8142732"/>
      <value value="18.1741454"/>
      <value value="18.5412875"/>
      <value value="18.9158463"/>
      <value value="19.2979718"/>
      <value value="19.6878166"/>
      <value value="20.0855369"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-b">
      <value value="&quot;single&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ess" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>ess-record</go>
    <timeLimit steps="300000"/>
    <metric>average-strategies</metric>
    <metric>average-spread</metric>
    <metric>average-age</metric>
    <metric>average-b</metric>
    <enumeratedValueSet variable="Demographic-Stochasticity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="function">
      <value value="&quot;MSS&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation-probability">
      <value value="0.1"/>
    </enumeratedValueSet>
    <steppedValueSet variable="disturbance" first="0" step="0.0020" last="0.5"/>
    <enumeratedValueSet variable="capacity">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-b">
      <value value="&quot;flat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="range">
      <value value="28.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="data-distance">
      <value value="5200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-b">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cutoff">
      <value value="0.17"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="species">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-2">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resolution">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction21">
      <value value="1"/>
    </enumeratedValueSet>
    <steppedValueSet variable="initial-rate-1" first="1" step="0.02" last="6"/>
    <enumeratedValueSet variable="mutate-alpha">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-rate">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-1">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-population">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction12">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualization">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="invasion-base" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>invasion-fitness-go</go>
    <timeLimit steps="1000"/>
    <metric>count invadors</metric>
    <metric>sum [growth] of invadors</metric>
    <metric>invasioncounter</metric>
    <enumeratedValueSet variable="Demographic-Stochasticity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-1">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-b">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resolution">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="function">
      <value value="&quot;MSS&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-b">
      <value value="&quot;invasion&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-1">
      <value value="0.1353353"/>
      <value value="0.1380692"/>
      <value value="0.1408584"/>
      <value value="0.1437039"/>
      <value value="0.146607"/>
      <value value="0.1495686"/>
      <value value="0.1525901"/>
      <value value="0.1556726"/>
      <value value="0.1588174"/>
      <value value="0.1620258"/>
      <value value="0.1652989"/>
      <value value="0.1686381"/>
      <value value="0.1720449"/>
      <value value="0.1755204"/>
      <value value="0.1790661"/>
      <value value="0.1826835"/>
      <value value="0.186374"/>
      <value value="0.190139"/>
      <value value="0.19398"/>
      <value value="0.1978987"/>
      <value value="0.2018965"/>
      <value value="0.2059751"/>
      <value value="0.2101361"/>
      <value value="0.2143811"/>
      <value value="0.2187119"/>
      <value value="0.2231302"/>
      <value value="0.2276377"/>
      <value value="0.2322363"/>
      <value value="0.2369278"/>
      <value value="0.241714"/>
      <value value="0.246597"/>
      <value value="0.2515786"/>
      <value value="0.2566608"/>
      <value value="0.2618457"/>
      <value value="0.2671353"/>
      <value value="0.2725318"/>
      <value value="0.2780373"/>
      <value value="0.283654"/>
      <value value="0.2893842"/>
      <value value="0.2952302"/>
      <value value="0.3011942"/>
      <value value="0.3072787"/>
      <value value="0.3134862"/>
      <value value="0.319819"/>
      <value value="0.3262798"/>
      <value value="0.3328711"/>
      <value value="0.3395955"/>
      <value value="0.3464558"/>
      <value value="0.3534547"/>
      <value value="0.3605949"/>
      <value value="0.3678794"/>
      <value value="0.3753111"/>
      <value value="0.3828929"/>
      <value value="0.3906278"/>
      <value value="0.398519"/>
      <value value="0.4065697"/>
      <value value="0.4147829"/>
      <value value="0.4231621"/>
      <value value="0.4317105"/>
      <value value="0.4404317"/>
      <value value="0.449329"/>
      <value value="0.458406"/>
      <value value="0.4676664"/>
      <value value="0.4771139"/>
      <value value="0.4867523"/>
      <value value="0.4965853"/>
      <value value="0.506617"/>
      <value value="0.5168513"/>
      <value value="0.5272924"/>
      <value value="0.5379444"/>
      <value value="0.5488116"/>
      <value value="0.5598984"/>
      <value value="0.5712091"/>
      <value value="0.5827483"/>
      <value value="0.5945205"/>
      <value value="0.6065307"/>
      <value value="0.6187834"/>
      <value value="0.6312836"/>
      <value value="0.6440364"/>
      <value value="0.6570468"/>
      <value value="0.67032"/>
      <value value="0.6838614"/>
      <value value="0.6976763"/>
      <value value="0.7117703"/>
      <value value="0.726149"/>
      <value value="0.7408182"/>
      <value value="0.7557837"/>
      <value value="0.7710516"/>
      <value value="0.7866279"/>
      <value value="0.8025188"/>
      <value value="0.8187308"/>
      <value value="0.8352702"/>
      <value value="0.8521438"/>
      <value value="0.8693582"/>
      <value value="0.8869204"/>
      <value value="0.9048374"/>
      <value value="0.9231163"/>
      <value value="0.9417645"/>
      <value value="0.9607894"/>
      <value value="0.9801987"/>
      <value value="1"/>
      <value value="1.0202013"/>
      <value value="1.0408108"/>
      <value value="1.0618365"/>
      <value value="1.0832871"/>
      <value value="1.1051709"/>
      <value value="1.1274969"/>
      <value value="1.1502738"/>
      <value value="1.1735109"/>
      <value value="1.1972174"/>
      <value value="1.2214028"/>
      <value value="1.2460767"/>
      <value value="1.2712492"/>
      <value value="1.2969301"/>
      <value value="1.3231298"/>
      <value value="1.3498588"/>
      <value value="1.3771278"/>
      <value value="1.4049476"/>
      <value value="1.4333294"/>
      <value value="1.4622846"/>
      <value value="1.4918247"/>
      <value value="1.5219616"/>
      <value value="1.5527072"/>
      <value value="1.584074"/>
      <value value="1.6160744"/>
      <value value="1.6487213"/>
      <value value="1.6820276"/>
      <value value="1.7160069"/>
      <value value="1.7506725"/>
      <value value="1.7860384"/>
      <value value="1.8221188"/>
      <value value="1.858928"/>
      <value value="1.8964809"/>
      <value value="1.9347923"/>
      <value value="1.9738777"/>
      <value value="2.0137527"/>
      <value value="2.0544332"/>
      <value value="2.0959355"/>
      <value value="2.1382762"/>
      <value value="2.1814723"/>
      <value value="2.2255409"/>
      <value value="2.2704998"/>
      <value value="2.316367"/>
      <value value="2.3631607"/>
      <value value="2.4108997"/>
      <value value="2.4596031"/>
      <value value="2.5092904"/>
      <value value="2.5599814"/>
      <value value="2.6116965"/>
      <value value="2.6644562"/>
      <value value="2.7182818"/>
      <value value="2.7731948"/>
      <value value="2.829217"/>
      <value value="2.886371"/>
      <value value="2.9446796"/>
      <value value="3.004166"/>
      <value value="3.0648542"/>
      <value value="3.1267684"/>
      <value value="3.1899333"/>
      <value value="3.2543742"/>
      <value value="3.3201169"/>
      <value value="3.3871877"/>
      <value value="3.4556135"/>
      <value value="3.5254215"/>
      <value value="3.5966397"/>
      <value value="3.6692967"/>
      <value value="3.7434214"/>
      <value value="3.8190435"/>
      <value value="3.8961933"/>
      <value value="3.9749016"/>
      <value value="4.0552"/>
      <value value="4.1371204"/>
      <value value="4.2206958"/>
      <value value="4.3059595"/>
      <value value="4.3929457"/>
      <value value="4.4816891"/>
      <value value="4.5722252"/>
      <value value="4.6645903"/>
      <value value="4.7588212"/>
      <value value="4.8549558"/>
      <value value="4.9530324"/>
      <value value="5.0530903"/>
      <value value="5.1551695"/>
      <value value="5.2593108"/>
      <value value="5.365556"/>
      <value value="5.4739474"/>
      <value value="5.5845285"/>
      <value value="5.6973434"/>
      <value value="5.8124374"/>
      <value value="5.9298564"/>
      <value value="6.0496475"/>
      <value value="6.1718584"/>
      <value value="6.2965383"/>
      <value value="6.4237368"/>
      <value value="6.5535049"/>
      <value value="6.6858944"/>
      <value value="6.8209585"/>
      <value value="6.958751"/>
      <value value="7.0993271"/>
      <value value="7.242743"/>
      <value value="7.3890561"/>
      <value value="7.5383249"/>
      <value value="7.6906092"/>
      <value value="7.8459698"/>
      <value value="8.0044689"/>
      <value value="8.1661699"/>
      <value value="8.3311375"/>
      <value value="8.4994376"/>
      <value value="8.6711377"/>
      <value value="8.8463063"/>
      <value value="9.0250135"/>
      <value value="9.2073309"/>
      <value value="9.3933313"/>
      <value value="9.5830892"/>
      <value value="9.7766804"/>
      <value value="9.9741825"/>
      <value value="10.1756743"/>
      <value value="10.3812366"/>
      <value value="10.5909515"/>
      <value value="10.8049029"/>
      <value value="11.0231764"/>
      <value value="11.2458593"/>
      <value value="11.4730407"/>
      <value value="11.7048115"/>
      <value value="11.9412644"/>
      <value value="12.182494"/>
      <value value="12.4285967"/>
      <value value="12.679671"/>
      <value value="12.9358173"/>
      <value value="13.1971382"/>
      <value value="13.463738"/>
      <value value="13.7357236"/>
      <value value="14.0132036"/>
      <value value="14.2962891"/>
      <value value="14.5850933"/>
      <value value="14.8797317"/>
      <value value="15.1803222"/>
      <value value="15.4869851"/>
      <value value="15.7998429"/>
      <value value="16.1190209"/>
      <value value="16.4446468"/>
      <value value="16.7768507"/>
      <value value="17.1157655"/>
      <value value="17.4615269"/>
      <value value="17.8142732"/>
      <value value="18.1741454"/>
      <value value="18.5412875"/>
      <value value="18.9158463"/>
      <value value="19.2979718"/>
      <value value="19.6878166"/>
      <value value="20.0855369"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-2">
      <value value="0.1353353"/>
      <value value="0.1380692"/>
      <value value="0.1408584"/>
      <value value="0.1437039"/>
      <value value="0.146607"/>
      <value value="0.1495686"/>
      <value value="0.1525901"/>
      <value value="0.1556726"/>
      <value value="0.1588174"/>
      <value value="0.1620258"/>
      <value value="0.1652989"/>
      <value value="0.1686381"/>
      <value value="0.1720449"/>
      <value value="0.1755204"/>
      <value value="0.1790661"/>
      <value value="0.1826835"/>
      <value value="0.186374"/>
      <value value="0.190139"/>
      <value value="0.19398"/>
      <value value="0.1978987"/>
      <value value="0.2018965"/>
      <value value="0.2059751"/>
      <value value="0.2101361"/>
      <value value="0.2143811"/>
      <value value="0.2187119"/>
      <value value="0.2231302"/>
      <value value="0.2276377"/>
      <value value="0.2322363"/>
      <value value="0.2369278"/>
      <value value="0.241714"/>
      <value value="0.246597"/>
      <value value="0.2515786"/>
      <value value="0.2566608"/>
      <value value="0.2618457"/>
      <value value="0.2671353"/>
      <value value="0.2725318"/>
      <value value="0.2780373"/>
      <value value="0.283654"/>
      <value value="0.2893842"/>
      <value value="0.2952302"/>
      <value value="0.3011942"/>
      <value value="0.3072787"/>
      <value value="0.3134862"/>
      <value value="0.319819"/>
      <value value="0.3262798"/>
      <value value="0.3328711"/>
      <value value="0.3395955"/>
      <value value="0.3464558"/>
      <value value="0.3534547"/>
      <value value="0.3605949"/>
      <value value="0.3678794"/>
      <value value="0.3753111"/>
      <value value="0.3828929"/>
      <value value="0.3906278"/>
      <value value="0.398519"/>
      <value value="0.4065697"/>
      <value value="0.4147829"/>
      <value value="0.4231621"/>
      <value value="0.4317105"/>
      <value value="0.4404317"/>
      <value value="0.449329"/>
      <value value="0.458406"/>
      <value value="0.4676664"/>
      <value value="0.4771139"/>
      <value value="0.4867523"/>
      <value value="0.4965853"/>
      <value value="0.506617"/>
      <value value="0.5168513"/>
      <value value="0.5272924"/>
      <value value="0.5379444"/>
      <value value="0.5488116"/>
      <value value="0.5598984"/>
      <value value="0.5712091"/>
      <value value="0.5827483"/>
      <value value="0.5945205"/>
      <value value="0.6065307"/>
      <value value="0.6187834"/>
      <value value="0.6312836"/>
      <value value="0.6440364"/>
      <value value="0.6570468"/>
      <value value="0.67032"/>
      <value value="0.6838614"/>
      <value value="0.6976763"/>
      <value value="0.7117703"/>
      <value value="0.726149"/>
      <value value="0.7408182"/>
      <value value="0.7557837"/>
      <value value="0.7710516"/>
      <value value="0.7866279"/>
      <value value="0.8025188"/>
      <value value="0.8187308"/>
      <value value="0.8352702"/>
      <value value="0.8521438"/>
      <value value="0.8693582"/>
      <value value="0.8869204"/>
      <value value="0.9048374"/>
      <value value="0.9231163"/>
      <value value="0.9417645"/>
      <value value="0.9607894"/>
      <value value="0.9801987"/>
      <value value="1"/>
      <value value="1.0202013"/>
      <value value="1.0408108"/>
      <value value="1.0618365"/>
      <value value="1.0832871"/>
      <value value="1.1051709"/>
      <value value="1.1274969"/>
      <value value="1.1502738"/>
      <value value="1.1735109"/>
      <value value="1.1972174"/>
      <value value="1.2214028"/>
      <value value="1.2460767"/>
      <value value="1.2712492"/>
      <value value="1.2969301"/>
      <value value="1.3231298"/>
      <value value="1.3498588"/>
      <value value="1.3771278"/>
      <value value="1.4049476"/>
      <value value="1.4333294"/>
      <value value="1.4622846"/>
      <value value="1.4918247"/>
      <value value="1.5219616"/>
      <value value="1.5527072"/>
      <value value="1.584074"/>
      <value value="1.6160744"/>
      <value value="1.6487213"/>
      <value value="1.6820276"/>
      <value value="1.7160069"/>
      <value value="1.7506725"/>
      <value value="1.7860384"/>
      <value value="1.8221188"/>
      <value value="1.858928"/>
      <value value="1.8964809"/>
      <value value="1.9347923"/>
      <value value="1.9738777"/>
      <value value="2.0137527"/>
      <value value="2.0544332"/>
      <value value="2.0959355"/>
      <value value="2.1382762"/>
      <value value="2.1814723"/>
      <value value="2.2255409"/>
      <value value="2.2704998"/>
      <value value="2.316367"/>
      <value value="2.3631607"/>
      <value value="2.4108997"/>
      <value value="2.4596031"/>
      <value value="2.5092904"/>
      <value value="2.5599814"/>
      <value value="2.6116965"/>
      <value value="2.6644562"/>
      <value value="2.7182818"/>
      <value value="2.7731948"/>
      <value value="2.829217"/>
      <value value="2.886371"/>
      <value value="2.9446796"/>
      <value value="3.004166"/>
      <value value="3.0648542"/>
      <value value="3.1267684"/>
      <value value="3.1899333"/>
      <value value="3.2543742"/>
      <value value="3.3201169"/>
      <value value="3.3871877"/>
      <value value="3.4556135"/>
      <value value="3.5254215"/>
      <value value="3.5966397"/>
      <value value="3.6692967"/>
      <value value="3.7434214"/>
      <value value="3.8190435"/>
      <value value="3.8961933"/>
      <value value="3.9749016"/>
      <value value="4.0552"/>
      <value value="4.1371204"/>
      <value value="4.2206958"/>
      <value value="4.3059595"/>
      <value value="4.3929457"/>
      <value value="4.4816891"/>
      <value value="4.5722252"/>
      <value value="4.6645903"/>
      <value value="4.7588212"/>
      <value value="4.8549558"/>
      <value value="4.9530324"/>
      <value value="5.0530903"/>
      <value value="5.1551695"/>
      <value value="5.2593108"/>
      <value value="5.365556"/>
      <value value="5.4739474"/>
      <value value="5.5845285"/>
      <value value="5.6973434"/>
      <value value="5.8124374"/>
      <value value="5.9298564"/>
      <value value="6.0496475"/>
      <value value="6.1718584"/>
      <value value="6.2965383"/>
      <value value="6.4237368"/>
      <value value="6.5535049"/>
      <value value="6.6858944"/>
      <value value="6.8209585"/>
      <value value="6.958751"/>
      <value value="7.0993271"/>
      <value value="7.242743"/>
      <value value="7.3890561"/>
      <value value="7.5383249"/>
      <value value="7.6906092"/>
      <value value="7.8459698"/>
      <value value="8.0044689"/>
      <value value="8.1661699"/>
      <value value="8.3311375"/>
      <value value="8.4994376"/>
      <value value="8.6711377"/>
      <value value="8.8463063"/>
      <value value="9.0250135"/>
      <value value="9.2073309"/>
      <value value="9.3933313"/>
      <value value="9.5830892"/>
      <value value="9.7766804"/>
      <value value="9.9741825"/>
      <value value="10.1756743"/>
      <value value="10.3812366"/>
      <value value="10.5909515"/>
      <value value="10.8049029"/>
      <value value="11.0231764"/>
      <value value="11.2458593"/>
      <value value="11.4730407"/>
      <value value="11.7048115"/>
      <value value="11.9412644"/>
      <value value="12.182494"/>
      <value value="12.4285967"/>
      <value value="12.679671"/>
      <value value="12.9358173"/>
      <value value="13.1971382"/>
      <value value="13.463738"/>
      <value value="13.7357236"/>
      <value value="14.0132036"/>
      <value value="14.2962891"/>
      <value value="14.5850933"/>
      <value value="14.8797317"/>
      <value value="15.1803222"/>
      <value value="15.4869851"/>
      <value value="15.7998429"/>
      <value value="16.1190209"/>
      <value value="16.4446468"/>
      <value value="16.7768507"/>
      <value value="17.1157655"/>
      <value value="17.4615269"/>
      <value value="17.8142732"/>
      <value value="18.1741454"/>
      <value value="18.5412875"/>
      <value value="18.9158463"/>
      <value value="19.2979718"/>
      <value value="19.6878166"/>
      <value value="20.0855369"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-alpha">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="capacity">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="species">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-population">
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="data-distance">
      <value value="5200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation-probability">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction21">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cutoff">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualization">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-rate">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="range">
      <value value="28.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disturbance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction12">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="invasion-dist" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>invasion-fitness-go</go>
    <timeLimit steps="1000"/>
    <metric>count invadors</metric>
    <metric>sum [growth] of invadors</metric>
    <metric>invasioncounter</metric>
    <enumeratedValueSet variable="Demographic-Stochasticity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-1">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-b">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resolution">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="function">
      <value value="&quot;MSS&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-b">
      <value value="&quot;invasion&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-1">
      <value value="0.1353353"/>
      <value value="0.1380692"/>
      <value value="0.1408584"/>
      <value value="0.1437039"/>
      <value value="0.146607"/>
      <value value="0.1495686"/>
      <value value="0.1525901"/>
      <value value="0.1556726"/>
      <value value="0.1588174"/>
      <value value="0.1620258"/>
      <value value="0.1652989"/>
      <value value="0.1686381"/>
      <value value="0.1720449"/>
      <value value="0.1755204"/>
      <value value="0.1790661"/>
      <value value="0.1826835"/>
      <value value="0.186374"/>
      <value value="0.190139"/>
      <value value="0.19398"/>
      <value value="0.1978987"/>
      <value value="0.2018965"/>
      <value value="0.2059751"/>
      <value value="0.2101361"/>
      <value value="0.2143811"/>
      <value value="0.2187119"/>
      <value value="0.2231302"/>
      <value value="0.2276377"/>
      <value value="0.2322363"/>
      <value value="0.2369278"/>
      <value value="0.241714"/>
      <value value="0.246597"/>
      <value value="0.2515786"/>
      <value value="0.2566608"/>
      <value value="0.2618457"/>
      <value value="0.2671353"/>
      <value value="0.2725318"/>
      <value value="0.2780373"/>
      <value value="0.283654"/>
      <value value="0.2893842"/>
      <value value="0.2952302"/>
      <value value="0.3011942"/>
      <value value="0.3072787"/>
      <value value="0.3134862"/>
      <value value="0.319819"/>
      <value value="0.3262798"/>
      <value value="0.3328711"/>
      <value value="0.3395955"/>
      <value value="0.3464558"/>
      <value value="0.3534547"/>
      <value value="0.3605949"/>
      <value value="0.3678794"/>
      <value value="0.3753111"/>
      <value value="0.3828929"/>
      <value value="0.3906278"/>
      <value value="0.398519"/>
      <value value="0.4065697"/>
      <value value="0.4147829"/>
      <value value="0.4231621"/>
      <value value="0.4317105"/>
      <value value="0.4404317"/>
      <value value="0.449329"/>
      <value value="0.458406"/>
      <value value="0.4676664"/>
      <value value="0.4771139"/>
      <value value="0.4867523"/>
      <value value="0.4965853"/>
      <value value="0.506617"/>
      <value value="0.5168513"/>
      <value value="0.5272924"/>
      <value value="0.5379444"/>
      <value value="0.5488116"/>
      <value value="0.5598984"/>
      <value value="0.5712091"/>
      <value value="0.5827483"/>
      <value value="0.5945205"/>
      <value value="0.6065307"/>
      <value value="0.6187834"/>
      <value value="0.6312836"/>
      <value value="0.6440364"/>
      <value value="0.6570468"/>
      <value value="0.67032"/>
      <value value="0.6838614"/>
      <value value="0.6976763"/>
      <value value="0.7117703"/>
      <value value="0.726149"/>
      <value value="0.7408182"/>
      <value value="0.7557837"/>
      <value value="0.7710516"/>
      <value value="0.7866279"/>
      <value value="0.8025188"/>
      <value value="0.8187308"/>
      <value value="0.8352702"/>
      <value value="0.8521438"/>
      <value value="0.8693582"/>
      <value value="0.8869204"/>
      <value value="0.9048374"/>
      <value value="0.9231163"/>
      <value value="0.9417645"/>
      <value value="0.9607894"/>
      <value value="0.9801987"/>
      <value value="1"/>
      <value value="1.0202013"/>
      <value value="1.0408108"/>
      <value value="1.0618365"/>
      <value value="1.0832871"/>
      <value value="1.1051709"/>
      <value value="1.1274969"/>
      <value value="1.1502738"/>
      <value value="1.1735109"/>
      <value value="1.1972174"/>
      <value value="1.2214028"/>
      <value value="1.2460767"/>
      <value value="1.2712492"/>
      <value value="1.2969301"/>
      <value value="1.3231298"/>
      <value value="1.3498588"/>
      <value value="1.3771278"/>
      <value value="1.4049476"/>
      <value value="1.4333294"/>
      <value value="1.4622846"/>
      <value value="1.4918247"/>
      <value value="1.5219616"/>
      <value value="1.5527072"/>
      <value value="1.584074"/>
      <value value="1.6160744"/>
      <value value="1.6487213"/>
      <value value="1.6820276"/>
      <value value="1.7160069"/>
      <value value="1.7506725"/>
      <value value="1.7860384"/>
      <value value="1.8221188"/>
      <value value="1.858928"/>
      <value value="1.8964809"/>
      <value value="1.9347923"/>
      <value value="1.9738777"/>
      <value value="2.0137527"/>
      <value value="2.0544332"/>
      <value value="2.0959355"/>
      <value value="2.1382762"/>
      <value value="2.1814723"/>
      <value value="2.2255409"/>
      <value value="2.2704998"/>
      <value value="2.316367"/>
      <value value="2.3631607"/>
      <value value="2.4108997"/>
      <value value="2.4596031"/>
      <value value="2.5092904"/>
      <value value="2.5599814"/>
      <value value="2.6116965"/>
      <value value="2.6644562"/>
      <value value="2.7182818"/>
      <value value="2.7731948"/>
      <value value="2.829217"/>
      <value value="2.886371"/>
      <value value="2.9446796"/>
      <value value="3.004166"/>
      <value value="3.0648542"/>
      <value value="3.1267684"/>
      <value value="3.1899333"/>
      <value value="3.2543742"/>
      <value value="3.3201169"/>
      <value value="3.3871877"/>
      <value value="3.4556135"/>
      <value value="3.5254215"/>
      <value value="3.5966397"/>
      <value value="3.6692967"/>
      <value value="3.7434214"/>
      <value value="3.8190435"/>
      <value value="3.8961933"/>
      <value value="3.9749016"/>
      <value value="4.0552"/>
      <value value="4.1371204"/>
      <value value="4.2206958"/>
      <value value="4.3059595"/>
      <value value="4.3929457"/>
      <value value="4.4816891"/>
      <value value="4.5722252"/>
      <value value="4.6645903"/>
      <value value="4.7588212"/>
      <value value="4.8549558"/>
      <value value="4.9530324"/>
      <value value="5.0530903"/>
      <value value="5.1551695"/>
      <value value="5.2593108"/>
      <value value="5.365556"/>
      <value value="5.4739474"/>
      <value value="5.5845285"/>
      <value value="5.6973434"/>
      <value value="5.8124374"/>
      <value value="5.9298564"/>
      <value value="6.0496475"/>
      <value value="6.1718584"/>
      <value value="6.2965383"/>
      <value value="6.4237368"/>
      <value value="6.5535049"/>
      <value value="6.6858944"/>
      <value value="6.8209585"/>
      <value value="6.958751"/>
      <value value="7.0993271"/>
      <value value="7.242743"/>
      <value value="7.3890561"/>
      <value value="7.5383249"/>
      <value value="7.6906092"/>
      <value value="7.8459698"/>
      <value value="8.0044689"/>
      <value value="8.1661699"/>
      <value value="8.3311375"/>
      <value value="8.4994376"/>
      <value value="8.6711377"/>
      <value value="8.8463063"/>
      <value value="9.0250135"/>
      <value value="9.2073309"/>
      <value value="9.3933313"/>
      <value value="9.5830892"/>
      <value value="9.7766804"/>
      <value value="9.9741825"/>
      <value value="10.1756743"/>
      <value value="10.3812366"/>
      <value value="10.5909515"/>
      <value value="10.8049029"/>
      <value value="11.0231764"/>
      <value value="11.2458593"/>
      <value value="11.4730407"/>
      <value value="11.7048115"/>
      <value value="11.9412644"/>
      <value value="12.182494"/>
      <value value="12.4285967"/>
      <value value="12.679671"/>
      <value value="12.9358173"/>
      <value value="13.1971382"/>
      <value value="13.463738"/>
      <value value="13.7357236"/>
      <value value="14.0132036"/>
      <value value="14.2962891"/>
      <value value="14.5850933"/>
      <value value="14.8797317"/>
      <value value="15.1803222"/>
      <value value="15.4869851"/>
      <value value="15.7998429"/>
      <value value="16.1190209"/>
      <value value="16.4446468"/>
      <value value="16.7768507"/>
      <value value="17.1157655"/>
      <value value="17.4615269"/>
      <value value="17.8142732"/>
      <value value="18.1741454"/>
      <value value="18.5412875"/>
      <value value="18.9158463"/>
      <value value="19.2979718"/>
      <value value="19.6878166"/>
      <value value="20.0855369"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-2">
      <value value="0.1353353"/>
      <value value="0.1380692"/>
      <value value="0.1408584"/>
      <value value="0.1437039"/>
      <value value="0.146607"/>
      <value value="0.1495686"/>
      <value value="0.1525901"/>
      <value value="0.1556726"/>
      <value value="0.1588174"/>
      <value value="0.1620258"/>
      <value value="0.1652989"/>
      <value value="0.1686381"/>
      <value value="0.1720449"/>
      <value value="0.1755204"/>
      <value value="0.1790661"/>
      <value value="0.1826835"/>
      <value value="0.186374"/>
      <value value="0.190139"/>
      <value value="0.19398"/>
      <value value="0.1978987"/>
      <value value="0.2018965"/>
      <value value="0.2059751"/>
      <value value="0.2101361"/>
      <value value="0.2143811"/>
      <value value="0.2187119"/>
      <value value="0.2231302"/>
      <value value="0.2276377"/>
      <value value="0.2322363"/>
      <value value="0.2369278"/>
      <value value="0.241714"/>
      <value value="0.246597"/>
      <value value="0.2515786"/>
      <value value="0.2566608"/>
      <value value="0.2618457"/>
      <value value="0.2671353"/>
      <value value="0.2725318"/>
      <value value="0.2780373"/>
      <value value="0.283654"/>
      <value value="0.2893842"/>
      <value value="0.2952302"/>
      <value value="0.3011942"/>
      <value value="0.3072787"/>
      <value value="0.3134862"/>
      <value value="0.319819"/>
      <value value="0.3262798"/>
      <value value="0.3328711"/>
      <value value="0.3395955"/>
      <value value="0.3464558"/>
      <value value="0.3534547"/>
      <value value="0.3605949"/>
      <value value="0.3678794"/>
      <value value="0.3753111"/>
      <value value="0.3828929"/>
      <value value="0.3906278"/>
      <value value="0.398519"/>
      <value value="0.4065697"/>
      <value value="0.4147829"/>
      <value value="0.4231621"/>
      <value value="0.4317105"/>
      <value value="0.4404317"/>
      <value value="0.449329"/>
      <value value="0.458406"/>
      <value value="0.4676664"/>
      <value value="0.4771139"/>
      <value value="0.4867523"/>
      <value value="0.4965853"/>
      <value value="0.506617"/>
      <value value="0.5168513"/>
      <value value="0.5272924"/>
      <value value="0.5379444"/>
      <value value="0.5488116"/>
      <value value="0.5598984"/>
      <value value="0.5712091"/>
      <value value="0.5827483"/>
      <value value="0.5945205"/>
      <value value="0.6065307"/>
      <value value="0.6187834"/>
      <value value="0.6312836"/>
      <value value="0.6440364"/>
      <value value="0.6570468"/>
      <value value="0.67032"/>
      <value value="0.6838614"/>
      <value value="0.6976763"/>
      <value value="0.7117703"/>
      <value value="0.726149"/>
      <value value="0.7408182"/>
      <value value="0.7557837"/>
      <value value="0.7710516"/>
      <value value="0.7866279"/>
      <value value="0.8025188"/>
      <value value="0.8187308"/>
      <value value="0.8352702"/>
      <value value="0.8521438"/>
      <value value="0.8693582"/>
      <value value="0.8869204"/>
      <value value="0.9048374"/>
      <value value="0.9231163"/>
      <value value="0.9417645"/>
      <value value="0.9607894"/>
      <value value="0.9801987"/>
      <value value="1"/>
      <value value="1.0202013"/>
      <value value="1.0408108"/>
      <value value="1.0618365"/>
      <value value="1.0832871"/>
      <value value="1.1051709"/>
      <value value="1.1274969"/>
      <value value="1.1502738"/>
      <value value="1.1735109"/>
      <value value="1.1972174"/>
      <value value="1.2214028"/>
      <value value="1.2460767"/>
      <value value="1.2712492"/>
      <value value="1.2969301"/>
      <value value="1.3231298"/>
      <value value="1.3498588"/>
      <value value="1.3771278"/>
      <value value="1.4049476"/>
      <value value="1.4333294"/>
      <value value="1.4622846"/>
      <value value="1.4918247"/>
      <value value="1.5219616"/>
      <value value="1.5527072"/>
      <value value="1.584074"/>
      <value value="1.6160744"/>
      <value value="1.6487213"/>
      <value value="1.6820276"/>
      <value value="1.7160069"/>
      <value value="1.7506725"/>
      <value value="1.7860384"/>
      <value value="1.8221188"/>
      <value value="1.858928"/>
      <value value="1.8964809"/>
      <value value="1.9347923"/>
      <value value="1.9738777"/>
      <value value="2.0137527"/>
      <value value="2.0544332"/>
      <value value="2.0959355"/>
      <value value="2.1382762"/>
      <value value="2.1814723"/>
      <value value="2.2255409"/>
      <value value="2.2704998"/>
      <value value="2.316367"/>
      <value value="2.3631607"/>
      <value value="2.4108997"/>
      <value value="2.4596031"/>
      <value value="2.5092904"/>
      <value value="2.5599814"/>
      <value value="2.6116965"/>
      <value value="2.6644562"/>
      <value value="2.7182818"/>
      <value value="2.7731948"/>
      <value value="2.829217"/>
      <value value="2.886371"/>
      <value value="2.9446796"/>
      <value value="3.004166"/>
      <value value="3.0648542"/>
      <value value="3.1267684"/>
      <value value="3.1899333"/>
      <value value="3.2543742"/>
      <value value="3.3201169"/>
      <value value="3.3871877"/>
      <value value="3.4556135"/>
      <value value="3.5254215"/>
      <value value="3.5966397"/>
      <value value="3.6692967"/>
      <value value="3.7434214"/>
      <value value="3.8190435"/>
      <value value="3.8961933"/>
      <value value="3.9749016"/>
      <value value="4.0552"/>
      <value value="4.1371204"/>
      <value value="4.2206958"/>
      <value value="4.3059595"/>
      <value value="4.3929457"/>
      <value value="4.4816891"/>
      <value value="4.5722252"/>
      <value value="4.6645903"/>
      <value value="4.7588212"/>
      <value value="4.8549558"/>
      <value value="4.9530324"/>
      <value value="5.0530903"/>
      <value value="5.1551695"/>
      <value value="5.2593108"/>
      <value value="5.365556"/>
      <value value="5.4739474"/>
      <value value="5.5845285"/>
      <value value="5.6973434"/>
      <value value="5.8124374"/>
      <value value="5.9298564"/>
      <value value="6.0496475"/>
      <value value="6.1718584"/>
      <value value="6.2965383"/>
      <value value="6.4237368"/>
      <value value="6.5535049"/>
      <value value="6.6858944"/>
      <value value="6.8209585"/>
      <value value="6.958751"/>
      <value value="7.0993271"/>
      <value value="7.242743"/>
      <value value="7.3890561"/>
      <value value="7.5383249"/>
      <value value="7.6906092"/>
      <value value="7.8459698"/>
      <value value="8.0044689"/>
      <value value="8.1661699"/>
      <value value="8.3311375"/>
      <value value="8.4994376"/>
      <value value="8.6711377"/>
      <value value="8.8463063"/>
      <value value="9.0250135"/>
      <value value="9.2073309"/>
      <value value="9.3933313"/>
      <value value="9.5830892"/>
      <value value="9.7766804"/>
      <value value="9.9741825"/>
      <value value="10.1756743"/>
      <value value="10.3812366"/>
      <value value="10.5909515"/>
      <value value="10.8049029"/>
      <value value="11.0231764"/>
      <value value="11.2458593"/>
      <value value="11.4730407"/>
      <value value="11.7048115"/>
      <value value="11.9412644"/>
      <value value="12.182494"/>
      <value value="12.4285967"/>
      <value value="12.679671"/>
      <value value="12.9358173"/>
      <value value="13.1971382"/>
      <value value="13.463738"/>
      <value value="13.7357236"/>
      <value value="14.0132036"/>
      <value value="14.2962891"/>
      <value value="14.5850933"/>
      <value value="14.8797317"/>
      <value value="15.1803222"/>
      <value value="15.4869851"/>
      <value value="15.7998429"/>
      <value value="16.1190209"/>
      <value value="16.4446468"/>
      <value value="16.7768507"/>
      <value value="17.1157655"/>
      <value value="17.4615269"/>
      <value value="17.8142732"/>
      <value value="18.1741454"/>
      <value value="18.5412875"/>
      <value value="18.9158463"/>
      <value value="19.2979718"/>
      <value value="19.6878166"/>
      <value value="20.0855369"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-alpha">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="capacity">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="species">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-population">
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="data-distance">
      <value value="5200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation-probability">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction21">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cutoff">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualization">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-rate">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="range">
      <value value="28.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disturbance">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction12">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="invasion-stochastic-all" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>invasion-fitness-go</go>
    <timeLimit steps="12200"/>
    <metric>count invadors</metric>
    <metric>sum [growth] of invadors</metric>
    <metric>invasioncounter</metric>
    <enumeratedValueSet variable="Demographic-Stochasticity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-1">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-b">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resolution">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="function">
      <value value="&quot;MSS&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-b">
      <value value="&quot;invasion&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-1">
      <value value="0.1353353"/>
      <value value="0.1380692"/>
      <value value="0.1408584"/>
      <value value="0.1437039"/>
      <value value="0.146607"/>
      <value value="0.1495686"/>
      <value value="0.1525901"/>
      <value value="0.1556726"/>
      <value value="0.1588174"/>
      <value value="0.1620258"/>
      <value value="0.1652989"/>
      <value value="0.1686381"/>
      <value value="0.1720449"/>
      <value value="0.1755204"/>
      <value value="0.1790661"/>
      <value value="0.1826835"/>
      <value value="0.186374"/>
      <value value="0.190139"/>
      <value value="0.19398"/>
      <value value="0.1978987"/>
      <value value="0.2018965"/>
      <value value="0.2059751"/>
      <value value="0.2101361"/>
      <value value="0.2143811"/>
      <value value="0.2187119"/>
      <value value="0.2231302"/>
      <value value="0.2276377"/>
      <value value="0.2322363"/>
      <value value="0.2369278"/>
      <value value="0.241714"/>
      <value value="0.246597"/>
      <value value="0.2515786"/>
      <value value="0.2566608"/>
      <value value="0.2618457"/>
      <value value="0.2671353"/>
      <value value="0.2725318"/>
      <value value="0.2780373"/>
      <value value="0.283654"/>
      <value value="0.2893842"/>
      <value value="0.2952302"/>
      <value value="0.3011942"/>
      <value value="0.3072787"/>
      <value value="0.3134862"/>
      <value value="0.319819"/>
      <value value="0.3262798"/>
      <value value="0.3328711"/>
      <value value="0.3395955"/>
      <value value="0.3464558"/>
      <value value="0.3534547"/>
      <value value="0.3605949"/>
      <value value="0.3678794"/>
      <value value="0.3753111"/>
      <value value="0.3828929"/>
      <value value="0.3906278"/>
      <value value="0.398519"/>
      <value value="0.4065697"/>
      <value value="0.4147829"/>
      <value value="0.4231621"/>
      <value value="0.4317105"/>
      <value value="0.4404317"/>
      <value value="0.449329"/>
      <value value="0.458406"/>
      <value value="0.4676664"/>
      <value value="0.4771139"/>
      <value value="0.4867523"/>
      <value value="0.4965853"/>
      <value value="0.506617"/>
      <value value="0.5168513"/>
      <value value="0.5272924"/>
      <value value="0.5379444"/>
      <value value="0.5488116"/>
      <value value="0.5598984"/>
      <value value="0.5712091"/>
      <value value="0.5827483"/>
      <value value="0.5945205"/>
      <value value="0.6065307"/>
      <value value="0.6187834"/>
      <value value="0.6312836"/>
      <value value="0.6440364"/>
      <value value="0.6570468"/>
      <value value="0.67032"/>
      <value value="0.6838614"/>
      <value value="0.6976763"/>
      <value value="0.7117703"/>
      <value value="0.726149"/>
      <value value="0.7408182"/>
      <value value="0.7557837"/>
      <value value="0.7710516"/>
      <value value="0.7866279"/>
      <value value="0.8025188"/>
      <value value="0.8187308"/>
      <value value="0.8352702"/>
      <value value="0.8521438"/>
      <value value="0.8693582"/>
      <value value="0.8869204"/>
      <value value="0.9048374"/>
      <value value="0.9231163"/>
      <value value="0.9417645"/>
      <value value="0.9607894"/>
      <value value="0.9801987"/>
      <value value="1"/>
      <value value="1.0202013"/>
      <value value="1.0408108"/>
      <value value="1.0618365"/>
      <value value="1.0832871"/>
      <value value="1.1051709"/>
      <value value="1.1274969"/>
      <value value="1.1502738"/>
      <value value="1.1735109"/>
      <value value="1.1972174"/>
      <value value="1.2214028"/>
      <value value="1.2460767"/>
      <value value="1.2712492"/>
      <value value="1.2969301"/>
      <value value="1.3231298"/>
      <value value="1.3498588"/>
      <value value="1.3771278"/>
      <value value="1.4049476"/>
      <value value="1.4333294"/>
      <value value="1.4622846"/>
      <value value="1.4918247"/>
      <value value="1.5219616"/>
      <value value="1.5527072"/>
      <value value="1.584074"/>
      <value value="1.6160744"/>
      <value value="1.6487213"/>
      <value value="1.6820276"/>
      <value value="1.7160069"/>
      <value value="1.7506725"/>
      <value value="1.7860384"/>
      <value value="1.8221188"/>
      <value value="1.858928"/>
      <value value="1.8964809"/>
      <value value="1.9347923"/>
      <value value="1.9738777"/>
      <value value="2.0137527"/>
      <value value="2.0544332"/>
      <value value="2.0959355"/>
      <value value="2.1382762"/>
      <value value="2.1814723"/>
      <value value="2.2255409"/>
      <value value="2.2704998"/>
      <value value="2.316367"/>
      <value value="2.3631607"/>
      <value value="2.4108997"/>
      <value value="2.4596031"/>
      <value value="2.5092904"/>
      <value value="2.5599814"/>
      <value value="2.6116965"/>
      <value value="2.6644562"/>
      <value value="2.7182818"/>
      <value value="2.7731948"/>
      <value value="2.829217"/>
      <value value="2.886371"/>
      <value value="2.9446796"/>
      <value value="3.004166"/>
      <value value="3.0648542"/>
      <value value="3.1267684"/>
      <value value="3.1899333"/>
      <value value="3.2543742"/>
      <value value="3.3201169"/>
      <value value="3.3871877"/>
      <value value="3.4556135"/>
      <value value="3.5254215"/>
      <value value="3.5966397"/>
      <value value="3.6692967"/>
      <value value="3.7434214"/>
      <value value="3.8190435"/>
      <value value="3.8961933"/>
      <value value="3.9749016"/>
      <value value="4.0552"/>
      <value value="4.1371204"/>
      <value value="4.2206958"/>
      <value value="4.3059595"/>
      <value value="4.3929457"/>
      <value value="4.4816891"/>
      <value value="4.5722252"/>
      <value value="4.6645903"/>
      <value value="4.7588212"/>
      <value value="4.8549558"/>
      <value value="4.9530324"/>
      <value value="5.0530903"/>
      <value value="5.1551695"/>
      <value value="5.2593108"/>
      <value value="5.365556"/>
      <value value="5.4739474"/>
      <value value="5.5845285"/>
      <value value="5.6973434"/>
      <value value="5.8124374"/>
      <value value="5.9298564"/>
      <value value="6.0496475"/>
      <value value="6.1718584"/>
      <value value="6.2965383"/>
      <value value="6.4237368"/>
      <value value="6.5535049"/>
      <value value="6.6858944"/>
      <value value="6.8209585"/>
      <value value="6.958751"/>
      <value value="7.0993271"/>
      <value value="7.242743"/>
      <value value="7.3890561"/>
      <value value="7.5383249"/>
      <value value="7.6906092"/>
      <value value="7.8459698"/>
      <value value="8.0044689"/>
      <value value="8.1661699"/>
      <value value="8.3311375"/>
      <value value="8.4994376"/>
      <value value="8.6711377"/>
      <value value="8.8463063"/>
      <value value="9.0250135"/>
      <value value="9.2073309"/>
      <value value="9.3933313"/>
      <value value="9.5830892"/>
      <value value="9.7766804"/>
      <value value="9.9741825"/>
      <value value="10.1756743"/>
      <value value="10.3812366"/>
      <value value="10.5909515"/>
      <value value="10.8049029"/>
      <value value="11.0231764"/>
      <value value="11.2458593"/>
      <value value="11.4730407"/>
      <value value="11.7048115"/>
      <value value="11.9412644"/>
      <value value="12.182494"/>
      <value value="12.4285967"/>
      <value value="12.679671"/>
      <value value="12.9358173"/>
      <value value="13.1971382"/>
      <value value="13.463738"/>
      <value value="13.7357236"/>
      <value value="14.0132036"/>
      <value value="14.2962891"/>
      <value value="14.5850933"/>
      <value value="14.8797317"/>
      <value value="15.1803222"/>
      <value value="15.4869851"/>
      <value value="15.7998429"/>
      <value value="16.1190209"/>
      <value value="16.4446468"/>
      <value value="16.7768507"/>
      <value value="17.1157655"/>
      <value value="17.4615269"/>
      <value value="17.8142732"/>
      <value value="18.1741454"/>
      <value value="18.5412875"/>
      <value value="18.9158463"/>
      <value value="19.2979718"/>
      <value value="19.6878166"/>
      <value value="20.0855369"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-2">
      <value value="0.1353353"/>
      <value value="0.1380692"/>
      <value value="0.1408584"/>
      <value value="0.1437039"/>
      <value value="0.146607"/>
      <value value="0.1495686"/>
      <value value="0.1525901"/>
      <value value="0.1556726"/>
      <value value="0.1588174"/>
      <value value="0.1620258"/>
      <value value="0.1652989"/>
      <value value="0.1686381"/>
      <value value="0.1720449"/>
      <value value="0.1755204"/>
      <value value="0.1790661"/>
      <value value="0.1826835"/>
      <value value="0.186374"/>
      <value value="0.190139"/>
      <value value="0.19398"/>
      <value value="0.1978987"/>
      <value value="0.2018965"/>
      <value value="0.2059751"/>
      <value value="0.2101361"/>
      <value value="0.2143811"/>
      <value value="0.2187119"/>
      <value value="0.2231302"/>
      <value value="0.2276377"/>
      <value value="0.2322363"/>
      <value value="0.2369278"/>
      <value value="0.241714"/>
      <value value="0.246597"/>
      <value value="0.2515786"/>
      <value value="0.2566608"/>
      <value value="0.2618457"/>
      <value value="0.2671353"/>
      <value value="0.2725318"/>
      <value value="0.2780373"/>
      <value value="0.283654"/>
      <value value="0.2893842"/>
      <value value="0.2952302"/>
      <value value="0.3011942"/>
      <value value="0.3072787"/>
      <value value="0.3134862"/>
      <value value="0.319819"/>
      <value value="0.3262798"/>
      <value value="0.3328711"/>
      <value value="0.3395955"/>
      <value value="0.3464558"/>
      <value value="0.3534547"/>
      <value value="0.3605949"/>
      <value value="0.3678794"/>
      <value value="0.3753111"/>
      <value value="0.3828929"/>
      <value value="0.3906278"/>
      <value value="0.398519"/>
      <value value="0.4065697"/>
      <value value="0.4147829"/>
      <value value="0.4231621"/>
      <value value="0.4317105"/>
      <value value="0.4404317"/>
      <value value="0.449329"/>
      <value value="0.458406"/>
      <value value="0.4676664"/>
      <value value="0.4771139"/>
      <value value="0.4867523"/>
      <value value="0.4965853"/>
      <value value="0.506617"/>
      <value value="0.5168513"/>
      <value value="0.5272924"/>
      <value value="0.5379444"/>
      <value value="0.5488116"/>
      <value value="0.5598984"/>
      <value value="0.5712091"/>
      <value value="0.5827483"/>
      <value value="0.5945205"/>
      <value value="0.6065307"/>
      <value value="0.6187834"/>
      <value value="0.6312836"/>
      <value value="0.6440364"/>
      <value value="0.6570468"/>
      <value value="0.67032"/>
      <value value="0.6838614"/>
      <value value="0.6976763"/>
      <value value="0.7117703"/>
      <value value="0.726149"/>
      <value value="0.7408182"/>
      <value value="0.7557837"/>
      <value value="0.7710516"/>
      <value value="0.7866279"/>
      <value value="0.8025188"/>
      <value value="0.8187308"/>
      <value value="0.8352702"/>
      <value value="0.8521438"/>
      <value value="0.8693582"/>
      <value value="0.8869204"/>
      <value value="0.9048374"/>
      <value value="0.9231163"/>
      <value value="0.9417645"/>
      <value value="0.9607894"/>
      <value value="0.9801987"/>
      <value value="1"/>
      <value value="1.0202013"/>
      <value value="1.0408108"/>
      <value value="1.0618365"/>
      <value value="1.0832871"/>
      <value value="1.1051709"/>
      <value value="1.1274969"/>
      <value value="1.1502738"/>
      <value value="1.1735109"/>
      <value value="1.1972174"/>
      <value value="1.2214028"/>
      <value value="1.2460767"/>
      <value value="1.2712492"/>
      <value value="1.2969301"/>
      <value value="1.3231298"/>
      <value value="1.3498588"/>
      <value value="1.3771278"/>
      <value value="1.4049476"/>
      <value value="1.4333294"/>
      <value value="1.4622846"/>
      <value value="1.4918247"/>
      <value value="1.5219616"/>
      <value value="1.5527072"/>
      <value value="1.584074"/>
      <value value="1.6160744"/>
      <value value="1.6487213"/>
      <value value="1.6820276"/>
      <value value="1.7160069"/>
      <value value="1.7506725"/>
      <value value="1.7860384"/>
      <value value="1.8221188"/>
      <value value="1.858928"/>
      <value value="1.8964809"/>
      <value value="1.9347923"/>
      <value value="1.9738777"/>
      <value value="2.0137527"/>
      <value value="2.0544332"/>
      <value value="2.0959355"/>
      <value value="2.1382762"/>
      <value value="2.1814723"/>
      <value value="2.2255409"/>
      <value value="2.2704998"/>
      <value value="2.316367"/>
      <value value="2.3631607"/>
      <value value="2.4108997"/>
      <value value="2.4596031"/>
      <value value="2.5092904"/>
      <value value="2.5599814"/>
      <value value="2.6116965"/>
      <value value="2.6644562"/>
      <value value="2.7182818"/>
      <value value="2.7731948"/>
      <value value="2.829217"/>
      <value value="2.886371"/>
      <value value="2.9446796"/>
      <value value="3.004166"/>
      <value value="3.0648542"/>
      <value value="3.1267684"/>
      <value value="3.1899333"/>
      <value value="3.2543742"/>
      <value value="3.3201169"/>
      <value value="3.3871877"/>
      <value value="3.4556135"/>
      <value value="3.5254215"/>
      <value value="3.5966397"/>
      <value value="3.6692967"/>
      <value value="3.7434214"/>
      <value value="3.8190435"/>
      <value value="3.8961933"/>
      <value value="3.9749016"/>
      <value value="4.0552"/>
      <value value="4.1371204"/>
      <value value="4.2206958"/>
      <value value="4.3059595"/>
      <value value="4.3929457"/>
      <value value="4.4816891"/>
      <value value="4.5722252"/>
      <value value="4.6645903"/>
      <value value="4.7588212"/>
      <value value="4.8549558"/>
      <value value="4.9530324"/>
      <value value="5.0530903"/>
      <value value="5.1551695"/>
      <value value="5.2593108"/>
      <value value="5.365556"/>
      <value value="5.4739474"/>
      <value value="5.5845285"/>
      <value value="5.6973434"/>
      <value value="5.8124374"/>
      <value value="5.9298564"/>
      <value value="6.0496475"/>
      <value value="6.1718584"/>
      <value value="6.2965383"/>
      <value value="6.4237368"/>
      <value value="6.5535049"/>
      <value value="6.6858944"/>
      <value value="6.8209585"/>
      <value value="6.958751"/>
      <value value="7.0993271"/>
      <value value="7.242743"/>
      <value value="7.3890561"/>
      <value value="7.5383249"/>
      <value value="7.6906092"/>
      <value value="7.8459698"/>
      <value value="8.0044689"/>
      <value value="8.1661699"/>
      <value value="8.3311375"/>
      <value value="8.4994376"/>
      <value value="8.6711377"/>
      <value value="8.8463063"/>
      <value value="9.0250135"/>
      <value value="9.2073309"/>
      <value value="9.3933313"/>
      <value value="9.5830892"/>
      <value value="9.7766804"/>
      <value value="9.9741825"/>
      <value value="10.1756743"/>
      <value value="10.3812366"/>
      <value value="10.5909515"/>
      <value value="10.8049029"/>
      <value value="11.0231764"/>
      <value value="11.2458593"/>
      <value value="11.4730407"/>
      <value value="11.7048115"/>
      <value value="11.9412644"/>
      <value value="12.182494"/>
      <value value="12.4285967"/>
      <value value="12.679671"/>
      <value value="12.9358173"/>
      <value value="13.1971382"/>
      <value value="13.463738"/>
      <value value="13.7357236"/>
      <value value="14.0132036"/>
      <value value="14.2962891"/>
      <value value="14.5850933"/>
      <value value="14.8797317"/>
      <value value="15.1803222"/>
      <value value="15.4869851"/>
      <value value="15.7998429"/>
      <value value="16.1190209"/>
      <value value="16.4446468"/>
      <value value="16.7768507"/>
      <value value="17.1157655"/>
      <value value="17.4615269"/>
      <value value="17.8142732"/>
      <value value="18.1741454"/>
      <value value="18.5412875"/>
      <value value="18.9158463"/>
      <value value="19.2979718"/>
      <value value="19.6878166"/>
      <value value="20.0855369"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-alpha">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="capacity">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="species">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-population">
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="data-distance">
      <value value="5200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation-probability">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction21">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cutoff">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualization">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-rate">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="range">
      <value value="28.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disturbance">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction12">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="simple run" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="500"/>
    <metric>sum [individuals] of turtles</metric>
    <enumeratedValueSet variable="data-distance">
      <value value="5200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-rate">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="range">
      <value value="28.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation-probability">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-population">
      <value value="0.47"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-b">
      <value value="&quot;single&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-1">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Demographic-Stochasticity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction21">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-1">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="species">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disturbance">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-b">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cutoff">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="capacity">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction12">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-alpha">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="function">
      <value value="&quot;MSS&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-2">
      <value value="14.8797317"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualization">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resolution">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="complexity" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2000"/>
    <metric>ticks</metric>
    <metric>sum [individuals] of species1</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="mutate-b">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Demographic-Stochasticity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction12">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction21">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="species">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-2">
      <value value="14.8797317"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="capacity">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation-probability">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="function">
      <value value="&quot;MSS&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-alpha">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-1">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="range">
      <value value="28.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disturbance">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cutoff">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resolution">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-population">
      <value value="0.47"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualization">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-1">
      <value value="0.1353353"/>
      <value value="0.1380692"/>
      <value value="0.1408584"/>
      <value value="0.1437039"/>
      <value value="0.146607"/>
      <value value="0.1495686"/>
      <value value="0.1525901"/>
      <value value="0.1556726"/>
      <value value="0.1588174"/>
      <value value="0.1620258"/>
      <value value="0.1652989"/>
      <value value="0.1686381"/>
      <value value="0.1720449"/>
      <value value="0.1755204"/>
      <value value="0.1790661"/>
      <value value="0.1826835"/>
      <value value="0.186374"/>
      <value value="0.190139"/>
      <value value="0.19398"/>
      <value value="0.1978987"/>
      <value value="0.2018965"/>
      <value value="0.2059751"/>
      <value value="0.2101361"/>
      <value value="0.2143811"/>
      <value value="0.2187119"/>
      <value value="0.2231302"/>
      <value value="0.2276377"/>
      <value value="0.2322363"/>
      <value value="0.2369278"/>
      <value value="0.241714"/>
      <value value="0.246597"/>
      <value value="0.2515786"/>
      <value value="0.2566608"/>
      <value value="0.2618457"/>
      <value value="0.2671353"/>
      <value value="0.2725318"/>
      <value value="0.2780373"/>
      <value value="0.283654"/>
      <value value="0.2893842"/>
      <value value="0.2952302"/>
      <value value="0.3011942"/>
      <value value="0.3072787"/>
      <value value="0.3134862"/>
      <value value="0.319819"/>
      <value value="0.3262798"/>
      <value value="0.3328711"/>
      <value value="0.3395955"/>
      <value value="0.3464558"/>
      <value value="0.3534547"/>
      <value value="0.3605949"/>
      <value value="0.3678794"/>
      <value value="0.3753111"/>
      <value value="0.3828929"/>
      <value value="0.3906278"/>
      <value value="0.398519"/>
      <value value="0.4065697"/>
      <value value="0.4147829"/>
      <value value="0.4231621"/>
      <value value="0.4317105"/>
      <value value="0.4404317"/>
      <value value="0.449329"/>
      <value value="0.458406"/>
      <value value="0.4676664"/>
      <value value="0.4771139"/>
      <value value="0.4867523"/>
      <value value="0.4965853"/>
      <value value="0.506617"/>
      <value value="0.5168513"/>
      <value value="0.5272924"/>
      <value value="0.5379444"/>
      <value value="0.5488116"/>
      <value value="0.5598984"/>
      <value value="0.5712091"/>
      <value value="0.5827483"/>
      <value value="0.5945205"/>
      <value value="0.6065307"/>
      <value value="0.6187834"/>
      <value value="0.6312836"/>
      <value value="0.6440364"/>
      <value value="0.6570468"/>
      <value value="0.67032"/>
      <value value="0.6838614"/>
      <value value="0.6976763"/>
      <value value="0.7117703"/>
      <value value="0.726149"/>
      <value value="0.7408182"/>
      <value value="0.7557837"/>
      <value value="0.7710516"/>
      <value value="0.7866279"/>
      <value value="0.8025188"/>
      <value value="0.8187308"/>
      <value value="0.8352702"/>
      <value value="0.8521438"/>
      <value value="0.8693582"/>
      <value value="0.8869204"/>
      <value value="0.9048374"/>
      <value value="0.9231163"/>
      <value value="0.9417645"/>
      <value value="0.9607894"/>
      <value value="0.9801987"/>
      <value value="1"/>
      <value value="1.0202013"/>
      <value value="1.0408108"/>
      <value value="1.0618365"/>
      <value value="1.0832871"/>
      <value value="1.1051709"/>
      <value value="1.1274969"/>
      <value value="1.1502738"/>
      <value value="1.1735109"/>
      <value value="1.1972174"/>
      <value value="1.2214028"/>
      <value value="1.2460767"/>
      <value value="1.2712492"/>
      <value value="1.2969301"/>
      <value value="1.3231298"/>
      <value value="1.3498588"/>
      <value value="1.3771278"/>
      <value value="1.4049476"/>
      <value value="1.4333294"/>
      <value value="1.4622846"/>
      <value value="1.4918247"/>
      <value value="1.5219616"/>
      <value value="1.5527072"/>
      <value value="1.584074"/>
      <value value="1.6160744"/>
      <value value="1.6487213"/>
      <value value="1.6820276"/>
      <value value="1.7160069"/>
      <value value="1.7506725"/>
      <value value="1.7860384"/>
      <value value="1.8221188"/>
      <value value="1.858928"/>
      <value value="1.8964809"/>
      <value value="1.9347923"/>
      <value value="1.9738777"/>
      <value value="2.0137527"/>
      <value value="2.0544332"/>
      <value value="2.0959355"/>
      <value value="2.1382762"/>
      <value value="2.1814723"/>
      <value value="2.2255409"/>
      <value value="2.2704998"/>
      <value value="2.316367"/>
      <value value="2.3631607"/>
      <value value="2.4108997"/>
      <value value="2.4596031"/>
      <value value="2.5092904"/>
      <value value="2.5599814"/>
      <value value="2.6116965"/>
      <value value="2.6644562"/>
      <value value="2.7182818"/>
      <value value="2.7731948"/>
      <value value="2.829217"/>
      <value value="2.886371"/>
      <value value="2.9446796"/>
      <value value="3.004166"/>
      <value value="3.0648542"/>
      <value value="3.1267684"/>
      <value value="3.1899333"/>
      <value value="3.2543742"/>
      <value value="3.3201169"/>
      <value value="3.3871877"/>
      <value value="3.4556135"/>
      <value value="3.5254215"/>
      <value value="3.5966397"/>
      <value value="3.6692967"/>
      <value value="3.7434214"/>
      <value value="3.8190435"/>
      <value value="3.8961933"/>
      <value value="3.9749016"/>
      <value value="4.0552"/>
      <value value="4.1371204"/>
      <value value="4.2206958"/>
      <value value="4.3059595"/>
      <value value="4.3929457"/>
      <value value="4.4816891"/>
      <value value="4.5722252"/>
      <value value="4.6645903"/>
      <value value="4.7588212"/>
      <value value="4.8549558"/>
      <value value="4.9530324"/>
      <value value="5.0530903"/>
      <value value="5.1551695"/>
      <value value="5.2593108"/>
      <value value="5.365556"/>
      <value value="5.4739474"/>
      <value value="5.5845285"/>
      <value value="5.6973434"/>
      <value value="5.8124374"/>
      <value value="5.9298564"/>
      <value value="6.0496475"/>
      <value value="6.1718584"/>
      <value value="6.2965383"/>
      <value value="6.4237368"/>
      <value value="6.5535049"/>
      <value value="6.6858944"/>
      <value value="6.8209585"/>
      <value value="6.958751"/>
      <value value="7.0993271"/>
      <value value="7.242743"/>
      <value value="7.3890561"/>
      <value value="7.5383249"/>
      <value value="7.6906092"/>
      <value value="7.8459698"/>
      <value value="8.0044689"/>
      <value value="8.1661699"/>
      <value value="8.3311375"/>
      <value value="8.4994376"/>
      <value value="8.6711377"/>
      <value value="8.8463063"/>
      <value value="9.0250135"/>
      <value value="9.2073309"/>
      <value value="9.3933313"/>
      <value value="9.5830892"/>
      <value value="9.7766804"/>
      <value value="9.9741825"/>
      <value value="10.1756743"/>
      <value value="10.3812366"/>
      <value value="10.5909515"/>
      <value value="10.8049029"/>
      <value value="11.0231764"/>
      <value value="11.2458593"/>
      <value value="11.4730407"/>
      <value value="11.7048115"/>
      <value value="11.9412644"/>
      <value value="12.182494"/>
      <value value="12.4285967"/>
      <value value="12.679671"/>
      <value value="12.9358173"/>
      <value value="13.1971382"/>
      <value value="13.463738"/>
      <value value="13.7357236"/>
      <value value="14.0132036"/>
      <value value="14.2962891"/>
      <value value="14.5850933"/>
      <value value="14.8797317"/>
      <value value="15.1803222"/>
      <value value="15.4869851"/>
      <value value="15.7998429"/>
      <value value="16.1190209"/>
      <value value="16.4446468"/>
      <value value="16.7768507"/>
      <value value="17.1157655"/>
      <value value="17.4615269"/>
      <value value="17.8142732"/>
      <value value="18.1741454"/>
      <value value="18.5412875"/>
      <value value="18.9158463"/>
      <value value="19.2979718"/>
      <value value="19.6878166"/>
      <value value="20.0855369"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-b">
      <value value="&quot;invasion&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="data-distance">
      <value value="5200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-rate">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="invasion-dist-Hassel" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>invasion-fitness-go</go>
    <timeLimit steps="1000"/>
    <metric>count invadors</metric>
    <metric>sum [growth] of invadors</metric>
    <metric>invasioncounter</metric>
    <enumeratedValueSet variable="Demographic-Stochasticity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-1">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-b">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resolution">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="function">
      <value value="&quot;Hassel&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-2">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-b">
      <value value="&quot;invasion&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-1">
      <value value="0.1353353"/>
      <value value="0.1380692"/>
      <value value="0.1408584"/>
      <value value="0.1437039"/>
      <value value="0.146607"/>
      <value value="0.1495686"/>
      <value value="0.1525901"/>
      <value value="0.1556726"/>
      <value value="0.1588174"/>
      <value value="0.1620258"/>
      <value value="0.1652989"/>
      <value value="0.1686381"/>
      <value value="0.1720449"/>
      <value value="0.1755204"/>
      <value value="0.1790661"/>
      <value value="0.1826835"/>
      <value value="0.186374"/>
      <value value="0.190139"/>
      <value value="0.19398"/>
      <value value="0.1978987"/>
      <value value="0.2018965"/>
      <value value="0.2059751"/>
      <value value="0.2101361"/>
      <value value="0.2143811"/>
      <value value="0.2187119"/>
      <value value="0.2231302"/>
      <value value="0.2276377"/>
      <value value="0.2322363"/>
      <value value="0.2369278"/>
      <value value="0.241714"/>
      <value value="0.246597"/>
      <value value="0.2515786"/>
      <value value="0.2566608"/>
      <value value="0.2618457"/>
      <value value="0.2671353"/>
      <value value="0.2725318"/>
      <value value="0.2780373"/>
      <value value="0.283654"/>
      <value value="0.2893842"/>
      <value value="0.2952302"/>
      <value value="0.3011942"/>
      <value value="0.3072787"/>
      <value value="0.3134862"/>
      <value value="0.319819"/>
      <value value="0.3262798"/>
      <value value="0.3328711"/>
      <value value="0.3395955"/>
      <value value="0.3464558"/>
      <value value="0.3534547"/>
      <value value="0.3605949"/>
      <value value="0.3678794"/>
      <value value="0.3753111"/>
      <value value="0.3828929"/>
      <value value="0.3906278"/>
      <value value="0.398519"/>
      <value value="0.4065697"/>
      <value value="0.4147829"/>
      <value value="0.4231621"/>
      <value value="0.4317105"/>
      <value value="0.4404317"/>
      <value value="0.449329"/>
      <value value="0.458406"/>
      <value value="0.4676664"/>
      <value value="0.4771139"/>
      <value value="0.4867523"/>
      <value value="0.4965853"/>
      <value value="0.506617"/>
      <value value="0.5168513"/>
      <value value="0.5272924"/>
      <value value="0.5379444"/>
      <value value="0.5488116"/>
      <value value="0.5598984"/>
      <value value="0.5712091"/>
      <value value="0.5827483"/>
      <value value="0.5945205"/>
      <value value="0.6065307"/>
      <value value="0.6187834"/>
      <value value="0.6312836"/>
      <value value="0.6440364"/>
      <value value="0.6570468"/>
      <value value="0.67032"/>
      <value value="0.6838614"/>
      <value value="0.6976763"/>
      <value value="0.7117703"/>
      <value value="0.726149"/>
      <value value="0.7408182"/>
      <value value="0.7557837"/>
      <value value="0.7710516"/>
      <value value="0.7866279"/>
      <value value="0.8025188"/>
      <value value="0.8187308"/>
      <value value="0.8352702"/>
      <value value="0.8521438"/>
      <value value="0.8693582"/>
      <value value="0.8869204"/>
      <value value="0.9048374"/>
      <value value="0.9231163"/>
      <value value="0.9417645"/>
      <value value="0.9607894"/>
      <value value="0.9801987"/>
      <value value="1"/>
      <value value="1.0202013"/>
      <value value="1.0408108"/>
      <value value="1.0618365"/>
      <value value="1.0832871"/>
      <value value="1.1051709"/>
      <value value="1.1274969"/>
      <value value="1.1502738"/>
      <value value="1.1735109"/>
      <value value="1.1972174"/>
      <value value="1.2214028"/>
      <value value="1.2460767"/>
      <value value="1.2712492"/>
      <value value="1.2969301"/>
      <value value="1.3231298"/>
      <value value="1.3498588"/>
      <value value="1.3771278"/>
      <value value="1.4049476"/>
      <value value="1.4333294"/>
      <value value="1.4622846"/>
      <value value="1.4918247"/>
      <value value="1.5219616"/>
      <value value="1.5527072"/>
      <value value="1.584074"/>
      <value value="1.6160744"/>
      <value value="1.6487213"/>
      <value value="1.6820276"/>
      <value value="1.7160069"/>
      <value value="1.7506725"/>
      <value value="1.7860384"/>
      <value value="1.8221188"/>
      <value value="1.858928"/>
      <value value="1.8964809"/>
      <value value="1.9347923"/>
      <value value="1.9738777"/>
      <value value="2.0137527"/>
      <value value="2.0544332"/>
      <value value="2.0959355"/>
      <value value="2.1382762"/>
      <value value="2.1814723"/>
      <value value="2.2255409"/>
      <value value="2.2704998"/>
      <value value="2.316367"/>
      <value value="2.3631607"/>
      <value value="2.4108997"/>
      <value value="2.4596031"/>
      <value value="2.5092904"/>
      <value value="2.5599814"/>
      <value value="2.6116965"/>
      <value value="2.6644562"/>
      <value value="2.7182818"/>
      <value value="2.7731948"/>
      <value value="2.829217"/>
      <value value="2.886371"/>
      <value value="2.9446796"/>
      <value value="3.004166"/>
      <value value="3.0648542"/>
      <value value="3.1267684"/>
      <value value="3.1899333"/>
      <value value="3.2543742"/>
      <value value="3.3201169"/>
      <value value="3.3871877"/>
      <value value="3.4556135"/>
      <value value="3.5254215"/>
      <value value="3.5966397"/>
      <value value="3.6692967"/>
      <value value="3.7434214"/>
      <value value="3.8190435"/>
      <value value="3.8961933"/>
      <value value="3.9749016"/>
      <value value="4.0552"/>
      <value value="4.1371204"/>
      <value value="4.2206958"/>
      <value value="4.3059595"/>
      <value value="4.3929457"/>
      <value value="4.4816891"/>
      <value value="4.5722252"/>
      <value value="4.6645903"/>
      <value value="4.7588212"/>
      <value value="4.8549558"/>
      <value value="4.9530324"/>
      <value value="5.0530903"/>
      <value value="5.1551695"/>
      <value value="5.2593108"/>
      <value value="5.365556"/>
      <value value="5.4739474"/>
      <value value="5.5845285"/>
      <value value="5.6973434"/>
      <value value="5.8124374"/>
      <value value="5.9298564"/>
      <value value="6.0496475"/>
      <value value="6.1718584"/>
      <value value="6.2965383"/>
      <value value="6.4237368"/>
      <value value="6.5535049"/>
      <value value="6.6858944"/>
      <value value="6.8209585"/>
      <value value="6.958751"/>
      <value value="7.0993271"/>
      <value value="7.242743"/>
      <value value="7.3890561"/>
      <value value="7.5383249"/>
      <value value="7.6906092"/>
      <value value="7.8459698"/>
      <value value="8.0044689"/>
      <value value="8.1661699"/>
      <value value="8.3311375"/>
      <value value="8.4994376"/>
      <value value="8.6711377"/>
      <value value="8.8463063"/>
      <value value="9.0250135"/>
      <value value="9.2073309"/>
      <value value="9.3933313"/>
      <value value="9.5830892"/>
      <value value="9.7766804"/>
      <value value="9.9741825"/>
      <value value="10.1756743"/>
      <value value="10.3812366"/>
      <value value="10.5909515"/>
      <value value="10.8049029"/>
      <value value="11.0231764"/>
      <value value="11.2458593"/>
      <value value="11.4730407"/>
      <value value="11.7048115"/>
      <value value="11.9412644"/>
      <value value="12.182494"/>
      <value value="12.4285967"/>
      <value value="12.679671"/>
      <value value="12.9358173"/>
      <value value="13.1971382"/>
      <value value="13.463738"/>
      <value value="13.7357236"/>
      <value value="14.0132036"/>
      <value value="14.2962891"/>
      <value value="14.5850933"/>
      <value value="14.8797317"/>
      <value value="15.1803222"/>
      <value value="15.4869851"/>
      <value value="15.7998429"/>
      <value value="16.1190209"/>
      <value value="16.4446468"/>
      <value value="16.7768507"/>
      <value value="17.1157655"/>
      <value value="17.4615269"/>
      <value value="17.8142732"/>
      <value value="18.1741454"/>
      <value value="18.5412875"/>
      <value value="18.9158463"/>
      <value value="19.2979718"/>
      <value value="19.6878166"/>
      <value value="20.0855369"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-2">
      <value value="0.1353353"/>
      <value value="0.1380692"/>
      <value value="0.1408584"/>
      <value value="0.1437039"/>
      <value value="0.146607"/>
      <value value="0.1495686"/>
      <value value="0.1525901"/>
      <value value="0.1556726"/>
      <value value="0.1588174"/>
      <value value="0.1620258"/>
      <value value="0.1652989"/>
      <value value="0.1686381"/>
      <value value="0.1720449"/>
      <value value="0.1755204"/>
      <value value="0.1790661"/>
      <value value="0.1826835"/>
      <value value="0.186374"/>
      <value value="0.190139"/>
      <value value="0.19398"/>
      <value value="0.1978987"/>
      <value value="0.2018965"/>
      <value value="0.2059751"/>
      <value value="0.2101361"/>
      <value value="0.2143811"/>
      <value value="0.2187119"/>
      <value value="0.2231302"/>
      <value value="0.2276377"/>
      <value value="0.2322363"/>
      <value value="0.2369278"/>
      <value value="0.241714"/>
      <value value="0.246597"/>
      <value value="0.2515786"/>
      <value value="0.2566608"/>
      <value value="0.2618457"/>
      <value value="0.2671353"/>
      <value value="0.2725318"/>
      <value value="0.2780373"/>
      <value value="0.283654"/>
      <value value="0.2893842"/>
      <value value="0.2952302"/>
      <value value="0.3011942"/>
      <value value="0.3072787"/>
      <value value="0.3134862"/>
      <value value="0.319819"/>
      <value value="0.3262798"/>
      <value value="0.3328711"/>
      <value value="0.3395955"/>
      <value value="0.3464558"/>
      <value value="0.3534547"/>
      <value value="0.3605949"/>
      <value value="0.3678794"/>
      <value value="0.3753111"/>
      <value value="0.3828929"/>
      <value value="0.3906278"/>
      <value value="0.398519"/>
      <value value="0.4065697"/>
      <value value="0.4147829"/>
      <value value="0.4231621"/>
      <value value="0.4317105"/>
      <value value="0.4404317"/>
      <value value="0.449329"/>
      <value value="0.458406"/>
      <value value="0.4676664"/>
      <value value="0.4771139"/>
      <value value="0.4867523"/>
      <value value="0.4965853"/>
      <value value="0.506617"/>
      <value value="0.5168513"/>
      <value value="0.5272924"/>
      <value value="0.5379444"/>
      <value value="0.5488116"/>
      <value value="0.5598984"/>
      <value value="0.5712091"/>
      <value value="0.5827483"/>
      <value value="0.5945205"/>
      <value value="0.6065307"/>
      <value value="0.6187834"/>
      <value value="0.6312836"/>
      <value value="0.6440364"/>
      <value value="0.6570468"/>
      <value value="0.67032"/>
      <value value="0.6838614"/>
      <value value="0.6976763"/>
      <value value="0.7117703"/>
      <value value="0.726149"/>
      <value value="0.7408182"/>
      <value value="0.7557837"/>
      <value value="0.7710516"/>
      <value value="0.7866279"/>
      <value value="0.8025188"/>
      <value value="0.8187308"/>
      <value value="0.8352702"/>
      <value value="0.8521438"/>
      <value value="0.8693582"/>
      <value value="0.8869204"/>
      <value value="0.9048374"/>
      <value value="0.9231163"/>
      <value value="0.9417645"/>
      <value value="0.9607894"/>
      <value value="0.9801987"/>
      <value value="1"/>
      <value value="1.0202013"/>
      <value value="1.0408108"/>
      <value value="1.0618365"/>
      <value value="1.0832871"/>
      <value value="1.1051709"/>
      <value value="1.1274969"/>
      <value value="1.1502738"/>
      <value value="1.1735109"/>
      <value value="1.1972174"/>
      <value value="1.2214028"/>
      <value value="1.2460767"/>
      <value value="1.2712492"/>
      <value value="1.2969301"/>
      <value value="1.3231298"/>
      <value value="1.3498588"/>
      <value value="1.3771278"/>
      <value value="1.4049476"/>
      <value value="1.4333294"/>
      <value value="1.4622846"/>
      <value value="1.4918247"/>
      <value value="1.5219616"/>
      <value value="1.5527072"/>
      <value value="1.584074"/>
      <value value="1.6160744"/>
      <value value="1.6487213"/>
      <value value="1.6820276"/>
      <value value="1.7160069"/>
      <value value="1.7506725"/>
      <value value="1.7860384"/>
      <value value="1.8221188"/>
      <value value="1.858928"/>
      <value value="1.8964809"/>
      <value value="1.9347923"/>
      <value value="1.9738777"/>
      <value value="2.0137527"/>
      <value value="2.0544332"/>
      <value value="2.0959355"/>
      <value value="2.1382762"/>
      <value value="2.1814723"/>
      <value value="2.2255409"/>
      <value value="2.2704998"/>
      <value value="2.316367"/>
      <value value="2.3631607"/>
      <value value="2.4108997"/>
      <value value="2.4596031"/>
      <value value="2.5092904"/>
      <value value="2.5599814"/>
      <value value="2.6116965"/>
      <value value="2.6644562"/>
      <value value="2.7182818"/>
      <value value="2.7731948"/>
      <value value="2.829217"/>
      <value value="2.886371"/>
      <value value="2.9446796"/>
      <value value="3.004166"/>
      <value value="3.0648542"/>
      <value value="3.1267684"/>
      <value value="3.1899333"/>
      <value value="3.2543742"/>
      <value value="3.3201169"/>
      <value value="3.3871877"/>
      <value value="3.4556135"/>
      <value value="3.5254215"/>
      <value value="3.5966397"/>
      <value value="3.6692967"/>
      <value value="3.7434214"/>
      <value value="3.8190435"/>
      <value value="3.8961933"/>
      <value value="3.9749016"/>
      <value value="4.0552"/>
      <value value="4.1371204"/>
      <value value="4.2206958"/>
      <value value="4.3059595"/>
      <value value="4.3929457"/>
      <value value="4.4816891"/>
      <value value="4.5722252"/>
      <value value="4.6645903"/>
      <value value="4.7588212"/>
      <value value="4.8549558"/>
      <value value="4.9530324"/>
      <value value="5.0530903"/>
      <value value="5.1551695"/>
      <value value="5.2593108"/>
      <value value="5.365556"/>
      <value value="5.4739474"/>
      <value value="5.5845285"/>
      <value value="5.6973434"/>
      <value value="5.8124374"/>
      <value value="5.9298564"/>
      <value value="6.0496475"/>
      <value value="6.1718584"/>
      <value value="6.2965383"/>
      <value value="6.4237368"/>
      <value value="6.5535049"/>
      <value value="6.6858944"/>
      <value value="6.8209585"/>
      <value value="6.958751"/>
      <value value="7.0993271"/>
      <value value="7.242743"/>
      <value value="7.3890561"/>
      <value value="7.5383249"/>
      <value value="7.6906092"/>
      <value value="7.8459698"/>
      <value value="8.0044689"/>
      <value value="8.1661699"/>
      <value value="8.3311375"/>
      <value value="8.4994376"/>
      <value value="8.6711377"/>
      <value value="8.8463063"/>
      <value value="9.0250135"/>
      <value value="9.2073309"/>
      <value value="9.3933313"/>
      <value value="9.5830892"/>
      <value value="9.7766804"/>
      <value value="9.9741825"/>
      <value value="10.1756743"/>
      <value value="10.3812366"/>
      <value value="10.5909515"/>
      <value value="10.8049029"/>
      <value value="11.0231764"/>
      <value value="11.2458593"/>
      <value value="11.4730407"/>
      <value value="11.7048115"/>
      <value value="11.9412644"/>
      <value value="12.182494"/>
      <value value="12.4285967"/>
      <value value="12.679671"/>
      <value value="12.9358173"/>
      <value value="13.1971382"/>
      <value value="13.463738"/>
      <value value="13.7357236"/>
      <value value="14.0132036"/>
      <value value="14.2962891"/>
      <value value="14.5850933"/>
      <value value="14.8797317"/>
      <value value="15.1803222"/>
      <value value="15.4869851"/>
      <value value="15.7998429"/>
      <value value="16.1190209"/>
      <value value="16.4446468"/>
      <value value="16.7768507"/>
      <value value="17.1157655"/>
      <value value="17.4615269"/>
      <value value="17.8142732"/>
      <value value="18.1741454"/>
      <value value="18.5412875"/>
      <value value="18.9158463"/>
      <value value="19.2979718"/>
      <value value="19.6878166"/>
      <value value="20.0855369"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-alpha">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="capacity">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="species">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-population">
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="data-distance">
      <value value="5200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation-probability">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction21">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cutoff">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualization">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-rate">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="range">
      <value value="28.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disturbance">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction12">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="invasion-dist-Ricker" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>invasion-fitness-go</go>
    <timeLimit steps="1000"/>
    <metric>count invadors</metric>
    <metric>sum [growth] of invadors</metric>
    <metric>invasioncounter</metric>
    <enumeratedValueSet variable="Demographic-Stochasticity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-1">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-b">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resolution">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="function">
      <value value="&quot;Ricker&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-2">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-b">
      <value value="&quot;invasion&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-1">
      <value value="0.1353353"/>
      <value value="0.1380692"/>
      <value value="0.1408584"/>
      <value value="0.1437039"/>
      <value value="0.146607"/>
      <value value="0.1495686"/>
      <value value="0.1525901"/>
      <value value="0.1556726"/>
      <value value="0.1588174"/>
      <value value="0.1620258"/>
      <value value="0.1652989"/>
      <value value="0.1686381"/>
      <value value="0.1720449"/>
      <value value="0.1755204"/>
      <value value="0.1790661"/>
      <value value="0.1826835"/>
      <value value="0.186374"/>
      <value value="0.190139"/>
      <value value="0.19398"/>
      <value value="0.1978987"/>
      <value value="0.2018965"/>
      <value value="0.2059751"/>
      <value value="0.2101361"/>
      <value value="0.2143811"/>
      <value value="0.2187119"/>
      <value value="0.2231302"/>
      <value value="0.2276377"/>
      <value value="0.2322363"/>
      <value value="0.2369278"/>
      <value value="0.241714"/>
      <value value="0.246597"/>
      <value value="0.2515786"/>
      <value value="0.2566608"/>
      <value value="0.2618457"/>
      <value value="0.2671353"/>
      <value value="0.2725318"/>
      <value value="0.2780373"/>
      <value value="0.283654"/>
      <value value="0.2893842"/>
      <value value="0.2952302"/>
      <value value="0.3011942"/>
      <value value="0.3072787"/>
      <value value="0.3134862"/>
      <value value="0.319819"/>
      <value value="0.3262798"/>
      <value value="0.3328711"/>
      <value value="0.3395955"/>
      <value value="0.3464558"/>
      <value value="0.3534547"/>
      <value value="0.3605949"/>
      <value value="0.3678794"/>
      <value value="0.3753111"/>
      <value value="0.3828929"/>
      <value value="0.3906278"/>
      <value value="0.398519"/>
      <value value="0.4065697"/>
      <value value="0.4147829"/>
      <value value="0.4231621"/>
      <value value="0.4317105"/>
      <value value="0.4404317"/>
      <value value="0.449329"/>
      <value value="0.458406"/>
      <value value="0.4676664"/>
      <value value="0.4771139"/>
      <value value="0.4867523"/>
      <value value="0.4965853"/>
      <value value="0.506617"/>
      <value value="0.5168513"/>
      <value value="0.5272924"/>
      <value value="0.5379444"/>
      <value value="0.5488116"/>
      <value value="0.5598984"/>
      <value value="0.5712091"/>
      <value value="0.5827483"/>
      <value value="0.5945205"/>
      <value value="0.6065307"/>
      <value value="0.6187834"/>
      <value value="0.6312836"/>
      <value value="0.6440364"/>
      <value value="0.6570468"/>
      <value value="0.67032"/>
      <value value="0.6838614"/>
      <value value="0.6976763"/>
      <value value="0.7117703"/>
      <value value="0.726149"/>
      <value value="0.7408182"/>
      <value value="0.7557837"/>
      <value value="0.7710516"/>
      <value value="0.7866279"/>
      <value value="0.8025188"/>
      <value value="0.8187308"/>
      <value value="0.8352702"/>
      <value value="0.8521438"/>
      <value value="0.8693582"/>
      <value value="0.8869204"/>
      <value value="0.9048374"/>
      <value value="0.9231163"/>
      <value value="0.9417645"/>
      <value value="0.9607894"/>
      <value value="0.9801987"/>
      <value value="1"/>
      <value value="1.0202013"/>
      <value value="1.0408108"/>
      <value value="1.0618365"/>
      <value value="1.0832871"/>
      <value value="1.1051709"/>
      <value value="1.1274969"/>
      <value value="1.1502738"/>
      <value value="1.1735109"/>
      <value value="1.1972174"/>
      <value value="1.2214028"/>
      <value value="1.2460767"/>
      <value value="1.2712492"/>
      <value value="1.2969301"/>
      <value value="1.3231298"/>
      <value value="1.3498588"/>
      <value value="1.3771278"/>
      <value value="1.4049476"/>
      <value value="1.4333294"/>
      <value value="1.4622846"/>
      <value value="1.4918247"/>
      <value value="1.5219616"/>
      <value value="1.5527072"/>
      <value value="1.584074"/>
      <value value="1.6160744"/>
      <value value="1.6487213"/>
      <value value="1.6820276"/>
      <value value="1.7160069"/>
      <value value="1.7506725"/>
      <value value="1.7860384"/>
      <value value="1.8221188"/>
      <value value="1.858928"/>
      <value value="1.8964809"/>
      <value value="1.9347923"/>
      <value value="1.9738777"/>
      <value value="2.0137527"/>
      <value value="2.0544332"/>
      <value value="2.0959355"/>
      <value value="2.1382762"/>
      <value value="2.1814723"/>
      <value value="2.2255409"/>
      <value value="2.2704998"/>
      <value value="2.316367"/>
      <value value="2.3631607"/>
      <value value="2.4108997"/>
      <value value="2.4596031"/>
      <value value="2.5092904"/>
      <value value="2.5599814"/>
      <value value="2.6116965"/>
      <value value="2.6644562"/>
      <value value="2.7182818"/>
      <value value="2.7731948"/>
      <value value="2.829217"/>
      <value value="2.886371"/>
      <value value="2.9446796"/>
      <value value="3.004166"/>
      <value value="3.0648542"/>
      <value value="3.1267684"/>
      <value value="3.1899333"/>
      <value value="3.2543742"/>
      <value value="3.3201169"/>
      <value value="3.3871877"/>
      <value value="3.4556135"/>
      <value value="3.5254215"/>
      <value value="3.5966397"/>
      <value value="3.6692967"/>
      <value value="3.7434214"/>
      <value value="3.8190435"/>
      <value value="3.8961933"/>
      <value value="3.9749016"/>
      <value value="4.0552"/>
      <value value="4.1371204"/>
      <value value="4.2206958"/>
      <value value="4.3059595"/>
      <value value="4.3929457"/>
      <value value="4.4816891"/>
      <value value="4.5722252"/>
      <value value="4.6645903"/>
      <value value="4.7588212"/>
      <value value="4.8549558"/>
      <value value="4.9530324"/>
      <value value="5.0530903"/>
      <value value="5.1551695"/>
      <value value="5.2593108"/>
      <value value="5.365556"/>
      <value value="5.4739474"/>
      <value value="5.5845285"/>
      <value value="5.6973434"/>
      <value value="5.8124374"/>
      <value value="5.9298564"/>
      <value value="6.0496475"/>
      <value value="6.1718584"/>
      <value value="6.2965383"/>
      <value value="6.4237368"/>
      <value value="6.5535049"/>
      <value value="6.6858944"/>
      <value value="6.8209585"/>
      <value value="6.958751"/>
      <value value="7.0993271"/>
      <value value="7.242743"/>
      <value value="7.3890561"/>
      <value value="7.5383249"/>
      <value value="7.6906092"/>
      <value value="7.8459698"/>
      <value value="8.0044689"/>
      <value value="8.1661699"/>
      <value value="8.3311375"/>
      <value value="8.4994376"/>
      <value value="8.6711377"/>
      <value value="8.8463063"/>
      <value value="9.0250135"/>
      <value value="9.2073309"/>
      <value value="9.3933313"/>
      <value value="9.5830892"/>
      <value value="9.7766804"/>
      <value value="9.9741825"/>
      <value value="10.1756743"/>
      <value value="10.3812366"/>
      <value value="10.5909515"/>
      <value value="10.8049029"/>
      <value value="11.0231764"/>
      <value value="11.2458593"/>
      <value value="11.4730407"/>
      <value value="11.7048115"/>
      <value value="11.9412644"/>
      <value value="12.182494"/>
      <value value="12.4285967"/>
      <value value="12.679671"/>
      <value value="12.9358173"/>
      <value value="13.1971382"/>
      <value value="13.463738"/>
      <value value="13.7357236"/>
      <value value="14.0132036"/>
      <value value="14.2962891"/>
      <value value="14.5850933"/>
      <value value="14.8797317"/>
      <value value="15.1803222"/>
      <value value="15.4869851"/>
      <value value="15.7998429"/>
      <value value="16.1190209"/>
      <value value="16.4446468"/>
      <value value="16.7768507"/>
      <value value="17.1157655"/>
      <value value="17.4615269"/>
      <value value="17.8142732"/>
      <value value="18.1741454"/>
      <value value="18.5412875"/>
      <value value="18.9158463"/>
      <value value="19.2979718"/>
      <value value="19.6878166"/>
      <value value="20.0855369"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-2">
      <value value="0.1353353"/>
      <value value="0.1380692"/>
      <value value="0.1408584"/>
      <value value="0.1437039"/>
      <value value="0.146607"/>
      <value value="0.1495686"/>
      <value value="0.1525901"/>
      <value value="0.1556726"/>
      <value value="0.1588174"/>
      <value value="0.1620258"/>
      <value value="0.1652989"/>
      <value value="0.1686381"/>
      <value value="0.1720449"/>
      <value value="0.1755204"/>
      <value value="0.1790661"/>
      <value value="0.1826835"/>
      <value value="0.186374"/>
      <value value="0.190139"/>
      <value value="0.19398"/>
      <value value="0.1978987"/>
      <value value="0.2018965"/>
      <value value="0.2059751"/>
      <value value="0.2101361"/>
      <value value="0.2143811"/>
      <value value="0.2187119"/>
      <value value="0.2231302"/>
      <value value="0.2276377"/>
      <value value="0.2322363"/>
      <value value="0.2369278"/>
      <value value="0.241714"/>
      <value value="0.246597"/>
      <value value="0.2515786"/>
      <value value="0.2566608"/>
      <value value="0.2618457"/>
      <value value="0.2671353"/>
      <value value="0.2725318"/>
      <value value="0.2780373"/>
      <value value="0.283654"/>
      <value value="0.2893842"/>
      <value value="0.2952302"/>
      <value value="0.3011942"/>
      <value value="0.3072787"/>
      <value value="0.3134862"/>
      <value value="0.319819"/>
      <value value="0.3262798"/>
      <value value="0.3328711"/>
      <value value="0.3395955"/>
      <value value="0.3464558"/>
      <value value="0.3534547"/>
      <value value="0.3605949"/>
      <value value="0.3678794"/>
      <value value="0.3753111"/>
      <value value="0.3828929"/>
      <value value="0.3906278"/>
      <value value="0.398519"/>
      <value value="0.4065697"/>
      <value value="0.4147829"/>
      <value value="0.4231621"/>
      <value value="0.4317105"/>
      <value value="0.4404317"/>
      <value value="0.449329"/>
      <value value="0.458406"/>
      <value value="0.4676664"/>
      <value value="0.4771139"/>
      <value value="0.4867523"/>
      <value value="0.4965853"/>
      <value value="0.506617"/>
      <value value="0.5168513"/>
      <value value="0.5272924"/>
      <value value="0.5379444"/>
      <value value="0.5488116"/>
      <value value="0.5598984"/>
      <value value="0.5712091"/>
      <value value="0.5827483"/>
      <value value="0.5945205"/>
      <value value="0.6065307"/>
      <value value="0.6187834"/>
      <value value="0.6312836"/>
      <value value="0.6440364"/>
      <value value="0.6570468"/>
      <value value="0.67032"/>
      <value value="0.6838614"/>
      <value value="0.6976763"/>
      <value value="0.7117703"/>
      <value value="0.726149"/>
      <value value="0.7408182"/>
      <value value="0.7557837"/>
      <value value="0.7710516"/>
      <value value="0.7866279"/>
      <value value="0.8025188"/>
      <value value="0.8187308"/>
      <value value="0.8352702"/>
      <value value="0.8521438"/>
      <value value="0.8693582"/>
      <value value="0.8869204"/>
      <value value="0.9048374"/>
      <value value="0.9231163"/>
      <value value="0.9417645"/>
      <value value="0.9607894"/>
      <value value="0.9801987"/>
      <value value="1"/>
      <value value="1.0202013"/>
      <value value="1.0408108"/>
      <value value="1.0618365"/>
      <value value="1.0832871"/>
      <value value="1.1051709"/>
      <value value="1.1274969"/>
      <value value="1.1502738"/>
      <value value="1.1735109"/>
      <value value="1.1972174"/>
      <value value="1.2214028"/>
      <value value="1.2460767"/>
      <value value="1.2712492"/>
      <value value="1.2969301"/>
      <value value="1.3231298"/>
      <value value="1.3498588"/>
      <value value="1.3771278"/>
      <value value="1.4049476"/>
      <value value="1.4333294"/>
      <value value="1.4622846"/>
      <value value="1.4918247"/>
      <value value="1.5219616"/>
      <value value="1.5527072"/>
      <value value="1.584074"/>
      <value value="1.6160744"/>
      <value value="1.6487213"/>
      <value value="1.6820276"/>
      <value value="1.7160069"/>
      <value value="1.7506725"/>
      <value value="1.7860384"/>
      <value value="1.8221188"/>
      <value value="1.858928"/>
      <value value="1.8964809"/>
      <value value="1.9347923"/>
      <value value="1.9738777"/>
      <value value="2.0137527"/>
      <value value="2.0544332"/>
      <value value="2.0959355"/>
      <value value="2.1382762"/>
      <value value="2.1814723"/>
      <value value="2.2255409"/>
      <value value="2.2704998"/>
      <value value="2.316367"/>
      <value value="2.3631607"/>
      <value value="2.4108997"/>
      <value value="2.4596031"/>
      <value value="2.5092904"/>
      <value value="2.5599814"/>
      <value value="2.6116965"/>
      <value value="2.6644562"/>
      <value value="2.7182818"/>
      <value value="2.7731948"/>
      <value value="2.829217"/>
      <value value="2.886371"/>
      <value value="2.9446796"/>
      <value value="3.004166"/>
      <value value="3.0648542"/>
      <value value="3.1267684"/>
      <value value="3.1899333"/>
      <value value="3.2543742"/>
      <value value="3.3201169"/>
      <value value="3.3871877"/>
      <value value="3.4556135"/>
      <value value="3.5254215"/>
      <value value="3.5966397"/>
      <value value="3.6692967"/>
      <value value="3.7434214"/>
      <value value="3.8190435"/>
      <value value="3.8961933"/>
      <value value="3.9749016"/>
      <value value="4.0552"/>
      <value value="4.1371204"/>
      <value value="4.2206958"/>
      <value value="4.3059595"/>
      <value value="4.3929457"/>
      <value value="4.4816891"/>
      <value value="4.5722252"/>
      <value value="4.6645903"/>
      <value value="4.7588212"/>
      <value value="4.8549558"/>
      <value value="4.9530324"/>
      <value value="5.0530903"/>
      <value value="5.1551695"/>
      <value value="5.2593108"/>
      <value value="5.365556"/>
      <value value="5.4739474"/>
      <value value="5.5845285"/>
      <value value="5.6973434"/>
      <value value="5.8124374"/>
      <value value="5.9298564"/>
      <value value="6.0496475"/>
      <value value="6.1718584"/>
      <value value="6.2965383"/>
      <value value="6.4237368"/>
      <value value="6.5535049"/>
      <value value="6.6858944"/>
      <value value="6.8209585"/>
      <value value="6.958751"/>
      <value value="7.0993271"/>
      <value value="7.242743"/>
      <value value="7.3890561"/>
      <value value="7.5383249"/>
      <value value="7.6906092"/>
      <value value="7.8459698"/>
      <value value="8.0044689"/>
      <value value="8.1661699"/>
      <value value="8.3311375"/>
      <value value="8.4994376"/>
      <value value="8.6711377"/>
      <value value="8.8463063"/>
      <value value="9.0250135"/>
      <value value="9.2073309"/>
      <value value="9.3933313"/>
      <value value="9.5830892"/>
      <value value="9.7766804"/>
      <value value="9.9741825"/>
      <value value="10.1756743"/>
      <value value="10.3812366"/>
      <value value="10.5909515"/>
      <value value="10.8049029"/>
      <value value="11.0231764"/>
      <value value="11.2458593"/>
      <value value="11.4730407"/>
      <value value="11.7048115"/>
      <value value="11.9412644"/>
      <value value="12.182494"/>
      <value value="12.4285967"/>
      <value value="12.679671"/>
      <value value="12.9358173"/>
      <value value="13.1971382"/>
      <value value="13.463738"/>
      <value value="13.7357236"/>
      <value value="14.0132036"/>
      <value value="14.2962891"/>
      <value value="14.5850933"/>
      <value value="14.8797317"/>
      <value value="15.1803222"/>
      <value value="15.4869851"/>
      <value value="15.7998429"/>
      <value value="16.1190209"/>
      <value value="16.4446468"/>
      <value value="16.7768507"/>
      <value value="17.1157655"/>
      <value value="17.4615269"/>
      <value value="17.8142732"/>
      <value value="18.1741454"/>
      <value value="18.5412875"/>
      <value value="18.9158463"/>
      <value value="19.2979718"/>
      <value value="19.6878166"/>
      <value value="20.0855369"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-alpha">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="capacity">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="species">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-population">
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="data-distance">
      <value value="5200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation-probability">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction21">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cutoff">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualization">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-rate">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="range">
      <value value="28.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disturbance">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction12">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="coexistence-dynamics" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="10000"/>
    <metric>population1</metric>
    <metric>population2</metric>
    <enumeratedValueSet variable="cutoff">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="data-distance">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-b">
      <value value="&quot;double&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction12">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-rate">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="capacity">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction21">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-1">
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="range">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Demographic-Stochasticity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation-probability">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-1">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resolution">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-alpha">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="species">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-population">
      <value value="0.47"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-b">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-2">
      <value value="5.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="function">
      <value value="&quot;MSS&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualization">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disturbance">
      <value value="0.05"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="competitive-exclusion" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <exitCondition>number-of-strategies &lt; 2</exitCondition>
    <metric>(count species1) &gt; 0</metric>
    <metric>(count species2) &gt; 0</metric>
    <metric>ticks</metric>
    <enumeratedValueSet variable="Demographic-Stochasticity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-rate">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="species">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-b">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="data-distance">
      <value value="5200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resolution">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="function">
      <value value="&quot;MSS&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disturbance">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction12">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation-probability">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="range">
      <value value="28.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cutoff">
      <value value="0.17"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="capacity">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-population">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-1">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-alpha">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualization">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction21">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-1">
      <value value="0.1353353"/>
      <value value="0.1380692"/>
      <value value="0.1408584"/>
      <value value="0.1437039"/>
      <value value="0.146607"/>
      <value value="0.1495686"/>
      <value value="0.1525901"/>
      <value value="0.1556726"/>
      <value value="0.1588174"/>
      <value value="0.1620258"/>
      <value value="0.1652989"/>
      <value value="0.1686381"/>
      <value value="0.1720449"/>
      <value value="0.1755204"/>
      <value value="0.1790661"/>
      <value value="0.1826835"/>
      <value value="0.186374"/>
      <value value="0.190139"/>
      <value value="0.19398"/>
      <value value="0.1978987"/>
      <value value="0.2018965"/>
      <value value="0.2059751"/>
      <value value="0.2101361"/>
      <value value="0.2143811"/>
      <value value="0.2187119"/>
      <value value="0.2231302"/>
      <value value="0.2276377"/>
      <value value="0.2322363"/>
      <value value="0.2369278"/>
      <value value="0.241714"/>
      <value value="0.246597"/>
      <value value="0.2515786"/>
      <value value="0.2566608"/>
      <value value="0.2618457"/>
      <value value="0.2671353"/>
      <value value="0.2725318"/>
      <value value="0.2780373"/>
      <value value="0.283654"/>
      <value value="0.2893842"/>
      <value value="0.2952302"/>
      <value value="0.3011942"/>
      <value value="0.3072787"/>
      <value value="0.3134862"/>
      <value value="0.319819"/>
      <value value="0.3262798"/>
      <value value="0.3328711"/>
      <value value="0.3395955"/>
      <value value="0.3464558"/>
      <value value="0.3534547"/>
      <value value="0.3605949"/>
      <value value="0.3678794"/>
      <value value="0.3753111"/>
      <value value="0.3828929"/>
      <value value="0.3906278"/>
      <value value="0.398519"/>
      <value value="0.4065697"/>
      <value value="0.4147829"/>
      <value value="0.4231621"/>
      <value value="0.4317105"/>
      <value value="0.4404317"/>
      <value value="0.449329"/>
      <value value="0.458406"/>
      <value value="0.4676664"/>
      <value value="0.4771139"/>
      <value value="0.4867523"/>
      <value value="0.4965853"/>
      <value value="0.506617"/>
      <value value="0.5168513"/>
      <value value="0.5272924"/>
      <value value="0.5379444"/>
      <value value="0.5488116"/>
      <value value="0.5598984"/>
      <value value="0.5712091"/>
      <value value="0.5827483"/>
      <value value="0.5945205"/>
      <value value="0.6065307"/>
      <value value="0.6187834"/>
      <value value="0.6312836"/>
      <value value="0.6440364"/>
      <value value="0.6570468"/>
      <value value="0.67032"/>
      <value value="0.6838614"/>
      <value value="0.6976763"/>
      <value value="0.7117703"/>
      <value value="0.726149"/>
      <value value="0.7408182"/>
      <value value="0.7557837"/>
      <value value="0.7710516"/>
      <value value="0.7866279"/>
      <value value="0.8025188"/>
      <value value="0.8187308"/>
      <value value="0.8352702"/>
      <value value="0.8521438"/>
      <value value="0.8693582"/>
      <value value="0.8869204"/>
      <value value="0.9048374"/>
      <value value="0.9231163"/>
      <value value="0.9417645"/>
      <value value="0.9607894"/>
      <value value="0.9801987"/>
      <value value="1"/>
      <value value="1.0202013"/>
      <value value="1.0408108"/>
      <value value="1.0618365"/>
      <value value="1.0832871"/>
      <value value="1.1051709"/>
      <value value="1.1274969"/>
      <value value="1.1502738"/>
      <value value="1.1735109"/>
      <value value="1.1972174"/>
      <value value="1.2214028"/>
      <value value="1.2460767"/>
      <value value="1.2712492"/>
      <value value="1.2969301"/>
      <value value="1.3231298"/>
      <value value="1.3498588"/>
      <value value="1.3771278"/>
      <value value="1.4049476"/>
      <value value="1.4333294"/>
      <value value="1.4622846"/>
      <value value="1.4918247"/>
      <value value="1.5219616"/>
      <value value="1.5527072"/>
      <value value="1.584074"/>
      <value value="1.6160744"/>
      <value value="1.6487213"/>
      <value value="1.6820276"/>
      <value value="1.7160069"/>
      <value value="1.7506725"/>
      <value value="1.7860384"/>
      <value value="1.8221188"/>
      <value value="1.858928"/>
      <value value="1.8964809"/>
      <value value="1.9347923"/>
      <value value="1.9738777"/>
      <value value="2.0137527"/>
      <value value="2.0544332"/>
      <value value="2.0959355"/>
      <value value="2.1382762"/>
      <value value="2.1814723"/>
      <value value="2.2255409"/>
      <value value="2.2704998"/>
      <value value="2.316367"/>
      <value value="2.3631607"/>
      <value value="2.4108997"/>
      <value value="2.4596031"/>
      <value value="2.5092904"/>
      <value value="2.5599814"/>
      <value value="2.6116965"/>
      <value value="2.6644562"/>
      <value value="2.7182818"/>
      <value value="2.7731948"/>
      <value value="2.829217"/>
      <value value="2.886371"/>
      <value value="2.9446796"/>
      <value value="3.004166"/>
      <value value="3.0648542"/>
      <value value="3.1267684"/>
      <value value="3.1899333"/>
      <value value="3.2543742"/>
      <value value="3.3201169"/>
      <value value="3.3871877"/>
      <value value="3.4556135"/>
      <value value="3.5254215"/>
      <value value="3.5966397"/>
      <value value="3.6692967"/>
      <value value="3.7434214"/>
      <value value="3.8190435"/>
      <value value="3.8961933"/>
      <value value="3.9749016"/>
      <value value="4.0552"/>
      <value value="4.1371204"/>
      <value value="4.2206958"/>
      <value value="4.3059595"/>
      <value value="4.3929457"/>
      <value value="4.4816891"/>
      <value value="4.5722252"/>
      <value value="4.6645903"/>
      <value value="4.7588212"/>
      <value value="4.8549558"/>
      <value value="4.9530324"/>
      <value value="5.0530903"/>
      <value value="5.1551695"/>
      <value value="5.2593108"/>
      <value value="5.365556"/>
      <value value="5.4739474"/>
      <value value="5.5845285"/>
      <value value="5.6973434"/>
      <value value="5.8124374"/>
      <value value="5.9298564"/>
      <value value="6.0496475"/>
      <value value="6.1718584"/>
      <value value="6.2965383"/>
      <value value="6.4237368"/>
      <value value="6.5535049"/>
      <value value="6.6858944"/>
      <value value="6.8209585"/>
      <value value="6.958751"/>
      <value value="7.0993271"/>
      <value value="7.242743"/>
      <value value="7.3890561"/>
      <value value="7.5383249"/>
      <value value="7.6906092"/>
      <value value="7.8459698"/>
      <value value="8.0044689"/>
      <value value="8.1661699"/>
      <value value="8.3311375"/>
      <value value="8.4994376"/>
      <value value="8.6711377"/>
      <value value="8.8463063"/>
      <value value="9.0250135"/>
      <value value="9.2073309"/>
      <value value="9.3933313"/>
      <value value="9.5830892"/>
      <value value="9.7766804"/>
      <value value="9.9741825"/>
      <value value="10.1756743"/>
      <value value="10.3812366"/>
      <value value="10.5909515"/>
      <value value="10.8049029"/>
      <value value="11.0231764"/>
      <value value="11.2458593"/>
      <value value="11.4730407"/>
      <value value="11.7048115"/>
      <value value="11.9412644"/>
      <value value="12.182494"/>
      <value value="12.4285967"/>
      <value value="12.679671"/>
      <value value="12.9358173"/>
      <value value="13.1971382"/>
      <value value="13.463738"/>
      <value value="13.7357236"/>
      <value value="14.0132036"/>
      <value value="14.2962891"/>
      <value value="14.5850933"/>
      <value value="14.8797317"/>
      <value value="15.1803222"/>
      <value value="15.4869851"/>
      <value value="15.7998429"/>
      <value value="16.1190209"/>
      <value value="16.4446468"/>
      <value value="16.7768507"/>
      <value value="17.1157655"/>
      <value value="17.4615269"/>
      <value value="17.8142732"/>
      <value value="18.1741454"/>
      <value value="18.5412875"/>
      <value value="18.9158463"/>
      <value value="19.2979718"/>
      <value value="19.6878166"/>
      <value value="20.0855369"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-2">
      <value value="0.1353353"/>
      <value value="0.1380692"/>
      <value value="0.1408584"/>
      <value value="0.1437039"/>
      <value value="0.146607"/>
      <value value="0.1495686"/>
      <value value="0.1525901"/>
      <value value="0.1556726"/>
      <value value="0.1588174"/>
      <value value="0.1620258"/>
      <value value="0.1652989"/>
      <value value="0.1686381"/>
      <value value="0.1720449"/>
      <value value="0.1755204"/>
      <value value="0.1790661"/>
      <value value="0.1826835"/>
      <value value="0.186374"/>
      <value value="0.190139"/>
      <value value="0.19398"/>
      <value value="0.1978987"/>
      <value value="0.2018965"/>
      <value value="0.2059751"/>
      <value value="0.2101361"/>
      <value value="0.2143811"/>
      <value value="0.2187119"/>
      <value value="0.2231302"/>
      <value value="0.2276377"/>
      <value value="0.2322363"/>
      <value value="0.2369278"/>
      <value value="0.241714"/>
      <value value="0.246597"/>
      <value value="0.2515786"/>
      <value value="0.2566608"/>
      <value value="0.2618457"/>
      <value value="0.2671353"/>
      <value value="0.2725318"/>
      <value value="0.2780373"/>
      <value value="0.283654"/>
      <value value="0.2893842"/>
      <value value="0.2952302"/>
      <value value="0.3011942"/>
      <value value="0.3072787"/>
      <value value="0.3134862"/>
      <value value="0.319819"/>
      <value value="0.3262798"/>
      <value value="0.3328711"/>
      <value value="0.3395955"/>
      <value value="0.3464558"/>
      <value value="0.3534547"/>
      <value value="0.3605949"/>
      <value value="0.3678794"/>
      <value value="0.3753111"/>
      <value value="0.3828929"/>
      <value value="0.3906278"/>
      <value value="0.398519"/>
      <value value="0.4065697"/>
      <value value="0.4147829"/>
      <value value="0.4231621"/>
      <value value="0.4317105"/>
      <value value="0.4404317"/>
      <value value="0.449329"/>
      <value value="0.458406"/>
      <value value="0.4676664"/>
      <value value="0.4771139"/>
      <value value="0.4867523"/>
      <value value="0.4965853"/>
      <value value="0.506617"/>
      <value value="0.5168513"/>
      <value value="0.5272924"/>
      <value value="0.5379444"/>
      <value value="0.5488116"/>
      <value value="0.5598984"/>
      <value value="0.5712091"/>
      <value value="0.5827483"/>
      <value value="0.5945205"/>
      <value value="0.6065307"/>
      <value value="0.6187834"/>
      <value value="0.6312836"/>
      <value value="0.6440364"/>
      <value value="0.6570468"/>
      <value value="0.67032"/>
      <value value="0.6838614"/>
      <value value="0.6976763"/>
      <value value="0.7117703"/>
      <value value="0.726149"/>
      <value value="0.7408182"/>
      <value value="0.7557837"/>
      <value value="0.7710516"/>
      <value value="0.7866279"/>
      <value value="0.8025188"/>
      <value value="0.8187308"/>
      <value value="0.8352702"/>
      <value value="0.8521438"/>
      <value value="0.8693582"/>
      <value value="0.8869204"/>
      <value value="0.9048374"/>
      <value value="0.9231163"/>
      <value value="0.9417645"/>
      <value value="0.9607894"/>
      <value value="0.9801987"/>
      <value value="1"/>
      <value value="1.0202013"/>
      <value value="1.0408108"/>
      <value value="1.0618365"/>
      <value value="1.0832871"/>
      <value value="1.1051709"/>
      <value value="1.1274969"/>
      <value value="1.1502738"/>
      <value value="1.1735109"/>
      <value value="1.1972174"/>
      <value value="1.2214028"/>
      <value value="1.2460767"/>
      <value value="1.2712492"/>
      <value value="1.2969301"/>
      <value value="1.3231298"/>
      <value value="1.3498588"/>
      <value value="1.3771278"/>
      <value value="1.4049476"/>
      <value value="1.4333294"/>
      <value value="1.4622846"/>
      <value value="1.4918247"/>
      <value value="1.5219616"/>
      <value value="1.5527072"/>
      <value value="1.584074"/>
      <value value="1.6160744"/>
      <value value="1.6487213"/>
      <value value="1.6820276"/>
      <value value="1.7160069"/>
      <value value="1.7506725"/>
      <value value="1.7860384"/>
      <value value="1.8221188"/>
      <value value="1.858928"/>
      <value value="1.8964809"/>
      <value value="1.9347923"/>
      <value value="1.9738777"/>
      <value value="2.0137527"/>
      <value value="2.0544332"/>
      <value value="2.0959355"/>
      <value value="2.1382762"/>
      <value value="2.1814723"/>
      <value value="2.2255409"/>
      <value value="2.2704998"/>
      <value value="2.316367"/>
      <value value="2.3631607"/>
      <value value="2.4108997"/>
      <value value="2.4596031"/>
      <value value="2.5092904"/>
      <value value="2.5599814"/>
      <value value="2.6116965"/>
      <value value="2.6644562"/>
      <value value="2.7182818"/>
      <value value="2.7731948"/>
      <value value="2.829217"/>
      <value value="2.886371"/>
      <value value="2.9446796"/>
      <value value="3.004166"/>
      <value value="3.0648542"/>
      <value value="3.1267684"/>
      <value value="3.1899333"/>
      <value value="3.2543742"/>
      <value value="3.3201169"/>
      <value value="3.3871877"/>
      <value value="3.4556135"/>
      <value value="3.5254215"/>
      <value value="3.5966397"/>
      <value value="3.6692967"/>
      <value value="3.7434214"/>
      <value value="3.8190435"/>
      <value value="3.8961933"/>
      <value value="3.9749016"/>
      <value value="4.0552"/>
      <value value="4.1371204"/>
      <value value="4.2206958"/>
      <value value="4.3059595"/>
      <value value="4.3929457"/>
      <value value="4.4816891"/>
      <value value="4.5722252"/>
      <value value="4.6645903"/>
      <value value="4.7588212"/>
      <value value="4.8549558"/>
      <value value="4.9530324"/>
      <value value="5.0530903"/>
      <value value="5.1551695"/>
      <value value="5.2593108"/>
      <value value="5.365556"/>
      <value value="5.4739474"/>
      <value value="5.5845285"/>
      <value value="5.6973434"/>
      <value value="5.8124374"/>
      <value value="5.9298564"/>
      <value value="6.0496475"/>
      <value value="6.1718584"/>
      <value value="6.2965383"/>
      <value value="6.4237368"/>
      <value value="6.5535049"/>
      <value value="6.6858944"/>
      <value value="6.8209585"/>
      <value value="6.958751"/>
      <value value="7.0993271"/>
      <value value="7.242743"/>
      <value value="7.3890561"/>
      <value value="7.5383249"/>
      <value value="7.6906092"/>
      <value value="7.8459698"/>
      <value value="8.0044689"/>
      <value value="8.1661699"/>
      <value value="8.3311375"/>
      <value value="8.4994376"/>
      <value value="8.6711377"/>
      <value value="8.8463063"/>
      <value value="9.0250135"/>
      <value value="9.2073309"/>
      <value value="9.3933313"/>
      <value value="9.5830892"/>
      <value value="9.7766804"/>
      <value value="9.9741825"/>
      <value value="10.1756743"/>
      <value value="10.3812366"/>
      <value value="10.5909515"/>
      <value value="10.8049029"/>
      <value value="11.0231764"/>
      <value value="11.2458593"/>
      <value value="11.4730407"/>
      <value value="11.7048115"/>
      <value value="11.9412644"/>
      <value value="12.182494"/>
      <value value="12.4285967"/>
      <value value="12.679671"/>
      <value value="12.9358173"/>
      <value value="13.1971382"/>
      <value value="13.463738"/>
      <value value="13.7357236"/>
      <value value="14.0132036"/>
      <value value="14.2962891"/>
      <value value="14.5850933"/>
      <value value="14.8797317"/>
      <value value="15.1803222"/>
      <value value="15.4869851"/>
      <value value="15.7998429"/>
      <value value="16.1190209"/>
      <value value="16.4446468"/>
      <value value="16.7768507"/>
      <value value="17.1157655"/>
      <value value="17.4615269"/>
      <value value="17.8142732"/>
      <value value="18.1741454"/>
      <value value="18.5412875"/>
      <value value="18.9158463"/>
      <value value="19.2979718"/>
      <value value="19.6878166"/>
      <value value="20.0855369"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-b">
      <value value="&quot;single&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="competitive-exclusion-differentStart" repetitions="1" runMetricsEveryStep="false">
    <setup>setup
ifelse (random 2 = 0)[
   ask species1 set inviduals - 50
   ask species2 set inviduals + 50
][
   ask species1 set inviduals + 50
   ask species2 set inviduals - 50
]</setup>
    <go>go</go>
    <exitCondition>number-of-strategies &lt; 2</exitCondition>
    <metric>(count species1) &gt; 0</metric>
    <metric>(count species2) &gt; 0</metric>
    <metric>ticks</metric>
    <enumeratedValueSet variable="Demographic-Stochasticity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-rate">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="species">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-b">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="data-distance">
      <value value="5200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resolution">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="function">
      <value value="&quot;MSS&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disturbance">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction12">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation-probability">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="range">
      <value value="28.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cutoff">
      <value value="0.17"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="capacity">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-population">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-1">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-alpha">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualization">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction21">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-1">
      <value value="0.1353353"/>
      <value value="0.1380692"/>
      <value value="0.1408584"/>
      <value value="0.1437039"/>
      <value value="0.146607"/>
      <value value="0.1495686"/>
      <value value="0.1525901"/>
      <value value="0.1556726"/>
      <value value="0.1588174"/>
      <value value="0.1620258"/>
      <value value="0.1652989"/>
      <value value="0.1686381"/>
      <value value="0.1720449"/>
      <value value="0.1755204"/>
      <value value="0.1790661"/>
      <value value="0.1826835"/>
      <value value="0.186374"/>
      <value value="0.190139"/>
      <value value="0.19398"/>
      <value value="0.1978987"/>
      <value value="0.2018965"/>
      <value value="0.2059751"/>
      <value value="0.2101361"/>
      <value value="0.2143811"/>
      <value value="0.2187119"/>
      <value value="0.2231302"/>
      <value value="0.2276377"/>
      <value value="0.2322363"/>
      <value value="0.2369278"/>
      <value value="0.241714"/>
      <value value="0.246597"/>
      <value value="0.2515786"/>
      <value value="0.2566608"/>
      <value value="0.2618457"/>
      <value value="0.2671353"/>
      <value value="0.2725318"/>
      <value value="0.2780373"/>
      <value value="0.283654"/>
      <value value="0.2893842"/>
      <value value="0.2952302"/>
      <value value="0.3011942"/>
      <value value="0.3072787"/>
      <value value="0.3134862"/>
      <value value="0.319819"/>
      <value value="0.3262798"/>
      <value value="0.3328711"/>
      <value value="0.3395955"/>
      <value value="0.3464558"/>
      <value value="0.3534547"/>
      <value value="0.3605949"/>
      <value value="0.3678794"/>
      <value value="0.3753111"/>
      <value value="0.3828929"/>
      <value value="0.3906278"/>
      <value value="0.398519"/>
      <value value="0.4065697"/>
      <value value="0.4147829"/>
      <value value="0.4231621"/>
      <value value="0.4317105"/>
      <value value="0.4404317"/>
      <value value="0.449329"/>
      <value value="0.458406"/>
      <value value="0.4676664"/>
      <value value="0.4771139"/>
      <value value="0.4867523"/>
      <value value="0.4965853"/>
      <value value="0.506617"/>
      <value value="0.5168513"/>
      <value value="0.5272924"/>
      <value value="0.5379444"/>
      <value value="0.5488116"/>
      <value value="0.5598984"/>
      <value value="0.5712091"/>
      <value value="0.5827483"/>
      <value value="0.5945205"/>
      <value value="0.6065307"/>
      <value value="0.6187834"/>
      <value value="0.6312836"/>
      <value value="0.6440364"/>
      <value value="0.6570468"/>
      <value value="0.67032"/>
      <value value="0.6838614"/>
      <value value="0.6976763"/>
      <value value="0.7117703"/>
      <value value="0.726149"/>
      <value value="0.7408182"/>
      <value value="0.7557837"/>
      <value value="0.7710516"/>
      <value value="0.7866279"/>
      <value value="0.8025188"/>
      <value value="0.8187308"/>
      <value value="0.8352702"/>
      <value value="0.8521438"/>
      <value value="0.8693582"/>
      <value value="0.8869204"/>
      <value value="0.9048374"/>
      <value value="0.9231163"/>
      <value value="0.9417645"/>
      <value value="0.9607894"/>
      <value value="0.9801987"/>
      <value value="1"/>
      <value value="1.0202013"/>
      <value value="1.0408108"/>
      <value value="1.0618365"/>
      <value value="1.0832871"/>
      <value value="1.1051709"/>
      <value value="1.1274969"/>
      <value value="1.1502738"/>
      <value value="1.1735109"/>
      <value value="1.1972174"/>
      <value value="1.2214028"/>
      <value value="1.2460767"/>
      <value value="1.2712492"/>
      <value value="1.2969301"/>
      <value value="1.3231298"/>
      <value value="1.3498588"/>
      <value value="1.3771278"/>
      <value value="1.4049476"/>
      <value value="1.4333294"/>
      <value value="1.4622846"/>
      <value value="1.4918247"/>
      <value value="1.5219616"/>
      <value value="1.5527072"/>
      <value value="1.584074"/>
      <value value="1.6160744"/>
      <value value="1.6487213"/>
      <value value="1.6820276"/>
      <value value="1.7160069"/>
      <value value="1.7506725"/>
      <value value="1.7860384"/>
      <value value="1.8221188"/>
      <value value="1.858928"/>
      <value value="1.8964809"/>
      <value value="1.9347923"/>
      <value value="1.9738777"/>
      <value value="2.0137527"/>
      <value value="2.0544332"/>
      <value value="2.0959355"/>
      <value value="2.1382762"/>
      <value value="2.1814723"/>
      <value value="2.2255409"/>
      <value value="2.2704998"/>
      <value value="2.316367"/>
      <value value="2.3631607"/>
      <value value="2.4108997"/>
      <value value="2.4596031"/>
      <value value="2.5092904"/>
      <value value="2.5599814"/>
      <value value="2.6116965"/>
      <value value="2.6644562"/>
      <value value="2.7182818"/>
      <value value="2.7731948"/>
      <value value="2.829217"/>
      <value value="2.886371"/>
      <value value="2.9446796"/>
      <value value="3.004166"/>
      <value value="3.0648542"/>
      <value value="3.1267684"/>
      <value value="3.1899333"/>
      <value value="3.2543742"/>
      <value value="3.3201169"/>
      <value value="3.3871877"/>
      <value value="3.4556135"/>
      <value value="3.5254215"/>
      <value value="3.5966397"/>
      <value value="3.6692967"/>
      <value value="3.7434214"/>
      <value value="3.8190435"/>
      <value value="3.8961933"/>
      <value value="3.9749016"/>
      <value value="4.0552"/>
      <value value="4.1371204"/>
      <value value="4.2206958"/>
      <value value="4.3059595"/>
      <value value="4.3929457"/>
      <value value="4.4816891"/>
      <value value="4.5722252"/>
      <value value="4.6645903"/>
      <value value="4.7588212"/>
      <value value="4.8549558"/>
      <value value="4.9530324"/>
      <value value="5.0530903"/>
      <value value="5.1551695"/>
      <value value="5.2593108"/>
      <value value="5.365556"/>
      <value value="5.4739474"/>
      <value value="5.5845285"/>
      <value value="5.6973434"/>
      <value value="5.8124374"/>
      <value value="5.9298564"/>
      <value value="6.0496475"/>
      <value value="6.1718584"/>
      <value value="6.2965383"/>
      <value value="6.4237368"/>
      <value value="6.5535049"/>
      <value value="6.6858944"/>
      <value value="6.8209585"/>
      <value value="6.958751"/>
      <value value="7.0993271"/>
      <value value="7.242743"/>
      <value value="7.3890561"/>
      <value value="7.5383249"/>
      <value value="7.6906092"/>
      <value value="7.8459698"/>
      <value value="8.0044689"/>
      <value value="8.1661699"/>
      <value value="8.3311375"/>
      <value value="8.4994376"/>
      <value value="8.6711377"/>
      <value value="8.8463063"/>
      <value value="9.0250135"/>
      <value value="9.2073309"/>
      <value value="9.3933313"/>
      <value value="9.5830892"/>
      <value value="9.7766804"/>
      <value value="9.9741825"/>
      <value value="10.1756743"/>
      <value value="10.3812366"/>
      <value value="10.5909515"/>
      <value value="10.8049029"/>
      <value value="11.0231764"/>
      <value value="11.2458593"/>
      <value value="11.4730407"/>
      <value value="11.7048115"/>
      <value value="11.9412644"/>
      <value value="12.182494"/>
      <value value="12.4285967"/>
      <value value="12.679671"/>
      <value value="12.9358173"/>
      <value value="13.1971382"/>
      <value value="13.463738"/>
      <value value="13.7357236"/>
      <value value="14.0132036"/>
      <value value="14.2962891"/>
      <value value="14.5850933"/>
      <value value="14.8797317"/>
      <value value="15.1803222"/>
      <value value="15.4869851"/>
      <value value="15.7998429"/>
      <value value="16.1190209"/>
      <value value="16.4446468"/>
      <value value="16.7768507"/>
      <value value="17.1157655"/>
      <value value="17.4615269"/>
      <value value="17.8142732"/>
      <value value="18.1741454"/>
      <value value="18.5412875"/>
      <value value="18.9158463"/>
      <value value="19.2979718"/>
      <value value="19.6878166"/>
      <value value="20.0855369"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-2">
      <value value="0.1353353"/>
      <value value="0.1380692"/>
      <value value="0.1408584"/>
      <value value="0.1437039"/>
      <value value="0.146607"/>
      <value value="0.1495686"/>
      <value value="0.1525901"/>
      <value value="0.1556726"/>
      <value value="0.1588174"/>
      <value value="0.1620258"/>
      <value value="0.1652989"/>
      <value value="0.1686381"/>
      <value value="0.1720449"/>
      <value value="0.1755204"/>
      <value value="0.1790661"/>
      <value value="0.1826835"/>
      <value value="0.186374"/>
      <value value="0.190139"/>
      <value value="0.19398"/>
      <value value="0.1978987"/>
      <value value="0.2018965"/>
      <value value="0.2059751"/>
      <value value="0.2101361"/>
      <value value="0.2143811"/>
      <value value="0.2187119"/>
      <value value="0.2231302"/>
      <value value="0.2276377"/>
      <value value="0.2322363"/>
      <value value="0.2369278"/>
      <value value="0.241714"/>
      <value value="0.246597"/>
      <value value="0.2515786"/>
      <value value="0.2566608"/>
      <value value="0.2618457"/>
      <value value="0.2671353"/>
      <value value="0.2725318"/>
      <value value="0.2780373"/>
      <value value="0.283654"/>
      <value value="0.2893842"/>
      <value value="0.2952302"/>
      <value value="0.3011942"/>
      <value value="0.3072787"/>
      <value value="0.3134862"/>
      <value value="0.319819"/>
      <value value="0.3262798"/>
      <value value="0.3328711"/>
      <value value="0.3395955"/>
      <value value="0.3464558"/>
      <value value="0.3534547"/>
      <value value="0.3605949"/>
      <value value="0.3678794"/>
      <value value="0.3753111"/>
      <value value="0.3828929"/>
      <value value="0.3906278"/>
      <value value="0.398519"/>
      <value value="0.4065697"/>
      <value value="0.4147829"/>
      <value value="0.4231621"/>
      <value value="0.4317105"/>
      <value value="0.4404317"/>
      <value value="0.449329"/>
      <value value="0.458406"/>
      <value value="0.4676664"/>
      <value value="0.4771139"/>
      <value value="0.4867523"/>
      <value value="0.4965853"/>
      <value value="0.506617"/>
      <value value="0.5168513"/>
      <value value="0.5272924"/>
      <value value="0.5379444"/>
      <value value="0.5488116"/>
      <value value="0.5598984"/>
      <value value="0.5712091"/>
      <value value="0.5827483"/>
      <value value="0.5945205"/>
      <value value="0.6065307"/>
      <value value="0.6187834"/>
      <value value="0.6312836"/>
      <value value="0.6440364"/>
      <value value="0.6570468"/>
      <value value="0.67032"/>
      <value value="0.6838614"/>
      <value value="0.6976763"/>
      <value value="0.7117703"/>
      <value value="0.726149"/>
      <value value="0.7408182"/>
      <value value="0.7557837"/>
      <value value="0.7710516"/>
      <value value="0.7866279"/>
      <value value="0.8025188"/>
      <value value="0.8187308"/>
      <value value="0.8352702"/>
      <value value="0.8521438"/>
      <value value="0.8693582"/>
      <value value="0.8869204"/>
      <value value="0.9048374"/>
      <value value="0.9231163"/>
      <value value="0.9417645"/>
      <value value="0.9607894"/>
      <value value="0.9801987"/>
      <value value="1"/>
      <value value="1.0202013"/>
      <value value="1.0408108"/>
      <value value="1.0618365"/>
      <value value="1.0832871"/>
      <value value="1.1051709"/>
      <value value="1.1274969"/>
      <value value="1.1502738"/>
      <value value="1.1735109"/>
      <value value="1.1972174"/>
      <value value="1.2214028"/>
      <value value="1.2460767"/>
      <value value="1.2712492"/>
      <value value="1.2969301"/>
      <value value="1.3231298"/>
      <value value="1.3498588"/>
      <value value="1.3771278"/>
      <value value="1.4049476"/>
      <value value="1.4333294"/>
      <value value="1.4622846"/>
      <value value="1.4918247"/>
      <value value="1.5219616"/>
      <value value="1.5527072"/>
      <value value="1.584074"/>
      <value value="1.6160744"/>
      <value value="1.6487213"/>
      <value value="1.6820276"/>
      <value value="1.7160069"/>
      <value value="1.7506725"/>
      <value value="1.7860384"/>
      <value value="1.8221188"/>
      <value value="1.858928"/>
      <value value="1.8964809"/>
      <value value="1.9347923"/>
      <value value="1.9738777"/>
      <value value="2.0137527"/>
      <value value="2.0544332"/>
      <value value="2.0959355"/>
      <value value="2.1382762"/>
      <value value="2.1814723"/>
      <value value="2.2255409"/>
      <value value="2.2704998"/>
      <value value="2.316367"/>
      <value value="2.3631607"/>
      <value value="2.4108997"/>
      <value value="2.4596031"/>
      <value value="2.5092904"/>
      <value value="2.5599814"/>
      <value value="2.6116965"/>
      <value value="2.6644562"/>
      <value value="2.7182818"/>
      <value value="2.7731948"/>
      <value value="2.829217"/>
      <value value="2.886371"/>
      <value value="2.9446796"/>
      <value value="3.004166"/>
      <value value="3.0648542"/>
      <value value="3.1267684"/>
      <value value="3.1899333"/>
      <value value="3.2543742"/>
      <value value="3.3201169"/>
      <value value="3.3871877"/>
      <value value="3.4556135"/>
      <value value="3.5254215"/>
      <value value="3.5966397"/>
      <value value="3.6692967"/>
      <value value="3.7434214"/>
      <value value="3.8190435"/>
      <value value="3.8961933"/>
      <value value="3.9749016"/>
      <value value="4.0552"/>
      <value value="4.1371204"/>
      <value value="4.2206958"/>
      <value value="4.3059595"/>
      <value value="4.3929457"/>
      <value value="4.4816891"/>
      <value value="4.5722252"/>
      <value value="4.6645903"/>
      <value value="4.7588212"/>
      <value value="4.8549558"/>
      <value value="4.9530324"/>
      <value value="5.0530903"/>
      <value value="5.1551695"/>
      <value value="5.2593108"/>
      <value value="5.365556"/>
      <value value="5.4739474"/>
      <value value="5.5845285"/>
      <value value="5.6973434"/>
      <value value="5.8124374"/>
      <value value="5.9298564"/>
      <value value="6.0496475"/>
      <value value="6.1718584"/>
      <value value="6.2965383"/>
      <value value="6.4237368"/>
      <value value="6.5535049"/>
      <value value="6.6858944"/>
      <value value="6.8209585"/>
      <value value="6.958751"/>
      <value value="7.0993271"/>
      <value value="7.242743"/>
      <value value="7.3890561"/>
      <value value="7.5383249"/>
      <value value="7.6906092"/>
      <value value="7.8459698"/>
      <value value="8.0044689"/>
      <value value="8.1661699"/>
      <value value="8.3311375"/>
      <value value="8.4994376"/>
      <value value="8.6711377"/>
      <value value="8.8463063"/>
      <value value="9.0250135"/>
      <value value="9.2073309"/>
      <value value="9.3933313"/>
      <value value="9.5830892"/>
      <value value="9.7766804"/>
      <value value="9.9741825"/>
      <value value="10.1756743"/>
      <value value="10.3812366"/>
      <value value="10.5909515"/>
      <value value="10.8049029"/>
      <value value="11.0231764"/>
      <value value="11.2458593"/>
      <value value="11.4730407"/>
      <value value="11.7048115"/>
      <value value="11.9412644"/>
      <value value="12.182494"/>
      <value value="12.4285967"/>
      <value value="12.679671"/>
      <value value="12.9358173"/>
      <value value="13.1971382"/>
      <value value="13.463738"/>
      <value value="13.7357236"/>
      <value value="14.0132036"/>
      <value value="14.2962891"/>
      <value value="14.5850933"/>
      <value value="14.8797317"/>
      <value value="15.1803222"/>
      <value value="15.4869851"/>
      <value value="15.7998429"/>
      <value value="16.1190209"/>
      <value value="16.4446468"/>
      <value value="16.7768507"/>
      <value value="17.1157655"/>
      <value value="17.4615269"/>
      <value value="17.8142732"/>
      <value value="18.1741454"/>
      <value value="18.5412875"/>
      <value value="18.9158463"/>
      <value value="19.2979718"/>
      <value value="19.6878166"/>
      <value value="20.0855369"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-b">
      <value value="&quot;single&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="competitive-exclusion-MSS-Mod" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <exitCondition>number-of-strategies &lt; 2</exitCondition>
    <metric>(count species1) &gt; 0</metric>
    <metric>(count species2) &gt; 0</metric>
    <metric>ticks</metric>
    <enumeratedValueSet variable="Demographic-Stochasticity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-rate">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="species">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-b">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="data-distance">
      <value value="5200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resolution">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="function">
      <value value="&quot;MSS-Mod&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disturbance">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction12">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation-probability">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="range">
      <value value="28.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cutoff">
      <value value="0.17"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="capacity">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-population">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-1">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-alpha">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualization">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction21">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-1">
      <value value="0.1353353"/>
      <value value="0.1380692"/>
      <value value="0.1408584"/>
      <value value="0.1437039"/>
      <value value="0.146607"/>
      <value value="0.1495686"/>
      <value value="0.1525901"/>
      <value value="0.1556726"/>
      <value value="0.1588174"/>
      <value value="0.1620258"/>
      <value value="0.1652989"/>
      <value value="0.1686381"/>
      <value value="0.1720449"/>
      <value value="0.1755204"/>
      <value value="0.1790661"/>
      <value value="0.1826835"/>
      <value value="0.186374"/>
      <value value="0.190139"/>
      <value value="0.19398"/>
      <value value="0.1978987"/>
      <value value="0.2018965"/>
      <value value="0.2059751"/>
      <value value="0.2101361"/>
      <value value="0.2143811"/>
      <value value="0.2187119"/>
      <value value="0.2231302"/>
      <value value="0.2276377"/>
      <value value="0.2322363"/>
      <value value="0.2369278"/>
      <value value="0.241714"/>
      <value value="0.246597"/>
      <value value="0.2515786"/>
      <value value="0.2566608"/>
      <value value="0.2618457"/>
      <value value="0.2671353"/>
      <value value="0.2725318"/>
      <value value="0.2780373"/>
      <value value="0.283654"/>
      <value value="0.2893842"/>
      <value value="0.2952302"/>
      <value value="0.3011942"/>
      <value value="0.3072787"/>
      <value value="0.3134862"/>
      <value value="0.319819"/>
      <value value="0.3262798"/>
      <value value="0.3328711"/>
      <value value="0.3395955"/>
      <value value="0.3464558"/>
      <value value="0.3534547"/>
      <value value="0.3605949"/>
      <value value="0.3678794"/>
      <value value="0.3753111"/>
      <value value="0.3828929"/>
      <value value="0.3906278"/>
      <value value="0.398519"/>
      <value value="0.4065697"/>
      <value value="0.4147829"/>
      <value value="0.4231621"/>
      <value value="0.4317105"/>
      <value value="0.4404317"/>
      <value value="0.449329"/>
      <value value="0.458406"/>
      <value value="0.4676664"/>
      <value value="0.4771139"/>
      <value value="0.4867523"/>
      <value value="0.4965853"/>
      <value value="0.506617"/>
      <value value="0.5168513"/>
      <value value="0.5272924"/>
      <value value="0.5379444"/>
      <value value="0.5488116"/>
      <value value="0.5598984"/>
      <value value="0.5712091"/>
      <value value="0.5827483"/>
      <value value="0.5945205"/>
      <value value="0.6065307"/>
      <value value="0.6187834"/>
      <value value="0.6312836"/>
      <value value="0.6440364"/>
      <value value="0.6570468"/>
      <value value="0.67032"/>
      <value value="0.6838614"/>
      <value value="0.6976763"/>
      <value value="0.7117703"/>
      <value value="0.726149"/>
      <value value="0.7408182"/>
      <value value="0.7557837"/>
      <value value="0.7710516"/>
      <value value="0.7866279"/>
      <value value="0.8025188"/>
      <value value="0.8187308"/>
      <value value="0.8352702"/>
      <value value="0.8521438"/>
      <value value="0.8693582"/>
      <value value="0.8869204"/>
      <value value="0.9048374"/>
      <value value="0.9231163"/>
      <value value="0.9417645"/>
      <value value="0.9607894"/>
      <value value="0.9801987"/>
      <value value="1"/>
      <value value="1.0202013"/>
      <value value="1.0408108"/>
      <value value="1.0618365"/>
      <value value="1.0832871"/>
      <value value="1.1051709"/>
      <value value="1.1274969"/>
      <value value="1.1502738"/>
      <value value="1.1735109"/>
      <value value="1.1972174"/>
      <value value="1.2214028"/>
      <value value="1.2460767"/>
      <value value="1.2712492"/>
      <value value="1.2969301"/>
      <value value="1.3231298"/>
      <value value="1.3498588"/>
      <value value="1.3771278"/>
      <value value="1.4049476"/>
      <value value="1.4333294"/>
      <value value="1.4622846"/>
      <value value="1.4918247"/>
      <value value="1.5219616"/>
      <value value="1.5527072"/>
      <value value="1.584074"/>
      <value value="1.6160744"/>
      <value value="1.6487213"/>
      <value value="1.6820276"/>
      <value value="1.7160069"/>
      <value value="1.7506725"/>
      <value value="1.7860384"/>
      <value value="1.8221188"/>
      <value value="1.858928"/>
      <value value="1.8964809"/>
      <value value="1.9347923"/>
      <value value="1.9738777"/>
      <value value="2.0137527"/>
      <value value="2.0544332"/>
      <value value="2.0959355"/>
      <value value="2.1382762"/>
      <value value="2.1814723"/>
      <value value="2.2255409"/>
      <value value="2.2704998"/>
      <value value="2.316367"/>
      <value value="2.3631607"/>
      <value value="2.4108997"/>
      <value value="2.4596031"/>
      <value value="2.5092904"/>
      <value value="2.5599814"/>
      <value value="2.6116965"/>
      <value value="2.6644562"/>
      <value value="2.7182818"/>
      <value value="2.7731948"/>
      <value value="2.829217"/>
      <value value="2.886371"/>
      <value value="2.9446796"/>
      <value value="3.004166"/>
      <value value="3.0648542"/>
      <value value="3.1267684"/>
      <value value="3.1899333"/>
      <value value="3.2543742"/>
      <value value="3.3201169"/>
      <value value="3.3871877"/>
      <value value="3.4556135"/>
      <value value="3.5254215"/>
      <value value="3.5966397"/>
      <value value="3.6692967"/>
      <value value="3.7434214"/>
      <value value="3.8190435"/>
      <value value="3.8961933"/>
      <value value="3.9749016"/>
      <value value="4.0552"/>
      <value value="4.1371204"/>
      <value value="4.2206958"/>
      <value value="4.3059595"/>
      <value value="4.3929457"/>
      <value value="4.4816891"/>
      <value value="4.5722252"/>
      <value value="4.6645903"/>
      <value value="4.7588212"/>
      <value value="4.8549558"/>
      <value value="4.9530324"/>
      <value value="5.0530903"/>
      <value value="5.1551695"/>
      <value value="5.2593108"/>
      <value value="5.365556"/>
      <value value="5.4739474"/>
      <value value="5.5845285"/>
      <value value="5.6973434"/>
      <value value="5.8124374"/>
      <value value="5.9298564"/>
      <value value="6.0496475"/>
      <value value="6.1718584"/>
      <value value="6.2965383"/>
      <value value="6.4237368"/>
      <value value="6.5535049"/>
      <value value="6.6858944"/>
      <value value="6.8209585"/>
      <value value="6.958751"/>
      <value value="7.0993271"/>
      <value value="7.242743"/>
      <value value="7.3890561"/>
      <value value="7.5383249"/>
      <value value="7.6906092"/>
      <value value="7.8459698"/>
      <value value="8.0044689"/>
      <value value="8.1661699"/>
      <value value="8.3311375"/>
      <value value="8.4994376"/>
      <value value="8.6711377"/>
      <value value="8.8463063"/>
      <value value="9.0250135"/>
      <value value="9.2073309"/>
      <value value="9.3933313"/>
      <value value="9.5830892"/>
      <value value="9.7766804"/>
      <value value="9.9741825"/>
      <value value="10.1756743"/>
      <value value="10.3812366"/>
      <value value="10.5909515"/>
      <value value="10.8049029"/>
      <value value="11.0231764"/>
      <value value="11.2458593"/>
      <value value="11.4730407"/>
      <value value="11.7048115"/>
      <value value="11.9412644"/>
      <value value="12.182494"/>
      <value value="12.4285967"/>
      <value value="12.679671"/>
      <value value="12.9358173"/>
      <value value="13.1971382"/>
      <value value="13.463738"/>
      <value value="13.7357236"/>
      <value value="14.0132036"/>
      <value value="14.2962891"/>
      <value value="14.5850933"/>
      <value value="14.8797317"/>
      <value value="15.1803222"/>
      <value value="15.4869851"/>
      <value value="15.7998429"/>
      <value value="16.1190209"/>
      <value value="16.4446468"/>
      <value value="16.7768507"/>
      <value value="17.1157655"/>
      <value value="17.4615269"/>
      <value value="17.8142732"/>
      <value value="18.1741454"/>
      <value value="18.5412875"/>
      <value value="18.9158463"/>
      <value value="19.2979718"/>
      <value value="19.6878166"/>
      <value value="20.0855369"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-2">
      <value value="0.1353353"/>
      <value value="0.1380692"/>
      <value value="0.1408584"/>
      <value value="0.1437039"/>
      <value value="0.146607"/>
      <value value="0.1495686"/>
      <value value="0.1525901"/>
      <value value="0.1556726"/>
      <value value="0.1588174"/>
      <value value="0.1620258"/>
      <value value="0.1652989"/>
      <value value="0.1686381"/>
      <value value="0.1720449"/>
      <value value="0.1755204"/>
      <value value="0.1790661"/>
      <value value="0.1826835"/>
      <value value="0.186374"/>
      <value value="0.190139"/>
      <value value="0.19398"/>
      <value value="0.1978987"/>
      <value value="0.2018965"/>
      <value value="0.2059751"/>
      <value value="0.2101361"/>
      <value value="0.2143811"/>
      <value value="0.2187119"/>
      <value value="0.2231302"/>
      <value value="0.2276377"/>
      <value value="0.2322363"/>
      <value value="0.2369278"/>
      <value value="0.241714"/>
      <value value="0.246597"/>
      <value value="0.2515786"/>
      <value value="0.2566608"/>
      <value value="0.2618457"/>
      <value value="0.2671353"/>
      <value value="0.2725318"/>
      <value value="0.2780373"/>
      <value value="0.283654"/>
      <value value="0.2893842"/>
      <value value="0.2952302"/>
      <value value="0.3011942"/>
      <value value="0.3072787"/>
      <value value="0.3134862"/>
      <value value="0.319819"/>
      <value value="0.3262798"/>
      <value value="0.3328711"/>
      <value value="0.3395955"/>
      <value value="0.3464558"/>
      <value value="0.3534547"/>
      <value value="0.3605949"/>
      <value value="0.3678794"/>
      <value value="0.3753111"/>
      <value value="0.3828929"/>
      <value value="0.3906278"/>
      <value value="0.398519"/>
      <value value="0.4065697"/>
      <value value="0.4147829"/>
      <value value="0.4231621"/>
      <value value="0.4317105"/>
      <value value="0.4404317"/>
      <value value="0.449329"/>
      <value value="0.458406"/>
      <value value="0.4676664"/>
      <value value="0.4771139"/>
      <value value="0.4867523"/>
      <value value="0.4965853"/>
      <value value="0.506617"/>
      <value value="0.5168513"/>
      <value value="0.5272924"/>
      <value value="0.5379444"/>
      <value value="0.5488116"/>
      <value value="0.5598984"/>
      <value value="0.5712091"/>
      <value value="0.5827483"/>
      <value value="0.5945205"/>
      <value value="0.6065307"/>
      <value value="0.6187834"/>
      <value value="0.6312836"/>
      <value value="0.6440364"/>
      <value value="0.6570468"/>
      <value value="0.67032"/>
      <value value="0.6838614"/>
      <value value="0.6976763"/>
      <value value="0.7117703"/>
      <value value="0.726149"/>
      <value value="0.7408182"/>
      <value value="0.7557837"/>
      <value value="0.7710516"/>
      <value value="0.7866279"/>
      <value value="0.8025188"/>
      <value value="0.8187308"/>
      <value value="0.8352702"/>
      <value value="0.8521438"/>
      <value value="0.8693582"/>
      <value value="0.8869204"/>
      <value value="0.9048374"/>
      <value value="0.9231163"/>
      <value value="0.9417645"/>
      <value value="0.9607894"/>
      <value value="0.9801987"/>
      <value value="1"/>
      <value value="1.0202013"/>
      <value value="1.0408108"/>
      <value value="1.0618365"/>
      <value value="1.0832871"/>
      <value value="1.1051709"/>
      <value value="1.1274969"/>
      <value value="1.1502738"/>
      <value value="1.1735109"/>
      <value value="1.1972174"/>
      <value value="1.2214028"/>
      <value value="1.2460767"/>
      <value value="1.2712492"/>
      <value value="1.2969301"/>
      <value value="1.3231298"/>
      <value value="1.3498588"/>
      <value value="1.3771278"/>
      <value value="1.4049476"/>
      <value value="1.4333294"/>
      <value value="1.4622846"/>
      <value value="1.4918247"/>
      <value value="1.5219616"/>
      <value value="1.5527072"/>
      <value value="1.584074"/>
      <value value="1.6160744"/>
      <value value="1.6487213"/>
      <value value="1.6820276"/>
      <value value="1.7160069"/>
      <value value="1.7506725"/>
      <value value="1.7860384"/>
      <value value="1.8221188"/>
      <value value="1.858928"/>
      <value value="1.8964809"/>
      <value value="1.9347923"/>
      <value value="1.9738777"/>
      <value value="2.0137527"/>
      <value value="2.0544332"/>
      <value value="2.0959355"/>
      <value value="2.1382762"/>
      <value value="2.1814723"/>
      <value value="2.2255409"/>
      <value value="2.2704998"/>
      <value value="2.316367"/>
      <value value="2.3631607"/>
      <value value="2.4108997"/>
      <value value="2.4596031"/>
      <value value="2.5092904"/>
      <value value="2.5599814"/>
      <value value="2.6116965"/>
      <value value="2.6644562"/>
      <value value="2.7182818"/>
      <value value="2.7731948"/>
      <value value="2.829217"/>
      <value value="2.886371"/>
      <value value="2.9446796"/>
      <value value="3.004166"/>
      <value value="3.0648542"/>
      <value value="3.1267684"/>
      <value value="3.1899333"/>
      <value value="3.2543742"/>
      <value value="3.3201169"/>
      <value value="3.3871877"/>
      <value value="3.4556135"/>
      <value value="3.5254215"/>
      <value value="3.5966397"/>
      <value value="3.6692967"/>
      <value value="3.7434214"/>
      <value value="3.8190435"/>
      <value value="3.8961933"/>
      <value value="3.9749016"/>
      <value value="4.0552"/>
      <value value="4.1371204"/>
      <value value="4.2206958"/>
      <value value="4.3059595"/>
      <value value="4.3929457"/>
      <value value="4.4816891"/>
      <value value="4.5722252"/>
      <value value="4.6645903"/>
      <value value="4.7588212"/>
      <value value="4.8549558"/>
      <value value="4.9530324"/>
      <value value="5.0530903"/>
      <value value="5.1551695"/>
      <value value="5.2593108"/>
      <value value="5.365556"/>
      <value value="5.4739474"/>
      <value value="5.5845285"/>
      <value value="5.6973434"/>
      <value value="5.8124374"/>
      <value value="5.9298564"/>
      <value value="6.0496475"/>
      <value value="6.1718584"/>
      <value value="6.2965383"/>
      <value value="6.4237368"/>
      <value value="6.5535049"/>
      <value value="6.6858944"/>
      <value value="6.8209585"/>
      <value value="6.958751"/>
      <value value="7.0993271"/>
      <value value="7.242743"/>
      <value value="7.3890561"/>
      <value value="7.5383249"/>
      <value value="7.6906092"/>
      <value value="7.8459698"/>
      <value value="8.0044689"/>
      <value value="8.1661699"/>
      <value value="8.3311375"/>
      <value value="8.4994376"/>
      <value value="8.6711377"/>
      <value value="8.8463063"/>
      <value value="9.0250135"/>
      <value value="9.2073309"/>
      <value value="9.3933313"/>
      <value value="9.5830892"/>
      <value value="9.7766804"/>
      <value value="9.9741825"/>
      <value value="10.1756743"/>
      <value value="10.3812366"/>
      <value value="10.5909515"/>
      <value value="10.8049029"/>
      <value value="11.0231764"/>
      <value value="11.2458593"/>
      <value value="11.4730407"/>
      <value value="11.7048115"/>
      <value value="11.9412644"/>
      <value value="12.182494"/>
      <value value="12.4285967"/>
      <value value="12.679671"/>
      <value value="12.9358173"/>
      <value value="13.1971382"/>
      <value value="13.463738"/>
      <value value="13.7357236"/>
      <value value="14.0132036"/>
      <value value="14.2962891"/>
      <value value="14.5850933"/>
      <value value="14.8797317"/>
      <value value="15.1803222"/>
      <value value="15.4869851"/>
      <value value="15.7998429"/>
      <value value="16.1190209"/>
      <value value="16.4446468"/>
      <value value="16.7768507"/>
      <value value="17.1157655"/>
      <value value="17.4615269"/>
      <value value="17.8142732"/>
      <value value="18.1741454"/>
      <value value="18.5412875"/>
      <value value="18.9158463"/>
      <value value="19.2979718"/>
      <value value="19.6878166"/>
      <value value="20.0855369"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-b">
      <value value="&quot;single&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Brancing-MSS-Mod" repetitions="1" runMetricsEveryStep="false">
    <setup>branching-setup</setup>
    <go>branching-go</go>
    <final>branching-end</final>
    <timeLimit steps="150000"/>
    <enumeratedValueSet variable="visualization">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-2">
      <value value="6.0496475"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="capacity">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Demographic-Stochasticity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resolution">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="range">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-rate-1">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-b">
      <value value="&quot;single&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cutoff">
      <value value="0.17"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-rate">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="data-distance">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-b">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="function">
      <value value="&quot;MSS-Mod&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutate-alpha">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction12">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-population">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation-probability">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disturbance">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="species">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="interaction21">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-strategy-1">
      <value value="4"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

@#$#@#$#@
1
@#$#@#$#@
