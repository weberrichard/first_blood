# CCC test modelling
CCC test refers to Common carotid Compression test, investigating the cerebral autoregulation. During this test, one of the carotid arteries is comressed, and a few physiological measurements are monitored in the meantime.

### Model
The "Abel_A62_63_control" model is used for the simulation. A few remarks on the model:
- Additional resistances are added to the circle of willis, since the original model did not account for the pressure drop caused by the sharp turns at edges A62, A63 and A77.
- Peripheral baroreflex (arterolar resistance changing induced by the baroreceptors as prescribed) and the myogenic response are the only active peripheral autoregulation processes in that model.
- The other form of baroreflex, the changing heart raet as a function of the systolic blood pressure is inactive, since it only has a small effect.
- This model is not closed, the peripheral models end in a ground.
- The transport of oxygen and carbon dioxide is not included.

![Alt text](AbelA63_62.png?raw=true "Title")

# How to use it for validation
- Obtain some measuremet data.
- Tune the model to match the homeostasis walues like pressure and velocity.
- Tune the autoregulaton parameters to match the CCC test itself. 

# Current simulation result
![Alt text](Abel_A62_63_control_velocity.png?raw=true "Title")
