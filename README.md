# insilico-perfusion

## Installation
```
python setup.py install
```

## Examples

### Create and optimize a vessel network
```
python examples/createVascularization.py
```

![alt text](pressure.gif "test")


### Run perfusion simulation on vessel network
```
python examples/calcPerfusion.py
```

![](concP.gif "intra-vascular marker concentration"){:height="200%" width="200%"}
![](concI.gif "interstitial marker concentration"){:height="200%" width="200%"}
![](conc.gif "marker concentration"){:height="200%" width="200%"}
