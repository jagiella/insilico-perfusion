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

![](concP.gif "intra-vascular marker concentration")
![](concI.gif "interstitial marker concentration")
![](conc.gif "marker concentration")
