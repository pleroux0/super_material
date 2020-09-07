# API

Only the interfaces and classes listed here are considered stable. Care will be taken to keep the public API constant, but other interfaces and classes can and will change between releases.

## Interfaces

The library consists of two public interfaces and a few implementations of those interfaces. The first interface,`GapEnergyInterface`, defines an interface with which to calculate the superconductor gap energy. The second interface, `SuperconductorConductivityInterface`defines an interface with which to calculate the superconductor conductivity.

### Gap Energy

```python
class SuperconductorConductivityInterface(ABC):
    @abstractmethod
    def evaluate(self, temperature: float, frequency: float) -> complex:
        """ Calculates the superconductor complex conductivity """
```

#### BCS

### Superconductor Conductivity

```python
class GapEnergyInterface(ABC):
    """ Interface to calculate the gap energy """

    @abstractmethod
    def evaluate(self, temperature: float) -> float:
        """ Evaluate the gap energy at the specific temperature """

    @abstractmethod
    def critical_temperature(self) -> float:
        """ Get the critical temperature of the gap energy """

    @abstractmethod
    def gap_energy_0(self) -> float:
        """ Get the gap energy at T = 0 K """

    def critical_frequency(self, temperature: float) -> float:
        return self.evaluate(temperature) / (h_bar * pi)
```

#### Zimmermann

#### Mattis-Bardeen

#### 

