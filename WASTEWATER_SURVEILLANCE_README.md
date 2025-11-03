# Wastewater Surveillance for Epidemic Simulation

## Overview

This document describes the wastewater-based epidemiological surveillance system added to the epidemic simulator. Wastewater surveillance provides an early warning system for disease outbreaks by detecting viral RNA in sewage before clinical symptoms appear.

## Three Key Features

### 1. RNA Viral Load in Wastewater
- **What**: Measures the concentration of viral RNA (copies/L) in wastewater from sewage catchment areas
- **How**: Aggregates viral shedding from all infected individuals in a sewershed
- **Output**: Daily viral concentration, total viral load, normalized concentration per 100k population

### 2. Detected Cases by Sewage Area
- **What**: Clinical surveillance data (detected cases, hospitalizations, deaths) aggregated by sewershed
- **How**: Maps individuals to their sewershed through home location and tracks clinical outcomes
- **Output**: Daily and cumulative counts of detected cases, hospitalizations, and deaths per sewershed

### 3. Detection Ratio Estimation
- **What**: Estimates the fraction of true infections detected by clinical surveillance
- **How**: Compares wastewater-estimated infections with clinically detected cases
- **Output**: `detection_ratio = detected_cases / estimated_true_infections`
- **Interpretation**:
  - Ratio of 0.8 = 80% of infections are detected (good surveillance)
  - Ratio of 0.3 = 70% of infections are missed (poor surveillance)
  - Underreporting factor = 1/ratio (e.g., 3.0 means 3 true infections per detected case)

---

## Files Added

### C++ Simulator
- **wastewater.h**: Header file with data structures and function declarations
- **wastewater.cc**: Implementation of viral shedding, RNA load calculation, and detection estimation
- **wastewater_config.json**: Example configuration file
- **Modified files**:
  - `models.h`: Added ENABLE_WASTEWATER_SURVEILLANCE flag
  - `simulator.cc`: Integrated wastewater surveillance into main loop
  - `makefile_base.mk`: Added wastewater.o to compilation
  - `outputs.h`: Added wastewater.h include

### JavaScript Simulator
- (See JavaScript Integration section below)

---

## Usage Instructions

### C++ Simulator

#### 1. Enable Wastewater Surveillance

Add to your input parameter file (e.g., `input_base/input_parameters.json`):

```json
{
  "ENABLE_WASTEWATER_SURVEILLANCE": true,
  "WASTEWATER_SAMPLING_FREQUENCY": 3.0,
  "WASTEWATER_NUM_SEWERSHEDS": 20
}
```

#### 2. Compile

```bash
cd cpp-simulator
make clean
make
```

#### 3. Run Simulation

```bash
./drive_simulator \
  --input_base ../simulator/input_files/ \
  --output_path ./output/ \
  # ... other parameters
```

#### 4. Output Files

The simulator will generate:
- **`output/wastewater_surveillance.csv`**: Complete wastewater surveillance data

CSV columns:
```
Day, Sewershed_ID, Population,
Viral_Concentration_copies_per_L, Normalized_Conc_per_100k, Total_Viral_Load,
Detected_Cases, Hospitalised, Deaths,
Cumulative_Detected, Cumulative_Hospitalised, Cumulative_Deaths,
Estimated_True_Infections, True_Active_Infections,
Detection_Ratio, Underreporting_Factor, Below_Detection_Limit
```

---

## Configuration Parameters

### Global Parameters (models.h)

```cpp
ENABLE_WASTEWATER_SURVEILLANCE = false   // Enable/disable feature
WASTEWATER_SAMPLING_FREQUENCY = 3.0      // Samples per week (typical: 3-7)
WASTEWATER_NUM_SEWERSHEDS = 20           // Number of sewershed catchments
```

### Viral Shedding Parameters (wastewater.h)

```cpp
peak_shedding_day = 6.0                  // Day when viral shedding peaks (5-7 typical)
peak_shedding_rate = 1e9                 // Peak RNA copies/person/day
shedding_curve_sigma = 3.0               // Width of shedding curve
shedding_duration_days = 30.0            // Total days of viral shedding
asymptomatic_shedding_fraction = 0.6     // Asymptomatic shed 60% of symptomatic
detection_limit_copies_per_L = 1000.0    // Lab detection threshold
```

### Sewershed Configuration (wastewater.h)

```cpp
wards_per_sewershed = 10                 // Wards grouped into each sewershed
sampling_frequency_per_week = 3.0        // How often samples are collected
collection_efficiency = 0.85             // Fraction of sewage collected (0.7-0.95)
dilution_factor = 1.0                    // Seasonal dilution (>1.0 in rainy season)
liters_per_person_per_day = 175.0        // Wastewater generation rate
```

### Detection Ratio Estimation

```cpp
viral_load_to_infection_factor = 1e-6    // Calibration: converts copies/L to infections/100k
minimum_infections_for_detection = 5     // Minimum infections for reliable WW signal
wastewater_lead_time_days = 5            // Days WW is ahead of clinical detection
```

---

## Scientific Background

### Viral Shedding Timeline

```
Infection Event (Day 0)
    ↓
Day 1-2: Shedding begins (EXPOSED phase)
    ↓
Day 3-5: Shedding increases (PRE-SYMPTOMATIC phase)
    ↓
Day 5-7: PEAK SHEDDING (early SYMPTOMATIC phase)
    ↓
Day 8-14: Shedding decreases (SYMPTOMATIC → RECOVERED)
    ↓
Day 15-30: Low-level shedding continues (RECOVERED)
```

**Key Point**: Wastewater detects infections 4-7 days before clinical symptoms appear, providing early warning.

### Wastewater vs. Clinical Surveillance

| Aspect | Clinical Surveillance | Wastewater Surveillance |
|--------|---------------------|------------------------|
| **Detection** | Only symptomatic + tested | All infected (symptomatic + asymptomatic) |
| **Timing** | After symptoms appear | 4-7 days earlier |
| **Coverage** | Depends on healthcare access | Entire population served by sewershed |
| **Cost** | High (individual testing) | Low (population-level sampling) |
| **Bias** | Testing bias, reporting delays | None (passive collection) |
| **Resolution** | Individual-level | Sewershed-level (population aggregate) |

### Detection Ratio Formula

```
detection_ratio = cumulative_detected_cases / estimated_true_infections

estimated_true_infections = f(viral_concentration, population)

where f() is the calibration function that converts RNA copies/L to infection estimates
```

---

## Interpreting Results

### Example Output

```csv
Day,Sewershed_ID,Population,Viral_Concentration_copies_per_L,Detected_Cases,
    True_Active_Infections,Detection_Ratio,Underreporting_Factor

30,5,50000,5000.5,12,40,0.30,3.33
60,5,50000,25000.8,85,250,0.34,2.94
90,5,50000,8000.2,45,120,0.375,2.67
```

**Interpretation**:
- **Day 30**: Only 30% of infections are detected (70% missed), underreporting factor of 3.3×
- **Day 60**: Outbreak peak - viral load 5× higher, but detection ratio still low
- **Day 90**: Post-peak - detection improving slightly (37.5% detected)

### Detection Ratio Benchmarks

| Detection Ratio | Interpretation | Surveillance Quality |
|----------------|----------------|---------------------|
| > 0.8 | Excellent surveillance - 80%+ infections detected | High |
| 0.5 - 0.8 | Good surveillance - half or more infections detected | Medium-High |
| 0.3 - 0.5 | Moderate surveillance - many infections missed | Medium |
| < 0.3 | Poor surveillance - majority of infections missed | Low |

---

## Calibration Guide

The most important parameter to calibrate is `viral_load_to_infection_factor`, which converts wastewater viral concentration to estimated infection counts.

### Method 1: Historical Data Calibration

If you have historical wastewater measurements and case counts:

```python
# Example calibration
historical_ww_concentration = [1000, 5000, 10000, 20000]  # copies/L
historical_true_cases_per_100k = [50, 250, 500, 1000]

# Linear regression or curve fitting
import numpy as np
factor = np.polyfit(historical_ww_concentration, historical_true_cases_per_100k, 1)[0]

# Typical range: 1e-7 to 1e-5
```

### Method 2: Literature-Based Estimation

From COVID-19 studies:
- Typical relationship: 1e4 - 1e5 copies/L ≈ 100-1000 infections per 100k population
- Suggested starting value: `viral_load_to_infection_factor = 1e-6`

### Method 3: Model Calibration

Run simulations with known infection rates and adjust `viral_load_to_infection_factor` until:

```
estimated_true_infections ≈ actual_infections (from ground truth)
```

---

## Advanced Features

### 1. Sewershed Configuration

Default: Groups wards automatically (10 wards per sewershed)

Custom configuration (future enhancement):
```json
{
  "sewersheds": [
    {
      "id": 0,
      "name": "North_Treatment_Plant",
      "ward_ids": [0, 1, 2, 3, 4],
      "flow_rate_L_per_day": 5000000,
      "sampling_frequency": 7
    }
  ]
}
```

### 2. Sampling Strategies

Configure sampling frequency to test different strategies:

```cpp
WASTEWATER_SAMPLING_FREQUENCY = 7.0   // Daily sampling (expensive)
WASTEWATER_SAMPLING_FREQUENCY = 3.0   // 3×/week (standard)
WASTEWATER_SAMPLING_FREQUENCY = 1.0   // Weekly (minimal cost)
```

**Research Question**: What is the minimum sampling frequency to detect outbreaks reliably?

### 3. Intervention Triggering

Future enhancement: Trigger interventions based on wastewater thresholds

```cpp
if (sewershed.viral_concentration > THRESHOLD) {
    trigger_ward_containment(sewershed.ward_ids);
}
```

---

## Visualization

### Python Plotting Example

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load wastewater data
ww = pd.read_csv('output/wastewater_surveillance.csv')

# Filter one sewershed
ss_data = ww[ww['Sewershed_ID'] == 5]

# Plot 1: Viral load over time
plt.figure(figsize=(12, 4))
plt.plot(ss_data['Day'], ss_data['Viral_Concentration_copies_per_L'])
plt.xlabel('Day')
plt.ylabel('Viral RNA (copies/L)')
plt.title('Wastewater Viral Load - Sewershed 5')
plt.yscale('log')
plt.grid(True)
plt.savefig('viral_load.png')

# Plot 2: Detection ratio over time
plt.figure(figsize=(12, 4))
plt.plot(ss_data['Day'], ss_data['Detection_Ratio'])
plt.xlabel('Day')
plt.ylabel('Detection Ratio')
plt.title('Clinical Surveillance Detection Ratio')
plt.ylim(0, 1)
plt.axhline(y=0.8, color='g', linestyle='--', label='Good surveillance (80%)')
plt.axhline(y=0.5, color='orange', linestyle='--', label='Moderate (50%)')
plt.axhline(y=0.3, color='r', linestyle='--', label='Poor (30%)')
plt.legend()
plt.grid(True)
plt.savefig('detection_ratio.png')

# Plot 3: Comparison - True vs Detected
plt.figure(figsize=(12, 4))
plt.plot(ss_data['Day'], ss_data['True_Active_Infections'], label='True Infections (Ground Truth)')
plt.plot(ss_data['Day'], ss_data['Cumulative_Detected'], label='Detected Cases (Clinical)')
plt.plot(ss_data['Day'], ss_data['Estimated_True_Infections'], label='Estimated from Wastewater', linestyle='--')
plt.xlabel('Day')
plt.ylabel('Number of Infections')
plt.title('Comparison: True vs Detected vs Wastewater-Estimated Infections')
plt.legend()
plt.grid(True)
plt.savefig('infection_comparison.png')
```

---

## JavaScript Simulator Integration

### Quick Integration (sim.js)

Add to global variables:

```javascript
// Wastewater surveillance parameters
var ENABLE_WASTEWATER = true;
var WASTEWATER_SAMPLING_FREQ = 3; // samples per week
var WASTEWATER_VIRAL_LOAD = [];
var WASTEWATER_DETECTED_CASES = [];
var WASTEWATER_DETECTION_RATIO = [];

// Viral shedding parameters
var PEAK_SHEDDING_DAY = 6;
var PEAK_SHEDDING_RATE = 1e9;
var SHEDDING_DURATION = 30;
```

Add to `runStep()` function:

```javascript
// Call once per day
if (steps % NUM_STEPS_PER_DAY === 0 && ENABLE_WASTEWATER) {
    updateWastewaterSurveillance(steps / NUM_STEPS_PER_DAY);
}

function updateWastewaterSurveillance(day) {
    let totalViralShedding = 0;
    let detectedCases = 0;
    let trueInfections = 0;

    // Calculate viral shedding from all infected individuals
    for (let i = 0; i < Individual.length; i++) {
        let state = Individual[i].infectionState;

        // Count true infections
        if (state > 0 && state < 7) { // Infected states
            trueInfections++;

            // Calculate days since infection
            let daysSinceInfection = day - Individual[i].timeOfInfection;

            // Calculate viral shedding
            if (daysSinceInfection >= 1 && daysSinceInfection <= SHEDDING_DURATION) {
                let shedding = PEAK_SHEDDING_RATE *
                              Math.exp(-0.5 * Math.pow((daysSinceInfection - PEAK_SHEDDING_DAY) / 3, 2));

                // Reduce for asymptomatic
                if (!Individual[i].symptomatic) {
                    shedding *= 0.6;
                }

                totalViralShedding += shedding;
            }
        }

        // Count detected cases (symptomatic with detection probability)
        if (Individual[i].symptomatic && Math.random() < CASE_INFECTION_RATIO) {
            detectedCases++;
        }
    }

    // Calculate concentration (simplified: assume 175L/person/day flow rate)
    let flowRate = Individual.length * 175; // L/day
    let concentration = totalViralShedding / flowRate;

    // Estimate infections from wastewater
    let estimatedInfections = concentration * 1e-6 * (Individual.length / 100000);

    // Calculate detection ratio
    let detectionRatio = detectedCases / (estimatedInfections || 1);

    // Store for plotting
    WASTEWATER_VIRAL_LOAD.push(concentration);
    WASTEWATER_DETECTED_CASES.push(detectedCases);
    WASTEWATER_DETECTION_RATIO.push(detectionRatio);
}
```

---

## Troubleshooting

### Issue: No wastewater output generated

**Solution**: Check that `ENABLE_WASTEWATER_SURVEILLANCE = true` in global parameters

### Issue: All detection ratios are 0

**Solution**: Increase `viral_load_to_infection_factor` (try 1e-5 or 1e-4)

### Issue: Detection ratios > 1.0

**Solution**: Decrease `viral_load_to_infection_factor` (try 1e-7 or 1e-8), or check that case_infection_ratio is reasonable

### Issue: Compilation error: "wastewater.h: No such file"

**Solution**: Ensure wastewater.h and wastewater.cc are in cpp-simulator directory, and makefile_base.mk includes wastewater.o

---

## Research Applications

### 1. Early Warning System Evaluation
**Question**: How much earlier does wastewater detect outbreaks compared to clinical surveillance?

**Method**: Compare wastewater signal peaks with clinical case peaks

### 2. Surveillance System Design
**Question**: What sampling frequency is cost-effective for outbreak detection?

**Method**: Run simulations with different `WASTEWATER_SAMPLING_FREQUENCY` values (1, 3, 7 per week)

### 3. Intervention Optimization
**Question**: Can wastewater-triggered interventions reduce deaths compared to case-triggered interventions?

**Method**: Implement wastewater threshold-based lockdowns

### 4. True Burden Estimation
**Question**: What fraction of infections are missed by clinical surveillance?

**Method**: Analyze `detection_ratio` and `underreporting_factor` across different intervention scenarios

---

## References

1. Wastewater-Based Epidemiology for COVID-19:
   - Peccia et al. (2020). "SARS-CoV-2 RNA concentrations in primary municipal sewage sludge as a leading indicator of COVID-19 outbreak dynamics"

2. Viral Shedding Dynamics:
   - Wölfel et al. (2020). "Virological assessment of hospitalized patients with COVID-2019"

3. Detection Ratio Estimation:
   - Wu et al. (2020). "Substantial underestimation of SARS-CoV-2 infection in the United States"

---

## Contact

For questions or issues with wastewater surveillance implementation:
- Create an issue on the GitHub repository
- Email the development team

---

## Version History

- **v1.0** (2025-11-03): Initial implementation with three core features
  - RNA viral load calculation
  - Detected cases by sewershed
  - Detection ratio estimation
