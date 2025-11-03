//Copyright [2020] [Indian Institute of Science, Bangalore & Tata Institute of Fundamental Research, Mumbai]
//SPDX-License-Identifier: Apache-2.0
#ifndef WASTEWATER_H_
#define WASTEWATER_H_

#include "models.h"
#include <vector>
#include <string>

// ============================================================================
// WASTEWATER SURVEILLANCE SYSTEM
// ============================================================================
// This module implements three key features:
// 1. RNA viral load in wastewater (copies/L)
// 2. Detected clinical cases aggregated by sewage area
// 3. Detection ratio estimation (detected cases / true cases from WW)
// ============================================================================

// ----------------------------------------------------------------------------
// VIRAL SHEDDING PARAMETERS
// ----------------------------------------------------------------------------
// Parameters controlling how infected individuals shed viral RNA
struct viral_shedding_params {
    // Peak shedding occurs around day 5-7 after infection
    double peak_shedding_day = 6.0;

    // Peak shedding rate in RNA copies per person per day
    double peak_shedding_rate = 1.0e9; // 1 billion copies/person/day at peak

    // Standard deviation for shedding curve (controls width)
    double shedding_curve_sigma = 3.0;

    // Total duration of viral shedding (days)
    double shedding_duration_days = 30.0;

    // Fraction of symptomatic shedding for asymptomatic cases (0.5-0.8)
    double asymptomatic_shedding_fraction = 0.6;

    // Minimum days before shedding starts (typically during exposed state)
    double shedding_start_day = 1.0;

    // Detection limit (minimum detectable concentration in copies/L)
    double detection_limit_copies_per_L = 1000.0;
};

// ----------------------------------------------------------------------------
// SEWERSHED STRUCTURE
// ----------------------------------------------------------------------------
// Represents a wastewater catchment area serving a population
struct sewershed {
    int sewershed_id;
    std::string sewershed_name;

    // Geographic coverage
    std::vector<int> ward_ids;          // Wards/communities draining to this sewershed
    std::vector<int> home_indices;      // Homes connected to this sewershed

    // Population and infrastructure
    count_type population = 0;           // Total population served
    double flow_rate_L_per_day = 0;      // Wastewater flow rate (liters/day)
    double collection_efficiency = 0.85; // Fraction of sewage collected (0.7-0.95)
    double dilution_factor = 1.0;        // Seasonal dilution (rainfall, etc.)

    // Sampling configuration
    double sampling_frequency_per_week = 3.0; // Number of samples per week
    count_type last_sample_day = 0;           // Last day when sample was collected
    bool sampling_active = true;              // Enable/disable sampling for this sewershed

    // ========================================================================
    // FEATURE 1: RNA VIRAL LOAD IN WASTEWATER
    // ========================================================================
    // Daily aggregated viral shedding from all infected individuals
    double daily_total_viral_shedding = 0;    // Total RNA copies shed by all infected
    double daily_wastewater_concentration = 0; // Concentration in copies/L
    double normalized_concentration = 0;       // Per 100,000 population

    // ========================================================================
    // FEATURE 2: DETECTED CASES IN SEWAGE AREA
    // ========================================================================
    // Clinical surveillance data aggregated for this sewershed
    count_type daily_detected_cases = 0;       // New detected cases today (clinical)
    count_type cumulative_detected_cases = 0;  // Total detected cases over time
    count_type daily_hospitalised = 0;         // Daily hospitalizations
    count_type cumulative_hospitalised = 0;    // Total hospitalizations
    count_type daily_deaths = 0;               // Daily deaths
    count_type cumulative_deaths = 0;          // Total deaths

    // ========================================================================
    // FEATURE 3: DETECTION RATIO ESTIMATION
    // ========================================================================
    // Estimate true infections from wastewater and compare with detected cases
    count_type estimated_true_infections = 0;  // Estimated from wastewater signal
    double detection_ratio = 0;                 // detected_cases / estimated_true_cases
    double underreporting_factor = 0;           // 1 / detection_ratio (multiplier)

    // Reset daily counters (called at start of each day)
    void reset_daily_counters() {
        daily_total_viral_shedding = 0;
        daily_wastewater_concentration = 0;
        daily_detected_cases = 0;
        daily_hospitalised = 0;
        daily_deaths = 0;
    }
};

// ----------------------------------------------------------------------------
// WASTEWATER SAMPLE DATA
// ----------------------------------------------------------------------------
// Individual sample record for output
struct wastewater_sample {
    int sewershed_id;
    int day;

    // Feature 1: RNA load
    double viral_concentration_copies_per_L;
    double normalized_concentration_per_100k;
    double total_viral_load_copies;

    // Feature 2: Detected cases
    count_type detected_cases_in_area;
    count_type hospitalised_in_area;
    count_type deaths_in_area;
    count_type cumulative_detected;
    count_type cumulative_hospitalised;
    count_type cumulative_deaths;

    // Feature 3: Detection ratio
    count_type estimated_true_infections;
    count_type true_active_infections;      // Ground truth from simulation
    double detection_ratio;
    double underreporting_factor;

    // Additional metrics
    count_type population;
    bool below_detection_limit;
};

// ----------------------------------------------------------------------------
// GLOBAL WASTEWATER CONFIGURATION
// ----------------------------------------------------------------------------
struct wastewater_config {
    bool wastewater_enabled = false;
    viral_shedding_params shedding_params;

    // Conversion factor: viral load to estimated infections
    // Calibrated based on: concentration (copies/L) -> infections per 100k
    double viral_load_to_infection_factor = 1.0e-6;

    // Minimum number of infections for reliable WW detection
    count_type minimum_infections_for_detection = 5;

    // Time lag for reporting (wastewater is 4-7 days ahead of clinical)
    int wastewater_lead_time_days = 5;
};

// ----------------------------------------------------------------------------
// FUNCTION DECLARATIONS
// ----------------------------------------------------------------------------

// Initialize sewersheds from configuration
void initialize_sewersheds(std::vector<sewershed>& sewersheds,
                          const std::vector<house>& homes,
                          const std::vector<community>& communities,
                          const wastewater_config& config);

// Calculate viral shedding for an individual based on days since infection
double calculate_viral_shedding(double days_since_infection,
                               bool is_symptomatic,
                               const viral_shedding_params& params);

// Update wastewater surveillance (called each simulation step)
void update_wastewater_surveillance(std::vector<sewershed>& sewersheds,
                                   const std::vector<agent>& agents,
                                   const std::vector<house>& homes,
                                   int current_time,
                                   const wastewater_config& config);

// Aggregate clinical cases by sewershed (Feature 2)
void aggregate_clinical_cases_by_sewershed(std::vector<sewershed>& sewersheds,
                                          const std::vector<agent>& agents,
                                          const std::vector<house>& homes,
                                          int current_time);

// Estimate true infections from wastewater and calculate detection ratio (Feature 3)
void estimate_detection_ratio(std::vector<sewershed>& sewersheds,
                             const wastewater_config& config);

// Collect wastewater sample (respects sampling frequency)
bool should_collect_sample(const sewershed& ss, int current_day);

// Add measurement noise to simulate real sampling
double add_measurement_noise(double concentration, double noise_std_fraction = 0.15);

// Get ground truth infection count in sewershed (for validation)
count_type get_true_active_infections_in_sewershed(const sewershed& ss,
                                                   const std::vector<agent>& agents,
                                                   const std::vector<house>& homes);

// Generate wastewater sample record for output
wastewater_sample generate_sample_record(const sewershed& ss,
                                        int current_day,
                                        const std::vector<agent>& agents,
                                        const std::vector<house>& homes);

// Output functions
void output_wastewater_samples(const std::string& output_path,
                              const std::vector<wastewater_sample>& samples);

void output_sewershed_timeseries(const std::string& output_path,
                                const std::vector<std::vector<wastewater_sample>>& timeseries);

#endif
